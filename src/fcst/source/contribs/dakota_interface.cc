//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: dakota_interface.cc
//    - Description: Used to read in an array of experimental data with column headers.
//                   For use with parameter estimation or generating polarization curves.
//    - Developers: Marc Secanell <secanell@ualberta.ca>
//    - $Id: dakota_interface.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <dakota_interface.h>

using namespace SIM;

//---------------------------------------------------------------------------
template <int dim>
DakotaInterface<dim>::DakotaInterface (const std::string input_file,
                                       ParameterHandler &param,
                                       const std::string dakota_parameters,
                                       const std::string dakota_results,
                                       FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>& app_lin,
                                       FuelCell::ApplicationCore::ApplicationWrapper& app)
:
input_file(input_file),
dakota_parameters(dakota_parameters),
dakota_results(dakota_results),
app_linear(&app_lin),
app(&app),
param(&param)
{
    FcstUtilities::log << "Creating standard DAKOTA Interface " << std::endl;
    FcstUtilities::log << "Dakota will read paramters in: "<< dakota_parameters << std::endl;
    FcstUtilities::log << "Analysis returns to: "<< dakota_results << std::endl;
    FcstUtilities::log.pop();
}

//---------------------------------------------------------------------------
template <int dim>
void
DakotaInterface<dim>::declare_parameters(ParameterHandler& param)
{
    param.enter_subsection("Grid generation");
    {
        param.declare_entry ("Number of Refinements",
                             "1",
                             Patterns::Integer());
    }
    param.leave_subsection();
    param.enter_subsection("Application");
    {
        param.declare_entry("Dakota version",
                            "5.1",
                            Patterns::Anything(),
                            "Version of Dakota used for the optimization.  v3.3 no longer supported");
    }
    param.leave_subsection();

    //-- Create entries on ParameterHandler object:
    app->declare_parameters(param);
}

//---------------------------------------------------------------------------
template <int dim>
void
DakotaInterface<dim>::initialize(ParameterHandler& param)
{
    param.enter_subsection("Grid generation");
    {
        n_ref = param.get_integer("Number of Refinements");
    }
    param.leave_subsection();
    param.enter_subsection("Application");
    {
        dakota_version = param.get("Dakota version");
    }
    param.leave_subsection();

}

//---------------------------------------------------------------------------
template <int dim>
void
DakotaInterface<dim>::run()
{
    // Create entries on ParameterHandler object:
    this->declare_parameters(*param);

    // Read data from file to ParamterHandler object:
    FcstUtilities::read_parameter_files(*param,
                                         input_file);

    // Initialize application:
    this->initialize(*param);

    // Read data from DAKOTA file and modify the necessary parameters
    //(i.e. the design variables)
    DakotaReadIn(dakota_parameters,
                 *app_linear,
                 *param);

    // Initialize application using ParameterHandler object:
    app->initialize(*param);

    // Vector to store solution:
    //BlockVector<double> sol;
    FuelCell::ApplicationCore::FEVector rhs;
    FuelCell::ApplicationCore::FEVector sol;
    FuelCell::ApplicationCore::FEVectors vectors;
    vectors.add_vector(sol,"Solution");
    // Tell app_linear that the solution will be transferred to refined meshes
    app_linear->add_vector_for_transfer(&sol);

    // -- Start adaptive refinement loop:
    for (unsigned int cycle=0; cycle<n_ref; ++cycle)
    {
        FcstUtilities::log<< "Cycle "<<cycle<<" : "<<std::endl;
        std::ostringstream streamOut;
        streamOut << cycle;
        std::string filename = "fuel_cell-grid-cycle"+streamOut.str();
        std::string sol_filename = "fuel_cell-sol-cycle"+streamOut.str();

        app_linear->get_data()->enter("Refinement", cycle); //.enter_data_scalar("Refinement", cycle);

        if (cycle == 0)
        {
            // initialize vectors with new grid
            app_linear->initialize_solution(sol);
            app->data_out("initial_sol",
            vectors);
        }
        else
        {
            app->estimate(vectors);
            app->remesh();
        }

        // -- Solve the nonlinear system of equations:
        try
        {
            app->solve(sol, vectors);
            app->grid_out(filename);
            FcstUtilities::log<< "Current density : "<<app->evaluate(vectors)<<"A/cm^2"<<std::endl;
            app->data_out(sol_filename,
                          vectors);
        }
        catch (const std::exception& e)
        {
            FcstUtilities::log << e.what() << std::endl;
        }
    }

    // -- Allocate space for responses and their derivatives:
    std::vector<double> f(app_linear->get_n_resp());
    std::vector<std::vector<double> > df_dl(app_linear->get_n_resp(),
                                            std::vector<double> (app_linear->get_n_dvar()));
    try
    {
        // -- Compute responses
        app_linear->responses (f,
                               vectors);
        bool gradient = false;
        for (unsigned int i=0; i<app_linear->get_n_resp(); ++i)
        {
            if (ASV[i] == 2 ||
                ASV[i] == 3)
                gradient = true;
        }
        if (gradient == true)
            // -- Compute the derivatives of the responses:
            app_linear->solve_direct(df_dl,
                                     vectors);//sol);
    }
    catch (const std::exception& e)
    {
        FcstUtilities::log << e.what() << std::endl;
    }

    DakotaWriteOut(dakota_results,
                   f,
                   df_dl);

}

//---------------------------------------------------------------------------
template <int dim>
void
DakotaInterface<dim>::DakotaReadIn (const std::string parameters,
                                    FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>& app,
                                    ParameterHandler &param)
{
    // This application program reads and writes parameter and response files directly
    //so no input/output filters are needed
    std::ifstream fin(parameters.c_str());
    unsigned int num_vars, num_fns;
    std::string vars_text, fns_text;

    // Get the parameter vector and ignore the labels
    //In Dakota 4.0:
    if (dakota_version == "3.3")
    {
        FcstUtilities::log << "DAKOTA version 3.3 is no longer supported by FCST,"
        << " please upgrade to version 4.0 or greater." << std::endl;
        ExcNotImplemented();
    }
    else // If 4.0 or greater, read the std parameter file type
        fin >> num_vars >> vars_text;

    std::vector<double> x(num_vars);
    std::vector<std::string> name_x(num_vars);
    for (unsigned int i=0; i<num_vars; i++)
    {
        fin >> x[i] >> name_x[i];
        fin.ignore(256, '\n');
    }

    // Get the ASV vector and ignore the labels
    fin >> num_fns >> fns_text;

    ASV.resize(num_fns);
    for (unsigned int i=0; i<num_fns; i++)
    {
        fin >> ASV[i];
        fin.ignore(256, '\n');
    }

    // -- CHECK:
    // -- Read data from the analysis code and make sure that they match:
    //read from parameter file all necessary data:
    unsigned int  n_dvar;
    unsigned int n_resp;
    std::vector<std::string> name_design_var(num_vars);
    std::ofstream fout;

    param.enter_subsection("Optimization Parameters");
    {
        //initialize design variables:
        param.enter_subsection("Design Variables");
        {
            n_dvar = param.get_integer("num_design_variables");
            // Modify size of the nome_response vector:
            name_design_var.clear();
            name_design_var.resize(n_dvar);

            //read data from file
            for (unsigned int i=0; i<n_dvar; ++i)
            {
                // obtain the name of the design variable:
                std::ostringstream streamOut;
                streamOut << i;
                // obtain the name of the design variables:
                std::string name = "DV_" + streamOut.str() + "_name";
                name_design_var[i] = param.get(name.c_str());
            }
        }
        param.leave_subsection();
        param.enter_subsection("Responses");
        {
            n_resp = param.get_integer("num_objectives")
            + param.get_integer("num_nl_constraints")
            + param.get_integer("num_eq_constraints");
        }
        param.leave_subsection();
    }
    param.leave_subsection();

    if (num_vars != n_dvar)
    {
        FcstUtilities::log << "Wrong number of variables for the fuel cell optimization\n";
        exit(-1);
    }
    for (unsigned int i=0; i<n_dvar; ++i)
    {
        if (name_x[i] != name_design_var[i])
        {

            FcstUtilities::log << "The design variable names on the Dakota and analysis file do not match.\n";
            FcstUtilities::log << "Make sure you modify both Dakota input file and analysis files.\n";
            FcstUtilities::log <<"From Dakota file : "<<name_x[i]<<". From input: "<<name_design_var[i];
            exit(-1);
        }
    }
    if (num_fns != n_resp)
    {
        FcstUtilities::log << "Wrong number of functions in fuel cell optimization\n";
        exit(-1);
    }

    //-- Set design variables to the values given by DAKOTA:
    for (unsigned int i = 0; i < n_dvar; i++)
    {
        FcstUtilities::modify_parameter_file(name_design_var[i], x[i], param); // parameter to change and its value (Note name_x does not mean anything)
    }
    // For debugging:
    /*
     * param.print_parameters(FcstUtilities::log,
     * ParameterHandler::Text);
     * abort();
     */
}

//---------------------------------------------------------------------------
template <int dim>
void
DakotaInterface<dim>::DakotaWriteOut (const std::string results,
				 const std::vector<double>& responses,
				 const std::vector<std::vector<double> >& dresponses_dl)
{

  std::ofstream fout(results.c_str());
  if (!fout)
    {
      FcstUtilities::log << "\nError: failure creating " << results << std::endl;
      exit(-1);
    }

  fout.precision(15); // 16 total digits
  fout.setf(std::ios::scientific);
  fout.setf(std::ios::right);

  for (unsigned int i=0; i<responses.size(); ++i)
    {
      if (ASV[i] == 1 ||
	  ASV[i] == 3)
	{
	  // **** Residual :
	  fout << responses[i]  << " f"<<i<<"\n";
	}
    }
  for (unsigned int i=0; i<responses.size(); ++i)
    {
      if (ASV[i] == 2 ||
	  ASV[i] == 3)
	{
	  // **** dR/dl:
	  fout << "[ ";
	  for (unsigned int j=0; j<dresponses_dl[i].size(); ++j)
	    {
	      fout<< dresponses_dl[i][j] <<"   ";
	    }
	  fout<<" ]\n ";
	}
    }

  fout.flush();
  fout.close();

}

//---------------------------------------------------------------------------
// Explicit instantations
template class SIM::DakotaInterface<deal_II_dimension>;
