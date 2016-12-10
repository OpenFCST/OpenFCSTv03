//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: dakota_direct_interface.cc
//    - Description: Class used to create an interface between FCST and DAKOTA
//    - Developers: Malte Krack, Peter Dobson, Marc Secanell <secanell@ualberta.ca>
//    - $Id: dakota_direct_interface.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <dakota_direct_interface.h>


#ifdef _WITH_DAKOTA

using namespace boost;

using namespace SIM;

//---------------------------------------------------------------------------
/*template <int dim>
DakotaDirectInterface<dim>::DakotaDirectInterface(boost::shared_ptr<const Dakota::ProblemDescDB > global_problem_db)
:
Dakota::DirectApplicInterface(*global_problem_db.get()),
data(shared_ptr<FuelCell::ApplicationCore::ApplicationData>(new FuelCell::ApplicationCore::ApplicationData()))
{
    problem_db = global_problem_db;
}
*/
//---------------------------------------------------------------------------

template <int dim>
DakotaDirectInterface<dim>::DakotaDirectInterface(DakotaApplication &fcst_interface,
                                                  boost::shared_ptr<const Dakota::ProblemDescDB > global_problem_db,
                                                  ParameterHandler& param,
                                                  boost::shared_ptr <FuelCell::ApplicationCore::ApplicationData> data,
                                                  boost::shared_ptr<SimulationSelector<dim> > sim_selector,
                                                  std::string &param_file)
:
fcst_interface(&fcst_interface),
Dakota::DirectApplicInterface(*global_problem_db.get()),
param(&param),
data(data),
sim_selector(sim_selector),
simulator_parameter_file_name(param_file)
{
    FcstUtilities::log << "Creating Direct DAKOTA Interface " << std::endl;
}


//---------------------------------------------------------------------------

template <int dim>
int DakotaDirectInterface<dim>::derived_map_ac(const Dakota::String& ac_name)
{
    FcstUtilities::log.push("OPT");
    // ac_name ought to be checked in the long run
    eval_timer.restart();

    shared_ptr<FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim> > app_linear;
    shared_ptr<FuelCell::ApplicationCore::ApplicationWrapper> newton;
    shared_ptr<FuelCell::ApplicationCore::AdaptiveRefinement<dim> > solver;

    app_linear = sim_selector->select_application();
    newton = sim_selector->select_solver(app_linear.get());
    solver = sim_selector->select_solver_method(app_linear.get(), newton.get());

    //solver->solve(simulator_parameter_file_name, *this->param);
    //this->param->clear();
    //FuelCell::ApplicationCore::declare_parameter_files(*this->param);
    solver->declare_parameters(*this->param);
    FcstUtilities::read_parameter_files(*this->param,
                                         simulator_parameter_file_name);

    // Adopt DAKOTA's design variable values, i.e. set_parameters()
    fcst_interface->synchronize_variables(n_dvar, n_resp, name_design_var, name_responses);
    // Read modified parameters:
    dakota_adopt_parameters(app_linear);
    // Initialize used to create objects, etc.
    solver->initialize(*this->param);


    // -- Allocate space for responses and their derivatives:
    std::vector<double> f(n_resp);
    std::vector<std::vector<double> > df_dl(n_resp,
                                            std::vector<double> (n_dvar));
    solver->run_app(f, df_dl);

    eval_timer.stop();
    FcstUtilities::log << "FuelCell function evaluation performed in " << eval_timer() << " seconds." << std::endl;

    dakota_assign_results(f, df_dl);
    FcstUtilities::log << "Results assigned to DAKOTA variables" << std::endl;
    FcstUtilities::log.pop();
    return 0;

}

//---------------------------------------------------------------------------
template <int dim>
void DakotaDirectInterface<dim>::dakota_adopt_parameters(shared_ptr<FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim> > app_linear)
{
    // Convert DAKOTA variable types into STL double and string vectors (easier way available??)
    std::vector<double> x(numACV);
    std::vector<std::string> name_x(numACV);
    for (unsigned int i = 0; i < numVars; i++)
    {
        x[i] = xC[i];
        name_x[i] = xCLabels[i];
    }

    app_linear->set_optimization_parameters(n_dvar, n_resp, name_design_var, name_responses);

    //-- Set design variables to the values given by DAKOTA:
    for (unsigned int i = 0; i < n_dvar; i++)
    {
        FcstUtilities::modify_parameter_file(name_design_var[i], x[i], *this->param); // parameter to change and its value (Note name_x does not mean anything)
    }

}

//---------------------------------------------------------------------------
template <int dim>
void
DakotaDirectInterface<dim>::dakota_assign_results(const std::vector<double>& responses,
                                                       const std::vector<std::vector<double> >& dresponses_dl)
{
    for (int i = 0; (unsigned int) i < numFns; ++i)
    {
        if (directFnASV[i] == 1 ||
            directFnASV[i] == 3)
        {
            fnVals[i] = responses[i];
        }
        if (directFnASV[i] == 2 ||
            directFnASV[i] == 3)
        {
            for (int j = 0; (unsigned int) j < dresponses_dl[i].size(); ++j)
            {
                fnGrads[i][j] = dresponses_dl[i][j];
            }
        }
    }
}

//---------------------------------------------------------------------------

template <int dim>
void DakotaDirectInterface<dim>::dakota_assign_results(const std::vector<double>& responses)
{
    for (int i = 0; (unsigned int) i < numFns; ++i)
    {
        if (directFnASV[i] == 1 || directFnASV[i] == 3)
        {
            fnVals[i] = responses[i];
        }
    }
}


//---------------------------------------------------------------------------
//----------------- DAKOTA LEAST-SQUARES INTERFACE --------------------------
//---------------------------------------------------------------------------

template <int dim>
DakotaLeastSquaresInterface<dim>::DakotaLeastSquaresInterface(DakotaApplication &fcst_interface,
                                                              boost::shared_ptr<const Dakota::ProblemDescDB > global_problem_db,
                                                              ParameterHandler& param,
                                                              boost::shared_ptr <FuelCell::ApplicationCore::ApplicationData> data,
                                                              boost::shared_ptr<SimulationSelector<dim> > sim_selector,
                                                              std::string &param_file)
:
DakotaDirectInterface<dim>::DakotaDirectInterface(fcst_interface, global_problem_db, param, data, sim_selector, param_file),
data(data)
{
    FcstUtilities::log << "\t for Parameter Estimation" << std::endl;
}

//---------------------------------------------------------------------------

template <int dim>
void DakotaLeastSquaresInterface<dim>::_initialize(ParameterHandler& param)
{
    param.enter_subsection("Optimization Parameters");
    {
        param.enter_subsection("Optimization Program Options");
        {
            NLS_data_file = param.get("NLS data file");
            NLS_residual_option = param.get("NLS residual option");
            NLS_residual_value = param.get("NLS residual value");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------

template <int dim>
int DakotaLeastSquaresInterface<dim>::derived_map_ac(const Dakota::String& ac_name)
{
    FcstUtilities::log.push("NLS");
    // ac_name ought to be checked in the long run
    nls_timer.restart();

    _initialize(*this->param);

    FcstUtilities::log << "Reading experimental data file : " << NLS_data_file << std::endl;
    ExperimentalData NLS_data(NLS_data_file);

    //get the experimental results vector
    NLS_data.extract_vector("Current", experimental_current);
    //output the remaining variable names and matrix of values to the screen
    NLS_data.print_data();

    //get the operating condition names and matrix
    NLS_data.get_experimental_values(OC_names, OC_values);

    NLS_residual.clear();
    NLS_residual.resize(experimental_current.size());

    //loop over number of data points
    for (unsigned int i = 0; i < experimental_current.size(); i++)
    {
        this->eval_timer.restart();
        FcstUtilities::log << "\n\nEvaluating Least-Square Term: " << i + 1 << " of " << experimental_current.size() << std::endl;
        //                for (unsigned int j = 0; j < OC_names.size(); j++)
        //                {
        //                     FcstUtilities::log << OC_names[j] << ": " << OC_values[i][j] << "\n";
        //                }

        param_names.clear();
        param_names.resize(OC_names.size());
        param_names = OC_names;

        param_values.clear();
        param_values.resize(OC_names.size());
        param_values = OC_values[i];

        shared_ptr<FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim> > app_linear;
        //shared_ptr<FuelCell::ApplicationCore::newtonBase> newton;
        shared_ptr<FuelCell::ApplicationCore::ApplicationWrapper> newton;
        shared_ptr<FuelCell::ApplicationCore::AdaptiveRefinement<dim> > solver;

        app_linear = this->sim_selector->select_application();
        newton = this->sim_selector->select_solver(app_linear.get());
        solver = this->sim_selector->select_solver_method(app_linear.get(), newton.get());

        //solver->solve(simulator_parameter_file_name, *this->param);
        //this->param->clear();
        //FuelCell::ApplicationCore::declare_parameter_files(*this->param);
        solver->declare_parameters(*this->param);
        FcstUtilities::read_parameter_files(*this->param,
                                             this->simulator_parameter_file_name);

        // Adopt DAKOTA's design variable values, i.e. set_parameters()
        this->fcst_interface->synchronize_variables(this->n_dvar, this->n_resp, this->name_design_var, this->name_responses);
        this->dakota_adopt_parameters(app_linear);

        solver->initialize(*this->param);

        // -- Allocate space for responses and their derivatives:
        std::vector<double> f(app_linear->get_n_resp());
        std::vector<std::vector<double> > df_dl(app_linear->get_n_resp(),
                                                std::vector<double> (app_linear->get_n_dvar()));

        solver->run_app(f,
                           df_dl);

        std::vector<std::string> name_resp = app_linear->get_name_responses();
        double numerical_current= 0.0;
        for (unsigned int r = 0; r<f.size(); ++r)
        {
            if (name_resp[r] == "current")
            {
                numerical_current = -f[r];
                break;
            }
            FcstUtilities::log << "FuelCell DakotaLeastSquaresInterface requires 'current' for a response." << std::endl;
            FcstUtilities::log << "Quitting Program";
            abort();
        }

        //evaluate the residual - weighted to percent deviation
        if (NLS_residual_option == "absolute")
        {
            NLS_residual[i] = (numerical_current - (experimental_current[i]));
            FcstUtilities::log << "Least-Square Residual: " << NLS_residual[i] << " A/cm^2" << std::endl;
        }
        else if (NLS_residual_option == "weighted")
        {
            // (I-I_e)/I_e  OR  I/I_e - 1
            NLS_residual[i] = (numerical_current / (experimental_current[i])) - 1.0;
            FcstUtilities::log << "Weighted Least-Square Residual: " << NLS_residual[i] << std::endl;
        }
        else if (NLS_residual_option == "sensitivity")
        {
            // (I-I_e)/I_e  OR  I/I_e - 1
            NLS_residual[i] = numerical_current;
            FcstUtilities::log << "Current Density for Sensitivity: " << NLS_residual[i] << std::endl;
        }

        this->eval_timer.stop();
        FcstUtilities::log << "FuelCell function evaluation performed in " << this->eval_timer() << " seconds." << std::endl;

    } //end of NLS 'for' loop

    nls_timer.stop();
    FcstUtilities::log << "FuelCell NLS Iteration performed in " << nls_timer() << " seconds." << std::endl;

    if (NLS_residual_value == "norm")
    {
        std::vector<double> NLS_residual_norm(1,0.0);
        for (unsigned int i=0; i<NLS_residual.size();++i)
            NLS_residual_norm[0] += pow(NLS_residual[i],2);
        NLS_residual_norm[0] = sqrt(NLS_residual_norm[0]);
        FcstUtilities::log << "Least-Square Residual Norm : " << NLS_residual_norm[0] << " A/cm^2" << std::endl;
        this->dakota_assign_results(NLS_residual_norm);
    }
    else
        this->dakota_assign_results(NLS_residual);

    FcstUtilities::log << "Results assigned to DAKOTA variables" << std::endl;
    NLS_residual.clear();
    FcstUtilities::log.pop();

    return 0;

}

//---------------------------------------------------------------------------
// Explicit instantations
template class SIM::DakotaDirectInterface<deal_II_dimension>;
template class SIM::DakotaLeastSquaresInterface<deal_II_dimension>;


#endif  //Dakota