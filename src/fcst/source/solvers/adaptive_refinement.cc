//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: adaptive_refinement.cc
//    - Description: Child of ApplicationWrapper used to implement adaptive refinement
//    - Developers: M. Secanell, Valentin N. Zingan, Aslan Kosakian
//
//---------------------------------------------------------------------------

#include "solvers/adaptive_refinement.h"

//---------------------------------------------------------------------------
template <int dim>
FuelCell::ApplicationCore::AdaptiveRefinement<dim>::AdaptiveRefinement
     (FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>& app_lin,
      FuelCell::ApplicationCore::ApplicationWrapper& app,
      const FuelCell::ApplicationCore::FEVector& solution)
  :
  FuelCell::ApplicationCore::ApplicationWrapper(app),
  app_linear(&app_lin),
  app(&app),
  solution(solution)

{
  FcstUtilities::log <<"->AdaptiveRefinement" << std::endl;
}

//---------------------------------------------------------------------------
template <int dim>
FuelCell::ApplicationCore::AdaptiveRefinement<dim>::~AdaptiveRefinement()
{}

//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::ApplicationCore::AdaptiveRefinement<dim>::declare_parameters(ParameterHandler& param) const
{
    // Initialize to parameters needed for evaluate (note: this is provisional)
    param.enter_subsection("Adaptive refinement");
    {
        param.declare_entry ("Number of Refinements",
                             "1",
                             Patterns::Integer(),
                             "This parameter is used to define the number of times the mesh will be adaptively refined. \n"
                             "The minimum value is one, i.e., only the original mesh is solved. At each adaptive refinement level, \n"
                             "either all the cell (global) or 30% of the cells with largest error (computed using an error estimator) are split into four. \n"
                             "The process is repeated at each refinement level. ");
        param.declare_entry("Output initial mesh",
                            "false",
                            Patterns::Bool(),
                            "Set flag to true if you want to output a EPS figure of the initial mesh "
                            "using the value in file initial mesh");
        param.declare_entry("Output initial mesh filename",
                            "initial_mesh",
                            Patterns::Anything(),
                            "File where the initial mesh will be output");
        param.declare_entry("Output intermediate solutions",
                            "false",
                            Patterns::Bool(),
                            "Set flag to true if you would like the solution at each grid refinement to be output. \n"
                            "Please note that outputting the solution is time consuming.");
        param.declare_entry("Output intermediate responses",
                            "false",
                            Patterns::Bool(),
                            "Compute the functionals in \texttt{Output variables} at each grid refinement. \n"
                            "Use this option if you want to perform a grid refinement study. Please note however that computing the functionals is time consuming.");
        param.declare_entry("Output final solution",
                            "true",
                            Patterns::Bool(),
                            "Output the final solution to a file.");
        param.declare_entry("Output final solution filename",
                            "fuel_cell_solution_DataFile",
                            Patterns::Anything(),
                            "File where the final solution will be output");
        param.declare_entry("Compute errors and convergence rates",
                            "false",
                             Patterns::Bool(),
                            "Internal option for developers. Set this value always to false.");
        param.declare_entry("Use nonlinear solver for linear problem",
                            "false",
                             Patterns::Bool(),
                            "Internal option for developers. Set this value always to false.");
        
    }
    param.leave_subsection();

//     param.enter_subsection("Application");
//     {
//         param.declare_entry ("Analytical gradients",
//                              "false",
//                              Patterns::Bool(),
//                              "(Only for optimization) Would you like to compute the analytical gradients of each output variable?");
//     }
//     param.leave_subsection();

    //--- Create entries on ParameterHandler object:
    app->declare_parameters(param);

}

//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::ApplicationCore::AdaptiveRefinement<dim>::initialize(ParameterHandler& param)
{
  // Initialize application using ParameterHandler object:
    param.enter_subsection("Adaptive refinement");
    {
        this->n_ref = param.get_integer("Number of Refinements");
        app_linear->get_data()->enter("n_ref",this->n_ref);
        this->data->enter("n_ref",this->n_ref); // for transient solver
        
        output_initial_mesh = param.get_bool("Output initial mesh");
        filename_initial_mesh = param.get("Output initial mesh filename");

        output_intermediate_sol = param.get_bool("Output intermediate solutions");
        output_final_sol = param.get_bool("Output final solution");
        filename_final_sol = param.get("Output final solution filename");

        output_intermediate_resp = param.get_bool("Output intermediate responses");

        L1_L2_error_and_convergence_rate    = param.get_bool("Compute errors and convergence rates");
        nonlinear_solver_for_linear_problem = param.get_bool("Use nonlinear solver for linear problem");

    }
    param.leave_subsection();
//     param.enter_subsection("Application");
//     {
//         gradients = param.get_bool ("Analytical gradients");
//     }
//     param.leave_subsection();
    gradients = false; // Analytical gradients are currently not working.
    
    app->initialize(param);
    
    this->data->enter_flag("perform_ar",true);
    this->data->enter_flag("cycle_ran_once",false);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
void
FuelCell::ApplicationCore::AdaptiveRefinement<dim>::solve(const std::string param_file,
                                                          ParameterHandler& param,
                                                          bool              solution_component_changes_between_data_files)
{
    // --- 1) we first declare parameters ---
    this->declare_parameters(param);
    
    // --- 2) we then read parameter file ---
    FcstUtilities::read_parameter_files(param,
                                        param_file);
    
    // --- 3) we initialize all related data ---
    this->initialize(param);
    
    // --- 4) we delegate the code flow to the next function ---
    this->run_app(solution_component_changes_between_data_files);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
void
FuelCell::ApplicationCore::AdaptiveRefinement<dim>::run_app(bool solution_component_changes_between_data_files)
{
    std::vector<double> resp(0);
    resp.clear();
    resp.resize(app_linear->get_n_resp());
    
    run_app (resp, solution_component_changes_between_data_files);
}

//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::ApplicationCore::AdaptiveRefinement<dim>::run_app( std::vector<double>& resp,
                                                             bool                 solution_component_changes_between_data_files )
{
    std::vector<std::vector<double> > dresp_dl(0);
    run_app (resp, dresp_dl, solution_component_changes_between_data_files);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
void
FuelCell::ApplicationCore::AdaptiveRefinement<dim>::run_app( std::vector<double>&               resp,
                                                             std::vector<std::vector<double> >& dresp_dl,
                                                             bool                               solution_component_changes_between_data_files )
{
    static unsigned int no_data_file = 1;
    
    // --- output grid ---
    if(output_initial_mesh)
        app_linear->grid_out(filename_initial_mesh);
    
    // --- vector to store solution and residual ---
    FuelCell::ApplicationCore::FEVectors vectors;
    vectors.add_vector(solution, "Solution");
    
    // --- solution will be transferred to refined meshes ---
    app_linear->add_vector_for_transfer(&solution);
        
    // --- a copy of the triangulation object to store the information on the coarse mesh ---
    app_linear->store_triangulation(coarse_triangulation);
    app_linear->initialize_solution(coarse_solution);
    
    /**
     * Variable @simulation_is_not_over is used to control the while loop below. It is only used
     * in transient simulations and is set to false as soon they reach the final time layer.
     */
    bool simulation_is_not_over=true;

    /**
     * Vector used to collect all responses (which are vectors as well). It is sent to the
     * data object so it can be accessed from the solver. For example, it is accessed from the
     * @TransientBase class for the responses output at the end of a transient simulation.
     */
    std::vector< std::vector<double> > all_responses;
    
    /*
     * This variable is sent to the data object to let the solver know if the vector of responses
     * @all_responses is ready to be collected. See @ThetaScheme for the usage example.
     */
    this->data->enter_flag("responses_ready",false);
    
    FEVector u_old;
    if (this->data->get_transient_solver() != FuelCell::ApplicationCore::TransientSolver::NONE)
    {
        u_old.reinit(solution);
        u_old=solution;
        app_linear->add_vector_for_transfer(&u_old); 
    }
    unsigned int n_ref_copy = this->n_ref;
    
    while (simulation_is_not_over)
    {     
        // Send n_ref to the data object
        this->data->enter("n_ref",this->n_ref);        
        
        // Variable @cycle_ran_once might be changed by the transient solver. See @ThetaScheme for usage.
        if (this->data->flag("cycle_ran_once"))
            this->n_ref=1;
        else
            this->n_ref=n_ref_copy;
            
        // --- adaptive refinement loop
        for(unsigned int cycle = 0; cycle < this->n_ref; ++cycle)
        {           
            FcstUtilities::log << "Cycle = " << cycle << " : " << std::endl;
            
            bool last_cycle = false;
            this->data->enter_flag("last_cycle",last_cycle);
            this->data->enter_flag("go_to_next_time_layer",false);
            
            std::ostringstream grid_streamOut;
            std::string        grid_filename;            
            grid_streamOut << cycle;
            grid_filename = "fuel_cell_grid_cycle_" + grid_streamOut.str();
            
            std::ostringstream sol_streamOut;
            std::string sol_filename;            
            sol_streamOut << std::setfill('0') << std::setw(5) << no_data_file;
            sol_filename = filename_final_sol+"_Cycle_" + grid_streamOut.str() + "_Sol_" + sol_streamOut.str();
            //sol_filename = "fuel_cell_solution_DataFile_" + sol_streamOut.str() + "_Cycle_" + grid_streamOut.str();
            
            app_linear->get_data()->enter("Refinement", cycle);
                      
            if (cycle == 0)
            {
                // --- initialize initial guess ---
                app_linear->initialize_solution(solution);
                if (this->data->get_transient_solver() != FuelCell::ApplicationCore::TransientSolver::NONE)
                {
                    u_old=solution;  // u_old is transferred to new meshes with solution, so there is no need to reinitialize it
                    this->data->enter_flag("go_to_next_time_layer",true);
                }                
            }
            else
            {
                // Estimate error per cell for the subsequent adaptive refinement.
                // The default implementation is supposed to be given in DoFApplicatioin<dim> class.
                // The custom implementation is usually given in the actual linear application class
                // using KellyErrorEstimator<dim>::estimate() function.
                FcstUtilities::log << "Entering estimate..."<<std::endl;
                app->estimate(vectors);                
                FcstUtilities::log << "Exited estimate"<<std::endl;
                
                if( nonlinear_solver_for_linear_problem )
                    app_linear->initialize_solution(solution);
                
                // Remesh and start with an initial "coarse" solution. Note that
                // in order to transfer the solution, add_vector_for_transfer needs
                // to be initialized (see before loop).
                if (this->data->flag("perform_ar"))
                    app->remesh();
                
                FcstUtilities::log << "Exited remesh"<<std::endl;  
            }
            
            // -- solve the nonlinear system of equations ---
            try
            {
                FcstUtilities::log << "Solving..." << std::endl;

                if( cycle == this->n_ref - 1 ) {
                    last_cycle = true;
                    this->data->enter_flag("last_cycle",true);
                }     
                
                if (this->data->get_transient_solver() != FuelCell::ApplicationCore::TransientSolver::NONE)
                    this->data->enter("old_transient_solution",u_old);
                app->solve(solution, vectors);
                
                FcstUtilities::log << "Solution: " << solution.n_blocks() << ':';
                for(unsigned int b = 0; b < solution.n_blocks(); ++b)
                    FcstUtilities::log << ' ' << solution.block(b).size();
                FcstUtilities::log << std::endl;
                
                if(output_intermediate_resp || last_cycle)
                {
                    //--- Resize the response vector and copy to it all responses
                    app->evaluate(vectors);
                    //Provisional while evaluate() is not changed
                    app_linear->responses (resp,vectors);
                    if (this->data->get_transient_solver() != FuelCell::ApplicationCore::TransientSolver::NONE)
                    {
                        // Add vector of responses to the end of all_responses
                        all_responses.push_back(resp);
                    }
                    //FcstUtilities::log<< "Current density : "<<-resp[0]<<"A/cm^2"<<std::endl;
                }
                
                if(output_intermediate_sol || (last_cycle && output_final_sol) )
                {
                    // Only for debug:
                    // app_linear->grid_out(grid_filename);
                    app->data_out(sol_filename,vectors);
                }
                
                if(last_cycle)
                {
                    app_linear->transfer_solution_to_coarse_mesh(coarse_triangulation,
                                                                 coarse_solution,
                                                                 solution);
                }
                
                
                if( L1_L2_error_and_convergence_rate )
                    app_linear->compute_L1_L2_error_and_convergence_rate(solution, cycle, convergence_tables);
                
                // -- Update convergence flag:
                app->get_data()->enter_flag("Newton convergence", true);
            }
            catch(const std::exception& e)
            {
                app->get_data()->enter_flag("Newton convergence", false);
                FcstUtilities::log << e.what() << std::endl;
            }
            
            FcstUtilities::log.pop();
            
        } // ARM LOOP
        
        if( L1_L2_error_and_convergence_rate )
            print_convergence_table();
        
        
        if (gradients == true)
        {
            try {
                dresp_dl.clear();
                dresp_dl.resize(app_linear->get_n_resp(), std::vector<double> (app_linear->get_n_dvar()));
                app_linear->solve_direct(dresp_dl,
                                         vectors);
            }
            catch (const std::exception& e)
            {
                FcstUtilities::log << e.what() << std::endl;
            }
        }
        
        std::vector< std::string > name_resp = app_linear->get_name_responses();
        
        for (unsigned int i=0; i<resp.size(); ++i)
            FcstUtilities::log<<"Response "<<name_resp[i].c_str()<<" is : "<<resp[i]<<std::endl;
        
        /**
         * If the model is steady-state, then set @simulation_is_not_over to false to stop the loop that is
         * designed for the transient problems. Otherwise, receive the value of that flag from the transient solver.
         */
        if (this->data->get_transient_solver() == FuelCell::ApplicationCore::TransientSolver::NONE)
            simulation_is_not_over=false;
        else
            simulation_is_not_over=this->data->flag("simulation_is_not_over");        
    
    no_data_file++;
    } // while simulation_is_not_over
    
    // Send all responses to the data object. Used in transient applications only.
    if (this->data->get_transient_solver() != FuelCell::ApplicationCore::TransientSolver::NONE)
    {
        this->data->enter("all_responses",all_responses);
        this->data->enter_flag("responses_ready",true);
        // Will not solve anything but will output total responses
        app->solve(solution, vectors);
    }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
void
FuelCell::ApplicationCore::AdaptiveRefinement<dim>::test_derivatives(const std::string input_file,
                                                    const std::string dvar,
                                                    const double value,
                                                    std::vector<double>& resp,
                                                    std::vector<std::vector<double> >& dresp,
                                                    const bool grad)
{
    //--- Declare parameters:
    ParameterHandler param;
    this->declare_parameters(param);

    //--- Read data from file to ParamterHandler object:
    FcstUtilities::read_parameter_files(param,
                                        input_file);
    //--- Modify data just read:
    std::vector<std::string> name_dvar(1);
    std::vector<double> value_dvar(1);
    name_dvar[0] = dvar;
    value_dvar[0] = value;

    //-- Set design variables to the values given by DAKOTA:
    for (unsigned int i = 0; i < value_dvar.size(); i++)
    {
        FcstUtilities::modify_parameter_file(name_dvar[i], value_dvar[i], param); // parameter to change and its value
    }

    //--- Initialize application:
    this->initialize(param);

    gradients = grad;

    // -- Solve the system of governing equations for the given set of parameters:
    run_app(resp, dresp);
}

//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::ApplicationCore::AdaptiveRefinement<dim>::print_parameters() const
{
    // Initialize Parameter handler
    ParameterHandler param;
    this->declare_parameters(param);

    std::ofstream marc("parameters_sample.prm");
    param.print_parameters(marc,
                           ParameterHandler::Text);
}

//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::ApplicationCore::AdaptiveRefinement<dim>::print_convergence_table()
{
    for( unsigned int i = 0; i < convergence_tables.size(); ++i )
    {
        // precision        
        convergence_tables[i].set_precision( "L1-error", 3);
        convergence_tables[i].set_scientific("L1-error", true);
        
        convergence_tables[i].set_precision( "L2-error", 3);
        convergence_tables[i].set_scientific("L2-error", true);
        
        // convergence rates        
        convergence_tables[i].evaluate_convergence_rates("L1-error", ConvergenceTable::reduction_rate_log2);
        convergence_tables[i].evaluate_convergence_rates("L2-error", ConvergenceTable::reduction_rate_log2);
        
        // output -> console        
        std::cout << std::endl;
        std::cout << "convergence table # " << i+1 << ":" << std::endl;
        
        std::cout << std::endl;
        convergence_tables[i].write_text(std::cout);
        
        // output -> LaTeX        
        convergence_tables[i].set_tex_caption("cycle"    , "cycle"      );
        convergence_tables[i].set_tex_caption("cells"    , "cells"      );
        convergence_tables[i].set_tex_caption("dofs"     , "dofs"       );
        convergence_tables[i].set_tex_caption("L1-error" , "$L_1$-error");
        convergence_tables[i].set_tex_caption("L2-error" , "$L_2$-error");
        
        convergence_tables[i].set_tex_format("cycle"    , "c");
        convergence_tables[i].set_tex_format("cells"    , "r");
        convergence_tables[i].set_tex_format("dofs"     , "r");
        convergence_tables[i].set_tex_format("L1-error" , "c");
        convergence_tables[i].set_tex_format("L2-error" , "c");
        
        std::ostringstream streamOut;
        streamOut << i+1;
        std::string   error_table_filename = "error_" + streamOut.str() + ".tex";
        std::ofstream error_table_file( error_table_filename.c_str() );
        convergence_tables[i].write_tex(error_table_file);
    }
    
}


//---------------------------------------------------------------------------
// Explicit instantations
template class FuelCell::ApplicationCore::AdaptiveRefinement<deal_II_dimension>;
