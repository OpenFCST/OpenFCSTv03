//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: dakota_application.cc
//    - Description: Used to read in an array of experimental data with column headers.
//                   For use with parameter estimation or generating polarization curves.
//    - Developers: Malte Krack, Peter Dobson <pdobson@ualberta.ca>, Marc Secanell <secanell@ualberta.ca>
//    - $Id: dakota_application.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include "contribs/dakota_application.h"

#ifdef _WITH_DAKOTA

using namespace Dakota;
using namespace SIM;


//---------------------------------------------------------------------------

DakotaApplication::DakotaApplication(boost::shared_ptr<Dakota::ProblemDescDB> global_problem_db,
                                     std::string parameter_file)
:
input_file(parameter_file),
problem_db(global_problem_db)
{

}

//---------------------------------------------------------------------------

void DakotaApplication::run()
{
     // lock the database
     problem_db->lock();

     // run the strategy
     selected_strategy.run_strategy();
}

//---------------------------------------------------------------------------

Dakota::RealVector DakotaApplication::vars_results()
{
    // Retrieve the final parameter values and returns the final strategy solution (variables)
    const Variables vars = selected_strategy.variables_results();
    return vars.all_continuous_variables();
}

//---------------------------------------------------------------------------

Dakota::RealVector DakotaApplication::resp_results()
{
    // Retrieve the final response values and returns the final strategy solution (responses)
    const Response resp = selected_strategy.response_results();
    return resp.function_values();
}

//---------------------------------------------------------------------------

void DakotaApplication::declare_parameters(ParameterHandler& param)
{
    param.enter_subsection("Application");
    {
        param.declare_entry("Dakota version",
                            "5.1",
                            Patterns::Anything(),
                            "Version of Dakota used for the optimization.  v3.3 no longer supported");
    }
    param.leave_subsection();

    param.enter_subsection("Grid generation");
    {
        param.declare_entry("Number of Refinements",
                            "1",
                            Patterns::Integer());
    }
    param.leave_subsection();

    param.enter_subsection("Optimization Parameters");
    {
        param.enter_subsection("Design Variables");
        {
            // design variables (allocate space for N design variable names and values):
            param.declare_entry("num_design_variables",
                                "1",
                                Patterns::Integer());
            param.declare_entry("design_variable_scaling_method",
                                "none",
                                Patterns::Selection("none|auto|value"));

            const unsigned int num_design_variables = 100;
            for (unsigned int index=0; index<num_design_variables; ++index)
            {
                std::ostringstream streamOut;
                streamOut << index;
                std::string name = "DV_" + streamOut.str() + "_name";
                std::string ip = "DV_" + streamOut.str() + "_ip";			// _ip = Initial Point
                std::string lb = "DV_" + streamOut.str() + "_lb";			// _lb = Lower Bound
                std::string ub = "DV_" + streamOut.str() + "_ub";			// _ub Upper Bound
                std::string part = "DV_" + streamOut.str() + "_partition";
                std::string scale = "DV_" + streamOut.str() + "_scale";
                std::string scale_method = "DV_" + streamOut.str() + "_scale_method";
                std::string step_size = "DV_" + streamOut.str() + "_step";

                param.declare_entry(name.c_str(),
                                    "no_name",
                                    Patterns::Anything());
                param.declare_entry(ip.c_str(),
                                    "1.0",
                                    Patterns::Double());
                param.declare_entry(lb.c_str(),
                                    "1.0",
                                    Patterns::Double());
                param.declare_entry(ub.c_str(),
                                    "1.0",
                                    Patterns::Double());
                param.declare_entry(part.c_str(),
                                    "2",
                                    Patterns::Integer());
                param.declare_entry(scale.c_str(),
                                    "1.0",
                                    Patterns::Double());
                param.declare_entry(scale_method.c_str(),
                                    "none",
                                    Patterns::Selection("none|auto|value"));
                                    param.declare_entry(step_size.c_str(),
                                                        "1.0e-4",
                                                        Patterns::Double());
                //std::string value = "DV_" + streamOut.str() + "_value";
                //param.declare_entry(value.c_str(),
                         //		      "0.0",
                      //      	      Patterns::Double());
            }
        }
        param.leave_subsection();

        param.enter_subsection("Responses");
        {
            //objective function and constraints
            param.declare_entry("num_objectives",
            "1",
            Patterns::Integer());
            //constraints (allocate space for N design variable names and values):
            param.declare_entry("num_constraints",
            "0",
            Patterns::Integer());
            param.declare_entry("num_nl_constraints",
            "0",
            Patterns::Integer());
            param.declare_entry("num_eq_constraints",
            "0",
            Patterns::Integer());
            const unsigned int num_resp = 100;
            for (unsigned int index=0; index<num_resp; ++index)
            {
                std::ostringstream streamOut;
                streamOut << index;
                std::string name = "RESP_" + streamOut.str() + "_name";
                param.declare_entry(name.c_str(),
                                    "no_name",
                                    Patterns::Anything());
                std::string lb = "RESP_" + streamOut.str() + "_lb";
                param.declare_entry(lb.c_str(),
                                    "-1e31",
                                    Patterns::Double());
                std::string ub = "RESP_" + streamOut.str() + "_ub";
                param.declare_entry(ub.c_str(),
                                    "1e31",
                                    Patterns::Double());
                std::string eq = "RESP_" + streamOut.str() + "_eq";
                param.declare_entry(eq.c_str(),
                                    "0.0",
                                    Patterns::Double());
            }
        }
        param.leave_subsection();

        // optimization solver parameters
        param.enter_subsection("Optimization Program Options");
        {
            param.declare_entry("optimization_program",
                                "DOT",
                                Patterns::Anything());
            param.declare_entry("Use dakota input file",
                                "false",
                                Patterns::Anything());
            param.declare_entry("Dakota_Input_File",
                                "dakota.in",
                                Patterns::Anything());
            param.declare_entry("Optimization strategy",
                                "single_method",
                                Patterns::Selection("single_method|pareto_set|multi_start|hybrid"));
            param.declare_entry("Optimization method",
                                "optpp_q_newton",
                                Patterns::Anything());
            param.declare_entry("Analytical gradients",
                                "false",
                                Patterns::Bool());
            param.declare_entry("Numerical gradients",
                                "false",
                                Patterns::Bool());
            param.declare_entry("Numerical gradient type",
                                "forward",
                                Patterns::Selection("forward|central"));
            param.declare_entry("Analytical hessians",
                                "false",
                                Patterns::Bool());
            param.declare_entry("Use NLS",
                                "false",
                                Patterns::Bool());
            param.declare_entry("NLS data file",
                                "experimental_data.txt",
                                Patterns::Anything());
            param.declare_entry("NLS residual option",
                                "weighted",
                                Patterns::Selection("absolute|weighted|sensitivity"));
            param.declare_entry("NLS residual value",
                                "vector",
                                Patterns::Selection("vector|norm"));
            param.declare_entry("Parametric Study",
                                "None",
                                Patterns::Selection("None|Constant volume fraction"));
            param.declare_entry("Maximum iterations",
                                "100",
                                Patterns::Integer());
            param.declare_entry("Maximum function evaluations",
                                "1000",
                                Patterns::Integer());
            param.declare_entry("Convergence tolerance",
                                "1.0e-4",
                                Patterns::Double());
            param.declare_entry("Constraint tolerance",
                                "1.0e-4",
                                Patterns::Double());
            param.declare_entry("Number of random pareto weights",
                                "0",
                                Patterns::Integer());
            param.declare_entry("Number of random start points",
                                "0",
                                Patterns::Integer());

            // method independent parameters
            param.declare_entry("Scaling",
                                "false",
                                Patterns::Bool());

            param.enter_subsection("Hybrid optimization");
            {
                param.declare_entry("Hybrid optimization strategy",
                                    "sequential",
                                    Patterns::Selection("sequential|embedded|collaborative"));
                param.declare_entry("Number of optimization methods",
                                    "3",
                                    Patterns::Integer());
                for (int i = 0; i < 3; i++)
                {
                    std::ostringstream streamOut;
                    streamOut << i;
                    std::string name = "Optimization method " + streamOut.str();
                    param.declare_entry(name.c_str(),
                                        "optpp_q_newton",
                                        Patterns::Anything());
                    std::string tol = "Convergence tolerance " + streamOut.str();
                    param.declare_entry(tol.c_str(),
                                        "1.0e-4",
                                        Patterns::Double());
                    std::string eval = "Maximum function evaluations " + streamOut.str();
                    param.declare_entry(eval.c_str(),
                                        "1000",
                                        Patterns::Double());
                    std::string iter = "Maximum iterations " + streamOut.str();
                    param.declare_entry(iter.c_str(),
                                        "100",
                                        Patterns::Double());
                }
            }//hybrid opt
            param.leave_subsection();

            param.enter_subsection("CONMIN");
            {

            }
            param.leave_subsection();

            param.enter_subsection("OPT++");
            {

                param.declare_entry("Steplength to boundary",
                "0.9",
                Patterns::Double());
                param.declare_entry("Gradient tolerance",
                "1.0e-4",
                Patterns::Double());
                param.declare_entry("Merit function",
                "argaez_tapia",
                Patterns::Anything());
                param.declare_entry("Centering parameter",
                "0.2",
                Patterns::Double());
                //				param.declare_entry("Search method",
                //								"value_based_line_search",
                //								Patterns::Selection("value_based_line_search|gradient_based_line_search|trust_region|tr_pds"));
            }//OPT++
            param.leave_subsection();

            param.enter_subsection("NL2SOL");
            {
                param.declare_entry("Initial trust radius",
                                    "-1",
                                    Patterns::Double());
                param.declare_entry("Absolute convergence tolerance",
                                    "-1",
                                    Patterns::Double());
                param.declare_entry("Parameter convergence tolerance",
                                    "-1",
                                    Patterns::Double());
                param.declare_entry("False convergence tolerance",
                                    "-1",
                                    Patterns::Double());
            }//NL2SOL
            param.leave_subsection();

            param.enter_subsection("NCSU");
            {
                param.declare_entry("solution target",
                "0",
                Patterns::Double());
                param.declare_entry("final solutions",
                "1",
                Patterns::Double());
                param.declare_entry("min boxsize limit",
                "1e-4",
                Patterns::Double());
                param.declare_entry("volume boxsize limit",
                "1e-8",
                Patterns::Double());
            }//NCSU
            param.leave_subsection();

            param.enter_subsection("coliny_cobyla");
            {
                param.declare_entry("Initial Delta",
                "1.0",
                Patterns::Double());
                param.declare_entry("Threshold Delta",
                "0.0001",
                Patterns::Double());
            }//coliny_pattern_search
            param.leave_subsection();

            param.enter_subsection("coliny_pattern_search");
            {
                param.declare_entry("Initial Delta",
                "1.0",
                Patterns::Double());
                param.declare_entry("Threshold Delta",
                "0.0001",
                Patterns::Double());
                param.declare_entry("Contraction Factor",
                "0.5",
                Patterns::Double());
            }//coliny_pattern_search
            param.leave_subsection();

            param.enter_subsection("asynch_pattern_search");
            {
                param.declare_entry("Initial Delta",
                "1.0",
                Patterns::Double());
                param.declare_entry("Threshold Delta",
                "0.0001",
                Patterns::Double());
                param.declare_entry("Contraction Factor",
                "0.5",
                Patterns::Double());
            }//asynch_pattern_search
            param.leave_subsection();

            param.enter_subsection("coliny_solis_wets");
            {
                param.declare_entry("Initial Delta",
                "1.0",
                Patterns::Double());
                param.declare_entry("Threshold Delta",
                "0.0001",
                Patterns::Double());
                param.declare_entry("Contraction Factor",
                "0.5",
                Patterns::Double());
            }//coliny_solis_wets
            param.leave_subsection();

            param.enter_subsection("moga");
            {
                param.declare_entry("Cross Over Type",
                                    "shuffle_random",
                                    Patterns::Selection("shuffle_random|multi_point_binary|muti_point_parameterized_binary|multi_point_real"));
                param.declare_entry("Initialization Type",
                                    "unique_random",
                                    Patterns::Selection("unique_random|simple_random|flat_file"));
                param.declare_entry("Mutation Type",
                                    "replace_uniform",
                                    Patterns::Selection("replace_uniform|bit_random|offset_cauchy|offset_uniform|offset_normal"));
                param.declare_entry("Number of Offspring",
                                    "2",
                                    Patterns::Double());
                param.declare_entry("Number of Parents",
                                    "2",
                                    Patterns::Double());
                param.declare_entry("Replacement Type",
                                    "below_limit = 6",
                                    Patterns::Anything()); // below_limit = 6 (change value if required)|roulette_wheel|unique_roulette_wheel|elitist
                param.declare_entry("Convergence Type",	// moga
                                    "metric_tracker",
                                    Patterns::Selection("metric_tracker"));
                param.declare_entry("Fitness Type",	// moga
                                    "domination_count",
                                    Patterns::Selection("domination_count|layer_rank"));

            }//moga
            param.leave_subsection();

            param.enter_subsection("soga");
            {
                param.declare_entry("Cross Over Type",
                                    "shuffle_random",
                                    Patterns::Selection("shuffle_random|multi_point_binary|muti_point_parameterized_binary|multi_point_real"));
                param.declare_entry("Initialization Type",
                                    "unique_random",
                                    Patterns::Selection("unique_random|simple_random|flat_file"));
                param.declare_entry("Mutation Type",
                                    "replace_uniform",
                                    Patterns::Selection("replace_uniform|bit_random|offset_cauchy|offset_uniform|offset_normal"));
                param.declare_entry("Number of Offspring",
                                    "2",
                                    Patterns::Double());
                param.declare_entry("Number of Parents",
                                    "2",
                                    Patterns::Double());
                param.declare_entry("Replacement Type",
                                    "elitist",
                                    Patterns::Anything()); // below_limit = 6 (change value if required)|roulette_wheel|unique_roulette_wheel|elitist
                param.declare_entry("Convergence Type",	// soga
                                    "average_fitness_tracker",
                                    Patterns::Selection("average_fitness_tracker|best_fitness_tracker"));
                param.declare_entry("Fitness Type",	// soga
                                    "merit_function",
                                    Patterns::Selection("merit_function"));

            }//soga
            param.leave_subsection();


        }//opt program options
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
void DakotaApplication::initialize(ParameterHandler& param)
{
    param.enter_subsection("Optimization Parameters");
    {
        //initialize design variables:
        param.enter_subsection("Design Variables");
        {
            n_dvar = param.get_integer("num_design_variables");
            name_design_var.clear();
            name_design_var.resize(n_dvar);
            dakota_name_design_var = Dakota::StringArray(n_dvar);
            ip_design_var = Dakota::RealVector(n_dvar);
            lb_design_var = Dakota::RealVector(n_dvar);
            ub_design_var = Dakota::RealVector(n_dvar);
            scale_design_var = Dakota::RealVector(n_dvar);
            scale_types = Dakota::StringArray(n_dvar);
            step_design_var = Dakota::RealVector(n_dvar);
            part_design_var.clear();
            part_design_var.resize(n_dvar);

            for (unsigned int i = 0; i < n_dvar; ++i)
            {
                // obtain the name of the design variable:
                std::ostringstream streamOut;
                streamOut << i;
                // obtain the name of the design variables, initial point, lower, and upper bounds:
                std::string name = "DV_" + streamOut.str() + "_name";
                std::string ip = "DV_" + streamOut.str() + "_ip";
                std::string lb = "DV_" + streamOut.str() + "_lb";
                std::string ub = "DV_" + streamOut.str() + "_ub";
                std::string part = "DV_" + streamOut.str() + "_partition";
                std::string scale = "DV_" + streamOut.str() + "_scale";
                std::string scale_method = "DV_" + streamOut.str() + "_scale_method";
                std::string step_size = "DV_" + streamOut.str() + "_step";
                // set design variable values from parameter file
                name_design_var[i] = param.get(name.c_str());
                dakota_name_design_var[i] = "DV_" + streamOut.str();
                //--
                ip_design_var[i] = param.get_double(ip.c_str());
                lb_design_var[i] = param.get_double(lb.c_str());
                ub_design_var[i] = param.get_double(ub.c_str());
                part_design_var[i] = param.get_double(part.c_str());
                scale_design_var[i] = param.get_double(scale.c_str());
                scale_types[i] = param.get(scale_method.c_str());
                step_design_var[i] = param.get_double(step_size.c_str());
            }
        }//design variables
        param.leave_subsection();

        param.enter_subsection("Responses");
        {
            n_obj = param.get_integer("num_objectives");
            n_con = param.get_integer("num_constraints");
            n_nl_con = param.get_integer("num_nl_constraints");
            n_eq_con = param.get_integer("num_eq_constraints");

            lb_nl_constraint = Dakota::RealVector(n_nl_con);
            ub_nl_constraint = Dakota::RealVector(n_nl_con);
            eq_constraint = Dakota::RealVector(n_eq_con);

            if ((n_nl_con + n_eq_con) > 0)  n_con = n_nl_con + n_eq_con;
            n_resp = n_obj + n_con;

            name_responses.clear();
            name_responses.resize(n_resp);
            //read data from file
            for (unsigned int i=0; i<n_resp; ++i)
            {
                // obtain the name of the design variable:
                std::ostringstream streamOut;
                streamOut << i;
                // obtain the name of the design variables:
                std::string name = "RESP_" + streamOut.str() + "_name";
                name_responses[i] = param.get(name.c_str());
            }

            for (int i = n_obj; i < (n_obj + n_nl_con + n_eq_con); i++)
            {
                // obtain the name of the constaint bounds
                std::ostringstream streamOut;
                streamOut << i;
                std::string lb = "RESP_" + streamOut.str() + "_lb";
                std::string ub = "RESP_" + streamOut.str() + "_ub";
                std::string eq = "RESP_" + streamOut.str() + "_eq";
                int j = i - n_obj;
                int k = j - n_nl_con;
                // set vector values from parameter file
                if (n_nl_con > j)
                {
                    lb_nl_constraint[j] = param.get_double(lb.c_str());
                    ub_nl_constraint[j] = param.get_double(ub.c_str());
                }
                if (k >= 0)
                {
                    eq_constraint[k] = param.get_double(eq.c_str());
                }
            }
        }//responses
        param.leave_subsection();

        param.enter_subsection("Optimization Program Options");
        {
            // openFCST specific parameters:
            use_dakota_input_file = param.get_bool("Use dakota input file");
            dakota_input = param.get("Dakota_Input_File");
            optimization_strategy = param.get("Optimization strategy");
            optimization_method = param.get("Optimization method");
            a_gradients = param.get_bool("Analytical gradients");
            n_gradients = param.get_bool("Numerical gradients");
            n_gradient_type = param.get("Numerical gradient type");
            hessians = param.get_bool("Analytical hessians");
            NLS_flag = param.get_bool("Use NLS");
            if (NLS_flag)
                optimization_method = "nl2sol"; //default: best algorithm for NLS method - override if specifiying optimization_method

            NLS_data_file = param.get("NLS data file");

            // Dakota Method Independent Controls:
            scaling = param.get_bool("Scaling");
            max_iter = param.get_integer("Maximum iterations");
            max_f_eval = param.get_integer("Maximum function evaluations");
            convergence_tol = param.get_double("Convergence tolerance");
            constraint_tol = param.get_double("Constraint tolerance");


            rand_pareto_weights = 0;
            rand_start_points = 0;


            //                if (optimization_strategy == "pareto_set")
            //                {
            rand_pareto_weights = param.get_integer("Number of random pareto weights");
            //                }
            //                else if (optimization_strategy == "multi_start")
            //                {
            rand_start_points = param.get_integer("Number of random start points");
            //                }
            //                else if (optimization_strategy == "hybrid")
            //                {
            param.enter_subsection("Hybrid optimization");
            {
                hybrid_strategy = param.get("Hybrid optimization strategy");
                n_hybrid_methods = param.get_integer("Number of optimization methods");

                if (optimization_strategy != "hybrid")
                    n_hybrid_methods = 1;

                model_ptr.resize(n_hybrid_methods);
                method_list.resize(n_hybrid_methods);
                hybrid_opt_method.resize(n_hybrid_methods);
                tol.resize(n_hybrid_methods);
                eval_max.resize(n_hybrid_methods);
                iter_max.resize(n_hybrid_methods);

                for (int i = 0; i < n_hybrid_methods; i++)
                {
                    std::ostringstream streamOut;
                    streamOut << i;
                    model_ptr[i] = "model_" + streamOut.str();
                    method_list[i] = "method_" + streamOut.str();

                    std::string name = "Optimization method " + streamOut.str();
                    hybrid_opt_method[i] = param.get(name.c_str());
                    std::string tol_name = "Convergence tolerance " + streamOut.str();
                    tol[i] = param.get_double(tol_name.c_str());
                    std::string eval_name = "Maximum function evaluations " + streamOut.str();
                    eval_max[i] = param.get_double(eval_name.c_str());
                    std::string iter_name = "Maximum iterations " + streamOut.str();
                    iter_max[i] = param.get_double(iter_name.c_str());
                }
            } //hybrid opt
            param.leave_subsection();


            param.enter_subsection("CONMIN");
            {

            }
            param.leave_subsection();

            //                if (optimization_method == "optpp_q_newton")
            //                {
            param.enter_subsection("OPT++");
            {
                step_to_boundary = param.get_double("Steplength to boundary");
                grad_tol = param.get_double("Gradient tolerance");
                merit_function = param.get("Merit function");
                centering_param = param.get_double("Centering parameter");
                //					search_method = param.get("Search method");
            }//opt++
            param.leave_subsection();
            //                }
            //                if (optimization_method == "nl2sol")
            //                {
            param.enter_subsection("NL2SOL");
            {
                trust_rad = param.get_double("Initial trust radius");
                abs_conv_tol = param.get_double("Absolute convergence tolerance");
                x_conv_tol = param.get_double("Parameter convergence tolerance");
                false_conv_tol = param.get_double("False convergence tolerance");
            }//nl2sol
            param.leave_subsection();
            //                }
            //                if (optimization_method == "ncsu_direct")
            //                {
            param.enter_subsection("NCSU");
            {
                solution_target = param.get_double("solution target");
                final_solutions = param.get_double("final solutions");
                min_boxsize_limit = param.get_double("min boxsize limit");
                volume_boxsize_limit = param.get_double("volume boxsize limit");
            }//ncsu_direct
            param.leave_subsection();

            //                }
            //                if (optimization_method == "coliny_cobyla")
            //                {
            param.enter_subsection("coliny_cobyla");
            {
                initial_delta = param.get_double("Initial Delta");
                threshold_delta = param.get_double("Threshold Delta");
            }//coliny_cobyla
            param.leave_subsection();

            //                }
            //                if (optimization_method == "coliny_pattern_search")
            //                {
            param.enter_subsection("coliny_pattern_search");
            {
                initial_delta = param.get_double("Initial Delta");
                threshold_delta = param.get_double("Threshold Delta");
                contraction_factor = param.get_double("Contraction Factor");
            }//coliny_pattern_search
            param.leave_subsection();

            //                }
            //                if (optimization_method == "asynch_pattern_search")
            //                {
            param.enter_subsection("asynch_pattern_search");
            {
                initial_delta = param.get_double("Initial Delta");
                threshold_delta = param.get_double("Threshold Delta");
                contraction_factor = param.get_double("Contraction Factor");
            }//asynch_pattern_search
            param.leave_subsection();

            //                }
            //                if (optimization_method == "coliny_solis_wets")
            //                {
            param.enter_subsection("asynch_pattern_search");
            {
                initial_delta = param.get_double("Initial Delta");
                threshold_delta = param.get_double("Threshold Delta");
                contraction_factor = param.get_double("Contraction Factor");
            }//coliny_solis_wets
            param.leave_subsection();


            //                }
            //                if (optimization_method == "moga")
            //                {
            param.enter_subsection("moga");
            {
                crossover_type = param.get("Cross Over Type");
                initialization_type = param.get("Initialization Type");
                mutation_type = param.get("Mutation Type");
                num_offspring = param.get_double("Number of Offspring");
                num_parents = param.get_double("Number of Parents");

                fitness_type_moga = param.get("Fitness Type");
                convergence_type_moga = param.get("Convergence Type");
                replacement_type_moga = param.get("Replacement Type");
            }//moga
            param.leave_subsection();

            //                }
            //                if (optimization_method == "soga")
            //                {
            param.enter_subsection("soga");
            {
                crossover_type = param.get("Cross Over Type");
                initialization_type = param.get("Initialization Type");
                mutation_type = param.get("Mutation Type");
                num_offspring = param.get_double("Number of Offspring");
                num_parents = param.get_double("Number of Parents");

                fitness_type = param.get("Fitness Type");
                convergence_type = param.get("Convergence Type");
                replacement_type = param.get("Replacement Type");
            }//soga
            param.leave_subsection();

        }
        param.leave_subsection(); //Optimization Program Options

    }
    param.leave_subsection(); //Optimization Parameters

}

//---------------------------------------------------------------------------
void DakotaApplication::manage_inputs(ParameterHandler& param)
{

    // Read data from input file
    FcstUtilities::read_parameter_files(param,
                                        input_file);

    // Retrieve parameter values from parameter handler
    initialize(param);

    if (!use_dakota_input_file)
    {
        if (optimization_strategy == "hybrid")
        {
            if (n_hybrid_methods > 1)
                write_hybrid_dakota_input_file();
            else
            {
                FcstUtilities::log << "Must have more than one optimization method for hybrid optimization.\n"
                << "Running single optimization strategy\n";
                optimization_strategy = "single_method";
            }
        }
        else if (optimization_strategy == "pareto_set")
            write_pareto_dakota_input_file();
        else if (optimization_strategy == "multi_start")
            write_multi_start_dakota_input_file();
        else if (optimization_strategy == "single_method")
            write_single_method_dakota_input_file();
        else
            FcstUtilities::log << "Optimization strategy not valid.\n";
    }

    // Parse augmented input file into data base
    problem_db->parse_inputs(dakota_input.c_str());

    for(int i = 0; i < n_hybrid_methods; i++)
    {
        problem_db->set_db_list_nodes(method_list[i]);

        // Set design variable information to the data base directly
        problem_db->set("variables.continuous_design.labels", dakota_name_design_var);
        problem_db->set("variables.continuous_design.lower_bounds", lb_design_var);
        problem_db->set("variables.continuous_design.upper_bounds", ub_design_var);
        problem_db->set("variables.continuous_design.scale_types", scale_types);
        problem_db->set("variables.continuous_design.scales", scale_design_var);
        if (optimization_strategy != "multi_start")
            problem_db->set("variables.continuous_design.initial_point", ip_design_var);

        if (n_nl_con > 0)
        {
            problem_db->set("responses.nonlinear_inequality_lower_bounds", lb_nl_constraint);
            problem_db->set("responses.nonlinear_inequality_upper_bounds", ub_nl_constraint);
        }
        if (n_eq_con > 0)
        {
            problem_db->set("responses.nonlinear_equality_targets", eq_constraint);
        }
    }

    // broadcast minimal DB specification to other processors (if needed)
    problem_db->broadcast();
    // Perform post-processing of minimal specification on all processors
    problem_db->post_process();

    // instantiate the strategy
    selected_strategy = Strategy(*problem_db.get());
}

//---------------------------------------------------------------------------

void DakotaApplication::write_hybrid_dakota_input_file()
{
    // Note: Best use for Hybrid method is without direct interface -- allows asynchronous evaluations.
    // TODO: needs development
    std::stringstream inputss;
    // interface
    inputss << "interface, \n\t id_interface = 'I1' \n\t direct,\n\t  analysis_driver = 'FuelCell Direct Interface'\n\n";
    // strategy
    inputss << "strategy, \n\t" << optimization_strategy << " " << hybrid_strategy << "\n\t";
    //TODO: Implement other methods (embedded, collaborative)
    if (hybrid_strategy == "sequential")
    {
        inputss << "method_list = ";
        for (int i = 0; i < n_hybrid_methods; i++)
            inputss << "'"<< method_list[i] <<"' ";
    }
    inputss << "\n\t graphics \n\t tabular_graphics_data \n\n";

    for (int i = 0; i < n_hybrid_methods; i++)
    {
        inputss << "model, \n\t"
        << "id_model = 'model_"<<i<<"'\n\t"
        << "single \n\t   "
        << "variables_pointer = 'V1' \n\t   "
        << "interface_pointer = 'I1' \n\t   "
        << "responses_pointer = 'R" << i << "'\n\n   ";
        write_responses(inputss, i, hybrid_opt_method[i]);
    }

    for (int i = 0; i < n_hybrid_methods; i++)
    {
        convergence_tol = tol[i];
        max_iter = iter_max[i];
        max_f_eval = eval_max[i];
        write_method(inputss, method_list[i], hybrid_opt_method[i], model_ptr[i]);
    }
    write_variables(inputss, 1);
    // Save augmented input file
    std::ofstream fout(dakota_input.c_str());
    fout << inputss.str().c_str();
    fout.close();
}

//---------------------------------------------------------------------------

void DakotaApplication::write_pareto_dakota_input_file()
{
    std::stringstream inputss;
    // interface
    inputss << "interface, \n\t id_interface = 'I1' \n\t direct,\n\t  analysis_driver = 'FuelCell Direct Interface'\n\n";
    // strategy
    inputss << "strategy, \n\t" << optimization_strategy << "\n\t"
    <<  "method_pointer = '" << method_list[0] << "'\n\t graphics \n\t tabular_graphics_data \n";
    if (rand_pareto_weights > 1)
        inputss << "\t random_weight_sets = " << rand_pareto_weights << "\n\n";
    else
    {
        FcstUtilities::log << "Warning: Using default pareto weights for 2 objective functions\n"
        << "See dakota_application.cc for weights, "
        << "or check options for using dakota input file\n";

        if(n_obj != 2)
        {
            FcstUtilities::log << "Number of objective functions not equal to 2!\n";
            abort();
        }

        inputss << "\t multi_objective_weight_sets = \n"
        << "\t\t 1.0 \t 0.0 \n"
        << "\t\t 0.9 \t 0.1 \n"
        << "\t\t 0.8 \t 0.2 \n"
        << "\t\t 0.7 \t 0.3 \n"
        << "\t\t 0.6 \t 0.4 \n"
        << "\t\t 0.5 \t 0.5 \n"
        << "\t\t 0.4 \t 0.6 \n"
        << "\t\t 0.3 \t 0.7 \n"
        << "\t\t 0.2 \t 0.8 \n"
        << "\t\t 0.1 \t 0.9 \n"
        << "\t\t 0.0 \t 1.0 \n"
        << std::endl;
    }

    inputss << "model, \n\t"
    << "id_model = '" << model_ptr[0] << "'\n\t"
    << "single \n\t   "
    << "variables_pointer = 'V1' \n\t"
    << "interface_pointer = 'I1' \n\t"
    << "responses_pointer = 'R1' \n\n";

    write_method(inputss, method_list[0], optimization_method, model_ptr[0]);
    write_responses(inputss, 1, optimization_method);
    write_variables(inputss,1);

    // Save input file
    std::ofstream fout(dakota_input.c_str());
    fout << inputss.str().c_str();
    fout.close();
}

//---------------------------------------------------------------------------

void DakotaApplication::write_multi_start_dakota_input_file()
{
    std::stringstream inputss;
    // interface
    inputss << "interface, \n\t id_interface = 'I1' \n\t direct,\n\t  analysis_driver = 'FuelCell Direct Interface'\n\n";
    // strategy
    inputss << "strategy, \n\t" << optimization_strategy << "\n\t"
    <<  "method_pointer = '" << method_list[0] << "'\n\t graphics \n\t tabular_graphics_data \n";
    if (rand_start_points > 1)
        inputss << "\t random_starts = " << rand_start_points << std::endl;
    else
    {
        FcstUtilities::log << "Warning: Using default starting points.\n"
        << "See dakota_application.cc\n"
        << "or check options for using dakota input file\n";
        inputss << "\t starting_points =  ";
        int num_starting_points = 4;
        for (int j = 0; j < num_starting_points; j++)
        {
            double k = double(j) / double(num_starting_points - 1);
            // 			FcstUtilities::log << "Random Start iteration: " << j << std::endl;
            // 			FcstUtilities::log << k << std::endl;
            for (unsigned int i = 0; i < n_dvar; i++)
            {
                inputss << lb_design_var[i] + k * ((ub_design_var[i] - lb_design_var[i])) << " ";
            }
            inputss << "\n\t\t\t ";
        }
    }
    inputss << "\n";

    inputss << "model, \n\t"
    << "id_model = '" << model_ptr[0] << "'\n\t"
    << "single \n\t   "
    << "variables_pointer = 'V1' \n\t"
    << "interface_pointer = 'I1' \n\t"
    << "responses_pointer = 'R1' \n\n";

    write_method(inputss, method_list[0], optimization_method, model_ptr[0]);
    write_responses(inputss,1, optimization_method);
    write_variables(inputss,1);

    // Save input file
    std::ofstream fout(dakota_input.c_str());
    fout << inputss.str().c_str();
    fout.close();
}

//---------------------------------------------------------------------------

void DakotaApplication::write_single_method_dakota_input_file()
{
    std::stringstream inputss;
    // interface
    inputss << "interface, \n\t id_interface = 'I1' \n\t direct,\n\t  analysis_driver = 'FuelCell Direct Interface'\n\n";
    // strategy
    inputss << "strategy, \n\t" << optimization_strategy << "\n\t"
    <<  "method_pointer = '" << method_list[0] << "'\n\t graphics \n\t tabular_graphics_data \n";

    inputss << "model, \n\t"
    << "id_model = '" << model_ptr[0] << "'\n\t"
    << "single \n\t   "
    << "variables_pointer = 'V1' \n\t"
    << "interface_pointer = 'I1' \n\t"
    << "responses_pointer = 'R1' \n\n";

    write_method(inputss, method_list[0], optimization_method, model_ptr[0]);
    write_responses(inputss,1, optimization_method);
    write_variables(inputss,1);

    // Save input file
    std::ofstream fout(dakota_input.c_str());
    fout << inputss.str().c_str();
    fout.close();
}

//---------------------------------------------------------------------------
void DakotaApplication::write_method(std::stringstream& inputss, std::string& method_id, std::string method_name, std::string model_ptr)
{
    inputss << "method, \n"
    << "\t id_method = '" << method_id << "'\n\t " << "model_pointer = '" << model_ptr << "'\n\t " << method_name << "\n";

    //method independent variables
    if (scaling)
        inputss << "\t scaling" << std::endl;

    /*if (optimization_method == "vector_parameter_study")
     *       {
     *             inputss << "\t  final_point = 0.8"<< std::endl;
     *             inputss << "\t  num_steps = "<< std::endl;
     *          for (int i=0;i<n_dvar;i++)
     *            inputss << part_design_var[i] << " ";
     }*/

    if (method_name == "multidim_parameter_study")
    {
        inputss << "\t  partitions = ";
        for (unsigned int i = 0; i < n_dvar; i++)
            inputss << part_design_var[i] << " ";
        inputss << "\n\n";
    }
    else
    {
        inputss << "\t convergence_tolerance = " << convergence_tol <<"\n";
        inputss << "\t constraint_tolerance = " << constraint_tol <<"\n";
        inputss << "\t max_iterations = " << max_iter << "\n";
        inputss << "\t max_function_evaluations = " << max_f_eval << "\n\n";

        if (method_name == "optpp_q_newton")
        {
            inputss << "\t steplength_to_boundary = " << step_to_boundary << "\n";
            inputss << "\t gradient_tolerance = " << grad_tol << "\n";
            inputss << "\t centering_parameter = " << centering_param << "\n";
            //inputss << "\t merit_function = '" << merit_function << "'\n\n";
            //	    inputss << "search_method,  \n \t" << search_method << "\n\n";
        }

        else if (method_name == "nl2sol")
        {
            ExperimentalData NLS_data(NLS_data_file);
            //NLS_data.print_data();
            m_data_points = NLS_data.get_num_rows();
            //FcstUtilities::log << "Number of data points: " << m_data_points << std::endl;

            inputss << "\t initial_trust_radius = " << trust_rad << "\n";
            inputss << "\t x_conv_tol = " << x_conv_tol << std::endl; //tolerance in parameter changes
            inputss << "\t absolute_conv_tol = " << abs_conv_tol << "\n"; //absolute tolerance in sum of squares term
            inputss << "\t false_conv_tol = " << false_conv_tol << "\n\n";
        }
        else if (method_name == "ncsu_direct")
        {
            inputss << "\t final_solutions = " << final_solutions << "\n";
            inputss << "\t solution_target = " << min_boxsize_limit << "\n";
            inputss << "\t min_boxsize_limit = " << min_boxsize_limit << "\n";
            inputss << "\t volume_boxsize_limit = " << volume_boxsize_limit << "\n\n";
        }

        else if (method_name == "coliny_cobyla")
        {
            inputss << "\t initial_delta = " << initial_delta << "\n";
            inputss << "\t threshold_delta = " << threshold_delta << "\n\n";
        }

        else if (method_name == "coliny_pattern_search")
        {
            inputss << "\t initial_delta = " << initial_delta << "\n";
            inputss << "\t contraction_factor = " << contraction_factor << "\n";
            inputss << "\t threshold_delta = " << threshold_delta << "\n\n";

        }

        else if (method_name == "asynch_pattern_search")
        {
            inputss << "\t initial_delta = " << initial_delta << "\n";
            inputss << "\t contraction_factor = " << contraction_factor << "\n";
            inputss << "\t threshold_delta = " << threshold_delta << "\n\n";

        }

        else if (method_name == "coliny_solis_wets")
        {
            inputss << "\t initial_delta = " << initial_delta << "\n";
            inputss << "\t contraction_factor = " << contraction_factor << "\n";
            inputss << "\t threshold_delta = " << threshold_delta << "\n\n";

        }

        else if (method_name == "moga")
        {
            inputss << "\t crossover_type " << crossover_type << "\n";
            inputss << "\t initialization_type " << initialization_type << "\n";
            inputss << "\t mutation_type " << mutation_type << "\n";
            inputss << "\t num_offspring " << num_offspring << "\n";
            inputss << "\t num_parents " << num_parents << "\n";


            inputss << "\t replacement_type " << replacement_type_moga << "\n";
            inputss << "\t convergence_type " << convergence_type_moga << "\n";
            inputss << "\t fitness_type " << fitness_type_moga << "\n\n";

        }


        else if (method_name == "soga")
        {
            inputss << "\t crossover_type " << crossover_type << "\n";
            inputss << "\t initialization_type " << initialization_type << "\n";
            inputss << "\t mutation_type " << mutation_type << "\n";
            inputss << "\t num_offspring " << num_offspring << "\n";
            inputss << "\t num_parents " << num_parents << "\n";


            inputss << "\t replacement_type " << replacement_type << "\n";
            inputss << "\t convergence_type " << convergence_type << "\n";
            inputss << "\t fitness_type " << fitness_type << "\n\n";

        }

    }
    // TODO: more methods should be included here

}

//---------------------------------------------------------------------------
void DakotaApplication::write_variables(std::stringstream& inputss, int id)
{
    inputss << "variables, \n\t id_variables = 'V" << id << "'\n\t continuous_design = " << n_dvar << "\n\n";
}

//---------------------------------------------------------------------------
void DakotaApplication::write_responses(std::stringstream& inputss, int id, std::string method_name)
{

    inputss << "responses, \n\t id_responses = 'R" << id << "'\n";
    if (method_name == "nl2sol")
        inputss << "\t num_least_squares_terms = " << m_data_points << "\n";
    else
        inputss << "\t num_objective_functions = " << n_obj << "\n";
    if (n_nl_con > 0)
        inputss << "\t num_nonlinear_inequality_constraints = " << n_nl_con << "\n";
    if (n_eq_con > 0)
        inputss << "\t num_nonlinear_equality_constraints = " << n_eq_con << "\n";
    if (a_gradients)
        inputss << "\t analytic_gradients " << "\n";
    else if (n_gradients && method_name != "ncsu_direct")
    {
        inputss << "\t numerical_gradients " << "\n";
        inputss << "\t \t method_source dakota " << "\n";
        inputss << "\t \t interval_type " << n_gradient_type << "\n";
        inputss << "\t \t fd_gradient_step_size = ";
        for (unsigned int i = 0; i < n_dvar; ++i)
            inputss << step_design_var[i] << " ";
        inputss << "\n";
    }
    else
        inputss << "\t no_gradients " << "\n";
    if (hessians)
        inputss << "\t analytic_hessians " << "\n";
    else
        inputss << "\t no_hessians " << "\n\n";
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void DakotaApplication::synchronize_variables(unsigned int &dvar, unsigned int &resp,
                                              std::vector<std::string> &name_dvar, std::vector<std::string> &name_resp)
{
    dvar = this->n_dvar;			// number of design variables
    resp = this->n_resp;			// number of responses
    name_dvar = this->name_design_var; // name of variables
    name_resp = this->name_responses;  // name of responses
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void DakotaApplication::assign_interface(Dakota::DirectApplicInterface* optimization_interface)
{
    // retrieve the list of Models from the Strategy
    ModelList& models = problem_db->model_list();
    // iterate over the Model list
    for (ModelLIter ml_iter = models.begin(); ml_iter != models.end(); ml_iter++)
    {
        Interface& interface = ml_iter->derived_interface();
        // set the correct list nodes within the DB prior to new instantiations
        problem_db->set_db_model_nodes(ml_iter->model_id());
        // plug in the new direct interface instance (DB does not need to be set)
        interface.assign_rep(optimization_interface, false);
    }
}


#endif  //Dakota