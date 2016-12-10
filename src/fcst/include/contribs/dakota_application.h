//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: dakota_application.h
//    - Description: Used to read in an array of experimental data with column headers.
//                   For use with parameter estimation or generating polarization curves.
//    - Developers: Malte Krack, Peter Dobson <pdobson@ualberta.ca>, Marc Secanell <secanell@ualberta.ca>
//    - $Id: dakota_application.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef dakota_application_h
#define dakota_application_h

#ifdef _WITH_DAKOTA

//deal.II:
#include <deal.II/base/parameter_handler.h>

////FCST classes
#include <contribs/experimental_data.h>
#include <utils/simulation_selector.h>
#include <utils/fcst_utilities.h>

//STL libraries:
#include<string>
#include<fstream>
#include<iostream>
#include<vector>

//Dakota:
#include <ParallelLibrary.hpp>
#include <ProblemDescDB.hpp>
#include <DakotaStrategy.hpp>
#include <DakotaModel.hpp>
#include <DirectApplicInterface.hpp>


namespace SIM
{
    /**
     * This class is used to interface with DAKOTA as an algorithm library. It reads the
     * optimization specifications from the parameter file specified in the main parameter file used to launch
     * FCST under subsection Simulator >> set optimization parameter file name.
     *
     * DakotaApplication is called by simulator_builder and it uses a direct application interface to the analysis code and runs the optimization.
     *
     * This member function defines all available design variables that can be used for optimization. Variables specified in the optimization file in:
     * @code
     * subsection Optimization Parameters
     *  subsection Design Variables
     *      set num_design_variables = 1
     *      set DV_0_name = subsection>>subsection>>subsection>>parameter
     *  end
     * end
     * @endcode
     * need to have been defined in this file. Note that the character ">>" is used to separate subsections. Several parameters
     * in the input file contain a list of comma separated parameters of the form material_id:value or boundary_id:value. If one of this
     * parameters needs to be used for optimization, it can be done using the following format to define the variable:
     * @code
     * set DV_0_name = subsection>>subsection>>subsection>>parameter:material_id
     * @endcode
     * The program will automatically replace the value in the list with the appropriate material_id
     *
     *
     * \author Peter Dobson,    2010
     * \author Malte Krack,    2010
     * \author Marc Secanell, 2012-14
     */
    class DakotaApplication
    {
    public:
        /**
         * Constructor for an object of this class. Optimization details are read from the
         * input_file, the optimization_interface is used for the function evaluations.
         */
        DakotaApplication(boost::shared_ptr<Dakota::ProblemDescDB> global_problem_db,
                          std::string parameter_file);

        /**
         * Declares & Initializes all parameters, creates Dakota input file.
         * Interfaces with Dakota to parse the input file, and populate the Dakota Problem Database
         * with variables and constraints.
         */
        void manage_inputs(ParameterHandler& param);

        /**
         * Assign the direct interface to the simulation code which is of type Dakota::DirectApplicInterface.
         * http://dakota.sandia.gov/docs/dakota/5.2/html-dev/DakLibrary.html
         * For an example,
         * see SIM::DakotaDirectInterface or SIM::DakotaLeastSquaresInterface
         */
        void assign_interface(Dakota::DirectApplicInterface* optimization_interface);

        /**
         * This function will first lock the problem_db preventing anychanges to the data.
         * It will then select the strategy for optimizating (E.g. SingleMethideStrategy,HybridStategy...) and
         * Run the optimization using the DakotaDirectInterface using the data which you have assigned by
         * using assign_interface()
         */
        void run();

        /**
         * Get boolean NLS flag. If this is set true in the text file the optimizating will use the
         * DakotaLeastSquaresInterface routine found in dakota_direct_interface.
         * NLS_flag is set using the opt_app_ .prm file.
         * @code
         *  subsection Optimization Program Options
         *      set Use NLS = true
         * @endcode
         * Additional inputs to consider in .prm file when using NLS.
         * @code
         *      set NLS data file = data.txt
         *      set NLS residual option = weighted
         *      set Numerical gradients = true
         *      set Numerical gradient type = central
         *      set Number of random start points = 4
         *      set Convergence tolerance = 1e-4
         *      set Constraint tolerance = 1e-4
         *      set NLS residual value = norm
         * @endcode
         */
        inline bool use_NLS()
        {
            return NLS_flag;
        }

        /**
         * Function retrieves and returns the final optimum design variables.
         */
        Dakota::RealVector vars_results();

        /**
         * Function retrieves and returns the optimum response values.
         */
        Dakota::RealVector resp_results();

        /**
         * Declare all parameters that are needed for optimization. These parameters are read in from the Opt_data file.
         */
        void declare_parameters(ParameterHandler& param);

        /**
         * This function is used to remotely axcess the number and name of design_variables, number and name of responses,
         * as they are private variables.
         */
        void synchronize_variables(unsigned int &n_dvar, unsigned int &n_resp,
                                   std::vector<std::string> &name_design_var, std::vector<std::string> &name_responses);

    private:
        /**
         * Read in parameters for the parameter handler in order to initialize data.These parameters are read in from the Opt_data file.
         */
        void initialize(ParameterHandler& param);
        /**
         * Write Dakota input file for hybrid optimization strategy
         */
        void write_hybrid_dakota_input_file();
        /**
         * Write Dakota input file for pareto optimization strategy
         */
        void write_pareto_dakota_input_file();
        /**
         * Write Dakota input file for multi-start optimization strategy
         */
        void write_multi_start_dakota_input_file();
        /**
         * Write Dakota input file for single optimization strategy
         */
        void write_single_method_dakota_input_file();
        /**
         * Member function to write Dakota input file based on method selection
         */
        void write_method(std::stringstream& inputss, std::string& method_id, std::string method_name, std::string model_ptr);
        /**
         * Member function to write Dakota input file based on specified response variables
         */
        void write_responses(std::stringstream& inputss, int id, std::string method_name);
        /**
         * Member function to write Dakota input file based on specified variables
         */
        void write_variables(std::stringstream& inputss, int id);

        /**
         * Exception thrown when a user defined design variables
         * is not found among available FCST design variables.
         */
        DeclException2(DesignVariableNotFoundInFCSTDesignVariables,
                       std::string,
                       std::string,
                       << "A "  << arg1 << " with name \"" << arg2 << "\" is not stored in available FCST design variables. Please define the variable in DakotaApplication class and in set_parameters in the appropriate class");

        /**
         * Analytical Gradients option for Dakota optimization
         */
        bool a_gradients;
        /**
         * Numerical Gradients option for Dakota optimization
         * @code
         * subsection Optimization Program Options
         *        set Numerical gradients = true
         * @endcode
         */
        bool n_gradients;
        /**
         * Numerical Gradients type for Dakota optimization
         * @code
         *  subsection Optimization Program Options
         *        set Numerical gradient type = central
         * @endcode
         * There are two types of "Numerical gradient types".
         * 		-forward
         * 		-central
         * These are declared in DakotaApplication::declare_parameters()
         * @code
         *        param.declare_entry("Numerical gradient type",
         *                            "forward",
         *                            Patterns::Selection("forward|central"));
         * @endcode
         */
        std::string n_gradient_type;
        /**
         * Analytical Hessian Option for Dakota optimization
         */
        bool hessians;

        /**
         * Number of design variables. This variable is set in the "opt_app_***_***.prm" file. Each design variable requires:
         * - Initial Point
         * - Lower Bound
         * - Upper Bound
         * - Scale (Method & Magnitude)
         * - Step size
         * An Example of the implemented prm file can be seen below.
         * @code
         *    subsection Design Variables
         * 		set num_design_variables = 1	 # 2
         * 		set DV_0_name = prc_N_c		 # L_CCL
         * 		set DV_1_name = T_cell
         *
         *     ####### Initial Point #######
         *     #######---------------#######
         * 		set DV_0_ip = 0.2
         *		set DV_1_ip = 280
         *
         *     #######	  	Lower Bound 	  #######
         *	  ####### lb < -1e30 for -inf   #######
         *     #-----------------------------------#
         *		set DV_0_lb = 0.2
         *		set DV_1_lb = 273.15
         *
         *      #######	 Upper Bound 	#######
         *      ####### ub > 1e30 for inf #######
         *      #-------------------------------#
         *		set DV_0_ub = 0.5						# 10e-4
         *		set DV_1_ub = 373.15
         *
         *      ####### Scales #######
         *      #######--------#######
         *		set DV_0_scale_method = value				# none | auto | value
         *		set DV_1_scale_method = value				# none | auto | value
         *
         *		set DV_0_scale = 0.1
         *		set DV_1_scale = 1
         *
         *      ####### Step size #######
         *      #######-----------#######
         *		set DV_0_step = 1e-4
         *		set DV_1_step = 1e-4
         *
         *  end
         * @endcode
         */
        unsigned int n_dvar;
        /**
         * Number of objective functions
         */
        unsigned int n_obj;
        /**
         * Number of non-linear constraints
         */
        unsigned int n_nl_con;
        /**
         * Number of equality constraints
         */
        unsigned int n_eq_con;
        /**
         * Total number of constraints = n_nl_con + n_eq_con
         */
        unsigned int n_con;
        /**
         * Number of responses = n_obj + n_con
         */
        unsigned int n_resp;
        /**
         * Member that stores the name of the design variables
         */
        std::vector<std::string> name_design_var;
        /**
         * Member that stores the name of all possible design variables that have been defined in FCST
         */
        std::vector<std::string> all_name_design_var;
        /**
         * Member that stores the name of the responses, i.e. objectives and constraints.
         * @note first we put the name of all objectives and then all constraints
         */
        std::vector<std::string> name_responses;

        /**
         * Name of the analysis input file
         */
        std::string input_file;
        /**
         * Pointer to the problem description data base
         */
        boost::shared_ptr<Dakota::ProblemDescDB> problem_db;

        //DAKOTA variables
        /**
         * Dakota version used
         */
        std::string dakota_version;

        /**
         * Flag to use an input file from Dakota
         * Specify interface, method, and strategy options, number of variables, number of responses/constrants.
         * Identifier strings must be specified as method_0, and sequentially (i.e. method_1, method_2, etc.) for hybrid method.
         * Do not specify variable string identifiers, bound values, constraint values, etc.
         * Variable bounds, responses and constraints set in parameter file
         */
        bool use_dakota_input_file;
        /**
         * Name of the Dakota input file
         */
        std::string dakota_input;
        /**
         * Optimization Strategy - Single, Pareto, Multi start point, Hybrid.
         * @code
         *   subsection Optimization Program Options
         *      set Optimization strategy = single_method	 		# single_method | multi_start | pareto_set | hybrid
         * @endcode
         */
        std::string optimization_strategy;
        /**
         * Optimization Method - OPT++, NL2SOL, multidim_parameter_study, etc.
         * @code
         *   subsection Optimization Program Options
         *      set Optimization method = optpp_q_newton			# multidim_parameter_study | optpp_q_newton | nl2sol | ncsu_direct
         * @endcode
         */
        std::string optimization_method;
        /**
         * Nonlinear Least-Squares Parameter Estimation Flag
         */
        bool NLS_flag;
        /**
         * Used to activate scaling. It will use the values in scale and scale_method in order to scale the design variables
         * (used by NLS and OPT+ methods).
         * Implementation of this can be seen above in number of design variables.
         */
        bool scaling;
        /**
         * Experimental data file for NLS parameter estimation - Operating Conditions & output
         */
        std::string NLS_data_file;
        /**
         * Maximum number of optimization iterations
         */
        int max_iter;
        /**
         * Maximum number of optimization function evaluations
         */
        int max_f_eval;
        /**
         * Optimization convergence tolerance
         */
        double convergence_tol;
        /**
         * Optimization constraint tolerance
         */
        double constraint_tol;
        /**
         * Number of random weights for generating objective function pareto set
         */
        int rand_pareto_weights;
        /**
         * Number of random weights for generating objective function pareto set
         */
        int rand_start_points;
        /**
         * Hybrid Optimization strategy - sequential, etc.
         */
        std::string hybrid_strategy;
        /**
         * Number of Hybrid Optimization methods
         */
        int n_hybrid_methods;
        /**
         * Hybrid Optimization method names
         */
        std::vector<std::string> hybrid_opt_method;
        /**
         * Hybrid Optimization method string identifiers
         */
        std::vector<std::string> method_list;
        /**
         * Hybrid Optimization model string identifiers
         */
        std::vector<std::string> model_ptr;
        /**
         * Hybrid Optimization method specific tolerances
         */
        std::vector<double> tol;
        /**
         * Hybrid Optimization method specific function evaluation maximums
         */
        std::vector<double> eval_max;
        /**
         * Hybrid Optimization method specific iteration maximums
         */
        std::vector<double> iter_max;


        /**
         * Parameter Study Partition vector
         * Set number of partitions for multidim_parameter_study
         * NOTE: Evaluated at n+1 points between lower and upper bound
         */
        std::vector<int> part_design_var;
        /**
         * OPT++ steplength
         */
        double step_to_boundary;
        /**
         * OPT++ gradient tolerance
         */
        double grad_tol;
        /**
         * OPT++ merit function
         */
        std::string merit_function;
        /**
         * NL2SOL initial trust radius
         */
        double trust_rad;
        /**
         * NL2SOL absolute convergence tolerance
         */
        double abs_conv_tol;
        /**
         * NL2SOL Parameter convergence tolerance
         */
        double x_conv_tol;
        /**
         * NL2SOL Parameter false convergence tolerance
         */
        double false_conv_tol;
        /**
         * NL2SOL Parameter storing the number of terms
         */
        int m_data_points;

        /**
         * NCSU parameter
         */
        double solution_target;
        /**
         * NCSU parameter
         */
        double final_solutions;
        /**
         * NCSU parameter
         */
        double min_boxsize_limit;
        /**
         * NCSU parameter
         */
        double volume_boxsize_limit;


        /**
         * Dakota array of design variable names
         */
        Dakota::StringArray dakota_name_design_var;
        /**
         * Dakota vector of design variable initial points
         */
        Dakota::RealVector ip_design_var;
        /**
         * Dakota vector of design variable lower bounds
         */
        Dakota::RealVector lb_design_var;
        /**
         * Dakota vector of design variable upper bounds
         */
        Dakota::RealVector ub_design_var;
        /**
         * Dakota vector of design variable scales
         */
        Dakota::RealVector scale_design_var;
        /**
         * Dakota vector of design variable scale methods
         */
        Dakota::StringArray scale_types;
        /**
         * Dakota vector of design variable step sizes for numerical gradients
         */
        Dakota::RealVector step_design_var;
        /**
         * Dakota vector of nonlinear constraint lower bounds
         */
        Dakota::RealVector lb_nl_constraint;
        /**
         * Dakota vector of nonlinear constraint upper bounds
         */
        Dakota::RealVector ub_nl_constraint;
        /**
         * Dakota vector of equality constraints
         */
        Dakota::RealVector eq_constraint;
        /**
         * DAKOTA strategy. From DAKOTA reference guide: "The Strategy class is the base class for the class hierarchy providing
         * the top level control in DAKOTA. The strategy is responsible for creating and managing iterators and models.
         * For memory efficiency and enhanced polymorphism, the strategy hierarchy employs the "letter/envelope idiom"
         * (see Coplien "Advanced C++", p. 133), for which the base class (Strategy) serves as the envelope and one of
         * the derived classes (selected in Strategy::get_strategy()) serves as the letter."
         *
         * For more information see: http://dakota.sandia.gov/licensing/votd/html-dev/
         */
        Dakota::Strategy selected_strategy;

        /**
         * Centering parameter for the opt++ strategy.
         */
        double centering_param;

        /**
         * Optimization SCOLIB (COLINY) - COBYLA (Constrained Optimization By Linear Appoximations) [coliny_cobyla], APPS [asynch_pattern_search], and coliny_pattern_search .
         * These are the variables required when running the COBYLA optimization methods. Below illustrates how these methods can be implemented into the "opt_app_***_***.prm" files.
         * @code
         * ######### Method Specific Parameters #########
         * #########	     SCOLIB (COLINY)      #########
         * #########----------------------------#########
         *
         * 	subsection coliny_cobyla
         * 		set Initial Delta = 1				# (default) 1
         * 		set Threshold Delta = 0.0001			# (default) 0.0001
         * 	end
         *
         * @endcode
         * If using "coliny_pattern_search", "coliny_solis_wets", OR  "asynch_pattern_search" these are the variables available.
         * @code
         * 	subsection coliny_pattern_search			# asynch_pattern_search | coliny_solis_wets
         * 	    set Initial Delta = 1				# (default) 1
         * 	    set Threshold Delta = 0.0001			# (default) 0.0001
         * 	    set Contraction Factor = 0.5			# (default) 0.5
         * 	end
         * @endcode
         */
        double initial_delta;
        double threshold_delta;
        double contraction_factor;
        //
        //    /**
        //      * Optimization - OPT++ Search method .
        //      */
        //     std::string search_method;
        //	bool changing_search_method;




        /**
         * Optimization JEGA (John Eddy Genetic Algorithm) - "moga" & "soga"
         * The JEGA library contains two global optimization methods:
         * MOGA - Multi-objective Genetic Algorithm
         * SOGA - Single-objective Genetic Algorithm
         *
         * An example of implementing these algorithms into the "opt_app_***_***.prm" file can be seen below.
         * When using either there are some variables that are common to both "moga" & "soga"
         *@code
         *
         * 	subsection soga							# moga
         *	    set Cross Over Type = shuffle_random		# (default) shuffle_random
         *	    set Initialization Type = unique_random		# (default) unique_random
         *	    set Mutation Type = replace_uniform			# (default) replace_uniform
         *	    set Number of Offspring = 2				# (default) 2
         *	    set Number of Parents = 2					# (default) 2
         *@endcode
         * Those which are uncommon to "moga" & "soga" are:
         *@code
         *
         * 	    set Replacement Type = elitist					# (default) elitist
         *	    set Convergence Type = average_fitness_tracker	# (default) average_fitness_tracker
         *	    set Fitness Type = merit_function
         * 	end
         *
         *@endcode
         *
         *
         / ***
         * JEGA - mutation_type
         */
        std::string mutation_type;
        /**
         * JEGA - crossover_type
         */
        std::string crossover_type;
        /**
         * JEGA - fitness_type
         */
        std::string fitness_type;
        std::string fitness_type_moga;
        /**
         * JEGA - replacement_type
         */
        std::string replacement_type;
        std::string replacement_type_moga;
        /**
         * JEGA - initialization_type
         */
        std::string initialization_type;
        /**
         * JEGA - convergence_type
         */
        std::string convergence_type;
        std::string convergence_type_moga;
        /**
         * JEGA - Number of Offspring and Parent
         */
        double num_offspring;
        double num_parents;



    };
}

#endif  //Dakota

#endif
