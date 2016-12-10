//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: parametric_study.h
//    - Description: Application used to perform parameteric studies on a given application.
//    - Developers: M. Secanell
//
//---------------------------------------------------------------------------

#ifndef _FUELCELL__PARAMETRIC_STUDY_H_
#define _FUELCELL__PARAMETRIC_STUDY_H_

// OpenFCST
#include <application_core/optimization_block_matrix_application.h>
#include <utils/simulation_selector.h>
#include <utils/fcst_utilities.h>

// STD include files
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>

namespace FuelCell
{
    /**
     * This application is used to perform a parametric study for a given fuel cell model.
     * This class reads the linear and a nonlinear applications to be used from the parameter file
     * and uses them in order to obtain a parameteric study for the given application.
     *
     * <h3> Usage </h3>
     *
     * This class is used in SimulationBuilder to run a complete parameteric study. All parameters
     * for this application should be defined in the main parameter file. There are two methods
     * for running a parametric study. The first involves specifying multiple parameters (up to 10) 
     * and a respective list that these parameters values will be set to. An example of this is: 
     * @code
     * subsection Simulator
     *     subsection Parametric study
     *         set Parameter 1 name   = subsection>>subsection>>name1
     *         set Parameter 1 values = 0.1, 0.2, 0.3  #Solves for value 0.1 first then sequentially moves through list
     *         
     *         set Parameter 2 name   = subsection>>subsection>>name2
     *         set Parameter 2 values = 1.0, 0.9, 0.8 
     *         ...
     *     end
     * end
     * @endcode
     * \note: On the first iteration it will solve the problem with name1 set to 0.1 and name2 set to 1.0. On the next iteration it will solve with name1 set to 0.2 and name2 set to 0.9, ...
     * \note: If you wish to change a value at a boundary or material ID then when specifying name1 add ":#" where # is the boundary/material ID number, e.g. Boundary data>>density_species_1 [g/cm^3]:2
     * \note: If you want to alter more than 10 parameters at a time, change max_num_parameters in parametric_study.h file and recompile OpenFCST. 
     *
     * The second method only allows for one parameter to be changed but it has the ease that you do not need to specify 
     * all the values that you wish to run the simulation at. Instead you enter an upper and lower bounds and the step size.
     * An example of this is:
     * @code
     * subsection Simulator
     *     subsection Parametric study
     *         set Parameter 1 name   = subsection>>subsection>>name
     *         set Initial value      = 0.5   #Enter the value you would like to start the parameteric study from.
     *         set Final value        = 1.0   #Enter the final value for the parameteric study.
     *         set Increment          = 0.1   #Spacing between points in the parameteric study
     *         set Adaptive Increment = false #Set to true if you would like to reduce the voltage increment adaptively if convergence could not be achieved with the larger value
     *         set Min. Increment     = 0.1   #If Adaptive Increment is set to true, this value controls the
     *                                        #minimum change in cell voltage before the parameter study gives up
     *                                        #and the voltage is updated again. Note that this value has to be more
     *                                        #than 0.01 V as a value of zero would lead to an infinite loop.
     *     end
     * end
     * @endcode
     *
     * To use the class, imply create an object, declare the parameters, read the file, initialize and run
     * @code
     * FuelCell::ParametricStudy<dim> curve;
     * ParameterHandler param;
     * curve.declare_parameters(param);
     * std::string simulator_parameter_file_name = "input_file.prm";
     * param.read_input(simulator_parameter_file_name,
     *                  true);
     * curve.initialize(param);
     * curve.run(param, simulator_parameter_file_name, sim_selector);
     * @endcode
     *
     * @author M. Secanell, 2014
     */
    template <int dim>
    class ParametricStudy
    {
    public:
        /**
         * Constructor
         */
        ParametricStudy();

        /**
         *
         */
        ~ParametricStudy();

        /**
         *
         */
        void declare_parameters (ParameterHandler& ) const;

        /**
         * Read parameters from file.
         */
        void initialize (ParameterHandler& );

        /**
         * Run a full polarization curve run
         */
        void run (ParameterHandler& param,
                  std::string simulator_parameter_file_name,
                  boost::shared_ptr<SimulationSelector<dim> > sim_selector);

    protected:
        /**
         * Modify cell voltage and any other parameters that you would like to modify
         * with respect to the inital input file.
         *
         */
        virtual void set_parameters(ParameterHandler& param,
                                    const shared_ptr<FuelCell::ApplicationCore::AdaptiveRefinement<dim> >& solver,
                                    const int iteration,
                                    const std::vector<std::string> parameter_name,
                                    const std::vector<double> param_value);

        /**
         * Run a single point in the polarization curve and return the current density and any other data
         */
        void run_point (ParameterHandler& param,
                        const std::string simulator_parameter_file_name,
                        const boost::shared_ptr<SimulationSelector<dim> > sim_selector,
                        const int iteration,
                        const std::vector<std::string> parameter_name,
                        const std::vector<double> param_value,
                        std::map<std::string, double>& functionals);
        
        /**
         * 
         */
        inline void print_iteration_info(const unsigned int iteration,
                                         const double param_value) const 
        {
            FcstUtilities::log<<"============================================================================="<<std::endl;
            FcstUtilities::log<<"Iteration : "<<iteration<<std::endl;
            FcstUtilities::log<<"Solving for..."<<std::endl;
            FcstUtilities::log<<"Parameter name: ";
            if(!parameter_name.empty())
                FcstUtilities::log<<parameter_name[0];
            FcstUtilities::log<<std::endl;
            FcstUtilities::log<<"Cell param_value: "<<param_value<<std::endl;
            FcstUtilities::log<<"======"<<std::endl;
        }
        
        inline void print_iteration_info(const unsigned int iteration,
                                         const std::vector<double> param_value) const 
        {
            FcstUtilities::log<<"============================================================================="<<std::endl;
            FcstUtilities::log<<"Iteration : "<<iteration<<std::endl;
            FcstUtilities::log<<"Solving for..."<<std::endl;
            for(unsigned int i = 0; i < parameter_name.size(); ++i)
            {
                FcstUtilities::log<<"Parameter " << i+1 << " name: "<<parameter_name[i]<<std::endl;
                FcstUtilities::log<<"Parameter " << i+1 << " value: "<<param_value[i]<<std::endl;
            }
            FcstUtilities::log<<"======"<<std::endl;
        }
        
        /**
         * Print parameters:
         */
        virtual void print_parameters() const;

        /**
         * Print parameteric study header into a file.
         */
        virtual void print_parameteric_study_header() ;

        /**
         * Store the value of the current density in the solution
         */
        void register_data(const double ,
                           std::map<std::string, double>& functionals);

        /**
         * Variable where the output file to store parameteric study results is stored
         */
        std::string parameter_filename;
        
        /**
         * String that defines the parameter that we would like to modify
         * The parameter should be of the from
         * - For normal parameter: Subsection_1>>Subsection_2>>Value
         * - For boundary value or graded: Subsection_1>>Subsection_2>>Material_id:Value
         */
        std::vector<std::string> parameter_name;

        /**
         * Initial value to run the parametric study, in the units specified for the value.
         *
         * @note This vaule is set via input file using the following section
         * @code
         * subsection Simulator
         *   subsection Parametric study
         *     set Initial Value =
         *   end
         * end
         * @endcode
         */
        double p_init;

        /**
         * Final parameter value used in the parameteric study
         *
         * @note This vaule is set via input file using the following section
         * @code
         * subsection Simulator
         *   subsection Parametric study
         *     set Final Value
         *   end
         * end
         * @endcode
         */
        double p_end;

        /**
         * Spacing between points in the parameteric study
         *
         * @note This vaule is set via input file using the following section
         * @code
         * subsection Simulator
         *   subsection Parameteric study
         *     set Increment
         *   end
         * end
         * @endcode
         */
        double dp;

        /**
         * Another way to define the sequence of parameter values
         * is to use the vector @p p_values instead of variables @p p_init, @p p_end, and @p dp.
         * This vector contains the discrete values
         * of a parameter of study.
         * @note This vector is set via input file using the following section
         * @code
         * subsection Simulator
         *   subsection Parameteric study
         *     set Parameter 1 values = v1, v2, v3, ..., vN
         *   end
         * end
         * @endcode
         */
        std::vector< std::vector<double> > p_values;

        /**
         * Set to true if you would like to reduce the voltage increment adaptively
         * if convergence could not be achieved with the larger value
         *
         * @note This vaule is set via input file using the following section
         * @code
         * subsection Simulator
         *   subsection Parametric study
         *     set Adaptive Increment?  = true (or false)
         *   end
         * end
         * @endcode
         */
        bool adaptive;

        /**
         *
         * @note This vaule is set via input file using the following section
         * @code
         * subsection Simulator
         *   subsection Parametric study
         *     set Min. Increment
         *   end
         * end
         * @endcode
         */
        double min_dp;

        int n_dpPts;

        /**
         * Convergence of the solution. True if the adaptive loop
         * was able to provide a solution.
         */
        bool convergence;

        /**
         * Map used to store the parametric study data
         */
        std::vector<std::vector<double> > curve;
        
        /**
         * File object for output logging
         */
        std::ofstream myfile;

        std::vector<std::string> name_responses;

        int n_resp;
        
        /**
         * Set the maximum number of parameters that you can alter in a Parametric Study
         */
        const unsigned int max_num_parameters = 10;

        /**
         * Vector used to store the solution from the previous iteration
         * if convergence has been achieved.
         */
        FuelCell::ApplicationCore::FEVector coarse_solution;
    };
}

#endif