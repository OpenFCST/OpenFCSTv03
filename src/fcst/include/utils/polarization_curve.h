//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: polarization_curve.h
//    - Description: Child of ApplicationWrapper used to implement adaptive refinement
//    - Developers: M. Secanell
//
//---------------------------------------------------------------------------

#ifndef _FUELCELL__POLARIZATION_CURVE_H_
#define _FUELCELL__POLARIZATION_CURVE_H_

// OpenFCST
#include <application_core/optimization_block_matrix_application.h>
#include <utils/simulation_selector.h>
#include <utils/fcst_utilities.h>
#include <utils/parametric_study.h>

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
     * This application is used to compute the polarization curve for a given fuel cell model. 
     * This class reads the linear and a nonlinear applications to be used from the parameter file
     * and uses them in order to obtain a polarization curve for the given application.
     * 
     * Note that this is a specific implementation of the #ParametricStudy class.
     * 
     * @warning
     * It expects that the application has the following capabilities:
     * a) The cell voltage should be set using the folloing flags in the input file:
     * @code
     * subsection Fuel cell data
     *     subsection Operating conditions
     *         set Adjust electronic potential boundary condition = true
     *         set Voltage cell = 0.85 #0.7
     *     end
     * end
     * @endcode
     * and the bipolar plate to GDL boundary ID must be known.
     * @code
     * subsection Grid generation
     *     subsection Internal mesh generator parameters 
     *         subsection Boundary ID
     *             set c_BPP/GDL = 2
     *         end
     *     end
     * end
     * @endcode
     * @note see app_cathode application for an example on voltage modification based on this idea.
     * 
     * b) The application can return an the current density. The current density is returned
     * under name "current"
     *
     * <h3> Usage </h3>
     * 
     * This class is used in SimulationBuilder to run a complete polarization curve. All parameters
     * for this application should be defined in the main parameter file under the following section:
     * @code
     * subsection Simulator
     *     subsection Polarization Curve
     *         set Initial voltage [V] = 1.1
     *         set Final voltage [V]   = 0.5
     *         set Increment [V]       = 0.1
     *         set Adaptive Increment  = false
     *         set Min. Increment [V]  = 0.1
     *     end
     * end
     * @endcode
     * 
     * To use the class, imply create an object, declare the parameters, read the file, initialize and run
     * @code
     * FuelCell::PolarizationCurve<dim> curve;
     * ParameterHandler param;
     * curve.declare_parameters(param);
     * std::string simulator_parameter_file_name = "input_file.prm";
     * param.read_input(simulator_parameter_file_name,
     *                  true);
     * curve.initialize(param);
     * curve.run(param, simulator_parameter_file_name, sim_selector);
     * @endcode
     *
     *
     * @author M. Secanell, 2014
     */
    template <int dim>
    class PolarizationCurve
    :
    public ParametricStudy<dim>
    {
    public:
        /**
         * Constructor
         */
        PolarizationCurve();
        
        /**
         * 
         */
        ~PolarizationCurve();
                /**
         * Declare all parameters that are needed for:
         * 
         * Currently the options implemented are:
         *
         * @code
         * subsection Simulator
         *   subsection Polarization curve
         *     set Initial voltage [V]  # initial voltage value
         *     set Final voltage [V]    # final voltage value
         *     set Increment [V]        # change in voltage between points
         *     set Adaptive Increment?  # allow value to change if convergence not achieved
         *     set Min. Increment [V]   # if voltage allowed to change, min. increment that should be used
         *   end
         * end
         * @endcode
         *
         */
        void declare_parameters (ParameterHandler& param) const;

        /**
         * Read parameters from file.
         */
        void initialize (ParameterHandler& param); 
                                
    private:
        /**
         * Modify cell voltage and any other parameters that you would like to modify
         * with respect to the inital input file.
         * 
         */
        void set_parameters(ParameterHandler& param,
                            const shared_ptr<FuelCell::ApplicationCore::AdaptiveRefinement<dim> >& solver,
                            const int iteration,
                            const std::vector<std::string> parameter_name,
                            const std::vector<double> param_value);
        
        /**
         * Print parameters:
         */
        void print_parameters() const;  
        
        /**
         * Print polarization curve header into a file.
         */
        void print_parameteric_study_header() ;                    
    };
}

#endif