// ------------------------------------------------------------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2016 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: scaling.h
// - Description: This class allows scaling of any of the vector valued equations by a constant to help with 
//                stability of finding a solution.
//
// - Developers: Chad Balen, University of Alberta
//
// ------------------------------------------------------------------------------------------------------------------------------

#ifndef _FUELCELL_SCALING__H
#define _FUELCELL_SCALING__H

//OpenFCST required header files:
#include <application_core/system_management.h>
#include <utils/fcst_utilities.h> // For converting string to map in initialize()

namespace FuelCell
{
    
    /**
     * Class used to store, read from file and define scaling factors to be applied to equations in equation matrix and solutions in solution residual.
     * It stores the following:
     * - bool to determine if scaling will be applied; \p applyScaling
     * - equations to scale and their respective scale factor in a map; \p ScalingMap
     * - pointer to the systemManagement object
     * 
     * \b NOTE: \p applyScaling does not necessarliy prevent use of function calls from Scaling class. That must be implemented in the application class. 
     *          This bool just provides the necessary information for implementation.
     * 
     * \b WARNING: you should have a thourough understanding of how the equations are applied before using this class. 
     *             As well, always confirm that you have applied them correctly and solutions match a case when no scaling occured.
     * 
     * In the input file, the following parameters can be specified (see declare_parameters ), with examples of how to use:
     * @code
     * subsection Equations
     * (...)
     *     set Apply scaling = true # true | false
     *     set Equation matrix scaling = Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - mass conservation - species 1:1.0e-4
     * end
     * @endcode 
     * 
     * 
     * <h3>Usage Details:</h3>
     * 
     * In order to use this class, first an object of the class needs to be created. Usually, one such objects exists in every application. To create the object,
     * include the .h file in the include application file and in the application data member section add the object. For example:
     * @code
     * #include "scaling.h"
     * 
     * // Then in the data member declaration (usually a private member)
     * FuelCell::Scaling Scale;
     * @endcode
     * 
     * Once the object is created, the section where the input data will be specified in the input file needs to be declared. To do so, in the declare_parameters section 
     * of your application call the following:
     * 
     * @code
     * //--------- IN DECLARE_PARAMETERS ------------------------------------------------------
     * template <int dim>
     * void 
     * NAME::AppCathodeKG<dim>::declare_parameters(ParameterHandler& param)
     * {
     *   (...)
     *   Scale.declare_parameters(param);
     *   (...)
     * }
     * @endcode
     *          
     * 
     * Finally, once the input file has been read by our application, your class needs to be initialized. You may want to use get_apply_scaling_bool() to set a bool
     * that will tell your applications if scaling will occur. This bool can then be checked before scaling functions are called. That way, if no scaling is required
     * then the functions do not need to be called. This is achieved using the function initialize()
     * @code
     * //--------- IN INITIALIZE ------------------------------------------------------
     * template <int dim>
     * void
     * NAME::AppCathodeKG<dim>::initialize(ParameterHandler& param)
     * {   
     *  (...) 
     *  Scale.initialize(param);
     *  bool applyScaling = Scale.get_apply_scaling_bool();
     * }
     * @endcode
     * 
     * You are now ready to use your Scaling object.
     * 
     * 
     * @author C.Balen, 2009-2016
     * 
     */
    class Scaling
    {
        public:
            /**
             * Constructor
             */
            Scaling(FuelCell::SystemManagement& system_management);
            
            /**
             * Destructor
             */
            ~Scaling();
            
            /**
             * Declare all necessary parameters in order to compute the coefficients
             * 
             * The parameters that can be specified in the input file are as follows:
             * 
             * @code
             * subsection Fuel cell data
             * (...)
             *   subsection Equations
             *     set Apply scaling = false
             *     set Equation matrix scaling = Electron Transport Equation:1e-4
             *   end
             * end
             * @endcode 
             */
            void declare_parameters (ParameterHandler &param) const;

            /**
             * Class used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            void initialize (ParameterHandler& param);
            
            /**
             * Function to return \p applyScaling bool. This way you can check bool before calling scaling functions to prevent wasted time running scaling functions.
             */
            const bool get_apply_scaling_bool() const
            {
                return applyScaling;
            }
            
            /**
             * Function scales cell_matrix, typically called in cell_matrix(). Call at then end of function.
             */
            void scale_cell_matrix(MatrixVector& cell_matrices);
            
            /**
             * Function scales cell_residual, typically called in cell_residual(). Call at then end of function.
             */
            void scale_cell_residual(FEVector& cell_res);
            
            /**
             * Function scales bdry_matrix, typically called in bdry_matrix(). Call at then end of function.
             * 
             * \b NOTE: at beginning of function create a local boundary residual like so,
             * 
             * FEVector local_bdry_res(this->block_info.local);
             * 
             * and use this FEVector in assemble_bdry_residual() functions. Finally call scale_bdry_residual()
             * and give the local boundary residual as first argument, and the actual boundary residual as second argument.
             * The correct actual boundary residual is then returned by reference.
             */
            void scale_bdry_matrix(MatrixVector& local_bdry_matrices, MatrixVector& bdry_matrices);
            
            /**
             * Function scales bdry_residual, typically called in bdry_residual(). Call at then end of function.
             * 
             * \b NOTE: at beginning of function create a local boundary matrix like so,
             * 
             * MatrixVector local_bdry_matrices; 
             * Scale.create_local_bdry_matrix(local_bdry_matrices, bdry_matrices); // See create_local_bdry_matrix() for more information
             * 
             * and use this MatrixVector in assemble_bdry_residual() functions. Finally call scale_bdry_matrix()
             * and give the local boundary matrix as first argument, and the actual boundary matrix as second argument.
             * The correct actual boundary matrix is then returned by reference.
             */
            void scale_bdry_residual(FEVector& local_bdry_res, FEVector& bdry_res);
            
            /**
             * Due to issues with MatrixVector this function, loops through all the blocks and correctly sizes \p local_bdry_matrices for you.
             * \p local_bdry_matrices is returned by reference.
             */
            void create_local_bdry_matrix(MatrixVector& local_bdry_matrices, MatrixVector& bdry_matrices);
            
        private:
            /**
             * Pointer to the external YourApplication<dim>::system_management object.
             */
            FuelCell::SystemManagement* system_management;
            
            /**
             * Bool to determine if user wants scalin to occur as based on parameter file. This way you can check bool before calling scaling functions to prevent wasted time running scaling functions.
             */
            bool applyScaling = false;
            
            /**
             * Map of the FuelCell::SystemManagement EquationNames (LHS) and their respective scaling factor.
             */
            std::map<std::string, double> ScalingMap;
            
    };
    
}

#endif // _FUELCELL_SCALING__H