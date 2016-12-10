//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: adaptive_refinement.h
//    - Description: Child of ApplicationWrapper used to implement adaptive refinement
//    - Developers: M. Secanell, Valentin N. Zingan
//
//---------------------------------------------------------------------------

#ifndef _FUELCELL__ADAPTIVEREFINEMENT_H_
#define _FUELCELL__ADAPTIVEREFINEMENT_H_

// Deal.II include files:
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/grid/tria.h>


// Fuel cell include files
#include <application_core/application_wrapper.h>
#include <application_core/optimization_block_matrix_application.h>
#include <utils/fcst_utilities.h>

// STD include files
#include <string>
#include <stdio.h>
#include <stdlib.h>

namespace FuelCell
{
    namespace ApplicationCore
    {
        /**
         * This class is initialized with an application that describes the linearization of the problem that
         * we would like to solve and the nonlinear solver that drives the process (usually a Newton loop).
         * Then, this class implements the adaptive refinement loop for the application.
         *
         * @author M. Secanell, 2006-13
         *
         * and the concept of sequential applications, simplified parameter study without Dakota, and convergence tables by
         *
         * @author Valentin N. Zingan
         */
        template <int dim>
        class AdaptiveRefinement
        :
        public ApplicationWrapper
        {
        public:
            /**
             * Constructor
             */
            AdaptiveRefinement ( FuelCell::ApplicationCore::ApplicationBase& application)
            :
            FuelCell::ApplicationCore::ApplicationWrapper(application)
            {
                solution = FuelCell::ApplicationCore::FEVector();
            }
            
            /**
             * Constructor
             */
            AdaptiveRefinement ( FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>& app_lin,
                                 FuelCell::ApplicationCore::ApplicationWrapper& app,
                                 const FuelCell::ApplicationCore::FEVector& solution = FuelCell::ApplicationCore::FEVector() );
            
            /**
             * 
             */
            ~AdaptiveRefinement();
            
            /**
             * Declare all parameters that are needed for:
             *   - control of adaptive refinement loop
             *   - control of output options during the loop, e.g. return the solution at each cycle
             *
             * Currently the options implemented are:
             *
             * @code
             * subsection Adaptive refinement
             *   set Number of Refinements = 1                   # Number of adaptive refinements used
             *   set Output initial mesh = false                 # Set flag to true if you want to output a EPS figure of the initial mesh using the value in file initial mesh
             *   set Output initial mesh filename = initial_mesh # File where the initial mesh will be output
             *   set Output initial solution = false             # Set flag to true if you want to output the initial solution to file
             *   set Initial solution filename = initial_sol     # File where the initial solution will be output
             *   set Output intermediate solutions = true        # Output the solution after each adaptive refinement level
             *   set Output intermediate responses = true        # Output any functional such as current density at each refinement level
             *   set Output final solution = true                # Output final solution to a file named fuel_cell-sol-cycleX where X is the cycle name.
             *   set Output solution for transfer = false        # Check whether a solution on a refined grid should be output to file on a coarse mesh
             *   set Read in initial solution from file = false  # Check whether a stored solution should be read from file and applied to the grid
             * end
             * @endcode
             *
             */
            void declare_parameters ( ParameterHandler& param ) const;
            
            /**
             * Set up how many equations are needed and
             * read in parameters for the parameter handler in order to initialize data
             */
            void initialize ( ParameterHandler& param );
            
            ///@name Solving functions
            //@{
            
            /**
             * Solve the nonlinear problem
             */
            void solve( const std::string param_file,
                        ParameterHandler& param,
                        bool              solution_component_changes_between_data_files = false );
            
            /**
             * Run application
             */
            void run_app(bool solution_component_changes_between_data_files = false);
            
            /**
             * Run application
             */
            void run_app ( std::vector<double>& resp,
                              bool solution_component_changes_between_data_files = false );
            
            /**
             * Run application
             */
            void run_app ( std::vector<double>& resp,
                              std::vector<std::vector<double> >& dresp_dl,
                              bool solution_component_changes_between_data_files = false );
            //@}
            
            /**
             * Member function used to test the derivatives
             */
            void test_derivatives ( const std::string input_file,
                                    const std::string dvar,
                                    const double value,
                                    std::vector<double>& resp,
                                    std::vector<std::vector<double> >& dresp,
                                    const bool gradient = true );
            
            /**
             * Print parameters:
             */
            void print_parameters() const;
            
            /**
             * This function returns
             * \p solution.
             */
            const FuelCell::ApplicationCore::FEVector& get_solution() const
            {
                return solution;
            }
            
            /**
             * This function returns
             * \p coarse_solution.
             */
            const FuelCell::ApplicationCore::FEVector& get_coarse_solution() const
            {
                return coarse_solution;
            }
            
            /**
             * This function returns
             * \p coarse_triangulation.
             */
            const Triangulation<dim>& get_coarse_triangulation() const
            {
                return coarse_triangulation;
            }
            
        private:
            /** */
            void print_convergence_table();
            
            /** 
             * Filename where to output the initial grid 
             */
            std::string filename_initial_mesh;
            
            /** Flag to specify if the initial grid should be saved to a file.
             * This flag is initialized via the input file using:
             * @code
             * subsection Adaptive refinement
             *   set Output initial mesh = false                 # Set flag to true if you want to output a EPS figure of the initial mesh using the value in file initial mesh
             * end
             * @endcode 
             */
            bool output_initial_mesh;
            

            
            /** 
             * Boolean flag that is set to true if intermediate solutions, i.e. solutions at
             * intermediate adaptive refinement levels should be saved to file. This flag is set
             * via input file with:
             * @code
             * subsection Adaptive refinement
             *  set Output intermediate solutions = false
             * end
             * @endcode
             */
            bool output_intermediate_sol;
            
            /** 
             * Boolean flag used to specify if the final solution should be output to a file.
             * This flag is set via the input file using:
             * @code
             * subsection Adaptive refinement
             *   set Output final solution = true                # Output final solution to a file named fuel_cell-sol-cycleX where X is the cycle name.
             * end
             * @endcode
             */
            bool output_final_sol;
            
            /**
             * Filename where the final solution will be output. The final format of the solution will be as follows:
             * filename_final_sol+"_Cycle_" + grid_streamOut.str() + "_Sol_" + sol_streamOut.str();
             * where Cycle refers to the adaptive refinement cycle and Sol refers to the point in the polarization curve.
             * 
             */
            std::string filename_final_sol;
            
            /** 
             * Bool flag used to specify if all responses/functionals should be evaluated at intermediate
             * refinement levels. This is useful when trying to estimate if the solution is grid independent
             * as eventually the functionals and solution should become independent of refinement level.
             * 
             * This flag is set via the input file using:
             * @code
             * subsection Adaptive refinement
             *   set Output intermediate responses = true        # Output any functional such as current density at each refinement level
             * end
             * @endcode
             */
            bool output_intermediate_resp;
                                  
            /**
             * Use \p true if you have the exact
             * or analytical solution to compare
             * your numerical results with.
             */
            bool L1_L2_error_and_convergence_rate;
            
            /**
             * Use \p true along with
             * zero initial guess defined
             * in your application to solve
             * a linear problem by means of
             * using any of the FCST nonlinear
             * solvers.
             */
            bool nonlinear_solver_for_linear_problem;
            
            /**
             * Number of initial refinements for the original mesh
             *
             * \todo3 Change the loop so that zero refinements means once through the loop
             */
            unsigned int n_ref;
            
            /** 
             * Compute the gradients?
             */
            bool gradients;    
            
            /**
             * Initial mesh.
             */
            Triangulation<dim> coarse_triangulation;
            
            /** 
             * Pointer to linear application 
             */
            FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim> *app_linear;
            
            /** 
             * Poiner to nonlinear application 
             */
            FuelCell::ApplicationCore::ApplicationWrapper *app;
            
            /**
             * Global FE solution at the support points of the most refined mesh.
             */
            FuelCell::ApplicationCore::FEVector solution;
            
            /**
             * Global FE solution at the support points of the initial mesh.
             */
            FuelCell::ApplicationCore::FEVector coarse_solution;
            
            /**
             * If the exact or analytical solution is available,
             * then those convergence tables collect all necessary information
             * on L1 and L2 norms of the error and convergence rates
             * for each particular component of the numerical solution.
             */
            std::vector< ConvergenceTable > convergence_tables;
        };
    }
}

#endif