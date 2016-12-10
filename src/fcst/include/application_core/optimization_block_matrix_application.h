// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: optimization_block_matrix_application.h
// - Description: This is a base class for all applications that provide optimization information
// - Developers: Marc Secanell, University of Alberta
// - $Id: optimization_block_matrix_application.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _FUEL_CELL_APPLICATION_CORE_OPTIMIZATION_BLOCK_MATRIX_APPLICATION_H_
#define _FUEL_CELL_APPLICATION_CORE_OPTIMIZATION_BLOCK_MATRIX_APPLICATION_H_

//-- dealII
#include <deal.II/base/convergence_table.h>
#include <deal.II/grid/persistent_tria.h>
#include <deal.II/dofs/dof_renumbering.h>

//-- OpenFCST
#include <utils/fcst_units.h>
#include <utils/fcst_constants.h>
#include <application_core/block_matrix_application.h>

namespace FuelCell
{
namespace ApplicationCore
{
        /**
         * Application handling matrices and assembling the linear system
         * to solve the sensitivity equations.
         *
         * This class inherits the information on Triangulation and DoFHandler
         * from its base class BlockMatrixApplication. It adds information on the
         * linear system for sensitivity analysis
         *
         * <h3> Desired changes </h3>
         * - Feb. 2014. All responses and design variables should be stored in a std::map<Enumeration, double> where
         * the enumeration corresponds to a string that is unique and used in Dakota.
         *
         * @author Marc Secanell, 2006
         */
        template <int dim>
        class OptimizationBlockMatrixApplication :
        public FuelCell::ApplicationCore::BlockMatrixApplication<dim>
        {
        public:
            /**
             * The Event set by OptimizationMatrixApplication if
             * if we want to compute the sensitivities
             * and a new matrix should be
             * assembled.
             */
            static const FuelCell::ApplicationCore::Event sensitivity_analysis;

            /**
             * Constructor for an object,
             * owning its own mesh and dof handler.
             *
             * Since this class will be used by MeshWorker::AssemblingIntegrator, it
             * needs the functions for cell, boundary and interior face integration specified
             * exactly as below.
             *
             * The base class Subscriptor (deal.ii class) is needed so that MeshWorker::AssemblingIntegrator
             * can store a SmartPointer to an object of this class
             *
             * @author Marc Secanell, 2009-13
             */
            OptimizationBlockMatrixApplication(FuelCell::ApplicationCore::DoFApplication<dim>&,
                                               bool triangulation_only);

            /**
             * Constructor for an object,
             * owning its own mesh and dof handler.
             *
             * If <tt>data</tt> is null,
             * a new handler for
             * application data is
             * generated.
             *
             * see constructor of DoFApplication.
             */
            OptimizationBlockMatrixApplication(boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data =  boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData>());


            /** Destructor */
            ~OptimizationBlockMatrixApplication(){};

            /**
             * Declare all parameters that are needed for:
             *   - the computation of the equation coefficients
             *   - the control of the linear system solution
             *   - ...
             */
            virtual void declare_parameters(ParameterHandler& param);

            /**
             */
            virtual void initialize(ParameterHandler& param);

            /**
             * Post-processing. Evaluate all responses, usually functionals, necessary
             * such as the objective function and the constraints of the
             * optimization problem
             * @note size of f is num_objectives + num_constraints
             */
            virtual void responses (std::vector<double>& f,
                                    const FuelCell::ApplicationCore::FEVectors& vectors);
            /**
             * This class is called by responses to make sure that all responses requested are
             * implemented in either cell_responses, global_responses or face_responses.
             */
            virtual void check_responses();

            /**
             * This class is used to evaluate all the functionals that require looping over cells.
             * This class is called by responses to evaluate the response at each cell.
             *
             */
            virtual void cell_responses (std::vector<double>& resp,
                                         const typename DoFApplication<dim>::CellInfo& info,
                                         const FuelCell::ApplicationCore::FEVector& sol);

            /**
             * This class is used to evaluate all the functionals that require looping over boundaries.
             * This class is called by responses to evaluate the response at each boundary.
             *
             */
            virtual void bdry_responses (std::vector<double>& resp,
                                         const typename DoFApplication<dim>::FaceInfo& info,
                                         const FuelCell::ApplicationCore::FEVector& sol);

            /**
             * This class is used to evaluate all responses that do not require looping over cells.
             * An example of one of this types of constraints is the solid volume fraction.
             */
            virtual void global_responses (std::vector<double>& resp,
                                           const FuelCell::ApplicationCore::FEVector& sol);

            /** This function is used to print the responses to screen (FcstUtilities::log) */
            virtual void print_responses (std::vector<double>& resp);


            /**
             * Post-processing. Evaluate the partial derivative
             * of all functionals necessary such as the objective function and the constraints of the
             * optimization problem. In this member class we obtain the term:
             * \f[
             * \frac{\partial f_m}{\partial \lambda_k}
             * \f]
             * @note size of df_dl is for f num_objectives + num_constraints
             * and for lambda = num_design_var
             *
             * TODO Needs to be tested again (2013).
             */
            virtual void dresponses_dl (std::vector<std::vector<double> >& df_dl,
                                        const FuelCell::ApplicationCore::FEVectors& src);
            /**
             * This class is used to evaluate the derivative of all the functionals that require looping over cells
             * with respect to the design variables.
             * This class is called by responses to evaluate the response at each cell.
             */
            virtual void cell_dresponses_dl(std::vector<std::vector<double> >& cell_df_dl,
                                            const typename DoFApplication<dim>::CellInfo& info,
                                            const FuelCell::ApplicationCore::FEVector& src);
            /**
             * This class is used to evaluate the sensitivities of all responses that do not require looping over cells
             * with respect of the design variables.
             * An example of one of this types of constraints is the solid volume fraction.
             */
            virtual void global_dresponses_dl(std::vector<std::vector<double> >& df_dl,
                                              const FuelCell::ApplicationCore::FEVector& src);
            /**
             * Loop over all cells
             * and assemble the vector
             * \f[
             * \frac{\partial f_m}{\partial u_i}
             * \f]
             * that is used to solve the sensitivity
             * equations by using the local
             * matrix integration functions
             * cell_dfunctional_du(), bdry_dfunctional_dlu()
             * and edge_dfunctinal_du() provided
             * by the derived class.
             * The vector index refers to m and the BlockVector
             * index refers to i.
             */
            virtual void dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector>& dst,
                                       const FuelCell::ApplicationCore::FEVectors& src);

            /**
             * Integration of local
             * bilinear form.
             */
            virtual void cell_dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector>& df_du,
                                            const typename DoFApplication<dim>::CellInfo& info,
                                            std::vector<std::vector<double> >& src);

            /**
             * This class is used to evaluate the sensitivities of all responses that do not require looping over cells
             * with respect of the design variables.
             * An example of one of this types of constraints is the solid volume fraction.
             */
            virtual void global_dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector>& df_du,
                                              const FuelCell::ApplicationCore::FEVector& src);

            /**
             * Loop over all cells
             * and assemble the vector
             * \f[
             * \frac{\partial R_n}{\partial \lambda_k}
             * \f]
             * that is used to solve the sensitivity
             * equations by using the local
             * matrix integration functions
             * cell_dresidual_dlambda(), bdry_dresidual_dlambda()
             * and edge_dresidual_dlambda() provided
             * by the derived class.
             * The vector index refers to k and the FEVector
             * index refers to n.
             */
            virtual void dresidual_dlambda(std::vector<FuelCell::ApplicationCore::FEVector>& dst,
                                           const FuelCell::ApplicationCore::FEVectors& src);
            /**
             * Integration of local
             * bilinear form.
             */
            virtual void cell_dresidual_dlambda(std::vector<FuelCell::ApplicationCore::FEVector>& cell_vector,
                                                const typename DoFApplication<dim>::CellInfo& cell,
                                                std::vector<std::vector<double> >& src);

            /**
             * Solver in order to obtain the
             * analytical sensitivities for the
             * given objective and constraints using
             * the direct method. This member function
             * calls the different assemble member functions and solves the
             * system of equations
             * @note The direct method is the most computationally efficient
             * method if the number of objectives and constraints
             * is larger than the number of design variables
             */
            void solve_direct (std::vector<std::vector<double> >& df_dl,
                               const FuelCell::ApplicationCore::FEVectors& sol);

            /**
             * Solver in order to obtain the
             * analytical sensitivities for the
             * given objective and constraints using
             * the adjoint method. This member function
             * calls the different assemble member functions and solves the
             * system of equations
             * @note The adjoint method is the most computationally efficient
             * method if the number of number of design variables
             * is larger than the number of objectives and constraints
             */
            void solve_adjoint (std::vector<std::vector<double> >& df_dl,
                                const FuelCell::ApplicationCore::FEVector& sol);
            /**
             * Member function that returns the number of responses
             */
            unsigned int get_n_resp() const;
            /**
             * Member function that returns the number of design variables
             */
            unsigned int get_n_dvar() const;
            /**
             * Member function that returns the name of the design variables
             */
            inline std::vector<std::string> get_name_dvar() const
            {
                return name_design_var;

            }

            /**
             * Member function that returns the name of the responses
             */
            inline std::vector<std::string> get_name_responses() const
            {
                return name_responses;

            }

            /**
             * Function used to return all possible response names in OpenFCST.
             */
            const std::vector<std::string> get_all_responses_names() const
            {
                return all_response_names;
            }
            /**
             *
             */
            /*void enter_data_scalar(std::string name,
             * const double& scalar)
             { *
                 this->get_data()->enter(name, scalar);
             }*/

            /**
             * Print the default parameter handler file
             */
            void print_default_parameter_file()
            {
                ParameterHandler prm;
                this->declare_parameters(prm);
                std::ofstream file;
                file.open("Default_parameter_file.prm");
                prm.print_parameters (file, ParameterHandler::Text);
                file.close();
            }

            /** Function to see if we are transferring a solution on a refined grid to the initial coarse grid */
            bool get_bool_transfer_solution()
            {return output_coarse_solution;}

            /**
             * This routine assigns the value of n_dvar, n_resp, name_design_var and name_responses to the
             * internal variables n_dvar, n_resp, name_design_var and name_response
             *
             */
            void set_optimization_parameters(unsigned int &n_dvar, unsigned int &n_resp,
                                             std::vector<std::string> &name_design_var, std::vector<std::string> &name_responses);

            /**
             *
             */
            void set_output_variables(std::vector<std::string> &dakota_name_responses)
            {
                name_responses = dakota_name_responses;
            };

            /**
             * If the exact or analytical solution is available,
             * then this function computes L1 and L2 norms of
             * the error and convergence rates
             * for each particular component of the numerical solution.
             *
             * This function must be overriden in the actual application classes,
             * see for example \p AppDarcyCTest class.
             *
             * @param solution           - numerical solution,
             * @param refinement_cycle   - refinement cycle used in \p AdaptiveRefinement class,
             * @param convergence_tables - convergence tables that collect all necessary information.
             */
            virtual void compute_L1_L2_error_and_convergence_rate(const FuelCell::ApplicationCore::FEVector&        solution,
                                                                  const unsigned int&              refinement_cycle,
                                                                  std::vector< ConvergenceTable >& convergence_tables) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__ << " called in Class " << info.name() << std::endl;

            }

        protected:
            /** Auxiliary routine to print the values df_dl
             * This routine should be called once df_df is assembled. This is done
             * on solve_direct
             */
            void print_dresponses_dl(std::vector<std::vector<double> > pdf_pdl)
            {
                for (unsigned int i=0; i<n_resp; ++i)
                    for (unsigned int j=0; j<n_dvar; ++j)
                        FcstUtilities::log<<"Df"<<i<<"/Dl"<<j<<" is equal to "<<pdf_pdl[i][j]<<std::endl;
            }

            /** Auxiliary routine to print the values of df_du
             * This routine should be called once df_du is assembled. This is done
             * on solve_direct
             */
            void print_dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector> df_du)
            {
                // Write partial df_du:

                for (unsigned int r=0; r<n_resp; ++r)
                {
                    // Output Du block by block
                    unsigned int n_blocks = df_du[r].n_blocks();
                    unsigned int size = df_du[r].size();
                    unsigned int comp_per_block = size/n_blocks;
                    for (unsigned int i=0; i<n_blocks; ++i)
                        for (unsigned int j=0; j<comp_per_block; ++j)
                        {
                            FcstUtilities::log<<"Element "<<j<<" in block "<<i<<" in STEP (Du) is "<<df_du[r].block(i)(j)<<std::endl;
                        }
                        FcstUtilities::log<<"NORM OF Df"<<r<<"/Du = "<<df_du[r].l2_norm()<<std::endl;
                }
            }

            /**
             * Member function used in order to extend the name of a file with the design variable name and value.
             * 
             * This is useful when doing visualization with ParaView.
             * 
             */
            const std::string extend_filename(const std::string& , const int precision = 3) const;
            
            /**
             * Definition of a pointer to a function that is called in a loop over
             * cells to compute the derivatives of either residual or responses with
             * respect to design variables.
             */
            typedef void (*CELL_Dvalues) (std::vector<FuelCell::ApplicationCore::FEVector>& ,
                                          const typename DoFApplication<dim>::CellInfo& ,
                                          std::vector<std::vector<double> >& );

            /**
             * This is an auxiliary function called by dresidual_dlambda and
             * dfunctional_du. Since both functions have to do the same thing,
             * loop over the cells, but they call a different member function for cell assembly,
             * I use a flag to change from one to the other and then I call this
             * generic function.
             * With this dresidual_dlambda and dfunctional_du have the same interface
             * and the changes happen to a protected class.
             */
            void dfunction (std::vector<FuelCell::ApplicationCore::FEVector>& dst,
                            const FuelCell::ApplicationCore::FEVectors& src,
                            bool dfunctional_du,
                            bool dresidual_dlambda);

            /**
             * Initialize the data member @p all_response_names with all the available Responses in OpenFCST
             * 
             */
            void set_all_response_names();

            /**
             * Number of design variables
             */
            unsigned int n_dvar;
            /**
             * Number of objective functions
             */
            unsigned int n_obj;
            /**
             * Number of responses = n_obj + n_con
             */
            unsigned int n_resp;
            /**
             * Member that stores the name of the design variables
             */
            std::vector<std::string> name_design_var;
            /**
             * Member that stores the name of the responses, i.e. objectives and constraints.
             * @note first we put the name of all objectives and then all constraints
             */
            std::vector<std::string> name_responses;

            /**
             * Member that stores the name of the output variables,
             * These names will be written to name_responses if optimization is not being used and ignored otherwise.
             */
            std::vector<std::string> name_output_var;

            /** Decision variable for whether the solution is to be output for transfer */
            bool output_coarse_solution;

            /** Decision variable for whether the application is being used for optimization */
            bool optimization;

            /**
             * @p true, if boundary responses are supposed to be computed.
             */
            bool boundary_responses;
            
            /**
             * Variable that holds a list of response names
             */
            std::vector <std::string> all_response_names;
            
            /**
             * Boundary id on which the boundary response is computed
             */
            std::vector<unsigned int> user_input_bdry;
        };

}
}

#endif