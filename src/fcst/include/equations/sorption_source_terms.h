//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: sorption_source_terms.h
//    - Description: This class is used to assemble cell matrix and cell residual
//      corresponding to sorption/desorption of water inside the catalyst layer.
//    - Developers: Marc Secanell, Madhur Bhaiya
//
//---------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_EQUATION_SORPTION_SOURCE_TERMS_H_
#define _FCST_FUELCELLSHOP_EQUATION_SORPTION_SOURCE_TERMS_H_

#include <equations/equation_base.h>
#include <layers/catalyst_layer.h>


namespace FuelCellShop
{
    namespace Equation
    {  
        ///@name Exceptions
        //@{
        /**
         * Exception thrown when a particular "Required" solution variable
         * is not found for sorption source terms,
         * in \p user-defined solution variables.
         */
        DeclException2(VariableNotFoundForSorption,
                       std::string,
                       std::string,
                       << arg1 << " should be one of the solution variables, in order to account for " << arg2 << " .");
        
        //@}    
        
        /**
         * This class assembles source terms corresponding to sorption/desorption of water inside the catalyst layer. Here, a vapor equilibriated membrane model is being used.
         * Under steady state operation, \f$ x_{H_2O} \f$, water vapor (gas) is in equilibrium with \f$ \lambda \f$, (sorbed water). This coupling is provided by the following term:
         * 
         * \f$ S_{\lambda} = k \frac{\rho_{dry}}{EW} ( \lambda_{eq} - \lambda ) \quad \in \quad \Omega \f$
         * - where this term is positive/negative when used in appropriate conjunction with <b>Membrane Water Content Transport Equation</b> and
         * <b>Ficks Transport Equation - water</b>, and can act as a source/sink [\p moles/\p(cm^3-s \p)].
         * - \f$ (\lambda_{eq} - \lambda) \f$, represents the difference between the equilibrium value for sorbed water and the actual value.
         * - \f$ \rho_{dry} \f$ is dry polymer electrolyte material density [\p gm/cm^3].
         * - \f$ EW = \frac{mass~of~dry~polymer~electrolyte~material~in~grams}{moles~of~SO_3^-} \f$ is equivalent weight.
         * - \f$ k \f$ = time constant [\p 1/s]. Its value can be provided using the parameter file. It is recommended to use the default
         * value, \em i.e. \b 10000.0, to guarantee effective coupling.
         * 
         * \remarks
         * - This is technically <b>not an Equation</b> class. It is only used to assemble source terms in conjunction with other equation classes.
         * - The coupling term is valid only for <b>Membrane Water Content Transport Equation</b> and <b>Ficks Transport Equation - water</b>, 
         * hence it is necessary to have both equations being solved for while using this class.
         * - Besides this, it also has the provision for accounting for heat relase during the sorption/desorption process. This can be enabled by setting the flag to \em TRUE in the parameter file.
         * It is necessary to have <b>Thermal Transport Equation</b> being solved for when the flag is set to \em TRUE, otherwise code will throw an \p error.
         * - Enthalpy of sorption [\p J/mol] is retrieved from Material::Nafion class, and can be modified using the parameter file entries corresponding to Nafion class.
         * - This class works with the following layer class only:
         *      - FuelCellShop::Layer::CatalystLayer<dim>
         * - <b>VERY IMPORTANT</b>: If this class is being considered, it's very important to use  #adjust_internal_cell_couplings
         * method of this class, before using  \p make_cell_couplings of \b SystemManagement at the application level.
         * 
         * <h3>Usage details</h3>
         * 
         * @code
         * // Creating source terms object (in Application Header file)
         * FuelCellShop::Equation::SorptionSourceTerms<dim> sorption_source;
         * 
         * // Declare parameters in application
         * sorption_source.declare_parameters(param);
         * 
         * // Initialize in application
         * sorption_source.initialize(param);
         * 
         * // Create a temporary vector in the application for storing couplings_map from all the equation used in the application.
         * std::vector<couplings_map> tmp;
         * 
         * ... // get_internal_cell_couplings from all the equations
         * 
         * sorption_source.adjust_internal_cell_couplings(tmp);
         * 
         * // Making cell couplings using SystemManagement object created in the application
         * system_management.make_cell_couplings(tmp);
         * 
         * // cell_matrix in application
         * // Do a check against layer and it should match with the layers currently working for this equation class.
         * // for eg: CCL is FuelCellShop::Layer::HomogeneousCL<dim> object.
         * sorption_source.assemble_cell_matrix(cell_matrices, cell_info, &CCL);
         * 
         * // cell_residual in application
         * sorption_source.assemble_cell_residual(cell_vector, cell_info, &CCL);
         * @endcode 
         * 
         * \author Madhur Bhaiya, 2013
         */
        template<int dim>
        class SorptionSourceTerms : public EquationBase<dim>
        {
        public:
            
            ///@name Constructors, destructor, and initalization
            //@{
            
            /**
             * Constructor.
             */
            SorptionSourceTerms(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = 
            boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >());
            
            /**
             * Destructor.
             */
            virtual ~SorptionSourceTerms();
            
            /**
             * Declare parameters.
             */
            virtual void declare_parameters(ParameterHandler& param) const;
            
            /**
             * Initialize parameters.
             */
            virtual void initialize(ParameterHandler& param);
            
            //@}
            
            ///@name Local CG FEM based assemblers
            //@{
            
            /**
             * Assemble local cell matrix.
             */
            virtual void assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector& cell_matrices,
                                              const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                              FuelCellShop::Layer::BaseLayer<dim>* const layer);
            
            /**
             * Assemble local cell residual.
             */
            virtual void assemble_cell_residual(FuelCell::ApplicationCore::FEVector& cell_residual,
                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                FuelCellShop::Layer::BaseLayer<dim>* const layer);
            
            //@}
            
            ///@name Accessors & Info
            //@{
            
            /**
             * This function is used to adjust \p std::vector \p < \p internal_cell_couplings \p >, which is generated after getting \p internal_cell_couplings
             * from all the equations being used in the application.
             * 
             * \note It's very important to use this function, if we are considering source terms due to water sorption/desorption in the catalyst layers. It's also noteworthy that
             * this function should be used after whole vector of \p internal_cell_couplings is created, and \b BEFORE \p make_cell_couplings of \b SystemManagement
             * is called.
             */
            virtual void adjust_internal_cell_couplings(std::vector< couplings_map >& ) const;
            
            /**
             * This function prints out
             * the \p info for this class.
             */
            virtual void print_equation_info() const;
            
            /**
             * Returns time constant \f$ k \f$ [\p 1/s].
             */
            double get_time_constant() const
            {
                return time_constant;
            }
            
            /**
             * Method to get whether the heat release/absorption due to sorption/desorption, inside the catalyst layers is \p ON or \p OFF (#flag_sorp_heat_cl).
             * The typical use of this method is for the \p \b post-processing routines in the application.
             */
            inline bool get_flag_sorp_heat_cl() const
            {
                return flag_sorp_heat_cl;
            }
            
            //@}
            
        protected:
            
            ///@name Boolean flags for Reaction heat source terms
            //@{
            
            /**
             * This boolean data member indicates that the heat release/absorption,
             * due to sorption/desorption inside the catalyst layers is \p ON or \p OFF.
             */
            bool flag_sorp_heat_cl;
            
            //@}
            
            ///@name Local CG FEM based assemblers - make_ functions
            //@{
            
            /**
             * This function computes Local CG FEM based
             * assemblers - constant data (generic).
             * 
             * \warning <b>For the DEVELOPERS:</b> \p block_index for VariableInfo are not filled here, but \p indices_exist flag is set to \b TRUE (as it
             * is needed at other places before \p cell_matrix or \p cell_residual assembly).
             * \p Developers need to be wary of this fact that the \p block_indices are still not filled yet, and they need to <b>fill them before</b> doing \p cell_matrix assembly.
             */
            virtual void make_assemblers_generic_constant_data();
            
            /**
             * This function computes
             * <p> Local CG FEM based assemblers - constant data (cell) </p>
             * and allocates the memory for
             * shape functions, shape function gradients, and 
             * \p JxW_cell in
             * <p> Local CG FEM based assemblers - variable data (cell) </p>.
             */
            virtual void make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info);
            
            /**
             * This function computes
             * <p> Local CG FEM based assemblers - variable data (cell) </p>.
             */
            virtual void make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                            FuelCellShop::Layer::BaseLayer<dim>* const layer);
            
            /**
             * This function is specifically created for assembly of cell matrices for the following equations, \em viz.
             * - <b>"Ficks Transport Equation - water"</b>
             * - <b>"Membrane Water Content Transport Equation"</b>
             * - <b>"Thermal Transport Equation</b>
             * 
             * Cell matrix terms for the abovementioned equations follow a similar pattern, using different factors; hence this function avoids
             * using repeated lines of code for each of these equations in  #assemble_cell_matrix method. The argument list is defined as follows:
             * - first corresponds to \p FuelCell::ApplicationCore::MatrixVector which we need to fill.
             * - second corresponds to typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo object; used to access fe elements for variables other than test function.
             * - second corresponds to \p std::string for the name of the equation for which we are doing the assembly. It should match with names mentioned above.
             * - third corresponds to \p FEValuesBase (fe discretization) used for the test function \em i.e. equation for which assembly is being done.
             * - fourth corresponds to shape functions for the test function, filled already by #make_assemblers_cell_variable_data.
             * - fifth corresponds to the factor used to multiply with the source term \f$~ \mathbf{\delta} \left( \frac{k \rho_{dry}}{EW} ( \lambda_{eq} - \lambda ) \right) \f$ 
             * for <b>Membrane Water Content Transport Equation</b> and <b>Ficks Transport Equation - water</b>.
             * 
             * \warning <b>For the DEVELOPERS:</b> This function should only be used after #make_assemblers_cell_variable_data has been called in #assemble_cell_matrix, otherwise it is bound to give wrong results.
             */            
            virtual void assemble_matrix_for_equation(FuelCell::ApplicationCore::MatrixVector& cell_matrices,
                                                      const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                      const std::string& eq_name,
                                                      const FEValuesBase<dim>& test_fe,
                                                      const std::vector< std::vector<double> >& test_shape_functions,
                                                      const double& sourceterm_factor);
            //@}
            
            ///@name Generic Constant Data
            //@{
            
            /**
             * VariableInfo structure corresponding to \p "water_molar_fraction".
             */
            VariableInfo x_water;
            
            /**
             * VariableInfo structure corresponding to \p "membrane_water_content".
             */
            VariableInfo lambda;
            
            /**
             * VariableInfo structure corresponding to \p "temperature_of_REV".
             */
            VariableInfo t_rev;
            
            /**
             * Time constant, \f$ k \f$ [\p 1/s]
             */
            double time_constant;
            
            //@}
            
            ///@name Local CG FEM based assemblers - variable data (cell)
            //@{
                        
            /**
             * Density [\p gm/cm^3] of the dry polymer electrolyte material in the cell.
             */
            double rho_dry_cell;
            
            /**
             * Equivalent weight of the polymer electrolyte material in the cell.
             */
            double EW_cell;
            
            /**
             * \f$ \mathbf{\lambda_{eq}} \f$ from sorption isotherm,
             * at all quadrature points in the cell.
             */
            std::vector<double> lambda_eq_cell;
            
            /**
             * Derivative of \f$ \mathbf{\lambda_{eq}} \f$ w.r.t \p "water_molar_fraction",
             * at all quadrature points in the cell.
             */
            std::vector<double> dlambdaEq_dxWater_cell;
            
            /**
             * Derivative of \f$ \mathbf{\lambda_{eq}} \f$ w.r.t \p "temperature_of_REV",
             * at all quadrature points in the cell.
             */
            std::vector<double> dlambdaEq_dT_cell;
            
            /**
             * Enthalpy of sorption of water, \f$ \Delta h_{sorp} \f$ [\p J/mol], 
             * at all quadrature points in the cell.
             */
            std::vector<double> h_sorp_cell;
            
            /**
             * Derivative of \f$ \Delta h_{sorp} \f$ w.r.t \p "temperature_of_REV", 
             * at all quadrature points in the cell.
             */
            std::vector<double> dhsorp_dT_cell;
            
            /**
             * \f$ \mathbf{T} \f$ shape functions.
             * 
             * \p phi_T_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{T} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_T_cell;
            
            /**
             * \f$ \mathbf{\lambda} \f$ shape functions.
             * 
             * \p phi_lambda_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{\lambda} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_lambda_cell;
            
            /**
             * \f$ \mathbf{x_{H_2O}} \f$ shape functions.
             * 
             * \p phi_xWater_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{x_{H_2O}} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_xWater_cell;
            
            //@}
            
            /**
             * Counter set to \em TRUE when \p cell_residual is being assembled.
             * This ensures that \b derivatives of source terms are computed only when \p cell_matrix is being assembled, in #make_assemblers_cell_variable_data.
             */
            bool cell_residual_counter;
            
            /**
             * Variable used to store the index in cell_info->global_data of the previous Newton solution
             * The solution at the previous iteration is used to compute cell_matrix and cell_residual
             */  
             unsigned int last_iter_cell;            
            
        };
             
    } // Equation
             
} // FuelCellShop

#endif