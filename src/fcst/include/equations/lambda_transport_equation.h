//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: lambda_transport_equation.h
//    - Description: Equation class for Lambda Transport (using Springer Model + Thermal Osmosis)
//    - Developers: Marc Secanell, Madhur Bhaiya, Valentin N. Zingan
//
//---------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_EQUATION_LAMBDA_TRANSPORT_EQUATION_H_
#define _FCST_FUELCELLSHOP_EQUATION_LAMBDA_TRANSPORT_EQUATION_H_

#include <equations/equation_base.h>

#include <layers/catalyst_layer.h>
#include <layers/membrane_layer.h>

// STD
#include <sstream>
#include <string>

namespace FuelCellShop
{
    namespace Equation
    {
        /**
         * This class deals with <b>Membrane Water Content Transport Equation</b>.
         *
         * This equation class solves for water transport inside the membrane using semi-empirical model proposed by Springer et al. According to Springer model,
         * sorbed water \f$( \lambda )\f$ inside the membrane is transported via electro-osmotic drag and fickian diffusion.
         * Besides this, under nonisothermal conditions, water flows from cold to hot side causing thermal osmosis.
         * These three transport modes are currently being considered in this equation class.
         *
         * It is to be noted that this model is suitable for water-vapor equilibriated membranes or isobaric membranes. Under such conditions, convective transport
         * can be neglected.
         *
         * It is solved with respect to:
         * - \f$ \mathbf{\lambda} \f$ \p(membrane_water_content \p)
         *
         * where, \f$ \lambda = \frac{moles~of~sorbed~water}{moles~of~SO_3^-} \f$.
         *
         * This equation can be written as:
         *
         *  \f$ \qquad \mathbf{\nabla} \cdot \left( \frac{n_d}{F} \hat{\sigma}_{m,eff} \mathbf{\nabla} \phi_m + \frac{\rho_{dry}}{EW} \hat{D}_{\lambda,eff} \mathbf{\nabla} \lambda + \frac{1}{M_{water}} \hat{D}_{T,eff} \mathbf{\nabla} T \right) = 0 \quad \in \quad \Omega \f$
         * - where, first term corresponds to <b>electro-osmotic drag</b>, second term corresponds to <b>fickian diffusion</b> and
         * third term corresponds to <b>thermo-osmotic diffusion</b>.
         * - \f$ \hat{\sigma}_{m,eff} \f$ is effective proton conductivity tensor [\p S/cm], which can be function of other variables,
         *  \em viz., \f$ \lambda \f$ \p(membrane_water_content \p) and \f$ T \f$ \p(temperature_of_REV \p).
         * - \f$ \hat{D}_{\lambda,eff} \f$ is effective lambda (sorbed water) diffusivity tensor [\p cm^2/s], which can be function of other variables,
         *  \em viz., \f$ \lambda \f$ \p(membrane_water_content \p) and \f$ T \f$ \p(temperature_of_REV \p).
         * - \f$ \hat{D}_{T,eff} \f$ is effective thermo-osmotic diffusivity tensor [\p gm/\p(cm-s-K \p)], which can be a function of \f$ T \f$ \p(temperature_of_REV \p).
         * - \f$ n_d \f$ is electro-osmotic drag coefficient.
         * - \f$ F \f$ is Faraday's constant.
         * - \f$ \rho_{dry} \f$ is dry polymer electrolyte material density [\p gm/cm^3].
         * - \f$ EW = \frac{mass~of~dry~polymer~electrolyte~material~in~grams}{moles~of~SO_3^-} \f$ is equivalent weight.
         * - \f$ M_{water} \f$ is molar weight of water [\p gm/moles]
         * - \f$ \mathbf{\phi_{m}} \f$ \p(protonic_electrical_potential \p)
         *
         * To be well-posed, these equations are equipped with the appropriate boundary conditions. All the boundary conditions can be described by
         * \p boundary_id \p(s \p) and \p boundary_type. Besides, this some boundary types require additional information, which can also be provided by the parameter file.
         * We consider following types of boundary conditions:
         *
         * - <b>No water flux / Symmetric:</b> A particular case of \p Neumann boundary condition.
         *
         * \remarks
         * - There is no provision to specify boundary indicators for \p No \p water \p flux or \p Symmetric boundary conditions, as FEM
         * formulation automatically implies a particular boundary is one of these cases, by default.
         * - This class works with the following layer classes only:
         *      - FuelCellShop::Layer::CatalystLayer<dim>
         *      - FuelCellShop::Layer::MembraneLayer<dim>
         * - Thermo-osmotic diffusion is set to \p OFF by default, but can be turned \p ON based on the flag set in the parameter file.
         * - There is no direct flag for setting electro-osmotic drag \p ON or \p OFF, it is defaulted to \p ON, if \f$ \phi_m \f$ is also
         * one of the solution variables. But even in this case, it can be indirectly set to \p OFF by setting \f$ n_d \f$, electro-osmotic drag coefficient
         * to \b 0 in the parameter file for \p Nafion material class.
         *
         * We solve the whole problem by linearizing the governing equation at each Newton iteration with subsequent
         * CG FEM discretization in space. The class contains the necessary class members to add the necessary contributions to cell_matrix
         * and cell_residual to the governing equations used to analyze lambda transport,
         *
         * <h3>Usage Details:</h3>
         *
         * @code
         * // Creating Equation object (in Application Header file)
         * FuelCellShop::Equation::LambdaTransportEquation<dim> lambda_transport;
         *
         * // Declare parameters in application
         * lambda_transport.declare_parameters(param);
         *
         * // Initialize in application
         * lambda_transport.initialize(param);
         *
         * // Create a temporary vector in the application for storing couplings_map from all the equation used in the application.
         * std::vector<couplings_map> tmp;
         * ... // other equations
         * tmp.push_back( lambda_transport.get_internal_cell_couplings() );
         *
         * // Look at SorptionSourceTerms class here, if source terms due to water sorption/desorption are to be considered.
         * // Making cell couplings using SystemManagement object created in the application
         * system_management.make_cell_couplings(tmp);
         *
         * // cell_matrix in application
         * // Do a check against layer and it should match with the layers currently working for this equation class.
         * // for eg: CCL is FuelCellShop::Layer::HomogeneousCL<dim> object.
         * lambda_transport.assemble_cell_matrix(cell_matrices, cell_info, &CCL);
         *
         * // cell_residual in application
         * lambda_transport.assemble_cell_residual(cell_vector, cell_info, &CCL);
         * @endcode
         *
         *
         * \note This class doesn't assemble for water sorption/desorption source terms; that is taken care off by \p SorptionSourceTerms class.
         * Please read the documentation of SorptionSourceTerms class, for additional methods to be implemented in the application.
         *
         * \warning If water sorption/desorption source terms are being considered, it's very important to use \p adjust_internal_cell_couplings
         * member function of SorptionSourceTerms class, before using  \p make_cell_couplings of \b SystemManagement at the application level.
         *
         * \author Madhur Bhaiya,      2013
         * \author Valentin N. Zingan, 2012-2014 - afterward improvements, optimization, checkings, CG FEM bug fixings
         */

        template<int dim>
        class LambdaTransportEquation : public EquationBase<dim>
        {
        public:

            ///@name Constructors, destructor, and initalization
            //@{

            /**
             * Constructor.
             */
            LambdaTransportEquation(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = 
            boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >());

            /**
             * Destructor.
             */
            virtual ~LambdaTransportEquation();

            /**
             * Declare parameters.
             */
            virtual void declare_parameters(ParameterHandler& param) const;

            /**
             * Initialize parameters. This class will call #make_internal_cell_couplings and #make_boundary_types.
             */
            virtual void initialize(ParameterHandler& param);

            //@}
            ///@name Local CG FEM based assemblers

            //@{
            /**
             * Assemble local cell matrix.
             */
            virtual void assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                              const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                              FuelCellShop::Layer::BaseLayer<dim>* const              layer);

            /**
             * Assemble local cell residual.
             */
            virtual void assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_rhs,
                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                FuelCellShop::Layer::BaseLayer<dim>* const              layer);

            /**
             * Assemble local boundary matrix.
             * Currently, <b>NOT IMPLEMENTED</b>.
             */
            virtual void assemble_bdry_matrix(FuelCell::ApplicationCore::MatrixVector&                                 bdry_matrices,
                                              const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                              FuelCellShop::Layer::BaseLayer<dim>* const              layer){};

            /**
             * Assemble local boundary residual.
             * Currently, <b>NOT IMPLEMENTED</b>.
             */
            virtual void assemble_bdry_residual(FuelCell::ApplicationCore::FEVector&                                     bdry_rhs,
                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                FuelCellShop::Layer::BaseLayer<dim>* const              layer){};

            //@}

            ///@name Accessors and info
            //@{

            /**
             * The function prints out the equation's info.
             */
            virtual void print_equation_info() const;

            //@}

        protected:

            ///@name Boolean flags for lambda transport modes
            //@{

            /**
             * Flag to indicate that lambda (sorbed water) transport by Thermo-osmosis is ON or OFF.
             */
            bool flag_thermoosmosis;

            //@}

            ///@name Local CG FEM based assemblers - make_ functions
            //@{

            /**
             * This function computes Local CG FEM based
             * assemblers - constant data (generic).
             */
            virtual void make_assemblers_generic_constant_data();

            /**
             * This function computes
             * <b> Local CG FEM based assemblers - constant data (cell) </b>
             * and allocates the memory for
             * \p shape \p functions, \p shape \p function \p gradients, and
             * \p JxW_cell in
             * <b> Local CG FEM based assemblers - variable data (cell) </b>.
             */
            virtual void make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info);

            /**
             * This function computes
             * <b> Local CG FEM based assemblers - constant data (boundary) </b>
             * and allocates the memory for
             * \p shape \p functions, #normal_vectors, and
             * \p JxW_bdry in
             * <b> Local CG FEM based assemblers - variable data (boundary) </b>.
             *
             * Currently, <b>NOT IMPLEMENTED</b>.
             */
            virtual void make_assemblers_bdry_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info){};

            /**
             * This function computes
             * <b> Local CG FEM based assemblers - variable data (cell) </b>.
             */
            virtual void make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                            FuelCellShop::Layer::BaseLayer<dim>* const layer);

            /**
             * This function computes
             * <b> Local CG FEM based assemblers - variable data (boundary) </b>.
             *
             * Currently, <b>NOT IMPLEMENTED</b>.
             */
            virtual void make_assemblers_bdry_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                            FuelCellShop::Layer::BaseLayer<dim>* const layer){};

            //@}

            ///@name Other make_ functions
            //@{

            /**
             * This function fills out
             * \p internal_cell_couplings.
             */
            virtual void make_internal_cell_couplings();

            /**
             * This function fills out
             * \p boundary_types.
             *
             * Currently, <b>NOT IMPLEMENTED</b>.
             */
            virtual void make_boundary_types(){};

            /**
             * This function fills out
             * \p output_types.
             *
             * Currently, <b>NOT IMPLEMENTED</b>.
             */
            virtual void make_output_types()
            {};

            //@}

            ///@name Generic Constant Data
            //@{

            /**
             * VariableInfo structure corresponding to \p "protonic_electrical_potential".
             */
            VariableInfo phi_m;

            /**
             * VariableInfo structure corresponding to \p "membrane_water_content".
             */
            VariableInfo lambda;

            /**
             * VariableInfo structure corresponding to \p "temperature_of_REV".
             */
            VariableInfo t_rev;

            /**
             * Universal Faraday's constant
             */
            double F;

            /**
             * Molar weight of water in grams/mole.
             */
            double M_water;

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
             * Electro-osmotic drag coefficient, at all quadrature points in the cell.
             */
            std::vector<double> nDrag_cell;

            /**
             * Derivative of electro-osmotic drag coefficient w.r.t. \p "membrane_water_content",
             * at all quadrature points in the cell.
             */
            std::vector<double> dnDrag_dlambda_cell;

            /**
             * Effective proton conductivity [\p S/cm],
             * at all quadrature points of the cell.
             */
            std::vector<double> sigmaMeff_cell;

            /**
             * Derivative of effective protonic conductivity
             * w.r.t \p "temperature_of_REV",
             * at all quadrature points in the cell.
             */
            std::vector<double> dsigmaMeff_dT_cell;

            /**
             * Derivative of effective protonic conductivity
             * w.r.t. \p "membrane_water_content",
             * at all quadrature points in the cell.
             */
            std::vector<double> dsigmaMeff_dlambda_cell;

            /**
             * Effecive lambda (sorbed water) diffusivity [\p cm^2/s],
             * at all quadrature points in the cell.
             */
            std::vector<double> Dlambdaeff_cell;

            /**
             * Derivative of effective lambda diffusivity
             * w.r.t. \p "temperature_of_REV",
             * at all quadrature points in the cell.
             */
            std::vector<double> dDlambdaeff_dT_cell;

            /**
             * Derivative of effective lambda diffusivity
             * w.r.t. \p "membrane_water_content",
             * at all quadrature points in the cell.
             */
            std::vector<double> dDlambdaeff_dlambda_cell;

            /**
             * Effecive thermo-osmotic diffusion coefficient [\p gm/\p(cm-s-K)],
             * at all quadrature points in the cell.
             */
            std::vector<double> DTeff_cell;

            /**
             * Derivative of effective thermo-osmotic diffusion coefficient
             * w.r.t. \p "temperature_of_REV",
             * at all quadrature points in the cell.
             */
            std::vector<double> dDTeff_dT_cell;

            /**
             * \f$ \mathbf{\lambda} \f$ shape functions.
             *
             * \p phi_lambda_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{\lambda} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_lambda_cell;

            /**
             * \f$ \mathbf{\lambda} \f$ shape function gradients.
             *
             * \p grad_phi_lambda_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{\lambda} \f$ shape function gradient
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< Tensor<1,dim> > > grad_phi_lambda_cell;

            /**
             * \f$ \mathbf{T} \f$ shape functions.
             *
             * \p phi_T_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{T} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< double > > phi_T_cell;

            /**
             * \f$ \mathbf{T} \f$ shape function gradients.
             *
             * \p grad_phi_T_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{T} \f$ shape function gradient
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< Tensor<1,dim> > > grad_phi_T_cell;

            /**
             * \f$ \mathbf{\phi_m} \f$ shape function gradients.
             *
             * \p grad_phi_phiM_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{\phi_m} \f$ shape function gradient
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< Tensor<1,dim> > > grad_phi_phiM_cell;

            //@}

            /**
             * Counter set to \em TRUE when \p cell_residual is being assembled.
             * This ensures that only effective transport properties are calculated, not their derivatives. (improves speed)
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