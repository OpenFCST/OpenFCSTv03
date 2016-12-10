//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: saturation_transport_equation.h
//    - Description: Equation class for liquid water saturation transport due to capillary diffusion
//    - Developers: Madhur Bhaiya
//
//---------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_EQUATION_SATURATION_TRANSPORT_EQUATION_H_
#define _FCST_FUELCELLSHOP_EQUATION_SATURATION_TRANSPORT_EQUATION_H_

// FCST includes
#include <utils/fcst_units.h>
#include <materials/PureLiquid.h>

#include <equations/equation_base.h>

#include <layers/gas_diffusion_layer.h>
#include <layers/micro_porous_layer.h>
#include <layers/catalyst_layer.h>

// STD
#include <sstream>
#include <string>

namespace FuelCellShop
{
    namespace Equation
    {
        /**
         * This class deals with <b>Liquid Water Saturation Transport Equation</b>.
         *
         * This equation class solves for capillary pressure driven liquid water saturation transport inside the porous layers of the PEMFC.
         * 
         * It is to be noted that this equation can only be solved alongside <b>Thermal Transport Equation</b> and 
         * <b>Ficks Transport Equation - water</b>. Thus these two equations are necessary while considering this equation class.
         * 
         * It is solved with respect to:
         * - \f$ \mathbf{s} \f$ \p(liquid_water_saturation \p)
         * 
         * where, \f$ s = \frac{V_{\text{liquid water}}}{V_{\text{pore}}} = \frac{\text{Volume of liquid water inside the pores}}{\text{Overall pore volume}}\f$.
         * 
         * This equation can be written as:
         * 
         *  \f$ \qquad -\mathbf{\nabla} \cdot \left[ \frac{\rho_l \hat{k}_l}{\mu_l} \left( \frac{\partial p_c}{\partial s}\right) \mathbf{\nabla} s \right] = K_{e/c} ~a_{lv} ~M_{H_2O} \left( p_{H_2O,v} - p_{sat} \right) \quad \in \quad \Omega \f$
         * 
         * - \f$ \rho_l \f$ is liquid water density \f$ \left[ \frac{g}{cm^3}\right] \f$.
         * - \f$ \hat{k}_l \f$ is liquid permeability tensor \f$ \left[ cm^2 \right] \f$, which can be a function of \p liquid_water_saturation, \f$ s \f$.
         * - \f$ \mu_l \f$ is liquid water viscosity \f$ \left[ \frac{g}{cm \cdot s} \right] \f$, which can be a function of \p(temperature_of_REV \p), \f$ T \f$.
         * - \f$ p_c \f$ is capillary pressure \f$ \left[ \frac{dyne}{cm^2} \right] \f$, which can be a function of \p liquid_water_saturation, \f$ s \f$ 
         * and \p(temperature_of_REV \p), \f$ T \f$.
         * - \f$ K_{e/c} \f$ is Evaporation-condensation (phase change) rate constant \f$ \left[ \frac{mol}{Pa \cdot cm^2 \cdot s} \right] \f$.
         * - \f$ a_{lv} \f$ is Liquid-vapour interfacial surface area per unit volume \f$ \left[ \frac{cm^2}{cm^3} \right] \f$.
         * - \f$ M_{H_2O} \f$ is molar weight of water \f$ \left[ \frac{g}{mol} \right] \f$.
         * - \f$ p_{H_2O,v} \f$ is partial pressure of the water vapour \f$ \left[ Pa \right] \f$, which is a function of \p water_molar_fraction, \f$ x_{H_2O} \f$.
         * - \f$ p_{sat} \f$ is saturation pressure \f$ \left[ Pa \right] \f$, which can be a function of \p(temperature_of_REV \p), \f$ T \f$.
         * 
         * To be well-posed, these equations are equipped with the appropriate boundary conditions. All the boundary conditions can be described by
         * \p boundary_id \p(s \p) and \p boundary_type. Besides, this some boundary types require additional information, which can also be provided by the parameter file.
         * We consider following types of boundary conditions:
         * 
         * - <b>No liquid water flux / Symmetric:</b> A particular case of \p Neumann boundary condition.
         * 
         * \remarks 
         * - There is no provision to specify boundary indicators for \p No \p liquid \p water \p flux or \p Symmetric boundary conditions, as FEM
         * formulation automatically implies a particular boundary is one of these cases, by default.
         * - This class works with the following layer classes only:
         *      - FuelCellShop::Layer::GasDiffusionLayer<dim>
         *      - FuelCellShop::Layer::MicroPorousLayer<dim>
         *      - FuelCellShop::Layer::CatalystLayer<dim>
         * 
         * We solve the whole problem by linearizing the governing equation at each Newton iteration with subsequent
         * CG FEM discretization in space. The class contains the necessary class members to add the necessary contributions to cell_matrix
         * and cell_residual to the governing equations used to analyze lambda transport,
         * 
         * <h3>Usage Details:</h3>
         * 
         * @code
         * // Creating Equation object (in Application Header file)
         * FuelCellShop::Equation::SaturationTransportEquation<dim> saturation_transport;
         * 
         * // Declare parameters in application
         * saturation_transport.declare_parameters(param);
         * 
         * // Initialize in application
         * saturation_transport.initialize(param);
         * 
         * // Create a temporary vector in the application for storing couplings_map from all the equation used in the application.
         * std::vector<couplings_map> tmp;
         * ... // other equations
         * tmp.push_back( saturation_transport.get_internal_cell_couplings() );
         * 
         * // Making cell couplings using SystemManagement object created in the application
         * system_management.make_cell_couplings(tmp);
         * 
         * // cell_matrix in application
         * // Do a check against layer and it should match with the layers currently working for this equation class.
         * // for eg: CCL is FuelCellShop::Layer::HomogeneousCL<dim> object.
         * saturation_transport.assemble_cell_matrix(cell_matrices, cell_info, &CCL);
         * 
         * // cell_residual in application
         * saturation_transport.assemble_cell_residual(cell_vector, cell_info, &CCL); 
         * @endcode
         * 
         * \note This class doesn't work if thermal transport equation and fick's diffusion equation for water vapour is not solved for in the application. This 
         * class also assembles the phase change source/sink term for the water vapour transport inside the layers. Also heat released/absorbed due to 
         * condensation/evaporation is also automatically assembled for the thermal transport equation.
         * 
         * \author Madhur Bhaiya, 2013
         */
        
        template<int dim>
        class SaturationTransportEquation : public EquationBase<dim>
        {
        public:
            
            ///@name Constructors, destructor, and initalization
            //@{
            
            /**
             * Constructor.
             */
            SaturationTransportEquation(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = 
            boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >());
            
            /**
             * Destructor.
             */
            virtual ~SaturationTransportEquation();
            
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
            
            /**
             * Returns evaporation rate constant, \f$ \left[ \frac{mol}{Pa \cdot cm^2 \cdot s}\right] \f$.
             */
            double get_kevap() const
            {
                return k_e;
            }
            
            /**
             * Returns condensation rate constant, \f$ \left[ \frac{mol}{Pa \cdot cm^2 \cdot s}\right] \f$.
             */
            double get_kcond() const
            {
                return k_c;
            }

            //@}
            
        protected:

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
             * VariableInfo structure corresponding to \p "liquid_water_saturation".
             */
            VariableInfo s_liquid_water;
            
            /**
             * VariableInfo structure corresponding to \p "water_molar_fraction".
             */
            VariableInfo x_water;
            
            /**
             * VariableInfo structure corresponding to \p "temperature_of_REV".
             */
            VariableInfo t_rev;
            
            /**
             * Block index for \p "water_molar_fraction" corresponding to <b>"Ficks Transport Equation - water"</b>.
             */
            unsigned int fickswater_blockindex_xwater;
            
            /**
             * Block index for \p "liquid_water_saturation" corresponding to <b>"Ficks Transport Equation - water"</b>.
             */
            unsigned int fickswater_blockindex_sliquid;
            
            /**
             * Block index for \p "temperature_of_REV" corresponding to <b>"Ficks Transport Equation - water"</b>.
             */
            unsigned int fickswater_blockindex_trev;
            
            /**
             * Block index for \p "water_molar_fraction" corresponding to <b>"Thermal Transport Equation"</b>.
             */
            unsigned int thermal_blockindex_xwater;
            
            /**
             * Block index for \p "liquid_water_saturation" corresponding to <b>"Thermal Transport Equation"</b>.
             */
            unsigned int thermal_blockindex_sliquid;
            
            /**
             * Block index for \p "temperature_of_REV" corresponding to <b>"Thermal Transport Equation"</b>.
             */
            unsigned int thermal_blockindex_trev;
            
            /**
             * Molar weight of water in grams/mole.
             */
            double M_water;
            
            /**
             * Evaporation rate constant, \f$ \left[ \frac{mol}{Pa \cdot cm^2 \cdot s}\right] \f$.
             */
            double k_e;
            /**
             * Condensation rate constant, \f$ \left[ \frac{mol}{Pa \cdot cm^2 \cdot s}\right] \f$.
             */
            double k_c;
            
            /**
             * Density of liquid water, \f$ \rho_l ~\left[ \frac{g}{cm^3} \right] \f$.
             */
            double rho_l;
            
            //@}
            
            ///@name Local CG FEM based assemblers - variable data (cell)
            //@{
            
            /**
             * Total pressure \f$ \left[ Pa \right] \f$ in the cell.
             */
            double p_cell;
            
            /**
             * Interfacial surface area per unit volume, \f$ a_{lv} ~\left[ \frac{cm^2}{cm^3} \right] \f$, at all quadrature points in the cell.
             */
            std::vector<double> area_lv_cell;
            
            /**
             * \f$ \frac{ \partial a_{lv}}{\partial s} \f$, at all quadrature points in the cell.
             */
            std::vector<double> darea_lv_ds_cell;
            
            /**
             * \f$ \frac{\rho_l \hat{k}_l}{\mu_l} \left(\frac{\partial p_c}{\partial s} \right) \f$, at all quadrature points in the cell.
             */
            std::vector< Tensor<2,dim> > rhok_mu_dpcds_cell;
            
            /**
             * \f$ \frac{\partial}{\partial s} \left[ \frac{\rho_l \hat{k}_l}{\mu_l} \left(\frac{\partial p_c}{\partial s} \right) \right] \f$, 
             * at all quadrature points in the cell.
             */
            std::vector< Tensor<2,dim> > ds_rhok_mu_dpcds_cell;
            
            /**
             * \f$ \frac{\partial}{\partial T} \left[ \frac{\rho_l \hat{k}_l}{\mu_l} \left(\frac{\partial p_c}{\partial s} \right) \right] \f$, 
             * at all quadrature points in the cell.
             */
            std::vector< Tensor<2,dim> > dT_rhok_mu_dpcds_cell;

            /**
             * \f$ \mathbf{s} \f$ shape functions.
             * 
             * \p phi_s_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{s} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_s_cell;
            
            /**
             * \f$ \mathbf{x_{H_2O}} \f$ shape functions.
             * 
             * \p phi_xwater_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{x_{H_2O}} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_xwater_cell;
            
            /**
             * \f$ \mathbf{T} \f$ shape functions.
             * 
             * \p phi_T_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{T} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_T_cell;
                   
            /**
             * \f$ \mathbf{s} \f$ shape function gradients.
             * 
             * \p grad_phi_s_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{s} \f$ shape function gradient
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< Tensor<1,dim> > > grad_phi_s_cell;
            
            /**
             * \f$ \mathbf{T} \f$ shape function gradients.
             * 
             * \p grad_phi_T_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{T} \f$ shape function gradient
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< Tensor<1,dim> > > grad_phi_T_cell;
            
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