// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: compressible_multi_component_KG_equations_coupled.h
// - Description: This class describes steady-state compressible and isothermal Kerkhof-Geboers
//   fluid transport equations for a single-phase multi-component case coupled with fuel cell
//   physics
// - Developers: Valentin N. Zingan, Chad Balen and Marc Secanell, University of Alberta
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_EQUATION_COMPRESSIBLE_MULTI_COMPONENT_KG_EQUATIONS_COUPLED_H_
#define _FCST_FUELCELLSHOP_EQUATION_COMPRESSIBLE_MULTI_COMPONENT_KG_EQUATIONS_COUPLED_H_

#include <equations/equation_base.h>
#include <layers/channel.h>
#include <layers/gas_diffusion_layer.h>
#include <layers/micro_porous_layer.h>
#include <layers/experimental_porous_layer.h>
#include <layers/catalyst_layer.h>
#include <utils/fcst_units.h>

namespace FuelCellShop
{
    namespace Equation
    {
        
        /**
         * \brief This class implements the multi-component mass transport equations proposed by Kerkhof-Geboers for fluid transport.
         *
         * These equations are
         *
         * - steady-state,
         * - compressible,
         * - isothermal,
         * - single-phase,
         * - multi-component,
         *
         * They take the following form:
         * - point-wise form in the channels of fuel cells,
         * - appropriately averaged over a so-called Representative Elementary Volume (REV) in the porous media regions,
         * - solved with respect to point-wise
         *   \f$ \{ \left( \rho_i, \mathbf{u}_i \right) \}_{i=1}^N \f$
         *   in channels and with respect to REV-wise
         *   \f$ \{ \left( < \rho_i >, < \mathbf{u}_i >^{\text{f}} \right) \}_{i=1}^N \f$
         *   in porous media,
         * - coupled with fuel cell physics in
         *   \p ReactionSourceTerms and \p SorptionSourceTerms classes.
         *
         * 
         * For any species \f$ i = 1,N \f$ \f$ \quad \f$ these equations can be written as
         *
         * \f$ \nabla \cdot \mathbf{F}_{ mass_i } = 0 \quad \text{in} \quad \Omega \f$ \f$ \quad - \quad \f$ \f$ i \f$-th mass conservation equation
         *
         * \f$ \nabla \cdot \hat{ \mathbf{F} }_{ mom_i } = \nabla \cdot \left( -p_i \hat{\boldsymbol{I}} + \hat{ \boldsymbol\sigma }_i \right) + \left( \mathbf{F}_i + \mathbf{D}_i + \rho_i \mathbf{g} \right)
         * \quad \text{in} \quad \Omega \f$ \f$ \quad - \quad \f$ \f$ i \f$-th momentum conservation equation
         *
         * where
         *
         * \f$ \mathbf{F}_{ mass_i } = \rho_i \mathbf{u}_i \f$ \f$ \quad - \quad \f$ \f$ i \f$-th mass flux
         *
         * \f$ \hat{ \mathbf{F} }_{ mom_i } = \rho_i \mathbf{u}_i \otimes \mathbf{u}_i \f$ \f$ \quad - \quad \f$ \f$ i \f$-th momentum flux
         *
         * \f$ p_i = \rho_i \frac{R}{M_i} T_{\text{mixture}}^{\text{const}} \f$ \f$ \quad - \quad \f$ \f$ i \f$-th partial pressure
         *
         * \f$ \hat{ \boldsymbol\sigma }_i = 2 \mu_i \nabla_s \epsilon \mathbf{u}_i + \lambda_i \left( \nabla \cdot \epsilon \mathbf{u}_i \right) \hat{ \mathbf{I} } \f$ \f$ \quad - \quad \f$ \f$ i \f$-th shear stress
         *
         * \f$
         *   \mathbf{F}_i =
         *   \begin{cases}
         *   \mathbf{0} \quad \text{in} \quad \Omega_c \\
         * - \mu_i \hat{ \mathbf{K} }^{-1} \epsilon \mathbf{u}_i - \hat{\boldsymbol\beta} \rho_i \big| \epsilon \mathbf{u}_i \big| \epsilon \mathbf{u}_i \quad \text{in} \quad \Omega_p
         *   \end{cases}
         * \f$ \f$ \quad - \quad \f$ \f$ i \f$-th drag force
         *
         * \f$ \mathbf{D}_i = \displaystyle \sum_{j=1}^N p_i p_j \hat{ \mathbf{D} }_{ij}^{-1} \left( \epsilon \mathbf{u}_j - \epsilon \mathbf{u}_i \right) \f$ \f$ \quad - \quad \f$ \f$ i \f$-th diffusion force
         *
         * \f$
         * \hat{ \mathbf{D} }_{ij}^{-1} =
         * \begin{cases}
         * \displaystyle \left( \sum_{l=1}^N p_l \cdot \mathscr{D}_{ij} \right)^{-1} \hat{ \mathbf{I} }          \quad \text{in} \quad \Omega_c \\
         * \displaystyle \left( \epsilon \sum_{l=1}^N p_l \cdot \mathscr{D}_{ij} \right)^{-1} \hat{ \mathbf{T} } \quad \text{in} \quad \Omega_p
         * \end{cases}
         * \f$ \f$ \quad - \quad \f$ the inverse of \f$ i \f$-th, \f$ j \f$-th Maxwell–Stefan isobaric diffusion coefficient
         * 
         * 
         * To be well-posed, these equations are equipped with appropriate boundary conditions.
         * Below we consider several types of boundaries and several types of boundary conditions assigned to those boundaries (here \f$ k \f$ denotes some portion of a boundary of a certain type).
         *
         * - Impermeable walls \f$ \left( \Gamma^{\text{walls},k}_i \right) \f$ \f$ \quad - \quad \f$ no-slip boundary condition
         *
         * \f$
         * \quad \mathbf{u}_i \Bigg|_{\Gamma^{\text{walls},k}_i} = \mathbf{0} \quad \text{OR}
         * \f$
         *
         * - Impermeable walls \f$ \left( \Gamma^{\text{walls},k}_i \right) \f$ \f$ \quad - \quad \f$ Navier slip boundary condition
         *
         * \f$
         * \mathbf{u}_i \cdot \mathbf{n} \Bigg|_{\Gamma^{\text{walls},k}_i} = 0 \quad \text{and} \quad
         * \left(1-\theta_i \right) \mathbf{u}_i \cdot \boldsymbol\tau_{\alpha} + \theta_i \left( -p_i \hat{\boldsymbol{I}} + \hat{ \boldsymbol\sigma }_i \right) \mathbf{n} \cdot \boldsymbol\tau_{\alpha}
         * \Bigg|_{\Gamma^{\text{walls},k}_i} = 0 \quad \text{with} \quad \alpha = 1, d-1 \quad \text{and} \quad 0 < \theta_i < 1 \quad \text{OR}
         * \f$
         *
         * - Impermeable walls \f$ \left( \Gamma^{\text{walls},k}_i \right) \f$ \f$ \quad - \quad \f$ Maxwell slip boundary condition
         *
         * \f$
         * \mathbf{u}_i \cdot \mathbf{n} \Bigg|_{\Gamma^{\text{walls},k}_i} = 0 \quad \text{and} \quad
         * \mathbf{u}_i \cdot \boldsymbol\tau_{\alpha} + \frac{2-\sigma_i}{\sigma_i \rho_i} \sqrt{\frac{\pi M_i}{2RT_{\text{mixture}}^{\text{const}}}}
         * \left( -p_i \hat{\boldsymbol{I}} + \hat{ \boldsymbol\sigma }_i \right) \mathbf{n} \cdot \boldsymbol\tau_{\alpha}
         * \Bigg|_{\Gamma^{\text{walls},k}_i} = 0 \quad \text{with} \quad \alpha = 1, d-1 \quad \text{OR}
         * \f$
         *
         * - Impermeable walls \f$ \left( \Gamma^{\text{walls},k}_i \right) \f$ \f$ \quad - \quad \f$ perfect slip boundary condition
         *
         * \f$
         * \mathbf{u}_i \cdot \mathbf{n} \Bigg|_{\Gamma^{\text{walls},k}_i} = 0 \quad \text{and} \quad
         * \left( -p_i \hat{\boldsymbol{I}} + \hat{ \boldsymbol\sigma }_i \right) \mathbf{n} \cdot \boldsymbol\tau_{\alpha} \Bigg|_{\Gamma^{\text{walls},k}_i} = 0 \quad
         * \text{with} \quad \alpha = 1, d-1
         * \f$
         *
         * - Symmetry line or plane \f$ \left( \Gamma^{\text{sym}} \right) \f$ \f$ \quad - \quad \f$ perfect slip boundary condition
         *
         * \f$
         * \forall i = 1,N: \quad \mathbf{u}_i \cdot \mathbf{n} \Bigg|_{\Gamma^{\text{sym}}} = 0 \quad \text{and} \quad
         * \left( -p_i \hat{\boldsymbol{I}} + \hat{ \boldsymbol\sigma }_i \right) \mathbf{n} \cdot \boldsymbol\tau_{\alpha} \Bigg|_{\Gamma^{\text{sym}}} = 0 \quad
         * \text{with} \quad \alpha = 1, d-1
         * \f$
         *
         * - Inlet-Outlet \f$ \left( \Gamma^{\text{in-out},k}_i \right) \f$ \f$ \quad - \quad \f$ Dirichlet boundary condition
         *
         * \f$
         * \rho_i \Bigg|_{\Gamma^{\text{in-out},k}_i} = \rho_i^{\text{prescribed},k} \quad \text{or} \quad \mathbf{u}_i \Bigg|_{\Gamma^{\text{in-out},k}_i} = \mathbf{u}_i^{\text{prescribed},k}
         * \quad \text{or} \quad \text{both} \quad \text{OR/AND}
         * \f$
         *
         * - Inlet-Outlet \f$ \left( \Gamma^{\text{in-out},k}_i \right) \f$ \f$ \quad - \quad \f$ normal stress free boundary condition
         *
         * \f$
         * \left( -p_i \hat{\boldsymbol{I}} + \hat{ \boldsymbol\sigma }_i \right) \mathbf{n} \Bigg|_{\Gamma^{\text{in-out},k}_i} = \mathbf{0} \quad \text{OR}
         * \f$
         *
         * - Inlet-Outlet \f$ \left( \Gamma^{\text{in-out},k}_i \right) \f$ \f$ \quad - \quad \f$ normal shear stress free boundary condition
         *
         * \f$
         * \hat{ \boldsymbol\sigma }_i \mathbf{n} \Bigg|_{\Gamma^{\text{in-out},k}_i} = \mathbf{0}
         * \f$
         * 
         * 
         * All these boundaries and boundary conditions can be specified in the parameters file as follows:
         *
         * @code
         * Impermeable walls [species] = [boundary_id]: [type]
         * @endcode
         *
         * where
         *
         * - [species]     = 1,2,3,4,5
         * - [boundary_id] = 0,1,2,3,...
         * - [type]        = no-slip | Navier slip | Maxwell slip | perfect slip
         *
         * @code
         * Symmetry line or plane [species] = [boundary_id]: [type]
         * @endcode
         *
         * where
         *
         * - [species]     = 1,2,3,4,5
         * - [boundary_id] = 0,1,2,3,...
         * - [type]        = perfect slip
         *
         * @code
         * Inlet-Outlet [species] = [boundary_id]: [type]
         * @endcode
         *
         * where
         *
         * - [species]     = 1,2,3,4,5
         * - [boundary_id] = 0,1,2,3,...
         * - [type]        = Dirichlet density | Dirichlet density and normal stress free | Dirichlet density and normal shear stress free | Dirichlet velocity | Dirichlet density and velocity | normal stress free |
         *                   normal shear stress free
         * 
         * 
         * These equations utilize
         *
         * - \p Channel,
         * - \p GasDiffusionLayer,
         * - \p MicroPorousLayer,
         * - \p CatalystLayer
         *
         * as layer classes and
         *
         * - \p GasMixture,
         * - \p PureGas
         *
         * as material classes and have the following output:
         *
         * - densities,  \f$ \rho_i       \f$ \f$ \quad \left[ \frac{\text{g}}{\text{cm}^3} \right] \f$
         * - velocities, \f$ \mathbf{u}_i \f$ \f$ \quad \left[ \frac{\text{cm}}{\text{sec}} \right] \f$
         *
         * We usually do not solve these equations in the membranes of fuel cells
         * which is fully justified under the normal conditions of operation.
         * 
         * 
         * The whole problem is solved by linearizing the governing equations
         * at each Newton iteration with subsequent CG FEM discretization in space.
         *
         * The implementation is performed by means of using only one weak formulation
         * for all \f$ N \left( d + 1 \right) \f$ scalar governing equations.
         * This approach guarantees the dim-independent programming of the equations at hand.
         * 
         * \note outer_product() is a deal.ii function
         *
         * \author Valentin N. Zingan, Chad Balen and Marc Secanell, 2013-15
         */
        
        template<int dim>
        class CompressibleMultiComponentKGEquationsCoupled : public EquationBase<dim>
        {
        public:
            
            ///@name Constructors, destructor, and initialization
            //@{
            
            /**
             * Constructor.
             */
            CompressibleMultiComponentKGEquationsCoupled(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = 
            boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >());
            
            /**
             * Destructor.
             */
            virtual ~CompressibleMultiComponentKGEquationsCoupled();
            
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
            virtual void assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                              const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                              FuelCellShop::Layer::BaseLayer<dim>* const                               layer);  
            /**
             * Assemble local cell residual.
             */
            virtual void assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                FuelCellShop::Layer::BaseLayer<dim>* const                               layer);
            
            /**
             * Assemble local boundary matrix.
             */
            virtual void assemble_bdry_matrix(FuelCell::ApplicationCore::MatrixVector&                                 bdry_matrices,
                                              const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                              FuelCellShop::Layer::BaseLayer<dim>* const                               layer);
            
            /**
             * Assemble local boundary residual.
             */
            virtual void assemble_bdry_residual(FuelCell::ApplicationCore::FEVector&                                     bdry_residual,
                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                FuelCellShop::Layer::BaseLayer<dim>* const                               layer);
            
            //@}
            
            ///@name Accessors and info
            //@{
            /**
             * This function returns
             * \p inlet_outlet_velocity_max.
             */
            const std::vector<double>& get_inlet_outlet_velocity_max() const
            {
                return inlet_outlet_velocity_max;
            }
            
            /**
             * This function returns
             * \p inlet_outlet_velocity_mixture_max.
             */
            const double get_inlet_outlet_velocity_mixture_max() const
            {
                return inlet_outlet_velocity_mixture_max;
            }
            
            /**
             * This function tells the application if the user wishes to 
             * use a parabolic profile for a boundary condition. If false
             * then user must specify equation to use with parameter
             * \p inlet-outlet \p velocity \p of \p the \p mixture \p equation
             */
            const bool get_use_parabolic_profile() const
            {
                return use_parabolic_profile;
            }
            
            /**
             * This function tells the application the height [cm] of the base of the
             * channel with respect to origin.
             */
            const double get_inlet_outlet_velocity_channel_base() const
            {
                return inlet_outlet_velocity_channel_base;
            }
            
            /**
             * This function tells the application the height [cm] of the roof of the
             * channel with respect to origin.
             */
            const double get_inlet_outlet_velocity_channel_roof() const
            {
                return inlet_outlet_velocity_channel_roof;
            }
            
            /**
             * This function tells the application the boundary id to apply velocity profile 
             * equation to. TODO: in the future it might be nice if OpenFCST could detect this based
             * on boundary conditions
             */
            const int get_inlet_outlet_boundary_ID() const
            {
                return inlet_outlet_boundary_ID;
            }
            
            /**
            * Outputs which component to apply variable boundary equation to (i.e. Pressure, VelocityX, VelocityY,
            * VelocityZ).
            */
            const std::string get_press_vel_comp_apply_to() const
            {
                return press_vel_comp_apply_to;
            }
            
            /**
             * If \p get_use_parabolic_profile() is false then user can specify an equation
             * to use for velocity profile at boundary. Ex. \p inlet-outlet \p velocity \p of 
             * \p the \p mixture \p equation = -625*y*y + 66.25*y - 0.755625
             */
            const std::string get_inlet_outlet_velocity_mixture_equation() const
            {
                return inlet_outlet_velocity_mixture_equation;
            }
            
            /**
             * \p get_gas_species_map() returns a map relating the different KG species 
             * equations to the desired gas species being simulated;
             * e.g 1:oxygen, 2:water vapour, 3:nitrogen means that the first 3 KG equations
             * (Mass Cons., Mom. Cons. X, Mom. Cons. Y) are for species oxygen, the next 
             * three are for water vapour, and the final three are for nitrogen.
             */
            const std::map<unsigned int, std::string> get_gas_species_map() const
            {
                return gas_species_map;
            }
            
            /**
             * If returns the number of species to use in the parameter file
             */
            const unsigned int get_num_of_species() const
            {
                return n_species;
            }
            
            /**
             * This function prints out
             * the equations info.
             */
            virtual void print_equation_info() const;
            
            //@}
            
        protected:
            
            ///@name Local CG FEM based assemblers - make_ functions
            //@{
            
            /**
             * This function computes
             * Computational constants
             * and
             * Local CG FEM based assemblers - constant data (generic).
             *
             * The template parameter INFO is either
             * typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo
             * or
             * typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo.
             */
            template<typename INFO>
            void make_assemblers_generic_constant_data(const INFO&                                InFo,
                                                       FuelCellShop::Layer::BaseLayer<dim>* const layer);
            
            /**
             * This function computes
             * Local CG FEM based assemblers - constant data (cell)
             * and
             * allocates the memory for
             *
             * - Local CG FEM based assemblers - variable data (previous Newton iteration - cell),
             * - Local CG FEM based assemblers - variable data (current Newton iteration - cell),
             * - Local CG FEM based assemblers - variable data (other - cell).
             */
            virtual void make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info);
            
            /**
             * This function computes
             * Local CG FEM based assemblers - constant data (boundary)
             * and
             * allocates the memory for
             *
             * - Local CG FEM based assemblers - variable data (previous Newton iteration - boundary),
             * - Local CG FEM based assemblers - variable data (current Newton iteration - boundary),
             * - Local CG FEM based assemblers - variable data (other - boundary).
             */
            virtual void make_assemblers_bdry_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info);
            
            /**
             * This function computes
             *
             * - Local CG FEM based assemblers - variable data (previous Newton iteration - cell),
             * - Local CG FEM based assemblers - variable data (current Newton iteration - cell),
             * - Local CG FEM based assemblers - variable data (other - cell).
             */
            virtual void make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                            FuelCellShop::Layer::BaseLayer<dim>* const                               layer);
            
            /**
             * This function computes
             *
             * - Local CG FEM based assemblers - variable data (previous Newton iteration - boundary),
             * - Local CG FEM based assemblers - variable data (current Newton iteration - boundary),
             * - Local CG FEM based assemblers - variable data (other - boundary).
             */
            virtual void make_assemblers_bdry_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                            FuelCellShop::Layer::BaseLayer<dim>* const                               layer);
            
            //@}
            
            ///@name Other - make_ functions
            //@{
            
            /**
             * This function fills out
             * \p internal_cell_couplings.
             */
            virtual void make_internal_cell_couplings();
            
            /**
             * This function fills out
             * \p matrix_block_indices.
             */
            virtual void make_matrix_block_indices();
            
            /**
             * This function fills out
             * \p residual_indices.
             */
            virtual void make_residual_indices();
            
            //@}
            
            //////////
            // DATA //
            //////////
            
            ///@name Boolean constants and form of the drag force in porous media
            //@{
            
            /**
             * This object indicates
             * whether the inertia term \f$ i \f$
             * is ON or OFF in channels.
             */
            std::vector<bool> inertia_in_channels;
            
            /**
             * This object indicates
             * whether the shear stress term \f$ i \f$
             * is ON or OFF in channels.
             */
            std::vector<bool> shear_stress_in_channels;
            
            /**
             * This object indicates
             * whether the gravity term \f$ i \f$
             * is ON or OFF in channels.
             */
            std::vector<bool> gravity_in_channels;
            
            /**
             * This object indicates
             * whether the inertia term \f$ i \f$
             * is ON or OFF in porous media.
             */
            std::vector<bool> inertia_in_porous_media;
            
            /**
             * This object indicates
             * whether the shear stress term \f$ i \f$
             * is ON or OFF in porous media.
             */
            std::vector<bool> shear_stress_in_porous_media;
            
            /**
             * This object indicates
             * whether the gravity term \f$ i \f$
             * is ON or OFF in porous media.
             */
            std::vector<bool> gravity_in_porous_media;
            
            /**
             * This object indicates
             * which form of the drag term \f$ i \f$
             * is supposed to be chosen.
             *
             * There are four options currently implemented in FCST.
             * These options are
             *
             * - none,
             * - Darcy,
             * - Forchheimer,
             * - Forchheimer modified.
             */
            std::vector<std::string> drag_in_porous_media;
            
            /**
             * Sometimes we implement
             * different types of slip
             * boundary conditions on
             * the curved impermeable walls
             * of the domain or the walls
             * that are not perfectly
             * aligned with the
             * coordinate axes.
             *
             * In this case, it might be
             * quite problematic (but still possible)
             * to strictly enforce the
             * normal component
             * of velocity \f$ i \f$
             * be equal
             * to 0.
             *
             * The alternative approach
             * is doing that by introducing
             * a penalization term \f$ \displaystyle \int_{ \Gamma^{\text{walls},k}_i } \frac{1}{\eta_i} \left( \boldsymbol\omega_i \cdot \mathbf{n} \right) \left( \mathbf{u}_i \cdot \mathbf{n} \right) dS \f$
             * with \f$ \eta_i \sim 10^{-10} - 10^{-12} \f$
             * into the weak formulation
             * of the equations
             * at hand.
             *
             * This object indicates
             * whether the penalty
             * method is used.
             */
            std::vector<bool> normal_velocity_is_suppressed_weakly;
            
            //@}
            
            ///@name Computational constants
            //@{
            
            /**
             * Number of species, \f$ N \f$.
             */
            unsigned int n_species;
            
            bool coupled_with_fuel_cell_physics;
            
            /**
             * Map of the gas species to their respective equation species number. Possible names for gasses include: oxygen, water, nitrogen. \f$ N \f$.
             */
            std::map<unsigned int, std::string> gas_species_map;
            
            /**
             * Universal gas constant, \f$ R = 8.314462175 \cdot 10^7 \quad \left[ \frac{\text{g } \text{cm}^2}{\text{mol K } \text{sec}^2} \right] \f$.
             */
            double R_universal;
            
            /**
             * Temperature of species mixture, \f$ T_{\text{mixture}}^{\text{const}} \quad \left[ \text{K} \right] \f$.
             */
            double T_mixture;
            
            /**
             * Pressure at inlet, \f$ P_{\text{in}}^{\text{const}} \quad \left[ \text{Pa} \right] \f$.
             */
            double P_in;
            
            /**
             * Gravitational acceleration, \f$ \mathbf{g} = \{ g_{\alpha} \}_{\alpha = 1}^d \quad \text{such that} \quad \forall \alpha \neq d : \quad g_{\alpha} = 0 \quad \left[ \frac{\text{cm}}{\text{sec}^2} \right] \quad
             * \text{and} \quad g_d = -9.81 \cdot 10^2 \quad \left[ \frac{\text{cm}}{\text{sec}^2} \right] \f$.
             */
            Tensor<1,dim> gravity_acceleration;
            
            /**
             * Unit tensor, \f$ \hat{ \mathbf{I} } = \{ \delta_{\alpha \beta} \}_{\alpha,\beta = 1}^d \f$.
             */
            SymmetricTensor<2,dim> unit;
            
            /**
             * Molar mass of pure gas, \f$ M_i \quad \left[ \frac{\text{g}}{\text{mol}} \right] \f$.
             */
            std::vector<double> molar_mass;
            
            /**
             * Dynamic viscosity of pure gas, \f$ \mu^o_i \quad \left[ \frac{\text{g}}{\text{cm sec}} \right] \f$.
             */
            std::vector<double> dynamic_viscosity;
            
            /**
             * Bulk viscosity of pure gas, \f$ \lambda_i \quad \left[ \frac{\text{g}}{\text{cm sec}} \right] \f$.
             */
            std::vector<double> bulk_viscosity;
            
            /**
             * Collision diameter of pure gas, \f$ \lambda_i \quad \left[ \frac{\text{g}}{\text{cm sec}} \right] \f$.
             */
            std::vector<double> collision_diameter;
            
            /**
             * Navier slip coefficient of pure gas, \f$ \theta_i \f$.
             */
            std::vector<double> theta;
            
            /**
             * Normal velocity suppression coefficient of pure gas, \f$ \eta_i \f$.
             */
            std::vector<double> eta;
            
            /**
             * Each entry of this structure defines
             * a Maxwell-Stefan isobaric diffusion coefficient of gas \f$ i \f$ in gas \f$ j \f$,
             * \f$ \displaystyle \sum_{l=1}^N p_l \cdot \mathscr{D}_{ij} \quad \left[ \frac{\text{g cm}}{\text{sec}^3} \right] \f$.
             */
            Table< 2, double > maxwell_stefan_isobaric_diffusion_coefficient;
            
            /**
             * Maximum inlet-outlet velocity of pure gas, \f$ \mathbf{u}_{i, \text{in,out}}^{\text{max}} \f$ \f$ \quad \left[ \frac{\text{cm}}{\text{sec}} \right] \f$.
             */
            std::vector<double> inlet_outlet_velocity_max;
            
            /**
             * Maximum inlet-outlet velocity of the mixture, \f$ \mathbf{u}_{\text{mix,in,out}}^{\text{max}} \f$ \f$ \quad \left[ \frac{\text{cm}}{\text{sec}} \right] \f$.
             */
            double      inlet_outlet_velocity_mixture_max;
            bool        use_parabolic_profile;
            double      inlet_outlet_velocity_channel_base;
            double      inlet_outlet_velocity_channel_roof;
            int         inlet_outlet_boundary_ID;
            std::string press_vel_comp_apply_to;
            std::string inlet_outlet_velocity_mixture_equation;
            
            /**
             * These constants \f$ \sigma_i \f$ are used in the Maxwell slip boundary condition.
             */
            std::vector<double> maxwell_constant;
            
            //@}
            
            ///@name Local CG FEM based assemblers - constant data (generic)
            //@{
            
            /**
             * Density extractors.
             */
            std::vector< FEValuesExtractors::Scalar > density;
            
            /**
             * Velocity extractors.
             */
            std::vector< FEValuesExtractors::Vector > velocity;
            
            //@}
            
            ///@name Local CG FEM based assemblers - constant data (cell)
            //@{
            
            /**
             * Implementation is
             * in the base class.
             */
            
            //@}
            
            ///@name Local CG FEM based assemblers - constant data (boundary)
            //@{
            
            /**
             * Implementation is
             * in the base class.
             */
            
            //@}
            
            ///@name Local CG FEM based assemblers - variable data (active mesh iterators)
            //@{
            
            /**
             * Implementation is
             * in the base class.
             */
            
            //@}
            
            ///@name Local CG FEM based assemblers - variable data (previous Newton iteration - cell)
            //@{
            
            /**
             * Density of each species
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector< std::vector<double> > density_cell_old;
            
            /**
             * Density gradient of each species
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector< std::vector< Tensor<1,dim> > > grad_density_cell_old;
            
            /**
             * Velocity of each species
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector< std::vector< Tensor<1,dim> > > velocity_cell_old;
            
            /**
             * Velocity divergence of each species
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector< std::vector<double> > div_velocity_cell_old;
            
            /**
             * Velocity symmetric gradient of each species
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector< std::vector< SymmetricTensor<2,dim> > > grads_velocity_cell_old;
            
            /**
             * Mass flux of each species
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector< std::vector< Tensor<1,dim> > > mass_flux_cell_old;
            
            /**
             * Momentum flux of each species
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector< std::vector< SymmetricTensor<2,dim> > > momentum_flux_cell_old;
            
            /**
             * Partial pressure of each species
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector< std::vector<double> > pressure_cell_old;
            
            /**
             * Shear stress of each species
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector< std::vector< SymmetricTensor<2,dim> > > shear_stress_cell_old;
            
            /**
             * Partial viscosity of mixture of each species, \f$ \eta_i \quad \left[ \frac{\text{g}}{\text{cm sec}} \right] \f$,
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             * NOTE: first index is the quadrature point and second index is the species.
             */
            std::vector< std::vector<double> > partial_viscosity_old;
            
            /**
             * Bulk viscosity of mixture of each species, \f$ \lambda_i \quad \left[ \frac{\text{g}}{\text{cm sec}} \right] \f$.
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             * NOTE: first index is the quadrature point and second index is the species.
             */
            std::vector< std::vector<double> > bulk_viscosity_old;
            
            /**
             * paramMatrix is used for calculating partial viscosity of mixture
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             * NOTE: first index is the quadrature point and second and third indices are the species.
             */
            std::vector< std::vector< std::vector<double> > > paramMatrix_old;
            
            /**
             * PInv_old is used for calculating the inverse of P when using the OmegaKG model for calculating partial viscosity
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector< FullMatrix<double> > PInv_old;
            
            /**
             * Drag force of each species
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector< std::vector< Tensor<1,dim> > > drag_cell_old;
            
            /**
             * Diffusion force of each species
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector< std::vector< Tensor<1,dim> > > diffusion_cell_old;
            
            /**
             * Gravity force of each species
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector< std::vector< Tensor<1,dim> > > gravity_cell_old;
            
            //@}
            
            ///@name Local CG FEM based assemblers - variable data (current Newton iteration - cell)
            //@{
            
            /**
             * Density shape functions.
             *
             * \p phi_density_cell \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th density shape function
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector<double> > > phi_density_cell;
            
            /**
             * Density shape function gradients.
             *
             * \p grad_phi_density_cell \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th density shape function gradient
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector< Tensor<1,dim> > > > grad_phi_density_cell;
            
            /**
             * Velocity shape functions.
             *
             * \p phi_velocity_cell \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th velocity shape function
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector< Tensor<1,dim> > > > phi_velocity_cell;
            
            /**
             * Velocity shape function divergences.
             *
             * \p div_phi_velocity_cell \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th velocity shape function divergence
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector<double> > > div_phi_velocity_cell;
            
            /**
             * Velocity shape function symmetric gradients.
             *
             * \p grads_phi_velocity_cell \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th velocity shape function symmetric gradient
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector< SymmetricTensor<2,dim> > > > grads_phi_velocity_cell;
            
            /**
             * Mass flux perturbations.
             *
             * \p delta_mass_flux_cell \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th mass flux perturbation
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector< Tensor<1,dim> > > > delta_mass_flux_cell;
            
            /**
             * Momentum flux perturbations.
             *
             * \p delta_momentum_flux_cell \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th momentum flux perturbation
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector< SymmetricTensor<2,dim> > > > delta_momentum_flux_cell;
            
            /**
             * Partial pressure perturbations.
             *
             * \p delta_pressure_cell \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th partial pressure perturbation
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector<double> > > delta_pressure_cell;
            
            /**
             * Partial dynamic viscosity perturbations.
             *
             * \p delta_partial_viscosity_cell \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th partial dynamic velocity perturbation
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector<double> > > delta_partial_viscosity_cell;
            
            /**
             * Partial bulk viscosity perturbations.
             *
             * \p delta_bulk_viscosity_cell \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th partial bulk velocity perturbation
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector<double> > > delta_bulk_viscosity_cell;
            
            /**
             * Shear stress perturbations.
             *
             * \p delta_shear_stress_cell \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th shear stress perturbation
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector< SymmetricTensor<2,dim> > > > delta_shear_stress_cell;
            
            /**
             * Drag force perturbations.
             *
             * \p delta_drag_cell \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th drag force perturbation
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector< Tensor<1,dim> > > > delta_drag_cell;
            
            /**
             * Diffusion force perturbations.
             *
             * \p delta_diffusion_cell \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th diffusion force perturbation
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector< Tensor<1,dim> > > > delta_diffusion_cell;
            
            /**
             * Gravity force perturbations.
             *
             * \p delta_gravity_cell \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th gravity force perturbation
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector< Tensor<1,dim> > > > delta_gravity_cell;
            
            //@}
            
            ///@name Local CG FEM based assemblers - variable data (other - cell)
            //@{
            
            /**
             * Implementation is
             * in the base class.
             */
            
            //@}
            
            ///@name Local CG FEM based assemblers - variable data (previous Newton iteration - boundary)
            //@{
            
            /**
             * Density of each species in the quadrature points of a boundary
             * at a previous Newton iteration.
             */
            std::vector< std::vector<double> > density_bdry_old;
            
            /**
             * Density gradient of each species in the quadrature points of a boundary
             * at a previous Newton iteration.
             */
            std::vector< std::vector< Tensor<1,dim> > > grad_density_bdry_old;
            
            /**
             * Velocity of each species in the quadrature points of a boundary
             * at a previous Newton iteration.
             */
            std::vector< std::vector< Tensor<1,dim> > > velocity_bdry_old;
            
            /**
             * Velocity divergence of each species in the quadrature points of a boundary
             * at a previous Newton iteration.
             */
            std::vector< std::vector<double> > div_velocity_bdry_old;
            
            /**
             * Velocity symmetric gradient of each species in the quadrature points of a boundary
             * at a previous Newton iteration.
             */
            std::vector< std::vector< SymmetricTensor<2,dim> > > grads_velocity_bdry_old;
            
            /**
             * Mass flux of each species in the quadrature points of a boundary
             * at a previous Newton iteration.
             */
            std::vector< std::vector< Tensor<1,dim> > > mass_flux_bdry_old;
            
            /**
             * Momentum flux of each species in the quadrature points of a boundary
             * at a previous Newton iteration.
             */
            std::vector< std::vector< SymmetricTensor<2,dim> > > momentum_flux_bdry_old;
            
            /**
             * Partial pressure of each species in the quadrature points of a boundary
             * at a previous Newton iteration.
             */
            std::vector< std::vector<double> > pressure_bdry_old;
            
            /**
             * Shear stress of each species in the quadrature points of a boundary
             * at a previous Newton iteration.
             */
            std::vector< std::vector< SymmetricTensor<2,dim> > > shear_stress_bdry_old;
            
            /**
             * Partial dynamic viscosity, \f$ \eta_i \quad \left[ \frac{\text{g}}{\text{cm sec}} \right] \f$,
             * for each species in the quadrature points of a boundary
             * at a previous Newton iteration.
             * NOTE: first index is the quadrature point and second index is the species.
             */
            std::vector< std::vector<double> > partial_viscosity_bdry_old;
            
            /**
             * Bulk viscosity of mixture of each species, \f$ \lambda_i \quad \left[ \frac{\text{g}}{\text{cm sec}} \right] \f$.
             * in the quadrature points of a boundary
             * at a previous Newton iteration.
             * NOTE: first index is the quadrature point and second index is the species.
             */
            std::vector< std::vector<double> > bulk_viscosity_bdry_old;
            
            /**
             * paramMatrix is used for calculating partial viscosity of mixture in the quadrature points of a boundary
             * at a previous Newton iteration.
             * NOTE: first index is the quadrature point and second and third indices are the species.
             */
            std::vector< std::vector< std::vector<double> > > paramMatrix_bdry_old;
            
            /**
             * PInv_bdry_old is used for calculating the inverse of P when using the OmegaKG model for calculating partial viscosity
             * in the quadrature points of a boundary
             * at a previous Newton iteration.
             */
            std::vector< FullMatrix<double> > PInv_bdry_old;
            
            //@}
            
            ///@name Local CG FEM based assemblers - variable data (current Newton iteration - boundary)
            //@{
            
            /**
             * Density shape functions.
             *
             * \p phi_density_bdry \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th density shape function
             * computed in \f$ q \f$-th quadrature point of a boundary
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector<double> > > phi_density_bdry;
            
            /**
             * Density shape function gradients.
             *
             * \p grad_phi_density_bdry \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th density shape function gradient
             * computed in \f$ q \f$-th quadrature point of a boundary
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector< Tensor<1,dim> > > > grad_phi_density_bdry;
            
            /**
             * Velocity shape functions.
             *
             * \p phi_velocity_bdry \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th velocity shape function
             * computed in \f$ q \f$-th quadrature point of a boundary
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector< Tensor<1,dim> > > > phi_velocity_bdry;
            
            /**
             * Velocity shape function divergences.
             *
             * \p div_phi_velocity_bdry \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th velocity shape function divergence
             * computed in \f$ q \f$-th quadrature point of a boundary
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector<double> > > div_phi_velocity_bdry;
            
            /**
             * Velocity shape function symmetric gradients.
             *
             * \p grads_phi_velocity_bdry \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th velocity shape function symmetric gradient
             * computed in \f$ q \f$-th quadrature point of a boundary
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector< SymmetricTensor<2,dim> > > > grads_phi_velocity_bdry;
            
            /**
             * Mass flux perturbations.
             *
             * \p delta_mass_flux_bdry \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th mass flux perturbation
             * computed in \f$ q \f$-th quadrature point of a boundary
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector< Tensor<1,dim> > > > delta_mass_flux_bdry;
            
            /**
             * Momentum flux perturbations.
             *
             * \p delta_momentum_flux_bdry \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th momentum flux perturbation
             * computed in \f$ q \f$-th quadrature point of a boundary
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector< SymmetricTensor<2,dim> > > > delta_momentum_flux_bdry;
            
            /**
             * Partial pressure perturbations.
             *
             * \p delta_pressure_bdry \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th partial pressure perturbation
             * computed in \f$ q \f$-th quadrature point of a boundary
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector<double> > > delta_pressure_bdry;
            
            /**
             * Shear stress perturbations.
             *
             * \p delta_shear_stress_bdry \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th shear stress perturbation
             * computed in \f$ q \f$-th quadrature point of a boundary
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector< SymmetricTensor<2,dim> > > > delta_shear_stress_bdry;
            
            /**
             * Partial dyanmic viscosity perturbations.
             *
             * \p delta_viscosity_mix_bdry \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th partial dyanmic viscosity perturbation
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector<double> > > delta_partial_viscosity_bdry;
            
            /**
             * Partial bulk viscosity perturbations.
             *
             * \p delta_viscosity_mix_bdry \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th partial bulk viscosity perturbation
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector<double> > > delta_bulk_viscosity_bdry;
            
            //@}
            
            ///@name Local CG FEM based assemblers - variable data (other - boundary)
            //@{
            
            /**
             * Symmetrized tensor product of
             * \p normal_vectors \p[ \p q \p]
             * and
             * \p phi_velocity_bdry \p[ \p s \p] \p[ \p q \p] \p[ \p k \p]
             * in the quadrature points of a boundary.
             */
            std::vector< std::vector< std::vector< SymmetricTensor<2,dim> > > > n_BY_phi_velocity_bdry_S;
            
            //@}
            ///@name Auxiliary classes
            //@{
            /**
             * Inline function to obtain update the porosity for a given layer class.
             * 
             * \note The pointer should be of type:
             * - FuelCellShop::Layer::GasDiffusionLayer<dim>
             * - FuelCellShop::Layer::MicroPorousLayer<dim>
             * - FuelCellShop::Layer::CatalystLayer<dim>
             * - FuelCellShop::Layer::ExperimentalPorousLayer<dim>
             */
            template <typename LAYER, typename INFO>
            inline void update_porosity(const LAYER* ptr, 
                                        const INFO info, 
                                        std::vector<double>& porosity)
            {
                if( ptr->get_porosity_is_constant() ) 
                    ptr->get_porosity(porosity);
                else
                    ptr->get_porosity(porosity, info.get_fe_val_unsplit().get_quadrature_points());
            }
            
            /**
             * Inline member function used to obtain the inverse permeability and the Forchheimer permeability
             * for a given layer.
             * 
             * \note The pointer should be of type:
             * - FuelCellShop::Layer::GasDiffusionLayer<dim>
             * - FuelCellShop::Layer::MicroPorousLayer<dim>
             * - FuelCellShop::Layer::CatalystLayer<dim>
             * - FuelCellShop::Layer::ExperimentalPorousLayer<dim>
             */
            template <typename LAYER, typename INFO>
            inline void update_permeability(const LAYER* ptr, 
                                            const INFO info, 
                                            std::vector< SymmetricTensor<2,dim> >& permeability_INV,
                                            std::vector< SymmetricTensor<2,dim> >& Forchheimer_permeability)
            {
                if( ptr->get_permeability_is_constant() )
                {
                    ptr->get_permeability_INV(permeability_INV);
                    ptr->get_Forchheimer_permeability(Forchheimer_permeability);
                }
                else
                {
                    ptr->get_permeability_INV(permeability_INV, info.get_fe_val_unsplit().get_quadrature_points());
                    ptr->get_Forchheimer_permeability(Forchheimer_permeability, info.get_fe_val_unsplit().get_quadrature_points());
                }
            }
            /**
             * Inline member function used to obtain the tortuosity for a given layer.
             * 
             * \note The pointer should be of type:
             * - FuelCellShop::Layer::GasDiffusionLayer<dim>
             * - FuelCellShop::Layer::MicroPorousLayer<dim>
             * - FuelCellShop::Layer::CatalystLayer<dim>
             * - FuelCellShop::Layer::ExperimentalPorousLayer<dim>
             */
            template <typename LAYER, typename INFO>
            inline void update_tortuosity(const LAYER* ptr, 
                                          const INFO info, 
                                          std::vector< SymmetricTensor<2,dim> >& tortuosity )
            {
                if( ptr->get_tortuosity_is_constant() )
                    ptr->get_tortuosity(tortuosity);
                else
                    ptr->get_tortuosity(tortuosity, info.get_fe_val_unsplit().get_quadrature_points());
            }
            
            /**
             * Inline member function used to update 
             * - temperature, 
             * - molar_mass,
             * - dynamic_viscosity, 
             * - maxwell_stefan_isobaric_diffusion_coefficient
             */
            template <typename LAYER>
            inline void update_gas_properties(const LAYER* ptr)
            {
                //TODO: these parameters should be coming from the same location not different classes
                T_mixture = ptr->get_gas_mixture()->get_temperature();
                P_in = ptr->get_gas_mixture()->get_total_pressure();
                
                for(unsigned int s = 0; s < n_species; ++s)
                {
                    molar_mass[s] = ptr->get_gas_mixture()->get_gases()[s]->get_molar_mass() * 1.0e3; // Convert kg/mol to 1000 g/mol
                    dynamic_viscosity[s] = ptr->get_gas_mixture()->get_gases()[s]->get_dynamic_viscosity(T_mixture) * 1.0e1; // Convert kg/m*s to 10 g/cm*s
                    collision_diameter[s] = ptr->get_gas_mixture()->get_gases()[s]->get_collision_diameter() * 1e-10; // Convert Angstroms to 1e-10 meters
                    maxwell_stefan_isobaric_diffusion_coefficient = ptr->get_gas_mixture()->get_ChapmanEnskog_isobaric_diffusion_coefficients(T_mixture); // [m^2/s]
                }
            }
            
            /**
             * This private member function is used in make_assemblers_cell_variable_data in order to compute the inverse of the
             * product of pressure and molecular diffusion coefficient. This value is only a function of temperature.
             * 
             * Parameters:
             * - prt: Pointer to the layer we would like to compute the coefficient for
             * - inv_pDeff: returned value for the coefficient. The vector is re-initialized in the class, so an empty vector can
             * be provided. It will be re-sized and filled as appropriate.
             * 
             * \warning This routine still needs work...
             * 
             */
            template <typename LAYER>
            inline void update_inv_diffusion_coefficients(const LAYER* ptr, 
                                                          std::vector< Table< 2, SymmetricTensor<2,dim> > >& inv_pDeff) const
            {
                //************************************
                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                // !!!!!!!!!!!!!! THIS NEEDS TO BE IMPROVED !!!!!!!!!!!!!!
                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                //************************************   
                Table< 2, Tensor<2,dim> > Deff;
                
                //I would like to replace this...
                ptr->effective_gas_diffusivity(Deff); // m^2/s
                // with
                /*
                ptr->get_T_and_p(T,p);
                SolutionVariable temperature = SolutionVariable(T, this->n_q_points_cell , temperature_of_REV);
                ptr->set_temperature( temperature );
                SolutionVariable pressure = SolutionVariable(1.0/101325.0, this->n_q_points_cell, temperature_of_REV); // Value of 1.0/101325.0 is equivalent to 1 Pa, i.e. we return D(T)
                const std::vector< PureGas *> gases = ptr->get_gases ();
                */
                /* If saturation:
                if (p_liquid_water.indices_exist)
                {
                    ptr->set_capillary_pressure( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][p_liquid_water.solution_index], capillary_pressure) );
                    deriv_flags.push_back(capillary_pressure);
                }
                */              
                
                unsigned int n_q_points_cell = inv_pDeff.size();
                
                for(unsigned int s = 0; s < n_species; ++s)
                    for(unsigned int s1 = 0; s1 < n_species; ++s1)
                    {
                        // std::vector< double > Deff_noniso;
                        // ptr->compute_gas_diffusion(gases[s], gases[s1]);
                        // ptr->effective_gas_diffusivity(Deff_noniso);            // m^2/s
                        
                        if( s1 != s )  
                        {
                            for(unsigned int q = 0; q < this->n_q_points_cell; ++q)   
                            {                                                               
                                Tensor<2, dim> pDeff_CGS = P_in*Deff(s,s1)*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);
                                // Units are Pa*cm^2/s, however Pa are NOT CGS, ( 1 Pa = 10 g/(cm s^2):
                                pDeff_CGS *= 10;
                                // Now, compute 1/(p*Deff):                                    
                                Tensor<2, dim> invTest;
                                invTest[0][0] = 1 / ( pDeff_CGS[0][0] );
                                invTest[1][1] = 1 / ( pDeff_CGS[1][1] );
                                #if deal_II_dimension == 3
                                    invTest[2][2] = 1 / ( pDeff_CGS[2][2] );
                                #endif
                                inv_pDeff[q](s,s1) = invTest;
                            }
                        }
                    }                    

                    //Print out Deff* for debugging
//                     bool printDeff = false;
//                     
//                     if( printDeff )
//                     {
//                         FcstUtilities::log << "-----------------------------------" <<std::endl;
//                         for( unsigned int s1 = 0; s1 < n_species; ++s1 )
//                         {
//                             for( unsigned int s2 = 0; s2 < n_species; ++s2 )
//                             {
//                                 if( s1 != s2 )
//                                 {
//                                     FcstUtilities::log << "Name layer is: "<<ptr->name_layer()<<std::endl;
//                                     FcstUtilities::log << "Species " << s1 << " and " << s2 << ": " <<std::endl;
//                                     FcstUtilities::log << "Deff[0][0]: " << Deff(s1,s2)[0][0]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2) << std::endl;
//                                     FcstUtilities::log << "Deff[1][1]: " << Deff(s1,s2)[1][1]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2) << std::endl;
//                                     
//                                     FcstUtilities::log << "inv_pDeff[0][0]: " << inv_pDeff[0](s1,s2)[0][0] << std::endl;
//                                     FcstUtilities::log << "inv_pDeff[1][1]: " << inv_pDeff[0](s1,s2)[1][1] << std::endl;
//                                 }
//                             }
//                         }
//                         FcstUtilities::log << "-----------------------------------" <<std::endl;
//                     }
            }
            
            /**
             * This private member function is used to calculate the partial viscosity of the mixture and bulk viscosity, from
             * dynamicViscosity, density, molar_mass. As well, it returns paramMatrix for possible later calculations (like when
             * calculating the variation). Since this function is used by make_assemblers_cell_variable_data and 
             * make_assemblers_bdry_variable_data which have different quadrature points one of the members is the number of 
             * quadrature points. \b NOTE: no need to size vectors paramMatrix, PInv, partialViscosity, and bulkViscosity
             * sizing is done by this function.
             * 
             * Parameters:
             * @param prt              - Pointer to the layer we would like to compute the coefficient for
             * @param quadraturePoints - number of quadrature points in cell/boundary
             * @param porosity         - porosity in cell, vector of quadrature points
             * @param density          - vector of vectors of density; first indice is species and second is quadrature point; [g/cm^3]
             * @param paramMatrix      - used internally; vector of vector of vectors used in calculation of partialViscosity; first indice is quadrature point, second is species 1 and second is species 2
             * @param PInv             - used internally; this parameter provides the inverse of P for the calculation of OmegaKG viscosity and variation of,
             * @param partialViscosity - vector of vectors of partial viscosity of mixture; first indice is quadrature point, second is species; [g/cm*s]
             * @param bulkViscosity:   - vector of vectors of bulk viscosities; first indice is quadrature point, second is species; [g/cm*s]
             * 
             */
            template <typename LAYER>
            inline void update_partial_viscosities(const LAYER*                                       ptr,
                                                   const unsigned int&                                quadraturePoints,
                                                   const std::vector<double>&                         porosity,
                                                   const std::vector< std::vector<double> >&          density,
                                                   std::vector< std::vector< std::vector<double> > >& paramMatrix,
                                                   std::vector< FullMatrix<double> >&                 PInv,
                                                   std::vector< std::vector<double> >&                partialViscosity,
                                                   std::vector< std::vector<double> >&                bulkViscosity) const
            {

                partialViscosity = ptr->get_gas_mixture()->get_isothermal_nonisobaric_partial_viscosity(T_mixture,          // [K]
                                                                                                        density,            // [g/cm^3]
                                                                                                        molar_mass,         // [g/mol]
                                                                                                        dynamic_viscosity,  // [g/cm*s]
                                                                                                        collision_diameter, // [m]
                                                                                                        porosity, 
                                                                                                        paramMatrix,
                                                                                                        PInv); // returns in [g/cm*s]
                
                bulkViscosity.clear(); //reset bulk viscosity vector
                bulkViscosity.resize(quadraturePoints, std::vector<double>(this->n_species, 0.0)); //size appropriately
                for(unsigned int q = 0; q < quadraturePoints; ++q)
                {
                    double porosityInv = 1 / porosity[q];
                    for(unsigned int s = 0; s < this->n_species; ++s)
                        bulkViscosity[q][s] = ptr->get_gas_mixture()->get_gases()[s]->get_bulk_viscosity(partialViscosity[q][s]); // [g/cm*s]
                }
            }
            
            /**
             * This private member function is used to calculate the variation in partial viscosity of the mixture.
             * 
             * Parameters:
             * @param prt                   - Pointer to the layer we would like to compute the coefficient for
             * @param partialViscosity      - vector of vectors of partial viscosity of mixture; first indice is quadrature point, second is species; [g/cm*s]
             * @param paramMatrix           - vector of vector of vectors used in calculation of partialViscosity; first indice is quadrature point, second is species 1 and second is species 2
             * @param porosity              - porosity in cell, vector of quadrature points
             * @param deltaDensity          - variation in density; same unit standard as density
             * @param density               - vector of vectors of density; first indice is species and second is quadrature point; [g/cm^3]
             * @param deltaPartialViscosity - vector returned by reference contains variation of partial viscosity; first indice is species, second is quadrature point, and third is dof
             */
            template <typename LAYER>
            inline void update_delta_partial_viscosities(const LAYER*                                             ptr,
                                                         const std::vector< std::vector<double> >&                partialViscosity,
                                                         const std::vector< std::vector< std::vector<double> > >& paramMatrix,
                                                               std::vector< FullMatrix<double> >&                 PInv,
                                                         const std::vector<double>&                               porosity,
                                                         const std::vector< std::vector< std::vector<double> > >& deltaDensity,
                                                         const std::vector< std::vector<double> >&                density,
                                                               std::vector< std::vector< std::vector<double> > >& deltaPartialViscosity) const
            {
                ptr->get_gas_mixture()->get_isothermal_nonisobaric_delta_partial_viscosity(T_mixture,          // [K]
                                                                                           paramMatrix,
                                                                                           PInv,
                                                                                           porosity,
                                                                                           molar_mass,         // [g/mol]
                                                                                           dynamic_viscosity,  // [g/cm*s]
                                                                                           collision_diameter, // [m]
                                                                                           deltaDensity,
                                                                                           density,            // [g/cm^3]
                                                                                           deltaPartialViscosity); // returns in cgs units
            }
            
            /**
             * Calls GasMixture function to loop through all variations of partial viscosity and calculates the varation in bulk viscosity according to
             * bulk_viscosity_mode. All vectors must be sized correctly before using this function.
             */
            template <typename LAYER>
            inline void update_delta_bulk_viscosities(const LAYER*                                             ptr,
                                                      const std::vector< std::vector< std::vector<double> > >& deltaPartialViscosity,
                                                            std::vector< std::vector< std::vector<double> > >& deltaBulkViscosity) const
            {
                ptr->get_gas_mixture()->get_delta_bulk_viscosity(ptr->get_gas_mixture()->get_gases(), 
                                                                 deltaPartialViscosity, 
                                                                 deltaBulkViscosity);
            }
            
        };
                 
    } // Equation
                 
} // FuelCellShop
             
#endif