// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: incompressible_single_component_NS_equations.h
// - Description: This class describes steady-state incompressible and isothermal Navier-Stokes
//   fluid transport equations for a single-phase single-component case
// - Developers: Valentin N. Zingan, University of Alberta
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_EQUATION_INCOMPRESSIBLE_SINGLE_COMPONENT_NS_EQUATIONS_H_
#define _FCST_FUELCELLSHOP_EQUATION_INCOMPRESSIBLE_SINGLE_COMPONENT_NS_EQUATIONS_H_

#include <equations/equation_base.h>
#include <layers/channel.h>
#include <layers/experimental_porous_layer.h>
#include <utils/fcst_units.h>

namespace FuelCellShop
{
    namespace Equation
    {

        /**
        * This class deals with
        * Navier-Stokes fluid transport equations.
        *
        * According to classification,
        * these equations are
        *
        * - steady-state,
        * - incompressible,
        * - isothermal,
        * - single-phase,
        * - single-component,
        *
        * - supposed to take
        *   the point-wise form
        *   in the channels of
        *   fuel cells,
        * - appropriately averaged
        *   over a so-called
        *   Representative Elementary Volume (REV)
        *   in the porous media regions,
        * - solved with respect to point-wise
        *   \f$ \left( p, \mathbf{u} \right) \f$
        *   in channels and with respect to REV-wise
        *   \f$ \left( <p>, <\mathbf{u}> \right) \f$
        *   in porous media.
        */

        /**
        *
        * These equations can be written as
        *
        * \f$ \nabla \cdot \mathbf{u} = 0 \quad \text{in} \quad \Omega \f$ \f$ \quad - \quad \f$ mass conservation equation
        *
        * \f$ \rho \nabla \cdot \left( \frac{1}{\epsilon} \mathbf{u} \otimes \mathbf{u} \right) = \nabla \cdot \left( -p \hat{\mathbf{I}} + \hat{\boldsymbol\sigma} \right) + \left( \mathbf{F} + \epsilon \rho \mathbf{g} \right)
        * \quad \text{in} \quad \Omega \f$ \f$ \quad - \quad \f$ momentum conservation equation
        *
        * where
        *
        * \f$ \hat{\boldsymbol\sigma} = 2 \mu \nabla_s \mathbf{u} \f$ \f$ \quad - \quad \f$ shear stress
        *
        * \f$
        *   \mathbf{F} =
        *   \begin{cases}
        *   \mathbf{0} \quad \text{in} \quad \Omega_c \\
        * - \mu \hat{ \mathbf{K} }^{-1} \mathbf{u} - \hat{\boldsymbol\beta} \rho \big| \mathbf{u} \big| \mathbf{u} \quad \text{in} \quad \Omega_p
        *   \end{cases}
        * \f$ \f$ \quad - \quad \f$ drag force
        */

        /**
        *
        * To be well-posed, these equations are equipped with appropriate boundary conditions.
        * Below we consider several types of boundaries and several types of boundary conditions assigned to those boundaries (here \f$ k \f$ denotes some portion of a boundary of a certain type).
        *
        * - Impermeable walls \f$ \left( \Gamma^{\text{walls},k} \right) \f$ \f$ \quad - \quad \f$ no-slip boundary condition
        *
        * \f$
        * \quad \mathbf{u} \Bigg|_{\Gamma^{\text{walls},k}} = \mathbf{0} \quad \text{OR}
        * \f$
        *
        * - Impermeable walls \f$ \left( \Gamma^{\text{walls},k} \right) \f$ \f$ \quad - \quad \f$ Navier slip boundary condition
        *
        * \f$
        * \mathbf{u} \cdot \mathbf{n} \Bigg|_{\Gamma^{\text{walls},k}} = 0 \quad \text{and} \quad
        * \left(1-\theta \right) \mathbf{u} \cdot \boldsymbol\tau_{\alpha} + \theta \left( -p \hat{\mathbf{I}} + \hat{\boldsymbol\sigma} \right) \mathbf{n} \cdot \boldsymbol\tau_{\alpha} \Bigg|_{\Gamma^{\text{walls},k}} = 0 \quad
        * \text{with} \quad \alpha = 1, d-1 \quad \text{and} \quad 0 < \theta < 1 \quad \text{OR}
        * \f$
        *
        * - Impermeable walls \f$ \left( \Gamma^{\text{walls},k} \right) \f$ \f$ \quad - \quad \f$ perfect slip boundary condition
        *
        * \f$
        * \mathbf{u} \cdot \mathbf{n} \Bigg|_{\Gamma^{\text{walls},k}} = 0 \quad \text{and} \quad
        * \left( -p \hat{\mathbf{I}} + \hat{\boldsymbol\sigma} \right) \mathbf{n} \cdot \boldsymbol\tau_{\alpha} \Bigg|_{\Gamma^{\text{walls},k}} = 0 \quad
        * \text{with} \quad \alpha = 1, d-1
        * \f$
        *
        * - Symmetry line or plane \f$ \left( \Gamma^{\text{sym}} \right) \f$ \f$ \quad - \quad \f$ perfect slip boundary condition
        *
        * \f$
        * \mathbf{u} \cdot \mathbf{n} \Bigg|_{\Gamma^{\text{sym}}} = 0 \quad \text{and} \quad
        * \left( -p \hat{\mathbf{I}} + \hat{\boldsymbol\sigma} \right) \mathbf{n} \cdot \boldsymbol\tau_{\alpha} \Bigg|_{\Gamma^{\text{sym}}} = 0 \quad
        * \text{with} \quad \alpha = 1, d-1
        * \f$
        *
        * - Inlet-Outlet \f$ \left( \Gamma^{\text{in-out},k} \right) \f$ \f$ \quad - \quad \f$ Dirichlet boundary condition
        *
        * \f$
        * p \Bigg|_{\Gamma^{\text{in-out},k}} = p^{\text{prescribed},k} \quad \text{or} \quad \mathbf{u} \Bigg|_{\Gamma^{\text{in-out},k}} = \mathbf{u}^{\text{prescribed},k}
        * \quad \text{or} \quad \text{both} \quad \text{OR/AND}
        * \f$
        *
        * - Inlet-Outlet \f$ \left( \Gamma^{\text{in-out},k} \right) \f$ \f$ \quad - \quad \f$ normal stress free boundary condition
        *
        * \f$
        * \left( -p \hat{\mathbf{I}} + \hat{\boldsymbol\sigma} \right) \mathbf{n} \Bigg|_{\Gamma^{\text{in-out},k}} = \mathbf{0}
        * \f$
        */

        /**
        *
        * All these boundaries and boundary conditions can be specified in the parameters file as follows:
        *
        * @code
        * Impermeable walls = [boundary_id]: [type]
        * @endcode
        *
        * where
        *
        * - [boundary_id] = 0,1,2,3,...
        * - [type]        = no-slip | Navier slip | perfect slip
        *
        * @code
        * Symmetry line or plane = [boundary_id]: [type]
        * @endcode
        *
        * where
        *
        * - [boundary_id] = 0,1,2,3,...
        * - [type]        = perfect slip
        *
        * @code
        * Inlet-Outlet = [boundary_id]: [type]
        * @endcode
        *
        * where
        *
        * - [boundary_id] = 0,1,2,3,...
        * - [type]        = Dirichlet pressure | Dirichlet pressure and normal stress free | Dirichlet velocity | Dirichlet pressure and velocity | normal stress free
        */

        /**
        *
        * These equations utilize
        *
        * - \p Channel,
        * - \p ExperimentalPorousLayer
        *
        * as layer classes
        * and
        *
        * - \p ExperimentalFluid,
        * - \p ExperimentalSolid
        *
        * as material classes
        * and
        * have the following
        * output:
        *
        * - pressure, \f$ p          \f$ \f$ \quad \left[ \text{Pa}                   \right] \f$
        * - velocity, \f$ \mathbf{u} \f$ \f$ \quad \left[ \frac{\text{m}}{\text{sec}} \right] \f$
        */

        /**
        *
        * The whole problem is solved
        * by linearizing the governing equations
        * at each Newton iteration with subsequent
        * CG FEM discretization in space.
        *
        * The implementation is
        * performed by means of using only
        * one weak formulation
        * for all \f$ d + 1 \f$
        * scalar governing equations.
        * This approach guarantees
        * the dim-independent
        * programming of
        * the equations
        * at hand.
        *
        * \author Valentin N. Zingan, 2012
        */

        template<int dim>
        class IncompressibleSingleComponentNSEquations : public EquationBase<dim>
        {
            public:

                ///@name Constructors, destructor, and initialization
                //@{

                /**
                * Constructor.
                */
                IncompressibleSingleComponentNSEquations(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = 
                            boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >());

                /**
                * Destructor.
                */
                virtual ~IncompressibleSingleComponentNSEquations();

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
                const double& get_inlet_outlet_velocity_max() const
                {
                    return inlet_outlet_velocity_max;
                }
                
                /**
                * User can specify an equation to use for density/velocity profile at boundary.
                * E. \p inlet-outlet \p equation = -625*y*y + 66.25*y - 0.755625
                */
                const std::string get_inlet_outlet_equation() const
                {
                    return inlet_outlet_equation;
                }
                
                /**
                * Outputs the Boundary ID to apply Dirichlet velocity or Dirichlet pressure to.
                */
                const unsigned int get_inlet_outlet_boundary_ID() const
                {
                    return inlet_outlet_boundary_ID;
                }
                
                /**
                * Outputs which equation to apply equation to (Pressure, VelocityX, VelocityY,
                * VelocityZ).
                */
                const std::string get_press_vel_comp_apply_to() const
                {
                    return press_vel_comp_apply_to;
                }

                /**
                * This function prints out
                * the equations info.
                */
                virtual void print_equation_info() const;

                //@}

            protected:

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
                * whether the inertia term
                * is ON or OFF in channels.
                */
                bool inertia_in_channels;

                /**
                * This object indicates
                * whether the shear stress term
                * is ON or OFF in channels.
                */
                bool shear_stress_in_channels;

                /**
                * This object indicates
                * whether the gravity term
                * is ON or OFF in channels.
                */
                bool gravity_in_channels;

                /**
                * This object indicates
                * whether the inertia term
                * is ON or OFF in porous media.
                */
                bool inertia_in_porous_media;

                /**
                * This object indicates
                * whether the shear stress term
                * is ON or OFF in porous media.
                */
                bool shear_stress_in_porous_media;

                /**
                * This object indicates
                * whether the gravity term
                * is ON or OFF in porous media.
                */
                bool gravity_in_porous_media;

                /**
                * This object indicates
                * which form of the drag term
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
                std::string drag_in_porous_media;

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
                * of velocity
                * be equal
                * to 0.
                *
                * The alternative approach
                * is doing that by introducing
                * a penalization term \f$ \displaystyle \int_{ \Gamma^{\text{walls},k} } \frac{1}{\eta} \left( \boldsymbol\omega \cdot \mathbf{n} \right) \left( \mathbf{u} \cdot \mathbf{n} \right) dS \f$
                * with \f$ \eta \sim 10^{-10} - 10^{-12} \f$
                * into the weak formulation
                * of the equations
                * at hand.
                *
                * This object indicates
                * whether the penalty
                * method is used.
                */
                bool normal_velocity_is_suppressed_weakly;

                //@}

                ///@name Computational constants
                //@{

                /**
                * Gravitational acceleration, \f$ \mathbf{g} = \{ g_{\alpha} \}_{\alpha = 1}^d \quad \text{such that} \quad \forall \alpha \neq d : \quad g_{\alpha} = 0 \quad \left[ \frac{\text{m}}{\text{sec}^2} \right] \quad
                * \text{and} \quad g_d = -9.81 \quad \left[ \frac{\text{m}}{\text{sec}^2} \right] \f$.
                */
                Tensor<1,dim> gravity_acceleration;

                /**
                * Unit tensor, \f$ \hat{ \mathbf{I} } = \{ \delta_{\alpha \beta} \}_{\alpha,\beta = 1}^d \f$.
                */
                SymmetricTensor<2,dim> unit;

                /**
                * Navier slip coefficient, \f$ \theta \f$.
                */
                double theta;

                /**
                * Normal velocity suppression coefficient, \f$ \eta \f$.
                */
                double eta;

                /**
                * Maximum inlet-outlet velocity, \f$ \mathbf{u}_{\text{in,out}}^{\text{max}} \f$ \f$ \quad \left[ \frac{\text{m}}{\text{sec}} \right] \f$.
                */
                double       inlet_outlet_velocity_max;
                std::string  inlet_outlet_equation;
                unsigned int inlet_outlet_boundary_ID;
                std::string  press_vel_comp_apply_to;

                //@}

        };

    } // Equation

} // FuelCellShop

#endif