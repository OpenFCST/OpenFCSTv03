// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: incompressible_single_component_NS_equations.cc
// - Description: This class describes steady-state incompressible and isothermal Navier-Stokes
//   fluid transport equations for a single-phase single-component case
// - Developers: Valentin N. Zingan, University of Alberta
//
// ----------------------------------------------------------------------------

#include "equations/incompressible_single_component_NS_equations.h"

namespace NAME = FuelCellShop::Equation;

//////////////////////////////////////////////////
//////////////////////////////////////////////////
// CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
//////////////////////////////////////////////////
//////////////////////////////////////////////////

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::IncompressibleSingleComponentNSEquations<dim>::IncompressibleSingleComponentNSEquations(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data)
:
NAME::EquationBase<dim>(system_management,data)
{
    FcstUtilities::log << "->FuelCellShop::Equation::IncompressibleSingleComponentNSEquations" << std::endl;

    gravity_acceleration = Constants::gravity_acceleration();
    unit                 = Constants::unit_tensor();

    this->make_internal_cell_couplings();

    this->counter.resize(3, true);
    this->counter[2] = false;
}

// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::IncompressibleSingleComponentNSEquations<dim>::~IncompressibleSingleComponentNSEquations()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::IncompressibleSingleComponentNSEquations<dim>::declare_parameters(ParameterHandler& param) const
{
    param.enter_subsection("Equations");
    {
        param.enter_subsection("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component");
        {
            param.enter_subsection("Boolean constants and form of the drag force in porous media");
            {
                param.declare_entry("Inertia in channels",
                                    "true",
                                    Patterns::Bool(),
                                    " ");
                param.declare_entry("Shear stress in channels",
                                    "true",
                                    Patterns::Bool(),
                                    " ");
                param.declare_entry("Gravity in channels",
                                    "false",
                                    Patterns::Bool(),
                                    " ");

                param.declare_entry("Inertia in porous media",
                                    "true",
                                    Patterns::Bool(),
                                    " ");
                param.declare_entry("Shear stress in porous media",
                                    "true",
                                    Patterns::Bool(),
                                    " ");
                param.declare_entry("Gravity in porous media",
                                    "false",
                                    Patterns::Bool(),
                                    " ");
                param.declare_entry("Drag in porous media",
                                    "none",
                                    Patterns::Selection("none | Darcy | Forchheimer | Forchheimer modified"),
                                    " ");

                param.declare_entry("Normal velocity is suppressed weakly",
                                    "false",
                                    Patterns::Bool(),
                                    " ");
            }
            param.leave_subsection();

            param.enter_subsection("Computational constants");
            {
                param.declare_entry("Navier slip coefficient",
                                    "0.0",
                                    Patterns::Double(0.0,1.0),
                                    " ");
                param.declare_entry("Normal velocity suppression coefficient",
                                    "1.0e-12",
                                    Patterns::Double(0.0,1.0),
                                    " ");
                param.declare_entry("Maximum inlet-outlet velocity [m/s]",
                                    "0.0",
                                    Patterns::Double(),
                                    " ");
                param.declare_entry("inlet-outlet equation",
                                    "",
                                    Patterns::Anything(),
                                    "For user defined density/velocity profile at inlet/outlet, ex -625*y*y + 66.25*y - 0.755625");
                param.declare_entry("pressure or velocity component to apply equation to",
                                    "",
                                    Patterns::Anything(),
                                    "Define which equation to apply equation to, ex Pressure, VelocityX, VelocityY, or VelocityZ");
            }
            param.leave_subsection();

            param.enter_subsection("Initial data");
            {
                param.declare_entry("Variable initial data",
                                    "false",
                                    Patterns::Bool(),
                                    " ");

                param.declare_entry("single_fluid_pressure",
                                    """",
                                    Patterns::Map( Patterns::Integer(0) , Patterns::Double() ),
                                    " ");

                param.declare_entry("single_fluid_velocity_X",
                                    """",
                                    Patterns::Map(   Patterns::Integer(0) , Patterns::Double() ),
                                    " ");

                param.declare_entry("single_fluid_velocity_Y",
                                    """",
                                    Patterns::Map( Patterns::Integer(0) , Patterns::Double() ),
                                    " ");

                param.declare_entry("single_fluid_velocity_Z",
                                    """",
                                    Patterns::Map( Patterns::Integer(0) , Patterns::Double() ),
                                    " ");
            }
            param.leave_subsection();

            param.enter_subsection("Boundary conditions");
            {
                param.declare_entry("Impermeable walls",
                                    """",
                                    Patterns::Map( Patterns::Integer(0) , Patterns::Selection("no-slip | Navier slip | perfect slip") ),
                                    " ");

                param.declare_entry("Symmetry line or plane",
                                    """",
                                    Patterns::Map( Patterns::Integer(0) , Patterns::Selection("perfect slip") ),
                                    " ");

                param.declare_entry("Inlet-Outlet",
                                    """",
                                    Patterns::Map( Patterns::Integer(0) , Patterns::Selection("Dirichlet pressure | Dirichlet pressure and normal stress free | Dirichlet velocity | Dirichlet pressure and velocity | normal stress free") ),
                                    " ");
            }
            param.leave_subsection();

            param.enter_subsection("Boundary data");
            {
                param.declare_entry("Variable boundary data",
                                    "false",
                                    Patterns::Bool(),
                                    " ");

                param.declare_entry("single_fluid_pressure",
                                    """",
                                    Patterns::Map( Patterns::Integer(0) , Patterns::Double() ),
                                    " ");

                param.declare_entry("single_fluid_velocity_X",
                                    """",
                                    Patterns::Map( Patterns::Integer(0) , Patterns::Double() ),
                                    " ");

                param.declare_entry("single_fluid_velocity_Y",
                                    """",
                                    Patterns::Map( Patterns::Integer(0) , Patterns::Double() ),
                                    " ");

                param.declare_entry("single_fluid_velocity_Z",
                                    """",
                                    Patterns::Map( Patterns::Integer(0) , Patterns::Double() ),
                                    " ");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

// ---            ---
// --- initialize ---
// ---            ---

template<int dim>
void
NAME::IncompressibleSingleComponentNSEquations<dim>::initialize(ParameterHandler& param)
{
    param.enter_subsection("Equations");
    {
        param.enter_subsection("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component");
        {
            param.enter_subsection("Boolean constants and form of the drag force in porous media");
            {
                inertia_in_channels          = param.get_bool("Inertia in channels");
                shear_stress_in_channels     = param.get_bool("Shear stress in channels");
                gravity_in_channels          = param.get_bool("Gravity in channels");
                inertia_in_porous_media      = param.get_bool("Inertia in porous media");
                shear_stress_in_porous_media = param.get_bool("Shear stress in porous media");
                gravity_in_porous_media      = param.get_bool("Gravity in porous media");
                drag_in_porous_media         = param.get("Drag in porous media");

                if( drag_in_porous_media != "none" )
                    this->counter[2] = true;

                normal_velocity_is_suppressed_weakly = param.get_bool("Normal velocity is suppressed weakly");
            }
            param.leave_subsection();

            param.enter_subsection("Computational constants");
            {
                theta                     = param.get_double("Navier slip coefficient");
                eta                       = param.get_double("Normal velocity suppression coefficient");
                inlet_outlet_velocity_max = param.get_double("Maximum inlet-outlet velocity [m/s]");
                inlet_outlet_equation     = param.get("inlet-outlet equation");
                press_vel_comp_apply_to   = param.get("pressure or velocity component to apply equation to");
            }
            param.leave_subsection();

            param.enter_subsection("Initial data");
            {
                this->variable_initial_data = param.get_bool("Variable initial data");

                if( !param.get("single_fluid_pressure").empty() )
                {
                    const std::map<types::material_id, double> tmp = FcstUtilities::string_to_map<types::material_id, double>( param.get("single_fluid_pressure") );
                    this->component_materialID_value["single_fluid_pressure"] = tmp;
                }

                if( !param.get("single_fluid_velocity_X").empty() )
                {
                    const std::map<types::material_id, double> tmp = FcstUtilities::string_to_map<types::material_id, double>( param.get("single_fluid_velocity_X") );
                    this->component_materialID_value["single_fluid_velocity_X"] = tmp;
                }

                if( !param.get("single_fluid_velocity_Y").empty() )
                {
                    const std::map<types::material_id, double> tmp = FcstUtilities::string_to_map<types::material_id, double>( param.get("single_fluid_velocity_Y") );
                    this->component_materialID_value["single_fluid_velocity_Y"] = tmp;
                }

                if( !param.get("single_fluid_velocity_Z").empty() )
                {
                    const std::map<types::material_id, double> tmp = FcstUtilities::string_to_map<types::material_id, double>( param.get("single_fluid_velocity_Z") );
                    this->component_materialID_value["single_fluid_velocity_Z"] = tmp;
                }
            }
            param.leave_subsection();

            param.enter_subsection("Boundary conditions");
            {
                if( !param.get("Impermeable walls").empty() )
                {
                    const std::map<unsigned int, std::string> tmp = FcstUtilities::string_to_map<unsigned int, std::string>( param.get("Impermeable walls") );

                    std::map<unsigned int, std::string>::const_iterator iter;

                    for( iter = tmp.begin(); iter != tmp.end(); ++iter )
                    {
                        BoundaryType ImpermeableWalls;
                        ImpermeableWalls.boundary_name      = "impermeable walls";
                        ImpermeableWalls.boundary_id        =  iter->first;
                        ImpermeableWalls.boundary_condition =  iter->second;

                        this->boundary_types.push_back(ImpermeableWalls);
                    }
                }

                if( !param.get("Symmetry line or plane").empty() )
                {
                    const std::map<unsigned int, std::string> tmp = FcstUtilities::string_to_map<unsigned int, std::string>( param.get("Symmetry line or plane") );

                    std::map<unsigned int, std::string>::const_iterator iter;

                    for( iter = tmp.begin(); iter != tmp.end(); ++iter )
                    {
                        BoundaryType SymmetryLineOrPlane;
                        SymmetryLineOrPlane.boundary_name      = "symmetry line or plane";
                        SymmetryLineOrPlane.boundary_id        =  iter->first;
                        SymmetryLineOrPlane.boundary_condition =  iter->second;

                        this->boundary_types.push_back(SymmetryLineOrPlane);
                    }
                }

                if( !param.get("Inlet-Outlet").empty() )
                {
                    const std::map<unsigned int, std::string> tmp = FcstUtilities::string_to_map<unsigned int, std::string>( param.get("Inlet-Outlet") );

                    std::map<unsigned int, std::string>::const_iterator iter;

                    for( iter = tmp.begin(); iter != tmp.end(); ++iter )
                    {
                        BoundaryType InletOutlet;
                        InletOutlet.boundary_name      = "inlet-outlet";
                        InletOutlet.boundary_id        =  iter->first;
                        InletOutlet.boundary_condition =  iter->second;

                        this->boundary_types.push_back(InletOutlet);
                        if( (InletOutlet.boundary_condition == "Dirichlet velocity") || (InletOutlet.boundary_condition == "Dirichlet pressure") )
                            inlet_outlet_boundary_ID = InletOutlet.boundary_id;
                    }
                }
            }
            param.leave_subsection();

            param.enter_subsection("Boundary data");
            {
                this->variable_boundary_data = param.get_bool("Variable boundary data");

                if( !param.get("single_fluid_pressure").empty() )
                {
                    const std::map<types::boundary_id, double> tmp = FcstUtilities::string_to_map<types::boundary_id, double>( param.get("single_fluid_pressure") );
                    this->component_boundaryID_value["single_fluid_pressure"] = tmp;
                }

                if( !param.get("single_fluid_velocity_X").empty() )
                {
                    const std::map<types::boundary_id, double> tmp = FcstUtilities::string_to_map<types::boundary_id, double>( param.get("single_fluid_velocity_X") );
                    this->component_boundaryID_value["single_fluid_velocity_X"] = tmp;
                }

                if( !param.get("single_fluid_velocity_Y").empty() )
                {
                    const std::map<types::boundary_id, double> tmp = FcstUtilities::string_to_map<types::boundary_id, double>( param.get("single_fluid_velocity_Y") );
                    this->component_boundaryID_value["single_fluid_velocity_Y"] = tmp;
                }

                if( !param.get("single_fluid_velocity_Z").empty() )
                {
                    const std::map<types::boundary_id, double> tmp = FcstUtilities::string_to_map<types::boundary_id, double>( param.get("single_fluid_velocity_Z") );
                    this->component_boundaryID_value["single_fluid_velocity_Z"] = tmp;
                }
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

///////////////////////////////////
///////////////////////////////////
// LOCAL CG FEM BASED ASSEMBLERS //
///////////////////////////////////
///////////////////////////////////

// ---                      ---
// --- assemble_cell_matrix ---
// ---                      ---

template<int dim>
void
NAME::IncompressibleSingleComponentNSEquations<dim>::assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                                                          const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                          FuelCellShop::Layer::BaseLayer<dim>* const                               layer)
{
    // --- types info ---

    const std::type_info& Channel                 = typeid(FuelCellShop::Layer::Channel<dim>);
    const std::type_info& ExperimentalPorousLayer = typeid(FuelCellShop::Layer::ExperimentalPorousLayer<dim>);

    const std::type_info& info = layer->get_base_type();

    // --- fluid properties ---

    double rho_f;
    double dvis_f;

    if     ( info == Channel )
    {
        try
        {
            FuelCellShop::Layer::Channel<dim>* ptr = dynamic_cast< FuelCellShop::Layer::Channel<dim>* >(layer);

            rho_f  = ptr->get_fluid()->get_density();
            dvis_f = ptr->get_fluid()->get_dynamic_viscosity();
        }
        catch(const std::bad_cast& e)
        {
            FcstUtilities::log << "Object is not of type FuelCellShop::Layer::Channel<dim>" << std::endl;
            FcstUtilities::log << e.what() << std::endl;
        }
    }
    else if( info == ExperimentalPorousLayer )
    {
        try
        {
            FuelCellShop::Layer::ExperimentalPorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::ExperimentalPorousLayer<dim>* >(layer);

            rho_f  = ptr->get_fluid()->get_density();
            dvis_f = ptr->get_fluid()->get_dynamic_viscosity();
        }
        catch(const std::bad_cast& e)
        {
            FcstUtilities::log << "Object is not of type FuelCellShop::Layer::ExperimentalPorousLayer<dim>" << std::endl;
            FcstUtilities::log << e.what() << std::endl;
        }
    }
    else
    {
        FcstUtilities::log << "Layer you specified is not implemented" << std::endl;
        AssertThrow( false , ExcNotImplemented() );
    }

    // --- porous media properties ---

    std::vector<double>                   porosity_list( cell_info.get_fe_val_unsplit().n_quadrature_points );
    std::vector< SymmetricTensor<2,dim> > permeability_inv_list( cell_info.get_fe_val_unsplit().n_quadrature_points );
    std::vector< SymmetricTensor<2,dim> > Forchheimer_permeability_list( cell_info.get_fe_val_unsplit().n_quadrature_points );

    if     ( info == Channel )
    {
        for(unsigned int q = 0; q < cell_info.get_fe_val_unsplit().n_quadrature_points; ++q)
        {
            porosity_list[q] = 1.0;
        }
    }
    else if( info == ExperimentalPorousLayer )
    {
        try
        {
            FuelCellShop::Layer::ExperimentalPorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::ExperimentalPorousLayer<dim>* >(layer);

            if( ptr->get_porosity_is_constant() )
            {
                ptr->get_porosity(porosity_list);
            }
            else
            {
                ptr->get_porosity(porosity_list, cell_info.get_fe_val_unsplit().get_quadrature_points());
            }

            if( ptr->get_permeability_is_constant() )
            {
                ptr->get_permeability_INV(permeability_inv_list);
                ptr->get_Forchheimer_permeability(Forchheimer_permeability_list);
            }
            else
            {
                ptr->get_permeability_INV(permeability_inv_list, cell_info.get_fe_val_unsplit().get_quadrature_points());
                ptr->get_Forchheimer_permeability(Forchheimer_permeability_list, cell_info.get_fe_val_unsplit().get_quadrature_points());
            }
        }
        catch(const std::bad_cast& e)
        {
            FcstUtilities::log << "Object is not of type FuelCellShop::Layer::ExperimentalPorousLayer<dim>" << std::endl;
            FcstUtilities::log << e.what() << std::endl;
        }
    }
    else
    {
        FcstUtilities::log << "Layer you specified is not implemented" << std::endl;
        AssertThrow( false , ExcNotImplemented() );
    }

    // --- STEP1 --- --- STEP1 --- --- STEP1 ---

    // --- extractors ---

    const unsigned int pressure_index       = this->system_management->solution_name_to_index("single_fluid_pressure");
    const unsigned int velocity_first_index = this->system_management->solution_name_to_index("single_fluid_velocity_X");

    const FEValuesExtractors::Scalar pressure(pressure_index);
    const FEValuesExtractors::Vector velocities(velocity_first_index);

    // --- fe values ---

    const FEValuesBase<dim>& fe_values_cell = cell_info.get_fe_val_unsplit();

    // --- shortcuts ---

    const unsigned int n_q_points_cell = fe_values_cell.n_quadrature_points;
    const unsigned int dofs_per_cell   = fe_values_cell.dofs_per_cell;

    // --- locals ---

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    
    for(unsigned int q = 0; q < n_q_points_cell; ++q)
    {
        permeability_inv_list[q]         *= Units::convert(1, Units::PER_UNIT2, Units::PER_C_UNIT2); // Use 1/m^2 instead of 1/cm^2
        Forchheimer_permeability_list[q] *= Units::convert(1, Units::PER_UNIT, Units::PER_C_UNIT); // Use 1/m instead of 1/cm
    }

    // --- STEP2 --- --- STEP2 --- --- STEP2 ---

    // --- olds ---

    std::vector< Tensor<1,dim> > u_old(n_q_points_cell);
    fe_values_cell[velocities].get_function_values( cell_info.global_data->vector( cell_info.global_data->n_vectors()-1 ),
                                                    u_old );

    // --- STEP3 --- --- STEP3 --- --- STEP3 ---

    // --- LOOP OVER QUADRATURE POINTS ---

    for(unsigned int q = 0; q < n_q_points_cell; ++q)
    {
        // --- step1 ---

        // --- JxW ---

        const double JxW = fe_values_cell.JxW(q);

        // --- shapes ---

        std::vector<double>          phi_p(dofs_per_cell);
        std::vector< Tensor<1,dim> > grad_phi_p(dofs_per_cell);

        std::vector< Tensor<1,dim> >          phi_u(dofs_per_cell);
        std::vector<double>                   div_phi_u(dofs_per_cell);
        std::vector< SymmetricTensor<2,dim> > grads_phi_u(dofs_per_cell);

        for(unsigned int k = 0; k < dofs_per_cell; ++k)
        {
            phi_p[k]      = fe_values_cell[pressure].value(k,q);
            grad_phi_p[k] = fe_values_cell[pressure].gradient(k,q);

            phi_u[k]       = fe_values_cell[velocities].value(k,q);
            div_phi_u[k]   = fe_values_cell[velocities].divergence(k,q);
            grads_phi_u[k] = fe_values_cell[velocities].symmetric_gradient(k,q);
        }

        // --- shape combinations ---

        std::vector< SymmetricTensor<2,dim> > delta_U(dofs_per_cell);
        std::vector< SymmetricTensor<2,dim> > delta_tau(dofs_per_cell);

        std::vector< Tensor<1,dim> > delta_F(dofs_per_cell);

        if     ( info == Channel )
        {
            if( inertia_in_channels )
            {
                for(unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                    Tensor<2,dim> delta_U_1;
                    outer_product(delta_U_1, phi_u[k], u_old[q]);

                    Tensor<2,dim> delta_U_2;
                    outer_product(delta_U_2, u_old[q], phi_u[k]);

                    delta_U[k] = ( 1.0 / porosity_list[q] ) * rho_f * ( delta_U_1 + delta_U_2 );
                }
            }

            if( shear_stress_in_channels )
            {
                for(unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                    delta_tau[k] = 2.0 * dvis_f * grads_phi_u[k];
                }
            }
        }
        else if( info == ExperimentalPorousLayer )
        {
            if( inertia_in_porous_media )
            {
                for(unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                    Tensor<2,dim> delta_U_1;
                    outer_product(delta_U_1, phi_u[k], u_old[q]);

                    Tensor<2,dim> delta_U_2;
                    outer_product(delta_U_2, u_old[q], phi_u[k]);

                    delta_U[k] = ( 1.0 / porosity_list[q] ) * rho_f * ( delta_U_1 + delta_U_2 );
                }
            }

            if( shear_stress_in_porous_media )
            {
                for(unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                    delta_tau[k] = 2.0 * dvis_f * grads_phi_u[k];
                }
            }

            if( drag_in_porous_media.compare("none") == 0 )
            {
                // do nothing
            }
            else if( drag_in_porous_media.compare("Darcy") == 0 )
            {
                for(unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                    Tensor<1,dim> delta_F_1;
                    delta_F_1 = - dvis_f * permeability_inv_list[q] * phi_u[k];

                    delta_F[k] = delta_F_1;
                }
            }
            else if( drag_in_porous_media.compare("Forchheimer") == 0 )
            {
                for(unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                    Tensor<1,dim> delta_F_1;
                    delta_F_1 = - dvis_f * permeability_inv_list[q] * phi_u[k];

                    double ratio;
                    ratio = ( u_old[q] * phi_u[k] ) / u_old[q].norm();

                    Tensor<1,dim> delta_F_2;
                    delta_F_2 = - rho_f * ratio * Forchheimer_permeability_list[q] * u_old[q];

                    Tensor<1,dim> delta_F_3;
                    delta_F_3 = - rho_f * u_old[q].norm() * Forchheimer_permeability_list[q] * phi_u[k];

                    delta_F[k] = delta_F_1 + delta_F_2 + delta_F_3;
                }
            }
            else if( drag_in_porous_media.compare("Forchheimer modified") == 0 )
            {
                for(unsigned int q = 0; q < n_q_points_cell; ++q)
                {

                }
            }
            else
            {
                FcstUtilities::log << "Form of the drag force you specified in porous media is not implemented" << std::endl;
                AssertThrow( false , ExcNotImplemented() );
            }
        }
        else
        {
            FcstUtilities::log << "Layer you specified is not implemented" << std::endl;
            AssertThrow( false , ExcNotImplemented() );
        }

        // --- step2 ---

        // --- LOOP OVER i ---

        for(unsigned int i = 0; i < dofs_per_cell; ++i)
        {
                // --- LOOP OVER j ---

                for(unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                    local_matrix(i,j) += phi_p[i]       * div_phi_u[j]* JxW
                                         -
                                         grads_phi_u[i] * delta_U[j]* JxW
                                         +
                                         grads_phi_u[i] * ( - phi_p[j]*unit + delta_tau[j] ) * JxW
                                         -
                                         phi_u[i]       * delta_F[j]* JxW;
                }
        }
    }

    // --- conversion ---

    if( this->counter[0] )
    {
        this->make_matrix_block_indices();
        this->counter[0] = false;
    }

    local_matrix *= -1.0;
    this->dealII_to_appframe(cell_matrices,
                             local_matrix,
                             this->matrix_block_indices);
}

// ---                        ---
// --- assemble_cell_residual ---
// ---                        ---

template<int dim>
void
NAME::IncompressibleSingleComponentNSEquations<dim>::assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
                                                                            const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                            FuelCellShop::Layer::BaseLayer<dim>* const                               layer)
{
    // --- types info ---

    const std::type_info& Channel                 = typeid(FuelCellShop::Layer::Channel<dim>);
    const std::type_info& ExperimentalPorousLayer = typeid(FuelCellShop::Layer::ExperimentalPorousLayer<dim>);

    const std::type_info& info = layer->get_base_type();

    // --- fluid properties ---

    double rho_f;
    double dvis_f;

    if( info == Channel )
    {
        try
        {
            FuelCellShop::Layer::Channel<dim>* ptr = dynamic_cast< FuelCellShop::Layer::Channel<dim>* >(layer);

            rho_f  = ptr->get_fluid()->get_density();
            dvis_f = ptr->get_fluid()->get_dynamic_viscosity();
        }
        catch(const std::bad_cast& e)
        {
            FcstUtilities::log << "Object is not of type FuelCellShop::Layer::Channel<dim>" << std::endl;
            FcstUtilities::log << e.what() << std::endl;
        }
    }
    else if( info == ExperimentalPorousLayer )
    {
        try
        {
            FuelCellShop::Layer::ExperimentalPorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::ExperimentalPorousLayer<dim>* >(layer);

            rho_f  = ptr->get_fluid()->get_density();
            dvis_f = ptr->get_fluid()->get_dynamic_viscosity();
        }
        catch(const std::bad_cast& e)
        {
            FcstUtilities::log << "Object is not of type FuelCellShop::Layer::ExperimentalPorousLayer<dim>" << std::endl;
            FcstUtilities::log << e.what() << std::endl;
        }
    }
    else
    {
        FcstUtilities::log << "Layer you specified is not implemented" << std::endl;
        AssertThrow( false , ExcNotImplemented() );
    }

    // --- porous media properties ---

    std::vector<double>                   porosity_list(                 cell_info.get_fe_val_unsplit().n_quadrature_points );
    std::vector< SymmetricTensor<2,dim> > permeability_inv_list(         cell_info.get_fe_val_unsplit().n_quadrature_points );
    std::vector< SymmetricTensor<2,dim> > Forchheimer_permeability_list( cell_info.get_fe_val_unsplit().n_quadrature_points );

    if     ( info == Channel )
    {
        for(unsigned int q = 0; q < cell_info.get_fe_val_unsplit().n_quadrature_points; ++q)
        {
            porosity_list[q] = 1.0;
        }
    }
    else if( info == ExperimentalPorousLayer )
    {
        try
        {
            FuelCellShop::Layer::ExperimentalPorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::ExperimentalPorousLayer<dim>* >(layer);

            if( ptr->get_porosity_is_constant() )
            {
                ptr->get_porosity(porosity_list);
            }
            else
            {
                ptr->get_porosity(porosity_list, cell_info.get_fe_val_unsplit().get_quadrature_points());
            }

            if( ptr->get_permeability_is_constant() )
            {
                ptr->get_permeability_INV(permeability_inv_list);
                ptr->get_Forchheimer_permeability(Forchheimer_permeability_list);
            }
            else
            {
                ptr->get_permeability_INV(permeability_inv_list, cell_info.get_fe_val_unsplit().get_quadrature_points());
                ptr->get_Forchheimer_permeability(Forchheimer_permeability_list, cell_info.get_fe_val_unsplit().get_quadrature_points());
            }
        }
        catch(const std::bad_cast& e)
        {
            FcstUtilities::log << "Object is not of type FuelCellShop::Layer::ExperimentalPorousLayer<dim>" << std::endl;
            FcstUtilities::log << e.what() << std::endl;
        }
    }
    else
    {
        FcstUtilities::log << "Layer you specified is not implemented" << std::endl;
        AssertThrow( false , ExcNotImplemented() );
    }
        
    // --- STEP1 --- --- STEP1 --- --- STEP1 ---

    // --- extractors ---

    const unsigned int pressure_index       = this->system_management->solution_name_to_index("single_fluid_pressure");
    const unsigned int velocity_first_index = this->system_management->solution_name_to_index("single_fluid_velocity_X");

    const FEValuesExtractors::Scalar pressure(pressure_index);
    const FEValuesExtractors::Vector velocities(velocity_first_index);

    // --- fe values ---

    const FEValuesBase<dim>& fe_values_cell = cell_info.get_fe_val_unsplit();

    // --- shortcuts ---

    const unsigned int n_q_points_cell = fe_values_cell.n_quadrature_points;
    const unsigned int dofs_per_cell   = fe_values_cell.dofs_per_cell;

    // --- locals ---

    Vector<double> local_residual(dofs_per_cell);
    
    for(unsigned int q = 0; q < n_q_points_cell; ++q)
    {
        permeability_inv_list[q]         *= Units::convert(1, Units::PER_UNIT2, Units::PER_C_UNIT2); // Use 1/m^2 instead of 1/cm^2
        Forchheimer_permeability_list[q] *= Units::convert(1, Units::PER_UNIT, Units::PER_C_UNIT); // Use 1/m instead of 1/cm
    }

    // --- STEP2 --- --- STEP2 --- --- STEP2 ---

    // --- olds ---

    std::vector<double>          p_old(n_q_points_cell);
    std::vector< Tensor<1,dim> > grad_p_old(n_q_points_cell);

    std::vector< Tensor<1,dim> >          u_old(n_q_points_cell);
    std::vector<double>                   div_u_old(n_q_points_cell);
    std::vector< SymmetricTensor<2,dim> > grads_u_old(n_q_points_cell);

    fe_values_cell[pressure].get_function_values(cell_info.global_data->vector( cell_info.global_data->n_vectors()-1 ), p_old);
    fe_values_cell[pressure].get_function_gradients(cell_info.global_data->vector( cell_info.global_data->n_vectors()-1 ), grad_p_old);

    fe_values_cell[velocities].get_function_values(cell_info.global_data->vector( cell_info.global_data->n_vectors()-1 ), u_old);
    fe_values_cell[velocities].get_function_divergences(cell_info.global_data->vector( cell_info.global_data->n_vectors()-1 ), div_u_old);
    fe_values_cell[velocities].get_function_symmetric_gradients(cell_info.global_data->vector( cell_info.global_data->n_vectors()-1 ), grads_u_old );

    // --- old combinations ---

    std::vector< SymmetricTensor<2,dim> > u_old_BY_u_old(n_q_points_cell);
    std::vector< SymmetricTensor<2,dim> > tau_old(n_q_points_cell);

    std::vector< Tensor<1,dim> > F_old(n_q_points_cell);
    std::vector< Tensor<1,dim> > G_old(n_q_points_cell);

    if( info == Channel )
    {
        if( inertia_in_channels )
        {
            for(unsigned int q = 0; q < n_q_points_cell; ++q)
            {
                Tensor<2,dim> tmp;
                outer_product(tmp, u_old[q], u_old[q]);

                u_old_BY_u_old[q] = ( 1.0 / porosity_list[q] ) * rho_f * tmp;
            }
        }

        if( shear_stress_in_channels )
        {
            for(unsigned int q = 0; q < n_q_points_cell; ++q)
            {
                tau_old[q] = 2.0 * dvis_f * grads_u_old[q];
            }
        }

        if( gravity_in_channels )
        {
            for(unsigned int q = 0; q < n_q_points_cell; ++q)
            {
                G_old[q] = porosity_list[q] * rho_f * gravity_acceleration;
            }
        }
    }
    else if( info == ExperimentalPorousLayer )
    {
        if( inertia_in_porous_media )
        {
            for(unsigned int q = 0; q < n_q_points_cell; ++q)
            {
                Tensor<2,dim> tmp;
                outer_product(tmp, u_old[q], u_old[q]);

                u_old_BY_u_old[q] = ( 1.0 / porosity_list[q] ) * rho_f * tmp;
            }
        }

        if( shear_stress_in_porous_media )
        {
            for(unsigned int q = 0; q < n_q_points_cell; ++q)
            {
                tau_old[q] = 2.0 * dvis_f * grads_u_old[q];
            }
        }

        if( gravity_in_porous_media )
        {
            for(unsigned int q = 0; q < n_q_points_cell; ++q)
            {
                G_old[q] = porosity_list[q] * rho_f * gravity_acceleration;
            }
        }

        if(      drag_in_porous_media.compare("none") == 0 )
        {
                // do nothing
        }
        else if( drag_in_porous_media.compare("Darcy") == 0 )
        {
            for(unsigned int q = 0; q < n_q_points_cell; ++q)
            {
                Tensor<1,dim> F_old_1;
                F_old_1 = - dvis_f * permeability_inv_list[q] * u_old[q];

                F_old[q] = F_old_1;
            }
        }
        else if( drag_in_porous_media.compare("Forchheimer") == 0 )
        {
            for(unsigned int q = 0; q < n_q_points_cell; ++q)
            {
                Tensor<1,dim> F_old_1;
                F_old_1 = - dvis_f * permeability_inv_list[q] * u_old[q];

                Tensor<1,dim> F_old_2;
                F_old_2 = - rho_f * u_old[q].norm() * Forchheimer_permeability_list[q] * u_old[q];

                F_old[q] = F_old_1 + F_old_2;
            }
        }
        else if( drag_in_porous_media.compare("Forchheimer modified") == 0 )
        {
            for(unsigned int q = 0; q < n_q_points_cell; ++q)
            {

            }
        }
        else
        {
            FcstUtilities::log << "Form of the drag force you specified in porous media is not implemented" << std::endl;
            AssertThrow( false , ExcNotImplemented() );
        }
    }
    else
    {
        FcstUtilities::log << "Layer you specified is not implemented" << std::endl;
        AssertThrow( false , ExcNotImplemented() );
    }

    // --- STEP3 --- --- STEP3 --- --- STEP3 ---

    // --- LOOP OVER QUADRATURE POINTS ---

    for(unsigned int q = 0; q < n_q_points_cell; ++q)
    {
        // --- step1 ---

        // --- JxW ---

        const double JxW = fe_values_cell.JxW(q);

        // --- shapes ---

        std::vector<double>          phi_p(dofs_per_cell);
        std::vector< Tensor<1,dim> > grad_phi_p(dofs_per_cell);

        std::vector< Tensor<1,dim> >          phi_u(dofs_per_cell);
        std::vector<double>                   div_phi_u(dofs_per_cell);
        std::vector< SymmetricTensor<2,dim> > grads_phi_u(dofs_per_cell);

        for(unsigned int k = 0; k < dofs_per_cell; ++k)
        {
            phi_p[k]      = fe_values_cell[pressure].value(k,q);
            grad_phi_p[k] = fe_values_cell[pressure].gradient(k,q);

            phi_u[k]       = fe_values_cell[velocities].value(k,q);
            div_phi_u[k]   = fe_values_cell[velocities].divergence(k,q);
            grads_phi_u[k] = fe_values_cell[velocities].symmetric_gradient(k,q);
        }

        // --- step2 ---

        // --- LOOP OVER i ---

        for(unsigned int i = 0; i < dofs_per_cell; ++i)
        {
            local_residual(i) += - phi_p[i]       *     div_u_old[q]                 * JxW
                                +
                                grads_phi_u[i] *     u_old_BY_u_old[q]            * JxW
                                -
                                grads_phi_u[i] * ( - p_old[q]*unit + tau_old[q] ) * JxW
                                +
                                phi_u[i]       *     F_old[q]                     * JxW
                                +
                                phi_u[i]       *     G_old[q]                     * JxW;
        }
    }

    // --- conversion ---

    if( this->counter[1] )
    {
        this->make_residual_indices();
        this->counter[1] = false;
    }

    this->dealII_to_appframe(cell_residual,
                             local_residual,
                             this->residual_indices);
}

// ---                      ---
// --- assemble_bdry_matrix ---
// ---                      ---

template<int dim>
void
NAME::IncompressibleSingleComponentNSEquations<dim>::assemble_bdry_matrix(FuelCell::ApplicationCore::MatrixVector&                                 bdry_matrices,
                                                                          const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                                          FuelCellShop::Layer::BaseLayer<dim>* const                               layer)
{
    // --- types info ---

    const std::type_info& Channel                 = typeid(FuelCellShop::Layer::Channel<dim>);
    const std::type_info& ExperimentalPorousLayer = typeid(FuelCellShop::Layer::ExperimentalPorousLayer<dim>);

    const std::type_info& info = layer->get_base_type();

    // --- fluid properties ---

    double rho_f;
    double dvis_f;

    if( info == Channel )
    {
        try
        {
            FuelCellShop::Layer::Channel<dim>* ptr = dynamic_cast< FuelCellShop::Layer::Channel<dim>* >(layer);

            rho_f  = ptr->get_fluid()->get_density();
            dvis_f = ptr->get_fluid()->get_dynamic_viscosity();
        }
        catch(const std::bad_cast& e)
        {
            FcstUtilities::log << "Object is not of type FuelCellShop::Layer::Channel<dim>" << std::endl;
            FcstUtilities::log << e.what() << std::endl;
        }
    }
    else if( info == ExperimentalPorousLayer )
    {
        try
        {
            FuelCellShop::Layer::ExperimentalPorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::ExperimentalPorousLayer<dim>* >(layer);

            rho_f  = ptr->get_fluid()->get_density();
            dvis_f = ptr->get_fluid()->get_dynamic_viscosity();
        }
        catch(const std::bad_cast& e)
        {
            FcstUtilities::log << "Object is not of type FuelCellShop::Layer::ExperimentalPorousLayer<dim>" << std::endl;
            FcstUtilities::log << e.what() << std::endl;
        }
    }
    else
    {
        FcstUtilities::log << "Layer you specified is not implemented" << std::endl;
        AssertThrow( false , ExcNotImplemented() );
    }

    // --- porous media properties ---

    std::vector<double> porosity_list( bdry_info.get_fe_val_unsplit().n_quadrature_points );

    if( info == Channel )
    {
        for(unsigned int q = 0; q < bdry_info.get_fe_val_unsplit().n_quadrature_points; ++q)
        {
            porosity_list[q] = 1.0;
        }
    }
    else if( info == ExperimentalPorousLayer )
    {
        try
        {
            FuelCellShop::Layer::ExperimentalPorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::ExperimentalPorousLayer<dim>* >(layer);

            if( ptr->get_porosity_is_constant() )
            {
                ptr->get_porosity(porosity_list);
            }
            else
            {
                ptr->get_porosity(porosity_list,
                                    bdry_info.get_fe_val_unsplit().get_quadrature_points());
            }
        }
        catch(const std::bad_cast& e)
        {
            FcstUtilities::log << "Object is not of type FuelCellShop::Layer::ExperimentalPorousLayer<dim>" << std::endl;
            FcstUtilities::log << e.what() << std::endl;
        }
    }
    else
    {
        FcstUtilities::log << "Layer you specified is not implemented" << std::endl;
        AssertThrow( false , ExcNotImplemented() );
    }

    // --- STEP1 --- --- STEP1 --- --- STEP1 ---

    // --- extractors ---

    const unsigned int pressure_index       = this->system_management->solution_name_to_index("single_fluid_pressure");
    const unsigned int velocity_first_index = this->system_management->solution_name_to_index("single_fluid_velocity_X");

    const FEValuesExtractors::Scalar pressure(pressure_index);
    const FEValuesExtractors::Vector velocities(velocity_first_index);

    // --- fe values ---

    const FEFaceValuesBase<dim>& fe_values_bdry = bdry_info.get_fe_val_unsplit();

    // --- shortcuts ---

    const unsigned int n_q_points_bdry = fe_values_bdry.n_quadrature_points;
    const unsigned int dofs_per_cell   = fe_values_bdry.dofs_per_cell;

    // --- normals ---

    std::vector< Point<dim> > normal_vectors(n_q_points_bdry);
    normal_vectors = fe_values_bdry.get_normal_vectors();

    // --- tangentials ---

    std::vector< std::vector< Point<dim> > > tangential_vectors;

    tangential_vectors.resize( dim-1, std::vector< Point<dim> >( n_q_points_bdry ));

    FemExtras::get_tangential_vectors(tangential_vectors, normal_vectors);

    // --- locals ---

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);

    // --- extras ---

    double normal_velocity_suppressed_weakly = 0.0;
    if( normal_velocity_is_suppressed_weakly )
        normal_velocity_suppressed_weakly = 1.0;

    // --- STEP2 --- --- STEP2 --- --- STEP2 ---

    // --- olds ---

    std::vector< Tensor<1,dim> > u_old(n_q_points_bdry);
    fe_values_bdry[velocities].get_function_values( bdry_info.global_data->vector( bdry_info.global_data->n_vectors()-1 ), u_old );

    // --- STEP3 --- --- STEP3 --- --- STEP3 ---

    // --- integral multipliers ---

    double int_1    = 1.0;
    double int_2    = 1.0;
    double p_factor = 1.0;

    // --- LOOP OVER BOUNDARIES ---

    std::vector< BoundaryType >::const_iterator iter;

    for( iter = this->boundary_types.begin(); iter != this->boundary_types.end(); ++iter )
    {
        if( this->belongs_to_boundary( bdry_info.dof_face->boundary_indicator(), iter->boundary_id ) )
        {
            // --- IMPERMEABLE WALLS ---

            if( iter->boundary_name.compare("impermeable walls") == 0 )
            {
                // --- NO-SLIP ---

                if( iter->boundary_condition.compare("no-slip") == 0 )
                {
                        // do nothing
                }

                // --- NAVIER SLIP ---

                else if( iter->boundary_condition.compare("Navier slip") == 0 )
                {
                    // --- LOOP OVER QUADRATURE POINTS ---

                    for(unsigned int q = 0; q < n_q_points_bdry; ++q)
                    {
                        // --- step1 ---

                        // --- JxW ---

                        const double JxW = fe_values_bdry.JxW(q);

                        // --- shapes ---

                        std::vector< Tensor<1,dim> > phi_u(dofs_per_cell);

                        for(unsigned int k = 0; k < dofs_per_cell; ++k)
                        {
                            phi_u[k] = fe_values_bdry[velocities].value(k,q);
                        }

                        // --- step2 ---

                        // --- LOOP OVER i ---

                        for(unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                            // --- LOOP OVER j ---

                            for(unsigned int j = 0; j < dofs_per_cell; ++j)
                            {
                                    double sum = 0.0;

                                    // --- LOOP OVER alpha ---

                                    for(unsigned int alpha = 0; alpha < dim - 1; ++alpha)
                                    {
                                            sum += ( phi_u[i] * tangential_vectors[alpha][q] ) * ( phi_u[j] * tangential_vectors[alpha][q] );
                                    }

                                    sum *= (1.0 - theta) / theta;

                                    local_matrix(i,j) += sum * JxW
                                                         +
                                                         normal_velocity_suppressed_weakly * (1.0/eta) * ( phi_u[i] * normal_vectors[q] ) * ( phi_u[j] * normal_vectors[q] ) * JxW;
                            }
                        }
                    }
                }

                // --- PERFECT SLIP ---

                else if( iter->boundary_condition.compare("perfect slip") == 0 )
                {
                    // --- LOOP OVER QUADRATURE POINTS ---

                    for(unsigned int q = 0; q < n_q_points_bdry; ++q)
                    {
                        // --- step1 ---

                        // --- JxW ---

                        const double JxW = fe_values_bdry.JxW(q);

                        // --- shapes ---

                        std::vector< Tensor<1,dim> > phi_u(dofs_per_cell);

                        for(unsigned int k = 0; k < dofs_per_cell; ++k)
                        {
                            phi_u[k] = fe_values_bdry[velocities].value(k,q);
                        }

                        // --- step2 ---

                        // --- LOOP OVER i ---

                        for(unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                            // --- LOOP OVER j ---

                            for(unsigned int j = 0; j < dofs_per_cell; ++j)
                            {
                                local_matrix(i,j) += normal_velocity_suppressed_weakly * (1.0/eta) * ( phi_u[i] * normal_vectors[q] ) * ( phi_u[j] * normal_vectors[q] ) * JxW;
                            }
                        }
                    }
                }

                // --- UNKNOWN ---

                else
                {
                        FcstUtilities::log << "Boundary condition you specified at \"impermeable walls\" is not implemented" << std::endl;
                        AssertThrow( false , ExcNotImplemented() );
                }
            }

            // --- SYMMETRY LINE OR PLANE ---

            else if( iter->boundary_name.compare("symmetry line or plane") == 0 )
            {
                // --- PERFECT SLIP ---

                if( iter->boundary_condition.compare("perfect slip") == 0 )
                {
                    // --- LOOP OVER QUADRATURE POINTS ---

                    for(unsigned int q = 0; q < n_q_points_bdry; ++q)
                    {
                        // --- step1 ---

                        // --- JxW ---

                        const double JxW = fe_values_bdry.JxW(q);

                        // --- shapes ---

                        std::vector< Tensor<1,dim> > phi_u(dofs_per_cell);

                        for(unsigned int k = 0; k < dofs_per_cell; ++k)
                        {
                            phi_u[k] = fe_values_bdry[velocities].value(k,q);
                        }

                        // --- step2 ---

                        // --- LOOP OVER i ---

                        for(unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                            // --- LOOP OVER j ---

                            for(unsigned int j = 0; j < dofs_per_cell; ++j)
                            {
                                    local_matrix(i,j) += normal_velocity_suppressed_weakly * (1.0/eta) * ( phi_u[i] * normal_vectors[q] ) * ( phi_u[j] * normal_vectors[q] ) * JxW;
                            }
                        }
                    }
                }

                // --- UNKNOWN ---

                else
                {
                        FcstUtilities::log << "Boundary condition you specified at \"symmetry line or plane\" is not implemented" << std::endl;
                        AssertThrow( false , ExcNotImplemented() );
                }
            }

            // --- INLET-OUTLET ---

            else if( iter->boundary_name.compare("inlet-outlet") == 0 )
            {
                // --- DIRICHLET PRESSURE ---

                if( iter->boundary_condition.compare("Dirichlet pressure") == 0 )
                {
                    p_factor = 0.0;
                }

                // --- DIRICHLET PRESSURE AND NORMAL STRESS FREE ---

                else if( iter->boundary_condition.compare("Dirichlet pressure and normal stress free") == 0 )
                {
                    int_2 = 0.0;
                }

                // --- DIRICHLET VELOCITY ---

                else if( iter->boundary_condition.compare("Dirichlet velocity") == 0 )
                {
                    int_1 = 0.0;
                    int_2 = 0.0;
                }

                // --- DIRICHLET PRESSURE AND VELOCITY ---

                else if( iter->boundary_condition.compare("Dirichlet pressure and velocity") == 0 )
                {
                    int_1 = 0.0;
                    int_2 = 0.0;
                }

                // --- NORMAL STRESS FREE ---

                else if( iter->boundary_condition.compare("normal stress free") == 0 )
                {
                    int_2 = 0.0;
                }

                // --- UNKNOWN ---

                else
                {
                    FcstUtilities::log << "Boundary condition you specified at \"inlet-outlet\" is not implemented" << std::endl;
                    AssertThrow( false , ExcNotImplemented() );
                }

                // --- LOOP OVER QUADRATURE POINTS ---

                for(unsigned int q = 0; q < n_q_points_bdry; ++q)
                {
                    // --- step1 ---

                    // --- JxW ---

                    const double JxW = fe_values_bdry.JxW(q);

                    // --- shapes ---

                    std::vector<double> phi_p(dofs_per_cell);

                    std::vector< Tensor<1,dim> >          phi_u(dofs_per_cell);
                    std::vector< SymmetricTensor<2,dim> > grads_phi_u(dofs_per_cell);

                    for(unsigned int k = 0; k < dofs_per_cell; ++k)
                    {
                        phi_p[k]       = fe_values_bdry[pressure].value(k,q);

                        phi_u[k]       = fe_values_bdry[velocities].value(k,q);
                        grads_phi_u[k] = fe_values_bdry[velocities].symmetric_gradient(k,q);
                    }

                    // --- shape combinations ---

                    std::vector< SymmetricTensor<2,dim> > n_BY_phi_u_S(dofs_per_cell);

                    std::vector< SymmetricTensor<2,dim> > delta_U(dofs_per_cell);
                    std::vector< SymmetricTensor<2,dim> > delta_tau(dofs_per_cell);

                    for(unsigned int k = 0; k < dofs_per_cell; ++k)
                    {
                        Tensor<2,dim> n_BY_phi_u;
                        outer_product(n_BY_phi_u, normal_vectors[q], phi_u[k]);
                        n_BY_phi_u_S[k] = 0.5 * ( n_BY_phi_u + transpose(n_BY_phi_u) );
                    }

                    if( info == Channel )
                    {
                        if( inertia_in_channels )
                        {
                            for(unsigned int k = 0; k < dofs_per_cell; ++k)
                            {
                                Tensor<2,dim> delta_U_1;
                                outer_product(delta_U_1, phi_u[k], u_old[q]);

                                Tensor<2,dim> delta_U_2;
                                outer_product(delta_U_2, u_old[q], phi_u[k]);

                                delta_U[k] = ( 1.0 / porosity_list[q] ) * rho_f * ( delta_U_1 + delta_U_2 );
                            }
                        }

                        if( shear_stress_in_channels )
                        {
                            for(unsigned int k = 0; k < dofs_per_cell; ++k)
                            {
                                delta_tau[k] = 2.0 * dvis_f * grads_phi_u[k];
                            }
                        }
                    }
                    else if( info == ExperimentalPorousLayer )
                    {
                        if( inertia_in_porous_media )
                        {
                            for(unsigned int k = 0; k < dofs_per_cell; ++k)
                            {
                                    Tensor<2,dim> delta_U_1;
                                    outer_product(delta_U_1, phi_u[k], u_old[q]);

                                    Tensor<2,dim> delta_U_2;
                                    outer_product(delta_U_2, u_old[q], phi_u[k]);

                                    delta_U[k] = ( 1.0 / porosity_list[q] ) * rho_f * ( delta_U_1 + delta_U_2 );
                            }
                        }

                        if( shear_stress_in_porous_media )
                        {
                            for(unsigned int k = 0; k < dofs_per_cell; ++k)
                            {
                                    delta_tau[k] = 2.0 * dvis_f * grads_phi_u[k];
                            }
                        }
                    }
                    else
                    {
                        FcstUtilities::log << "Layer you specified is not implemented" << std::endl;
                        AssertThrow( false , ExcNotImplemented() );
                    }

                    // --- step2 ---

                    // --- LOOP OVER i ---

                    for(unsigned int i = 0; i < dofs_per_cell; ++i)
                    {
                        // --- LOOP OVER j ---

                        for(unsigned int j = 0; j < dofs_per_cell; ++j)
                        {
                            local_matrix(i,j) += int_1 * n_BY_phi_u_S[i] * delta_U[j] * JxW
                                                    -
                                                    int_2 * n_BY_phi_u_S[i] * ( - p_factor*phi_p[j]*unit + delta_tau[j] ) * JxW;
                        }
                    }
                }
            }
        }
    }

    // --- conversion ---

    if( this->counter[0] )
    {
        this->make_matrix_block_indices();
        this->counter[0] = false;
    }

    local_matrix *= -1.0;
    this->dealII_to_appframe(bdry_matrices,
                             local_matrix,
                             this->matrix_block_indices);
}

// ---                        ---
// --- assemble_bdry_residual ---
// ---                        ---

template<int dim>
void
NAME::IncompressibleSingleComponentNSEquations<dim>::assemble_bdry_residual(FuelCell::ApplicationCore::FEVector&                                     bdry_residual,
                                                                            const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                                            FuelCellShop::Layer::BaseLayer<dim>* const                               layer)
{
    // --- types info ---

    const std::type_info& Channel                 = typeid(FuelCellShop::Layer::Channel<dim>);
    const std::type_info& ExperimentalPorousLayer = typeid(FuelCellShop::Layer::ExperimentalPorousLayer<dim>);

    const std::type_info& info = layer->get_base_type();

    // --- fluid properties ---

    double rho_f;
    double dvis_f;

    if     ( info == Channel )
    {
        try
        {
            FuelCellShop::Layer::Channel<dim>* ptr = dynamic_cast< FuelCellShop::Layer::Channel<dim>* >(layer);

            rho_f  = ptr->get_fluid()->get_density();
            dvis_f = ptr->get_fluid()->get_dynamic_viscosity();
        }
        catch(const std::bad_cast& e)
        {
            FcstUtilities::log << "Object is not of type FuelCellShop::Layer::Channel<dim>" << std::endl;
            FcstUtilities::log << e.what() << std::endl;
        }
    }
    else if( info == ExperimentalPorousLayer )
    {
        try
        {
            FuelCellShop::Layer::ExperimentalPorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::ExperimentalPorousLayer<dim>* >(layer);

            rho_f  = ptr->get_fluid()->get_density();
            dvis_f = ptr->get_fluid()->get_dynamic_viscosity();
        }
        catch(const std::bad_cast& e)
        {
            FcstUtilities::log << "Object is not of type FuelCellShop::Layer::ExperimentalPorousLayer<dim>" << std::endl;
            FcstUtilities::log << e.what() << std::endl;
        }
    }
    else
    {
        FcstUtilities::log << "Layer you specified is not implemented" << std::endl;
        AssertThrow( false , ExcNotImplemented() );
    }

    // --- porous media properties ---

    std::vector<double> porosity_list( bdry_info.get_fe_val_unsplit().n_quadrature_points );

    if( info == Channel )
    {
        for(unsigned int q = 0; q < bdry_info.get_fe_val_unsplit().n_quadrature_points; ++q)
        {
            porosity_list[q] = 1.0;
        }
    }
    else if( info == ExperimentalPorousLayer )
    {
        try
        {
            FuelCellShop::Layer::ExperimentalPorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::ExperimentalPorousLayer<dim>* >(layer);

            if( ptr->get_porosity_is_constant() )
            {
                ptr->get_porosity(porosity_list);
            }
            else
            {
                ptr->get_porosity(porosity_list, bdry_info.get_fe_val_unsplit().get_quadrature_points());
            }
        }
        catch(const std::bad_cast& e)
        {
            FcstUtilities::log << "Object is not of type FuelCellShop::Layer::ExperimentalPorousLayer<dim>" << std::endl;
            FcstUtilities::log << e.what() << std::endl;
        }
    }
    else
    {
        FcstUtilities::log << "Layer you specified is not implemented" << std::endl;
        AssertThrow( false , ExcNotImplemented() );
    }

    // --- STEP1 --- --- STEP1 --- --- STEP1 ---

    // --- extractors ---

    const unsigned int pressure_index       = this->system_management->solution_name_to_index("single_fluid_pressure");
    const unsigned int velocity_first_index = this->system_management->solution_name_to_index("single_fluid_velocity_X");

    const FEValuesExtractors::Scalar pressure(pressure_index);
    const FEValuesExtractors::Vector velocities(velocity_first_index);

    // --- fe values ---

    const FEFaceValuesBase<dim>& fe_values_bdry = bdry_info.get_fe_val_unsplit();

    // --- shortcuts ---

    const unsigned int n_q_points_bdry = fe_values_bdry.n_quadrature_points;
    const unsigned int dofs_per_cell   = fe_values_bdry.dofs_per_cell;

    // --- normals ---

    std::vector< Point<dim> > normal_vectors(n_q_points_bdry);
    normal_vectors = fe_values_bdry.get_normal_vectors();

    // --- tangentials ---

    std::vector< std::vector< Point<dim> > > tangential_vectors;

    tangential_vectors.resize( dim-1, std::vector< Point<dim> >( n_q_points_bdry ));

    FemExtras::get_tangential_vectors(tangential_vectors, normal_vectors);

    // --- locals ---

    Vector<double> local_residual(dofs_per_cell);

    // --- extras ---

    double normal_velocity_suppressed_weakly = 0.0;
    if( normal_velocity_is_suppressed_weakly )
        normal_velocity_suppressed_weakly = 1.0;

    // --- STEP2 --- --- STEP2 --- --- STEP2 ---

    // --- olds ---

    std::vector<double> p_old(n_q_points_bdry);

    std::vector< Tensor<1,dim> >          u_old(n_q_points_bdry);
    std::vector< SymmetricTensor<2,dim> > grads_u_old(n_q_points_bdry);

    fe_values_bdry[pressure].get_function_values(bdry_info.global_data->vector( bdry_info.global_data->n_vectors()-1 ), p_old);

    fe_values_bdry[velocities].get_function_values(bdry_info.global_data->vector( bdry_info.global_data->n_vectors()-1 ), u_old);
    fe_values_bdry[velocities].get_function_symmetric_gradients(bdry_info.global_data->vector( bdry_info.global_data->n_vectors()-1 ), grads_u_old);

    // --- old combinations ---

    std::vector< SymmetricTensor<2,dim> > u_old_BY_u_old(n_q_points_bdry);
    std::vector< SymmetricTensor<2,dim> > tau_old(n_q_points_bdry);

    if( info == Channel )
    {
        if( inertia_in_channels )
        {
            for(unsigned int q = 0; q < n_q_points_bdry; ++q)
            {
                Tensor<2,dim> tmp;
                outer_product(tmp, u_old[q], u_old[q]);

                u_old_BY_u_old[q] = ( 1.0 / porosity_list[q] ) * rho_f * tmp;
            }
        }

        if( shear_stress_in_channels )
        {
            for(unsigned int q = 0; q < n_q_points_bdry; ++q)
            {
                tau_old[q] = 2.0 * dvis_f * grads_u_old[q];
            }
        }
    }
    else if( info == ExperimentalPorousLayer )
    {
        if( inertia_in_porous_media )
        {
            for(unsigned int q = 0; q < n_q_points_bdry; ++q)
            {
                Tensor<2,dim> tmp;
                outer_product(tmp, u_old[q], u_old[q]);

                u_old_BY_u_old[q] = ( 1.0 / porosity_list[q] ) * rho_f * tmp;
            }
        }

        if( shear_stress_in_porous_media )
        {
            for(unsigned int q = 0; q < n_q_points_bdry; ++q)
            {
                tau_old[q] = 2.0 * dvis_f * grads_u_old[q];
            }
        }
    }
    else
    {
        FcstUtilities::log << "Layer you specified is not implemented" << std::endl;
        AssertThrow( false , ExcNotImplemented() );
    }

    // --- STEP3 --- --- STEP3 --- --- STEP3 ---

    // --- integral multipliers ---

    double int_1 = 1.0;
    double int_2 = 1.0;

    // --- LOOP OVER BOUNDARIES ---

    std::vector< BoundaryType >::const_iterator iter;

    for( iter = this->boundary_types.begin(); iter != this->boundary_types.end(); ++iter )
    {
        if( this->belongs_to_boundary( bdry_info.dof_face->boundary_indicator(),
                                        iter->boundary_id ) )
        {
            // --- IMPERMEABLE WALLS ---

            if( iter->boundary_name.compare("impermeable walls") == 0 )
            {
                // --- NO-SLIP ---

                if( iter->boundary_condition.compare("no-slip") == 0 )
                {
                        // do nothing
                }

                // --- NAVIER SLIP ---

                else if( iter->boundary_condition.compare("Navier slip") == 0 )
                {
                    // --- LOOP OVER QUADRATURE POINTS ---

                    for(unsigned int q = 0; q < n_q_points_bdry; ++q)
                    {
                        // --- step1 ---

                        // --- JxW ---

                        const double JxW = fe_values_bdry.JxW(q);

                        // --- shapes ---

                        std::vector< Tensor<1,dim> > phi_u(dofs_per_cell);

                        for(unsigned int k = 0; k < dofs_per_cell; ++k)
                        {
                            phi_u[k] = fe_values_bdry[velocities].value(k,q);
                        }

                        // --- step2 ---

                        // --- LOOP OVER i ---

                        for(unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                            double sum = 0.0;

                            // --- LOOP OVER alpha ---

                            for(unsigned int alpha = 0; alpha < dim - 1; ++alpha)
                            {
                                sum += ( phi_u[i] * tangential_vectors[alpha][q] ) * ( u_old[q] * tangential_vectors[alpha][q] );
                            }

                            sum *= (1.0 - theta) / theta;

                            local_residual(i) += - sum * JxW
                                                 -
                                                 normal_velocity_suppressed_weakly * (1.0/eta) * ( phi_u[i] * normal_vectors[q] ) * ( u_old[q] * normal_vectors[q] ) * JxW;
                        }
                    }
                }

                // --- PERFECT SLIP ---

                else if( iter->boundary_condition.compare("perfect slip") == 0 )
                {
                    // --- LOOP OVER QUADRATURE POINTS ---

                    for(unsigned int q = 0; q < n_q_points_bdry; ++q)
                    {
                        // --- step1 ---

                        // --- JxW ---

                        const double JxW = fe_values_bdry.JxW(q);

                        // --- shapes ---

                        std::vector< Tensor<1,dim> > phi_u(dofs_per_cell);

                        for(unsigned int k = 0; k < dofs_per_cell; ++k)
                        {
                            phi_u[k] = fe_values_bdry[velocities].value(k,q);
                        }

                        // --- step2 ---

                        // --- LOOP OVER i ---

                        for(unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                            local_residual(i) += - normal_velocity_suppressed_weakly * (1.0/eta) * ( phi_u[i] * normal_vectors[q] ) * ( u_old[q] * normal_vectors[q] ) * JxW;
                        }
                    }
                }

                // --- UNKNOWN ---

                else
                {
                    FcstUtilities::log << "Boundary condition you specified at \"impermeable walls\" is not implemented" << std::endl;
                    AssertThrow( false , ExcNotImplemented() );
                }
            }

            // --- SYMMETRY LINE OR PLANE ---

            else if( iter->boundary_name.compare("symmetry line or plane") == 0 )
            {
                // --- PERFECT SLIP ---

                if( iter->boundary_condition.compare("perfect slip") == 0 )
                {
                    // --- LOOP OVER QUADRATURE POINTS ---

                    for(unsigned int q = 0; q < n_q_points_bdry; ++q)
                    {
                        // --- step1 ---

                        // --- JxW ---

                        const double JxW = fe_values_bdry.JxW(q);

                        // --- shapes ---

                        std::vector< Tensor<1,dim> > phi_u(dofs_per_cell);

                        for(unsigned int k = 0; k < dofs_per_cell; ++k)
                        {
                            phi_u[k] = fe_values_bdry[velocities].value(k,q);
                        }

                        // --- step2 ---

                        // --- LOOP OVER i ---

                        for(unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                            local_residual(i) += - normal_velocity_suppressed_weakly * (1.0/eta) * ( phi_u[i] * normal_vectors[q] ) * ( u_old[q] * normal_vectors[q] ) * JxW;
                        }
                    }
                }

                // --- UNKNOWN ---

                else
                {
                    FcstUtilities::log << "Boundary condition you specified at \"symmetry line or plane\" is not implemented" << std::endl;
                    AssertThrow( false , ExcNotImplemented() );
                }
            }

            // --- INLET-OUTLET ---

            else if( iter->boundary_name.compare("inlet-outlet") == 0 )
            {
                // --- DIRICHLET PRESSURE ---

                if( iter->boundary_condition.compare("Dirichlet pressure") == 0 )
                {
                        // do nothing
                }

                // --- DIRICHLET PRESSURE AND NORMAL STRESS FREE ---

                else if( iter->boundary_condition.compare("Dirichlet pressure and normal stress free") == 0 )
                {
                    int_2 = 0.0;
                }

                // --- DIRICHLET VELOCITY ---

                else if( iter->boundary_condition.compare("Dirichlet velocity") == 0 )
                {
                    int_1 = 0.0;
                    int_2 = 0.0;
                }

                // --- DIRICHLET PRESSURE AND VELOCITY ---

                else if( iter->boundary_condition.compare("Dirichlet pressure and velocity") == 0 )
                {
                    int_1 = 0.0;
                    int_2 = 0.0;
                }

                // --- NORMAL STRESS FREE ---

                else if( iter->boundary_condition.compare("normal stress free") == 0 )
                {
                    int_2 = 0.0;
                }

                // --- UNKNOWN ---

                else
                {
                    FcstUtilities::log << "Boundary condition you specified at \"inlet-outlet\" is not implemented" << std::endl;
                    AssertThrow( false , ExcNotImplemented() );
                }

                // --- LOOP OVER QUADRATURE POINTS ---

                for(unsigned int q = 0; q < n_q_points_bdry; ++q)
                {
                    // --- step1 ---

                    // --- JxW ---

                    const double JxW = fe_values_bdry.JxW(q);

                    // --- shapes ---

                    std::vector< Tensor<1,dim> > phi_u(dofs_per_cell);

                    for(unsigned int k = 0; k < dofs_per_cell; ++k)
                    {
                        phi_u[k] = fe_values_bdry[velocities].value(k,q);
                    }

                    // --- shape combinations ---

                    std::vector< SymmetricTensor<2,dim> > n_BY_phi_u_S(dofs_per_cell);

                    for(unsigned int k = 0; k < dofs_per_cell; ++k)
                    {
                        Tensor<2,dim> n_BY_phi_u;
                        outer_product(n_BY_phi_u, normal_vectors[q], phi_u[k]);
                        n_BY_phi_u_S[k] = 0.5 * ( n_BY_phi_u + transpose(n_BY_phi_u) );
                    }

                    // --- step2 ---

                    // --- LOOP OVER i ---

                    for(unsigned int i = 0; i < dofs_per_cell; ++i)
                    {
                        local_residual(i) += - int_1 * n_BY_phi_u_S[i] * u_old_BY_u_old[q] * JxW
                                             +
                                             int_2 * n_BY_phi_u_S[i] * ( - p_old[q]*unit + tau_old[q] ) * JxW;
                    }
                }
            }
        }
    }

    // --- conversion ---

    if( this->counter[1] )
    {
            this->make_residual_indices();
            this->counter[1] = false;
    }

    this->dealII_to_appframe(bdry_residual,
                             local_residual,
                             this->residual_indices);
}

/////////////////////////////
/////////////////////////////
// OTHER - make_ FUNCTIONS //
/////////////////////////////
/////////////////////////////

// ---                              ---
// --- make_internal_cell_couplings ---
// ---                              ---

template<int dim>
void
NAME::IncompressibleSingleComponentNSEquations<dim>::make_internal_cell_couplings()
{
    std::map< std::string, DoFTools::Coupling > tmp;

    #ifdef _1D_
        tmp["single_fluid_pressure"]   = DoFTools::none;
        tmp["single_fluid_velocity_X"] = DoFTools::always;
        this->internal_cell_couplings["Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - mass conservation"]       = tmp;
        tmp.clear();

        tmp["single_fluid_pressure"]   = DoFTools::always;
        tmp["single_fluid_velocity_X"] = DoFTools::always;
        this->internal_cell_couplings["Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X"] = tmp;
        tmp.clear();

    #endif

    #ifdef _2D_
        tmp["single_fluid_pressure"]   = DoFTools::none;
        tmp["single_fluid_velocity_X"] = DoFTools::always;
        tmp["single_fluid_velocity_Y"] = DoFTools::always;
        this->internal_cell_couplings["Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - mass conservation"]       = tmp;
        tmp.clear();

        tmp["single_fluid_pressure"]   = DoFTools::always;
        tmp["single_fluid_velocity_X"] = DoFTools::always;
        tmp["single_fluid_velocity_Y"] = DoFTools::always;
        this->internal_cell_couplings["Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X"] = tmp;
        tmp.clear();

        tmp["single_fluid_pressure"]   = DoFTools::always;
        tmp["single_fluid_velocity_X"] = DoFTools::always;
        tmp["single_fluid_velocity_Y"] = DoFTools::always;
        this->internal_cell_couplings["Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Y"] = tmp;
        tmp.clear();

    #endif

    #ifdef _3D_

        tmp["single_fluid_pressure"]   = DoFTools::none;
        tmp["single_fluid_velocity_X"] = DoFTools::always;
        tmp["single_fluid_velocity_Y"] = DoFTools::always;
        tmp["single_fluid_velocity_Z"] = DoFTools::always;
        this->internal_cell_couplings["Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - mass conservation"]       = tmp;

        tmp.clear();

        tmp["single_fluid_pressure"]   = DoFTools::always;
        tmp["single_fluid_velocity_X"] = DoFTools::always;
        tmp["single_fluid_velocity_Y"] = DoFTools::always;
        tmp["single_fluid_velocity_Z"] = DoFTools::always;
        this->internal_cell_couplings["Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X"] = tmp;

        tmp.clear();

        tmp["single_fluid_pressure"]   = DoFTools::always;
        tmp["single_fluid_velocity_X"] = DoFTools::always;
        tmp["single_fluid_velocity_Y"] = DoFTools::always;
        tmp["single_fluid_velocity_Z"] = DoFTools::always;
        this->internal_cell_couplings["Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Y"] = tmp;

        tmp.clear();

        tmp["single_fluid_pressure"]   = DoFTools::always;
        tmp["single_fluid_velocity_X"] = DoFTools::always;
        tmp["single_fluid_velocity_Y"] = DoFTools::always;
        tmp["single_fluid_velocity_Z"] = DoFTools::always;
        this->internal_cell_couplings["Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Z"] = tmp;

        tmp.clear();

    #endif
}

// ---                           ---
// --- make_matrix_block_indices ---
// ---                           ---

template<int dim>
void
NAME::IncompressibleSingleComponentNSEquations<dim>::make_matrix_block_indices()
{
    unsigned int index;

    // --- 1D --- --- 1D --- --- 1D ---
    #ifdef _1D_
        // --- mass conservation ---
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - mass conservation",
                                                            "single_fluid_velocity_X");
        this->matrix_block_indices.push_back(index);
        // --- momentum conservation X ---
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X",
                                                            "single_fluid_pressure");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X",
                                                            "single_fluid_velocity_X");
        this->matrix_block_indices.push_back(index);
    #endif

    // --- 2D --- --- 2D --- --- 2D ---
    #ifdef _2D_
        // --- mass conservation ---
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - mass conservation",
                                                            "single_fluid_velocity_X");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - mass conservation",
                                                            "single_fluid_velocity_Y");
        this->matrix_block_indices.push_back(index);
        // --- momentum conservation X ---
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X",
                                                            "single_fluid_pressure");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X",
                                                            "single_fluid_velocity_X");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X",
                                                            "single_fluid_velocity_Y");
        this->matrix_block_indices.push_back(index);

        // --- momentum conservation Y ---
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Y",
                                                            "single_fluid_pressure");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Y",
                                                            "single_fluid_velocity_X");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Y",
                                                            "single_fluid_velocity_Y");
        this->matrix_block_indices.push_back(index);
    #endif

    // --- 3D --- --- 3D --- --- 3D ---
    #ifdef _3D_

    // --- mass conservation ---
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - mass conservation",
                                                            "single_fluid_velocity_X");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - mass conservation",
                                                            "single_fluid_velocity_Y");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - mass conservation",
                                                            "single_fluid_velocity_Z");
        this->matrix_block_indices.push_back(index);
        // --- momentum conservation X ---
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X",
                                                            "single_fluid_pressure");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X",
                                                            "single_fluid_velocity_X");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X",
                                                            "single_fluid_velocity_Y");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X",
                                                            "single_fluid_velocity_Z");
        this->matrix_block_indices.push_back(index);

        // --- momentum conservation Y ---
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Y",
                                                            "single_fluid_pressure");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Y",
                                                            "single_fluid_velocity_X");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Y",
                                                            "single_fluid_velocity_Y");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Y",
                                                            "single_fluid_velocity_Z");
        this->matrix_block_indices.push_back(index);

        // --- momentum conservation Z ---

        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Z",
                                                            "single_fluid_pressure");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Z",
                                                            "single_fluid_velocity_X");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Z",
                                                            "single_fluid_velocity_Y");
        this->matrix_block_indices.push_back(index);
        index = this->system_management->matrix_block_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Z",
                                                            "single_fluid_velocity_Z");
        this->matrix_block_indices.push_back(index);

    #endif
}

// ---                       ---
// --- make_residual_indices ---
// ---                       ---

template<int dim>
void
NAME::IncompressibleSingleComponentNSEquations<dim>::make_residual_indices()
{
    unsigned int index;

    // --- 1D --- --- 1D --- --- 1D ---
    #ifdef _1D_
        // --- mass conservation ---
        index = this->system_management->equation_name_to_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - mass conservation");
        this->residual_indices.push_back(index);
        // --- momentum conservation X ---
        index = this->system_management->equation_name_to_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X");
        this->residual_indices.push_back(index);
    #endif

    // --- 2D --- --- 2D --- --- 2D ---
    #ifdef _2D_
        // --- mass conservation ---
        index = this->system_management->equation_name_to_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - mass conservation");
        this->residual_indices.push_back(index);
        // --- momentum conservation X ---
        index = this->system_management->equation_name_to_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X");
        this->residual_indices.push_back(index);
        // --- momentum conservation Y ---
        index = this->system_management->equation_name_to_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Y");
        this->residual_indices.push_back(index);
    #endif

    // --- 3D --- --- 3D --- --- 3D ---
    #ifdef _3D_
        // --- mass conservation ---
        index = this->system_management->equation_name_to_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - mass conservation");
        this->residual_indices.push_back(index);
        // --- momentum conservation X ---
        index = this->system_management->equation_name_to_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation X");
        this->residual_indices.push_back(index);
        // --- momentum conservation Y ---
        index = this->system_management->equation_name_to_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Y");
        this->residual_indices.push_back(index);
        // --- momentum conservation Z ---
        index = this->system_management->equation_name_to_index("Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component - momentum conservation Z");
        this->residual_indices.push_back(index);
    #endif
}

////////////////////////
////////////////////////
// ACCESSORS AND INFO //
////////////////////////
////////////////////////

// ---                     ---
// --- print_equation_info ---
// ---                     ---

template<int dim>
void
NAME::IncompressibleSingleComponentNSEquations<dim>::print_equation_info() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "------------------------------";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "Parameters for Navier-Stokes Fluid Transport Equations - steady-state - incompressible - isothermal - single-phase - single-component:";
    FcstUtilities::log << std::endl;

    /////////////
    // CHANNEL //
    /////////////
    if( !this->counter[2] )
    {
        FcstUtilities::log << "Boolean constants:";
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Inertia in channels                  = " << inertia_in_channels;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Shear stress in channels             = " << shear_stress_in_channels;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Gravity in channels                  = " << gravity_in_channels;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Normal velocity is suppressed weakly = " << normal_velocity_is_suppressed_weakly;
        FcstUtilities::log << std::endl;
        FcstUtilities::log << std::endl;
    }

    //////////////////////////////
    // CHANNEL AND POROUS LAYER //
    //////////////////////////////
    else
    {
        FcstUtilities::log << "Boolean constants and form of the drag force in porous media:";
        FcstUtilities::log << std::endl;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Inertia in channels                    = " << inertia_in_channels;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Shear stress in channels               = " << shear_stress_in_channels;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Gravity in channels                    = " << gravity_in_channels;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Inertia in porous media                = " << inertia_in_porous_media;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Shear stress in porous media           = " << shear_stress_in_porous_media;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Gravity in porous media                = " << gravity_in_porous_media;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Form of the drag force in porous media = " << drag_in_porous_media;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Normal velocity is suppressed weakly   = " << normal_velocity_is_suppressed_weakly;
        FcstUtilities::log << std::endl;
        FcstUtilities::log << std::endl;
    }

    FcstUtilities::log << "Computational constants:";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;

    if( normal_velocity_is_suppressed_weakly )
    {
        FcstUtilities::log << "Normal velocity suppression coefficient = " << eta;
        FcstUtilities::log << std::endl;
        FcstUtilities::log << std::endl;
    }

    if( inlet_outlet_velocity_max != 0.0 )
    {
        FcstUtilities::log << "Maximum inlet-outlet velocity [m/s] = " << inlet_outlet_velocity_max;
        FcstUtilities::log << std::endl;
        FcstUtilities::log << std::endl;
    }

    FcstUtilities::log << "Initial data:";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;

    for(component_materialID_value_map::const_iterator iter  = this->component_materialID_value.begin();
                                                        iter != this->component_materialID_value.end(); ++iter)
    {
        std::map<types::material_id, double> tmp = iter->second;

        for(std::map<types::material_id, double>::const_iterator iter2  = tmp.begin(); iter2 != tmp.end(); ++iter2)
        {
            FcstUtilities::log << "Name of the solution component  = " << iter->first   << std::endl;
            FcstUtilities::log << "Material id                     = " << iter2->first  << std::endl;
            FcstUtilities::log << "Value of the solution component = " << iter2->second << std::endl;
            FcstUtilities::log << std::endl;
        }
    }

    FcstUtilities::log << "Boundary conditions:";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;

    std::vector< BoundaryType >::const_iterator iter;

    for( iter = this->boundary_types.begin(); iter != this->boundary_types.end(); ++iter )
    {
        FcstUtilities::log << iter->boundary_name      << std::endl;
        FcstUtilities::log << iter->boundary_id        << std::endl;
        FcstUtilities::log << iter->boundary_condition << std::endl;
        if( iter->boundary_condition == "Navier slip" )
                FcstUtilities::log << "Navier slip coefficient = " << theta << std::endl;
        FcstUtilities::log << std::endl;
    }

    FcstUtilities::log << "Boundary data:";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;

    for(component_boundaryID_value_map::const_iterator iter  = this->component_boundaryID_value.begin();
                                                        iter != this->component_boundaryID_value.end(); ++iter)
    {
        std::map<types::boundary_id, double> tmp = iter->second;

        for(std::map<types::boundary_id, double>::const_iterator iter2  = tmp.begin(); iter2 != tmp.end(); ++iter2)
        {
            FcstUtilities::log << "Name of the solution component  = " << iter->first   << std::endl;
            FcstUtilities::log << "Boundary id                     = " << iter2->first  << std::endl;
            FcstUtilities::log << "Value of the solution component = " << iter2->second << std::endl;
            FcstUtilities::log << std::endl;
        }
    }

    FcstUtilities::log << "------------------------------";
    FcstUtilities::log << std::endl;
}

/////////////////////////////
/////////////////////////////
// EXPLICIT INSTANTIATIONS //
/////////////////////////////
/////////////////////////////

// ---                                          ---
// --- IncompressibleSingleComponentNSEquations ---
// ---                                          ---

template class NAME::IncompressibleSingleComponentNSEquations<deal_II_dimension>;