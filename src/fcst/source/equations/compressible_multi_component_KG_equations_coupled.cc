// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: compressible_multi_component_KG_equations_coupled.cc
// - Description: This class describes steady-state compressible and isothermal Kerkhof-Geboers
//   fluid transport equations for a single-phase multi-component case coupled with fuel cell
//   physics
// - Developers: Valentin N. Zingan, Chad Balen and Marc Secanell, University of Alberta
//
// ----------------------------------------------------------------------------

#include "equations/compressible_multi_component_KG_equations_coupled.h"

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
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::CompressibleMultiComponentKGEquationsCoupled
    (FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data)
:
NAME::EquationBase<dim>(system_management,data)
{
    FcstUtilities::log << "->FuelCellShop::Equation::CompressibleMultiComponentKGEquationsCoupled" << std::endl;

    R_universal           = Constants::R()*1.0e7; // Use CGS units, i.e., (g cm2)/(s2 K mol)
    gravity_acceleration  = Constants::gravity_acceleration()*Units::convert(1, Units::C_UNIT, Units::UNIT); // Use cm/s^2 instead of m/s^2
    unit                  = Constants::unit_tensor();
    
    this->counter.resize(4, true);
    this->counter[3] = false;
}

// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::~CompressibleMultiComponentKGEquationsCoupled()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::declare_parameters(ParameterHandler& param) const
{
    const unsigned int n_species_max = 5;

    param.enter_subsection("Equations");
    {
        param.enter_subsection("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component");
        {
            param.declare_entry("Coupled with fuel cell physics",
                                "false",
                                Patterns::Bool(),
                                " ");
            param.declare_entry("Number of species",
                                "1",
                                Patterns::Integer(1,5),
                                " ");
            param.declare_entry("Gas species",
                                """",
                                Patterns::Map( Patterns::Integer(), Patterns::Anything() ),
                                "Format id1:value1, id2:value2, ... where id# is an integer specifying the equation number for the gas species "
                                "and value# is a string of any of the following names: oxygen, water, nitrogen"); // TODO: fix the comment later

            param.enter_subsection("Boolean constants and form of the drag force in porous media");
            {
                for(unsigned int index = 1; index <= n_species_max; ++index)
                {
                    std::ostringstream streamOut;
                    streamOut << index;
                    
                    //Declare parameter for each species
                    std::string name = "Inertia in channels " + streamOut.str();
                    param.declare_entry(name.c_str(),
                                        "true",
                                        Patterns::Bool(),
                                        " ");
                    
                    name = "Shear stress in channels " + streamOut.str();
                    param.declare_entry(name.c_str(),
                                        "true",
                                        Patterns::Bool(),
                                        " ");
                    
                    name = "Gravity in channels " + streamOut.str();
                    param.declare_entry(name.c_str(),
                                        "false",
                                        Patterns::Bool(),
                                        " ");
                    
                    name = "Inertia in porous media " + streamOut.str();
                    param.declare_entry(name.c_str(),
                                        "true",
                                        Patterns::Bool(),
                                        " ");
                    
                    name = "Shear stress in porous media " + streamOut.str();
                    param.declare_entry(name.c_str(),
                                        "true",
                                        Patterns::Bool(),
                                        " ");
                    
                    name = "Gravity in porous media " + streamOut.str();
                    param.declare_entry(name.c_str(),
                                        "false",
                                        Patterns::Bool(),
                                        " ");
                    
                    name = "Drag in porous media " + streamOut.str();
                    param.declare_entry(name.c_str(),
                                        "none",
                                        Patterns::Selection("none | Darcy | Forchheimer | Forchheimer modified"),
                                        " ");
                    
                    name = "Normal velocity is suppressed weakly " + streamOut.str();
                    param.declare_entry(name.c_str(),
                                        "false",
                                        Patterns::Bool(),
                                        " ");  
                }
            }
            param.leave_subsection();

            param.enter_subsection("Computational constants");
            {
                for(unsigned int index = 1; index <= n_species_max; ++index)
                {
                    std::ostringstream streamOut;
                    streamOut << index;
                    
                    //Declare parameter for each species
                    std::string name = "Navier slip coefficient " + streamOut.str();
                    param.declare_entry(name.c_str(),
                                        "0.0",
                                        Patterns::Double(0.0,1.0),
                                        " ");
                    
                    name = "Normal velocity suppression coefficient " + streamOut.str();
                    param.declare_entry(name.c_str(),
                                        "1.0e-12",
                                        Patterns::Double(0.0,1.0),
                                        " ");
                    
                    name = "Maximum inlet-outlet velocity " + streamOut.str() + " [cm/s]";
                    param.declare_entry(name.c_str(),
                                        "0.0",
                                        Patterns::Double(),
                                        " ");
                    
                    name = "Maxwell slip coefficient " + streamOut.str();
                    param.declare_entry(name.c_str(),
                                        "1.0",
                                        Patterns::Double(),
                                        " ");
                }

                param.declare_entry("Maximum inlet-outlet velocity of the mixture [cm/s]",
                                    "1.0", //Coupled
                                    Patterns::Double(),
                                    " ");
                
                param.declare_entry("Use parabolic profile",
                                    "true",
                                    Patterns::Bool(),
                                    "Bool for whether OpenFCST will calculate parabolic profile at inlet-outlet or use user defined equation");
                
                param.declare_entry("inlet-outlet velocity channel base [cm]",
                                    "0.0",
                                    Patterns::Double(),
                                    "If user wants a parabolic velocity profile, OpenFCST needs the base of the channel to be specified");
                
                param.declare_entry("inlet-outlet velocity channel roof [cm]",
                                    "0.0",
                                    Patterns::Double(),
                                    "If user wants a parabolic velocity profile, OpenFCST needs the roof of the channel to be specified");
                
                param.declare_entry("pressure or velocity component that variable boundary condition will be applied to",
                                    "",
                                    Patterns::Anything(),
                                    "Define which solution variable to apply variable boundary condition equation to, ex Pressure, VelocityX, VelocityY, or VelocityZ");
                                
                param.declare_entry("inlet-outlet velocity of the mixture equation",
                                    "",
                                    Patterns::Anything(),
                                    "For user defined velocity profile at inlet/outlet, ex -625*y*y + 66.25*y - 0.755625");
            }
            param.leave_subsection();

            param.enter_subsection("Initial data");
            {
                param.declare_entry("Variable initial data",
                                    "false",
                                    Patterns::Bool(),
                                    " ");

                for(unsigned int index = 1; index <= n_species_max; ++index)
                {
                    std::ostringstream streamOut;
                    streamOut << index;
                    
                    //Declare parameter for each species
                    std::string name = "density_species_" + streamOut.str() + " [g/cm^3]";
                    param.declare_entry(name.c_str(),
                                        "",
                                        Patterns::Map(   Patterns::Integer(0,255) , Patterns::Double()   ),
                                        " ");
                    
                    name = "velocity_species_" + streamOut.str() + "_X [cm/s]";
                    param.declare_entry(name.c_str(),
                                        "",
                                        Patterns::Map(   Patterns::Integer(0,255) , Patterns::Double()   ),
                                        " ");
                    
                    name = "velocity_species_" + streamOut.str() + "_Y [cm/s]";
                    param.declare_entry(name.c_str(),
                                        "",
                                        Patterns::Map(   Patterns::Integer(0,255) , Patterns::Double()   ),
                                        " ");
                    
                    name = "velocity_species_" + streamOut.str() + "_Z [cm/s]";
                    param.declare_entry(name.c_str(),
                                        "",
                                        Patterns::Map(   Patterns::Integer(0,255) , Patterns::Double()   ),
                                        " ");
                }
            }
            param.leave_subsection();

            param.enter_subsection("Boundary conditions");
            {
                for(unsigned int index = 1; index <= n_species_max; ++index)
                {
                    std::ostringstream streamOut;
                    streamOut << index;
                    
                    //Declare parameter for each species
                    std::string name = "Impermeable walls " + streamOut.str();
                    param.declare_entry(name.c_str(),
                                        "",
                                        Patterns::Map(   Patterns::Integer(0,255) , Patterns::Selection("no-slip | Navier slip | Maxwell slip | perfect slip")   ),
                                        " ");
                    
                    name = "Symmetry line or plane " + streamOut.str();
                    param.declare_entry(name.c_str(),
                                        "",
                                        Patterns::Map(   Patterns::Integer(0,255) , Patterns::Selection("perfect slip")   ),
                                        " ");
                    
                    name = "Inlet-Outlet " + streamOut.str();
                    param.declare_entry(name.c_str(),
                                        "",
                                        Patterns::Map(   Patterns::Integer(0,255) , Patterns::Selection("Dirichlet density | Dirichlet density and normal stress free | Dirichlet density and normal shear stress free | Dirichlet velocity | Dirichlet density and velocity | normal stress free | normal shear stress free")   ),
                                        " ");
                }
            }
            param.leave_subsection();

            param.enter_subsection("Boundary data");
            {
                param.declare_entry("Variable boundary data",
                                    "false",
                                    Patterns::Bool(),
                                    " ");

                for(unsigned int index = 1; index <= n_species_max; ++index)
                {
                    std::ostringstream streamOut;
                    streamOut << index;
                    
                    //Declare parameter for each species
                    std::string name = "density_species_" + streamOut.str() + " [g/cm^3]";
                    param.declare_entry(name.c_str(),
                                        "",
                                        Patterns::Map(   Patterns::Integer(0,255) , Patterns::Double()   ),
                                        " ");
                    
                    name = "velocity_species_" + streamOut.str() + "_X [cm/s]";
                    param.declare_entry(name.c_str(),
                                        "",
                                        Patterns::Map(   Patterns::Integer(0,255) , Patterns::Double()   ),
                                        " ");
                    
                    name = "velocity_species_" + streamOut.str() + "_Y [cm/s]";
                    param.declare_entry(name.c_str(),
                                        "",
                                        Patterns::Map(   Patterns::Integer(0,255) , Patterns::Double()   ),
                                        " ");
                
                    name = "velocity_species_" + streamOut.str() + "_Z [cm/s]";
                    param.declare_entry(name.c_str(),
                                        "",
                                        Patterns::Map(   Patterns::Integer(0,255) , Patterns::Double()   ),
                                        " ");
                }
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
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::initialize(ParameterHandler& param)
{
    param.enter_subsection("Equations");
    {
        param.enter_subsection("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component");
        {
            coupled_with_fuel_cell_physics = param.get_bool("Coupled with fuel cell physics");
            n_species = param.get_integer("Number of species");
            gas_species_map = FcstUtilities::string_to_map<unsigned int, std::string>( param.get("Gas species") );
            
            // resize species info containers
            molar_mass.resize(n_species);
            dynamic_viscosity.resize(n_species);
            collision_diameter.resize(n_species);
            
            param.enter_subsection("Boolean constants and form of the drag force in porous media");
            {
                for(unsigned int index = 1; index <= n_species; ++index)
                {
                    std::ostringstream streamOut;
                    streamOut << index;
                    //Store parameter values for each species
                    std::string name = "Inertia in channels " + streamOut.str();
                    inertia_in_channels.push_back( param.get_bool(name.c_str()) );
                    
                    name = "Shear stress in channels " + streamOut.str();
                    shear_stress_in_channels.push_back( param.get_bool(name.c_str()) );
                    
                    name = "Gravity in channels " + streamOut.str();
                    gravity_in_channels.push_back( param.get_bool(name.c_str()) );
                    
                    name = "Inertia in porous media " + streamOut.str();
                    inertia_in_porous_media.push_back( param.get_bool(name.c_str()) );
                    
                    name = "Shear stress in porous media " + streamOut.str();
                    shear_stress_in_porous_media.push_back( param.get_bool(name.c_str()) );
                    
                    name = "Gravity in porous media " + streamOut.str();
                    gravity_in_porous_media.push_back( param.get_bool(name.c_str()) );
                    
                    name = "Drag in porous media " + streamOut.str();
                    drag_in_porous_media.push_back( param.get(name.c_str()) );
                    
                    if( drag_in_porous_media[index-1] != "none" )
                        this->counter[3] = true;
                    
                    name = "Normal velocity is suppressed weakly " + streamOut.str();
                    normal_velocity_is_suppressed_weakly.push_back( param.get_bool(name.c_str()) );
                }
            }
            param.leave_subsection();

            param.enter_subsection("Computational constants");
            {
                for(unsigned int index = 1; index <= n_species; ++index)
                {
                    std::ostringstream streamOut;
                    streamOut << index;
                    
                    //Store parameter values for each species
                    std::string name = "Navier slip coefficient " + streamOut.str();
                    theta.push_back( param.get_double(name.c_str()) );
                    
                    name = "Normal velocity suppression coefficient " + streamOut.str();
                    eta.push_back( param.get_double(name.c_str()) );
                    
                    name = "Maximum inlet-outlet velocity " + streamOut.str() + " [cm/s]";
                    inlet_outlet_velocity_max.push_back( param.get_double(name.c_str()) );
                    
                    name = "Maxwell slip coefficient " + streamOut.str();
                    maxwell_constant.push_back( param.get_double(name.c_str()) );
                }

                inlet_outlet_velocity_mixture_max      = param.get_double("Maximum inlet-outlet velocity of the mixture [cm/s]");
                use_parabolic_profile                  = param.get_bool("Use parabolic profile");
                inlet_outlet_velocity_channel_base     = param.get_double("inlet-outlet velocity channel base [cm]");
                inlet_outlet_velocity_channel_roof     = param.get_double("inlet-outlet velocity channel roof [cm]");
                press_vel_comp_apply_to                = param.get("pressure or velocity component that variable boundary condition will be applied to");
                inlet_outlet_velocity_mixture_equation = param.get("inlet-outlet velocity of the mixture equation");
            }
            param.leave_subsection();

            param.enter_subsection("Initial data");
            {
                this->variable_initial_data = param.get_bool("Variable initial data");

                for(unsigned int index = 1; index <= n_species; ++index)
                {
                    std::ostringstream streamOut;
                    streamOut << index;
                    
                    //Store parameter values for each species
                    std::string name = "density_species_" + streamOut.str();
                    std::string paramName = name + " [g/cm^3]";
                    if( !param.get(paramName.c_str()).empty() )
                    {
                        const std::map<types::material_id, double> tmp = FcstUtilities::string_to_map<types::material_id, double>( param.get(paramName.c_str()) );
                        this->component_materialID_value[name.c_str()] = tmp;
                    }
                    
                    name = "velocity_species_" + streamOut.str() + "_X";
                    paramName = name + " [cm/s]";
                    if( !param.get(paramName.c_str()).empty() )
                    {
                        const std::map<types::material_id, double> tmp = FcstUtilities::string_to_map<types::material_id, double>( param.get(paramName.c_str()) );
                        this->component_materialID_value[name.c_str()] = tmp;
                    }
                    
                    name = "velocity_species_" + streamOut.str() + "_Y";
                    paramName = name + " [cm/s]";
                    if( !param.get(paramName.c_str()).empty() )
                    {
                        const std::map<types::material_id, double> tmp = FcstUtilities::string_to_map<types::material_id, double>( param.get(paramName.c_str()) );
                        this->component_materialID_value[name.c_str()] = tmp;
                    }
                    
                    name = "velocity_species_" + streamOut.str() + "_Z";
                    paramName = name + " [cm/s]";
                    if( !param.get(paramName.c_str()).empty() )
                    {
                        const std::map<types::material_id, double> tmp = FcstUtilities::string_to_map<types::material_id, double>( param.get(paramName.c_str()) );
                        this->component_materialID_value[name.c_str()] = tmp;
                    }
                }
            }
            param.leave_subsection();

            param.enter_subsection("Boundary conditions");
            {
                this->multi_boundary_types.resize(n_species);

                for(unsigned int index = 1; index <= n_species; ++index)
                {
                    std::ostringstream streamOut;
                    streamOut << index;
                    
                    //Store parameter values for each species
                    std::string name = "Impermeable walls " + streamOut.str();
                    if( !param.get(name.c_str()).empty() )
                    {
                        const std::map<unsigned int, std::string> tmp = FcstUtilities::string_to_map<unsigned int, std::string>( param.get(name.c_str()) );

                        std::map<unsigned int, std::string>::const_iterator iter;

                        for( iter = tmp.begin(); iter != tmp.end(); ++iter )
                        {
                            BoundaryType ImpermeableWalls;
                            ImpermeableWalls.boundary_name      = "impermeable walls";
                            ImpermeableWalls.boundary_id        =  iter->first;
                            ImpermeableWalls.boundary_condition =  iter->second;

                            this->multi_boundary_types[index-1].push_back(ImpermeableWalls);
                        }
                    }
                    
                    name = "Symmetry line or plane " + streamOut.str();
                    if( !param.get(name.c_str()).empty() )
                    {
                        const std::map<unsigned int, std::string> tmp = FcstUtilities::string_to_map<unsigned int, std::string>( param.get(name.c_str()) );

                        std::map<unsigned int, std::string>::const_iterator iter;

                        for( iter = tmp.begin(); iter != tmp.end(); ++iter )
                        {
                            BoundaryType SymmetryLineOrPlane;
                            SymmetryLineOrPlane.boundary_name      = "symmetry line or plane";
                            SymmetryLineOrPlane.boundary_id        =  iter->first;
                            SymmetryLineOrPlane.boundary_condition =  iter->second;

                            this->multi_boundary_types[index-1].push_back(SymmetryLineOrPlane);
                        }
                    }
                    
                    name = "Inlet-Outlet " + streamOut.str();
                    if( !param.get(name.c_str()).empty() )
                    {
                        const std::map<unsigned int, std::string> tmp = FcstUtilities::string_to_map<unsigned int, std::string>( param.get(name.c_str()) );

                        std::map<unsigned int, std::string>::const_iterator iter;

                        for( iter = tmp.begin(); iter != tmp.end(); ++iter )
                        {
                            BoundaryType InletOutlet;
                            InletOutlet.boundary_name      = "inlet-outlet";
                            InletOutlet.boundary_id        =  iter->first;
                            InletOutlet.boundary_condition =  iter->second;

                            this->multi_boundary_types[index-1].push_back(InletOutlet);
                            if( ( InletOutlet.boundary_condition == "Dirichlet velocity") )
                                // Should I add this?
                                /*
                                || (InletOutlet.boundary_condition == "Dirichlet density") 
                                || (InletOutlet.boundary_condition == "Dirichlet density and normal stress free") 
                                || (InletOutlet.boundary_condition == "Dirichlet density and normal shear stress free") )
                                */
                                inlet_outlet_boundary_ID = InletOutlet.boundary_id;
                        }
                    }
                }
            }
            param.leave_subsection();

            param.enter_subsection("Boundary data");
            {
                this->variable_boundary_data = param.get_bool("Variable boundary data");

                for(unsigned int index = 1; index <= n_species; ++index)
                {
                    std::ostringstream streamOut;
                    streamOut << index;
                    
                    //Store parameter values for each species
                    std::string name = "density_species_" + streamOut.str();
                    std::string paramName = name + " [g/cm^3]";
                    if( !param.get(paramName.c_str()).empty() )
                    {
                        const std::map<types::boundary_id, double> tmp = FcstUtilities::string_to_map<types::boundary_id, double>( param.get(paramName.c_str()) );
                        this->component_boundaryID_value[name.c_str()] = tmp;
                    }
                    
                    name = "velocity_species_" + streamOut.str() + "_X";
                    paramName = name + " [cm/s]";
                    if( !param.get(paramName.c_str()).empty() )
                    {
                        const std::map<types::boundary_id, double> tmp = FcstUtilities::string_to_map<types::boundary_id, double>( param.get(paramName.c_str()) );
                        this->component_boundaryID_value[name.c_str()] = tmp;
                    }
                    
                    name = "velocity_species_" + streamOut.str() + "_Y";
                    paramName = name + " [cm/s]";
                    if( !param.get(paramName.c_str()).empty() )
                    {
                        const std::map<types::boundary_id, double> tmp = FcstUtilities::string_to_map<types::boundary_id, double>( param.get(paramName.c_str()) );
                        this->component_boundaryID_value[name.c_str()] = tmp;
                    }
                    
                    name = "velocity_species_" + streamOut.str() + "_Z";
                    paramName = name + " [cm/s]";
                    if( !param.get(paramName.c_str()).empty() )
                    {
                        const std::map<types::boundary_id, double> tmp = FcstUtilities::string_to_map<types::boundary_id, double>( param.get(paramName.c_str()) );
                        this->component_boundaryID_value[name.c_str()] = tmp;
                    }
                }
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();

    this->make_internal_cell_couplings();
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
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                                                              const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                              FuelCellShop::Layer::BaseLayer<dim>* const                               layer)
{
    //Generic Constant Data
    if( this->counter[0] )
    {
        this->make_assemblers_generic_constant_data(cell_info, layer);
        this->counter[0] = false;
    }

    //Cell Constant Data
    if( this->counter[1] )
    {
        this->make_assemblers_cell_constant_data(cell_info);
        this->counter[1] = false;
    }

    //Cell Variable Data
    this->make_assemblers_cell_variable_data(cell_info, layer);

    //Local Cell Matrix
    FullMatrix<double> local_matrix(this->dofs_per_cell, this->dofs_per_cell);

    //Loop Over Species
    for(unsigned int s = 0; s < n_species; ++s)
    {
        //Loop Over Quadrature Points
        for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
        {
            //Loop Over i
            for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
            {
                //Loop Over j
                for(unsigned int j = 0; j < this->dofs_per_cell; ++j)
                {
                    local_matrix(i,j) += this->JxW_cell[q] * (
                                           grad_phi_density_cell[s][q][i]   * delta_mass_flux_cell[s][q][j]
                                           +
                                           grads_phi_velocity_cell[s][q][i] * delta_momentum_flux_cell[s][q][j]
                                           +
                                           grads_phi_velocity_cell[s][q][i] * ( delta_pressure_cell[s][q][j] * unit - delta_shear_stress_cell[s][q][j] )
                                           +
                                           phi_velocity_cell[s][q][i]       * ( delta_drag_cell[s][q][j] + delta_diffusion_cell[s][q][j] + delta_gravity_cell[s][q][j] ) ) ;
                }
            }
        }
    }

    this->dealII_to_appframe(cell_matrices,
                             local_matrix,
                             this->matrix_block_indices);
}

// ---                        ---
// --- assemble_cell_residual ---
// ---                        ---

template<int dim>
void
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
                                                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                                FuelCellShop::Layer::BaseLayer<dim>* const                               layer)
{
    //Generic Constant Data
    if( this->counter[0] )
    {
        this->make_assemblers_generic_constant_data(cell_info, layer);
        this->counter[0] = false;
    }

    //Cell Constant Data
    if( this->counter[1] )
    {
        this->make_assemblers_cell_constant_data(cell_info);
        this->counter[1] = false;
    }

    //Cell Variable Data
    this->make_assemblers_cell_variable_data(cell_info, layer);

    //Local Cell Residual
    Vector<double> local_residual(this->dofs_per_cell);

    //Loop Over Species
    for(unsigned int s = 0; s < n_species; ++s)
    {
        //Loop Over Quadrature Points
        for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
        {
            //Loop Over i
            for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
            {
                local_residual(i) += this->JxW_cell[q] * ( 
                                       grad_phi_density_cell[s][q][i]   *   mass_flux_cell_old[s][q] 
                                        +
                                       grads_phi_velocity_cell[s][q][i] *   momentum_flux_cell_old[s][q]
                                        +
                                       grads_phi_velocity_cell[s][q][i] * ( pressure_cell_old[s][q] * unit - shear_stress_cell_old[s][q] )
                                        +
                                       phi_velocity_cell[s][q][i]       * ( drag_cell_old[s][q] + diffusion_cell_old[s][q] + gravity_cell_old[s][q] ) );
            }
        }
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
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::assemble_bdry_matrix(FuelCell::ApplicationCore::MatrixVector&                                 bdry_matrices,
                                                                              const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                                              FuelCellShop::Layer::BaseLayer<dim>* const                               layer)
{
    //Generic Constant Data
    if( this->counter[0] )
    {
        this->make_assemblers_generic_constant_data(bdry_info, layer);
        this->counter[0] = false;
    }
    
    //Boundary Constant Data
    if( this->counter[2] )
    {
        this->make_assemblers_bdry_constant_data(bdry_info);
        this->counter[2] = false;
    }
    
    //Boundary Variable Data
    this->make_assemblers_bdry_variable_data(bdry_info, layer);

    //Types Info
    const std::type_info& Channel                 = typeid(FuelCellShop::Layer::Channel<dim>);
    const std::type_info& GasDiffusionLayer       = typeid(FuelCellShop::Layer::GasDiffusionLayer<dim>);
    const std::type_info& MicroPorousLayer        = typeid(FuelCellShop::Layer::MicroPorousLayer<dim>);
    const std::type_info& CatalystLayer           = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    const std::type_info& info                    = layer->get_base_type();

    //Local Boundary Matrix
    FullMatrix<double> local_matrix(this->dofs_per_cell, this->dofs_per_cell);
    
    //Loop Over Species
    for(unsigned int s = 0; s < n_species; ++s)
    {
        //Initialize Integral Multipliers
        double int_1 = 1.0; double int_2 = 1.0;
        double int_3 = 1.0; double int_4 = 1.0;
        double rho_factor = 1.0; double u_factor = 1.0;
        
        //Extras
        double inertia = 0.0;
        double shear_stress = 0.0;
        double normal_velocity_suppressed_weakly = 0.0;
        
        if ( (info == Channel) && inertia_in_channels[s] )
            inertia = 1.0;
        else if ( (info == Channel) && shear_stress_in_channels[s] )
            shear_stress = 1.0;
        else if ( ( info == GasDiffusionLayer || info == MicroPorousLayer || info == CatalystLayer ) && inertia_in_porous_media[s] )
            inertia = 1.0;
        else if ( ( info == GasDiffusionLayer || info == MicroPorousLayer || info == CatalystLayer ) && shear_stress_in_porous_media[s] )
            shear_stress = 1.0;        
        else if ( normal_velocity_is_suppressed_weakly[s] )
            normal_velocity_suppressed_weakly = 1.0;
        /*
        else{
            std::cerr << "Layer you specified is not implemented in CompressibleMultiComponentKGEquationsCoupled<dim>::assemble_bdry_matrix" << std::endl;
            AssertThrow( false , ExcNotImplemented() );
        }
        */
        //Boundaries
        const std::vector< BoundaryType > tmp = this->multi_boundary_types[s];

        std::vector< BoundaryType >::const_iterator iter;
        
        //Loop Over Boundaries
        for( iter = tmp.begin(); iter != tmp.end(); ++iter )
        {
            if( this->belongs_to_boundary( bdry_info.dof_face->boundary_indicator(), iter->boundary_id ) )
            {
                if( iter->boundary_name.compare("impermeable walls") == 0 )
                {
                    if( iter->boundary_condition.compare("no-slip") == 0 )
                    {
                            // do nothing
                    }
                    else if( iter->boundary_condition.compare("Navier slip") == 0 )
                    {
                        //Loop Over Quadrature Points and DOFs
                        for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
                            for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                for(unsigned int j = 0; j < this->dofs_per_cell; ++j)
                                {
                                    double sum = 0.0;
                                    
                                    //Loop Over Alpha
                                    for(unsigned int alpha = 0; alpha < dim - 1; ++alpha)
                                    {
                                            sum += ( phi_velocity_bdry[s][q][i] * this->tangential_vectors[alpha][q] ) * ( phi_velocity_bdry[s][q][j] * this->tangential_vectors[alpha][q] );
                                    }

                                    sum *= (1.0 - theta[s]) / theta[s];

                                    local_matrix(i,j) += sum * this->JxW_bdry[q]
                                                         +
                                                         normal_velocity_suppressed_weakly * (1.0/eta[s]) * ( phi_velocity_bdry[s][q][i] * this->normal_vectors[q] ) * ( phi_velocity_bdry[s][q][j] * this->normal_vectors[q] ) * this->JxW_bdry[q];
                                }
                    }
                    else if( iter->boundary_condition.compare("Maxwell slip") == 0 )
                    {
                        //Loop Over Quadrature Points and DOFs
                        for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
                            for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                for(unsigned int j = 0; j < this->dofs_per_cell; ++j)
                                {
                                    double sum = 0.0;
                                    
                                    //Loop Over Alpha
                                    for(unsigned int alpha = 0; alpha < dim - 1; ++alpha)
                                    {
                                        sum += ( phi_velocity_bdry[s][q][i] * this->tangential_vectors[alpha][q] ) * ( delta_mass_flux_bdry[s][q][j] * this->tangential_vectors[alpha][q] );
                                    }

                                    sum *= 1.0/maxwell_constant[s];

                                    local_matrix(i,j) += sum * this->JxW_bdry[q]
                                                         +
                                                         normal_velocity_suppressed_weakly * (1.0/eta[s]) * ( phi_velocity_bdry[s][q][i] * this->normal_vectors[q] ) * ( phi_velocity_bdry[s][q][j] * this->normal_vectors[q] ) * this->JxW_bdry[q];
                                }
                    }
                    else if( iter->boundary_condition.compare("perfect slip") == 0 )
                    {
                        //Loop Over Quadrature Points and dofs
                        for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
                            for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                for(unsigned int j = 0; j < this->dofs_per_cell; ++j)
                                {
                                    local_matrix(i,j) += normal_velocity_suppressed_weakly * (1.0/eta[s]) * 
                                                         ( phi_velocity_bdry[s][q][i] * this->normal_vectors[q] ) * 
                                                         ( phi_velocity_bdry[s][q][j] * this->normal_vectors[q] ) * this->JxW_bdry[q];
                                }
                    }
                    else //Unknown, throw error
                    {
                        std::cerr << "Boundary condition you specified at \"impermeable walls\" is not implemented" << std::endl;
                        AssertThrow( false , ExcNotImplemented() );
                    }
                }
                else if( iter->boundary_name.compare("symmetry line or plane") == 0 )
                {
                    if( iter->boundary_condition.compare("perfect slip") == 0 )
                    {
                        //Loop Over Quadrature Points and DOFs
                        for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
                            for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                for(unsigned int j = 0; j < this->dofs_per_cell; ++j)
                                {
                                    local_matrix(i,j) += normal_velocity_suppressed_weakly * (1.0/eta[s]) * 
                                                         ( phi_velocity_bdry[s][q][i] * this->normal_vectors[q] ) * 
                                                         ( phi_velocity_bdry[s][q][j] * this->normal_vectors[q] ) * this->JxW_bdry[q];
                                }
                    }
                    else //Unknown, throw error
                    {
                        std::cerr << "Boundary condition you specified at \"symmetry line or plane\" is not implemented" << std::endl;
                        AssertThrow( false , ExcNotImplemented() );
                    }
                }
                else if( iter->boundary_name.compare("inlet-outlet") == 0 )
                {
                    if( iter->boundary_condition.compare("Dirichlet density") == 0 ) {
                        int_1      = 0.0; int_2      = 1.0;
                        int_3      = 1.0; int_4      = 1.0;
                        rho_factor = 0.0; u_factor   = 1.0;
                    }
                    else if( iter->boundary_condition.compare("Dirichlet density and normal stress free") == 0 )  {
                        int_1      = 0.0; int_2      = 1.0;
                        int_3      = 0.0; int_4      = 0.0;
                        rho_factor = 0.0; u_factor   = 1.0;
                    }
                    else if( iter->boundary_condition.compare("Dirichlet density and normal shear stress free") == 0 )   {
                        int_1      = 0.0; int_2      = 1.0;
                        int_3      = 1.0; int_4      = 0.0;
                        rho_factor = 0.0; u_factor   = 1.0;
                    }
                    else if( iter->boundary_condition.compare("Dirichlet velocity") == 0 )   {
                        int_1      = 1.0; int_2      = 0.0;
                        int_3      = 0.0; int_4      = 0.0;
                        rho_factor = 1.0; u_factor   = 0.0;
                    }
                    else if( iter->boundary_condition.compare("Dirichlet density and velocity") == 0 )    {
                        int_1      = 0.0; int_2      = 0.0;
                        int_3      = 0.0; int_4      = 0.0;
                        rho_factor = 0.0; u_factor   = 0.0;
                    }
                    else if( iter->boundary_condition.compare("normal stress free") == 0 )  {
                        int_1      = 1.0; int_2      = 1.0;
                        int_3      = 0.0; int_4      = 0.0;
                        rho_factor = 1.0; u_factor   = 1.0;
                    }
                    else if( iter->boundary_condition.compare("normal shear stress free") == 0 )   {
                        int_1      = 1.0; int_2      = 1.0;
                        int_3      = 1.0; int_4      = 0.0;
                        rho_factor = 1.0; u_factor   = 1.0;
                    }
                    else //Unknown, throw error
                    {
                        std::cerr << "Boundary condition you specified at \"inlet-outlet\" is not implemented" << std::endl;
                        AssertThrow( false , ExcNotImplemented() );
                    }
                    
                    //Loop Over Quadrature Points and DOFs
                    for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
                        for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                            for(unsigned int j = 0; j < this->dofs_per_cell; ++j)
                            {
                                Tensor<1,dim> tmp1, tmp2;
                                Tensor<2,dim> tmp3, tmp4, tmp5;                                
                                
                                tmp1 = phi_density_bdry[s][q][j] * velocity_bdry_old[s][q];
                                tmp2 = density_bdry_old[s][q]    * phi_velocity_bdry[s][q][j];
                                
                                outer_product(tmp3, velocity_bdry_old[s][q], velocity_bdry_old[s][q]);                                
                                outer_product(tmp4, phi_velocity_bdry[s][q][j], velocity_bdry_old[s][q]);                                
                                outer_product(tmp5, velocity_bdry_old[s][q], phi_velocity_bdry[s][q][j]);
                                
                                tmp3 *= phi_density_bdry[s][q][j];
                                tmp4 *= density_bdry_old[s][q];
                                tmp5 *= density_bdry_old[s][q];
                                
                                SymmetricTensor<2,dim> tmp;
                                tmp = rho_factor * tmp3 + tmp4 + tmp5;
                                
                                local_matrix(i,j) += this->JxW_bdry[q] * ( 
                                                       int_1 * phi_density_bdry[s][q][i] * ( tmp1 + u_factor * tmp2 ) * this->normal_vectors[q]
                                                       +
                                                       int_2 * inertia * n_BY_phi_velocity_bdry_S[s][q][i] * tmp
                                                       -
                                                       int_3 * n_BY_phi_velocity_bdry_S[s][q][i] * 
                                                       ( - rho_factor * delta_pressure_bdry[s][q][j] * unit + int_4 * shear_stress * delta_shear_stress_bdry[s][q][j] ) );
                            }
                }
            }
        }
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
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::assemble_bdry_residual(FuelCell::ApplicationCore::FEVector&                                     bdry_residual,
                                                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                                                FuelCellShop::Layer::BaseLayer<dim>* const                               layer)
{
    //Generic Constant Data
    if( this->counter[0] )
    {
        this->make_assemblers_generic_constant_data(bdry_info, layer);
        this->counter[0] = false;
    }
    
    //Boundary Constant Data
    if( this->counter[2] )
    {
        this->make_assemblers_bdry_constant_data(bdry_info);
        this->counter[2] = false;
    }
    
    //Boundary Variable Data
    this->make_assemblers_bdry_variable_data(bdry_info, layer);
    
    //Local Boundary Residual
    Vector<double> local_residual(this->dofs_per_cell);
    
    //Loop Over Species
    for(unsigned int s = 0; s < n_species; ++s)
    {
        //Integral Multipliers
        double int_1 = 1.0;
        double int_2 = 1.0;
        double int_3 = 1.0;
        double int_4 = 1.0;

        //Extras
        double normal_velocity_suppressed_weakly = 0.0;
        if( normal_velocity_is_suppressed_weakly[s] )
            normal_velocity_suppressed_weakly = 1.0;
        
        //Boundaries
        const std::vector< BoundaryType > tmp = this->multi_boundary_types[s];

        std::vector< BoundaryType >::const_iterator iter;
        
        //Loop Over Boundaryies
        for( iter = tmp.begin(); iter != tmp.end(); ++iter )
        {
            if( this->belongs_to_boundary( bdry_info.dof_face->boundary_indicator(),
                                           iter->boundary_id ) )
            {
                //Check boundary conditions to determine how to setup boundary residual
                if( iter->boundary_name.compare("impermeable walls") == 0 )
                {
                    if( iter->boundary_condition.compare("no-slip") == 0 )
                    {
                        // do nothing
                    }
                    else if( iter->boundary_condition.compare("Navier slip") == 0 )
                    {
                        //Loop Over Quadrature Points
                        for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
                        {
                            //Loop Over i
                            for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                            {
                                double sum = 0.0;
                                
                                //Loop Over Alpha
                                for(unsigned int alpha = 0; alpha < dim - 1; ++alpha)
                                {
                                    sum += ( phi_velocity_bdry[s][q][i] * this->tangential_vectors[alpha][q] ) * ( velocity_bdry_old[s][q] * this->tangential_vectors[alpha][q] );
                                }

                                sum *= (1.0 - theta[s]) / theta[s];

                                local_residual(i) += - sum * this->JxW_bdry[q]
                                                        -
                                                        normal_velocity_suppressed_weakly * (1.0/eta[s]) * ( phi_velocity_bdry[s][q][i] * this->normal_vectors[q] ) * ( velocity_bdry_old[s][q] * this->normal_vectors[q] ) * this->JxW_bdry[q];
                            }
                        }
                    }
                    else if( iter->boundary_condition.compare("Maxwell slip") == 0 )
                    {
                        //Loop Over Quadrature Points
                        for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
                        {
                            //Loop Over i
                            for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                            {
                                double sum = 0.0;
                                
                                //Loop Over Alpha
                                for(unsigned int alpha = 0; alpha < dim - 1; ++alpha)
                                {
                                    sum += ( phi_velocity_bdry[s][q][i] * this->tangential_vectors[alpha][q] ) * ( mass_flux_bdry_old[s][q] * this->tangential_vectors[alpha][q] );
                                }

                                sum *= 1.0/maxwell_constant[s];

                                local_residual(i) += - sum * this->JxW_bdry[q]
                                                        -
                                                        normal_velocity_suppressed_weakly * (1.0/eta[s]) * ( phi_velocity_bdry[s][q][i] * this->normal_vectors[q] ) * ( velocity_bdry_old[s][q] * this->normal_vectors[q] ) * this->JxW_bdry[q];
                            }
                        }
                    }
                    else if( iter->boundary_condition.compare("perfect slip") == 0 )
                    {
                        //Loop Over Quadrature Points
                        for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
                        {
                            //Loop Over i
                            for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                            {
                                local_residual(i) += - normal_velocity_suppressed_weakly * (1.0/eta[s]) * ( phi_velocity_bdry[s][q][i] * this->normal_vectors[q] ) * ( velocity_bdry_old[s][q] * this->normal_vectors[q] ) * this->JxW_bdry[q];
                            }
                        }
                    }
                    else //Unknown, throw error
                    {
                        std::cerr << "Boundary condition you specified at \"impermeable walls\" is not implemented" << std::endl;
                        AssertThrow( false , ExcNotImplemented() );
                    }
                }
                else if( iter->boundary_name.compare("symmetry line or plane") == 0 )
                {
                        if( iter->boundary_condition.compare("perfect slip") == 0 )
                        {
                            //Loop Over Quadrature Points
                            for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
                            {
                                //Loop Over i
                                for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                {
                                    local_residual(i) += - normal_velocity_suppressed_weakly * (1.0/eta[s]) * ( phi_velocity_bdry[s][q][i] * this->normal_vectors[q] ) * ( velocity_bdry_old[s][q] * this->normal_vectors[q] ) * this->JxW_bdry[q];
                                }
                            }
                        }
                        else //Unknown, throw error
                        {
                            std::cerr << "Boundary condition you specified at \"symmetry line or plane\" is not implemented" << std::endl;
                            AssertThrow( false , ExcNotImplemented() );
                        }
                }
                else if( iter->boundary_name.compare("inlet-outlet") == 0 )
                {
                    if( iter->boundary_condition.compare("Dirichlet density") == 0 )
                    {
                        int_1 = 0.0;
                        int_2 = 1.0;
                        int_3 = 1.0;
                        int_4 = 1.0;
                    }
                    else if( iter->boundary_condition.compare("Dirichlet density and normal stress free") == 0 )
                    {
                        int_1 = 0.0;
                        int_2 = 1.0;
                        int_3 = 0.0;
                        int_4 = 0.0;
                    }
                    else if( iter->boundary_condition.compare("Dirichlet density and normal shear stress free") == 0 )
                    {
                        int_1 = 0.0;
                        int_2 = 1.0;
                        int_3 = 1.0;
                        int_4 = 0.0;
                    }
                    else if( iter->boundary_condition.compare("Dirichlet velocity") == 0 )
                    {
                        int_1 = 1.0;
                        int_2 = 0.0;
                        int_3 = 0.0;
                        int_4 = 0.0;
                    }
                    else if( iter->boundary_condition.compare("Dirichlet density and velocity") == 0 )
                    {
                        int_1 = 0.0;
                        int_2 = 0.0;
                        int_3 = 0.0;
                        int_4 = 0.0;
                    }
                    else if( iter->boundary_condition.compare("normal stress free") == 0 )
                    {
                        int_1 = 1.0;
                        int_2 = 1.0;
                        int_3 = 0.0;
                        int_4 = 0.0;
                    }
                    else if( iter->boundary_condition.compare("normal shear stress free") == 0 )
                    {
                        int_1 = 1.0;
                        int_2 = 1.0;
                        int_3 = 1.0;
                        int_4 = 0.0;
                    }
                    else //Unknown, throw error
                    {
                        std::cerr << "Boundary condition you specified at \"inlet-outlet\" is not implemented" << std::endl;
                        AssertThrow( false , ExcNotImplemented() );
                    }
                    //Looping Over Quadrature Points
                    for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
                    {
                        //Looping Over i
                        for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                        {
                            local_residual(i) += - int_1 * phi_density_bdry[s][q][i]         *     mass_flux_bdry_old[s][q] * this->normal_vectors[q]                   * this->JxW_bdry[q]
                                                 - int_2 * n_BY_phi_velocity_bdry_S[s][q][i] *     momentum_flux_bdry_old[s][q]                                         * this->JxW_bdry[q]
                                                 + int_3 * n_BY_phi_velocity_bdry_S[s][q][i] * ( - pressure_bdry_old[s][q]*unit + int_4 * shear_stress_bdry_old[s][q] ) * this->JxW_bdry[q];
                        }
                    }
                }
            }
        }
    }
    
     this->dealII_to_appframe(bdry_residual,
                              local_residual,
                              this->residual_indices);
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
// LOCAL CG FEM BASED ASSEMBLERS - make_ FUNCTIONS //
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

// ---                                       ---
// --- make_assemblers_generic_constant_data ---
// ---                                       ---

template<int dim>
template<typename INFO>
void
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::make_assemblers_generic_constant_data(const INFO&                                InFo,
                                                                                               FuelCellShop::Layer::BaseLayer<dim>* const layer)
{
    try
    {
        //Types Info
        const std::type_info& Channel                 = typeid(FuelCellShop::Layer::Channel<dim>);
        const std::type_info& GasDiffusionLayer       = typeid(FuelCellShop::Layer::GasDiffusionLayer<dim>);
        const std::type_info& MicroPorousLayer        = typeid(FuelCellShop::Layer::MicroPorousLayer<dim>);
        const std::type_info& CatalystLayer           = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
        const std::type_info& info                    = layer->get_base_type();

        //Check what material in, get necessary properties
        if( info == Channel )
        {
            FuelCellShop::Layer::Channel<dim>* ptr = dynamic_cast< FuelCellShop::Layer::Channel<dim>* >(layer);
            update_gas_properties(ptr);
        }
        else if( info == GasDiffusionLayer || info == MicroPorousLayer || info == CatalystLayer )
        {
            FuelCellShop::Layer::PorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::PorousLayer<dim>* >(layer);
            update_gas_properties(ptr);
        }
        else //Unknown, throw error
        {
            std::cerr << "Layer you specified is not implemented" << std::endl;
            AssertThrow( false , ExcNotImplemented() );
        }

        AssertThrow( n_species == molar_mass.size()       , ExcDimensionMismatch(n_species, molar_mass.size())        );
        AssertThrow( n_species == dynamic_viscosity.size(), ExcDimensionMismatch(n_species, dynamic_viscosity.size()) );
            
        //Change units on various properties
        for(unsigned int s1 = 0; s1 < n_species; ++s1)
        {
            for(unsigned int s2 = 0; s2 < n_species; ++s2)
                maxwell_stefan_isobaric_diffusion_coefficient(s1,s2) *= 1.0e5; //Multiply by 10^5 * P where P = 1 Pascal
        }
    }
    catch(const std::bad_cast& e)
    {
        std::cerr << "Object is not of type FuelCellShop::Layer::CatalystLayer<dim>" << std::endl;
        std::cerr << e.what() << std::endl;
    }
    
    for(unsigned int s = 0; s < n_species; ++s)
    {
        double square_root;
        square_root = std::sqrt( (Constants::Pi()*molar_mass[s])/(2.0*R_universal*T_mixture) );
        
        maxwell_constant[s] = (2.0-maxwell_constant[s])*square_root/maxwell_constant[s];
    }

    //Setup density and velocity indice vectors, etc.
    std::vector<unsigned int> density_indices(n_species);
    std::vector<unsigned int> velocity_first_indices(n_species);

    for(unsigned int index = 0; index < n_species; ++index)
    {
        density_indices[index]        = index*(dim+1)          + this->system_management->solution_name_to_index("density_species_1");
        velocity_first_indices[index] = density_indices[index] + 1;
        
        density.push_back(  FEValuesExtractors::Scalar(density_indices[index]       ));
        velocity.push_back( FEValuesExtractors::Vector(velocity_first_indices[index]));
    }

    this->dofs_per_cell = InFo.get_fe_val_unsplit().dofs_per_cell;

    this->make_matrix_block_indices();
    this->make_residual_indices();
}

// ---                                    ---
// --- make_assemblers_cell_constant_data ---
// ---                                    ---

template<int dim>
void
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
    this->n_q_points_cell = cell_info.get_fe_val_unsplit().n_quadrature_points;

    //Resize vectors
    density_cell_old.resize       (n_species, std::vector<double>                  (this->n_q_points_cell));
    grad_density_cell_old.resize  (n_species, std::vector< Tensor<1,dim> >         (this->n_q_points_cell));
    velocity_cell_old.resize      (n_species, std::vector< Tensor<1,dim> >         (this->n_q_points_cell));
    div_velocity_cell_old.resize  (n_species, std::vector<double>                  (this->n_q_points_cell));
    grads_velocity_cell_old.resize(n_species, std::vector< SymmetricTensor<2,dim> >(this->n_q_points_cell));
    mass_flux_cell_old.resize     (n_species, std::vector< Tensor<1,dim> >         (this->n_q_points_cell));
    momentum_flux_cell_old.resize (n_species, std::vector< SymmetricTensor<2,dim> >(this->n_q_points_cell));
    pressure_cell_old.resize      (n_species, std::vector<double>                  (this->n_q_points_cell));
    shear_stress_cell_old.resize  (n_species, std::vector< SymmetricTensor<2,dim> >(this->n_q_points_cell));
    drag_cell_old.resize          (n_species, std::vector< Tensor<1,dim> >         (this->n_q_points_cell));
    diffusion_cell_old.resize     (n_species, std::vector< Tensor<1,dim> >         (this->n_q_points_cell));
    gravity_cell_old.resize       (n_species, std::vector< Tensor<1,dim> >         (this->n_q_points_cell));

    phi_density_cell.resize            (n_species, std::vector< std::vector<double> >                  (this->n_q_points_cell, std::vector<double>                  (this->dofs_per_cell)));
    grad_phi_density_cell.resize       (n_species, std::vector< std::vector< Tensor<1,dim> > >         (this->n_q_points_cell, std::vector< Tensor<1,dim> >         (this->dofs_per_cell)));
    phi_velocity_cell.resize           (n_species, std::vector< std::vector< Tensor<1,dim> > >         (this->n_q_points_cell, std::vector< Tensor<1,dim> >         (this->dofs_per_cell)));
    div_phi_velocity_cell.resize       (n_species, std::vector< std::vector<double> >                  (this->n_q_points_cell, std::vector<double>                  (this->dofs_per_cell)));
    grads_phi_velocity_cell.resize     (n_species, std::vector< std::vector< SymmetricTensor<2,dim> > >(this->n_q_points_cell, std::vector< SymmetricTensor<2,dim> >(this->dofs_per_cell)));
    delta_mass_flux_cell.resize        (n_species, std::vector< std::vector< Tensor<1,dim> > >         (this->n_q_points_cell, std::vector< Tensor<1,dim> >         (this->dofs_per_cell)));
    delta_momentum_flux_cell.resize    (n_species, std::vector< std::vector< SymmetricTensor<2,dim> > >(this->n_q_points_cell, std::vector< SymmetricTensor<2,dim> >(this->dofs_per_cell)));
    delta_pressure_cell.resize         (n_species, std::vector< std::vector<double> >                  (this->n_q_points_cell, std::vector<double>                  (this->dofs_per_cell)));
    delta_shear_stress_cell.resize     (n_species, std::vector< std::vector< SymmetricTensor<2,dim> > >(this->n_q_points_cell, std::vector< SymmetricTensor<2,dim> >(this->dofs_per_cell)));
    delta_partial_viscosity_cell.resize(n_species, std::vector< std::vector<double> >                  (this->n_q_points_cell, std::vector<double>                  (this->dofs_per_cell)));
    delta_bulk_viscosity_cell.resize   (n_species, std::vector< std::vector<double> >                  (this->n_q_points_cell, std::vector<double>                  (this->dofs_per_cell)));
    delta_drag_cell.resize             (n_species, std::vector< std::vector< Tensor<1,dim> > >         (this->n_q_points_cell, std::vector< Tensor<1,dim> >         (this->dofs_per_cell)));
    delta_diffusion_cell.resize        (n_species, std::vector< std::vector< Tensor<1,dim> > >         (this->n_q_points_cell, std::vector< Tensor<1,dim> >         (this->dofs_per_cell)));
    delta_gravity_cell.resize          (n_species, std::vector< std::vector< Tensor<1,dim> > >         (this->n_q_points_cell, std::vector< Tensor<1,dim> >         (this->dofs_per_cell)));

    this->JxW_cell.resize(this->n_q_points_cell);
}

// ---                                    ---
// --- make_assemblers_bdry_constant_data ---
// ---                                    ---

template<int dim>
void
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::make_assemblers_bdry_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info)
{
    this->n_q_points_bdry = bdry_info.get_fe_val_unsplit().n_quadrature_points;

    //Resize vectors
    density_bdry_old.resize       (n_species, std::vector<double>                  (this->n_q_points_bdry));
    grad_density_bdry_old.resize  (n_species, std::vector< Tensor<1,dim> >         (this->n_q_points_bdry));
    velocity_bdry_old.resize      (n_species, std::vector< Tensor<1,dim> >         (this->n_q_points_bdry));
    div_velocity_bdry_old.resize  (n_species, std::vector<double>                  (this->n_q_points_bdry));
    grads_velocity_bdry_old.resize(n_species, std::vector< SymmetricTensor<2,dim> >(this->n_q_points_bdry));
    mass_flux_bdry_old.resize     (n_species, std::vector< Tensor<1,dim> >         (this->n_q_points_bdry));
    momentum_flux_bdry_old.resize (n_species, std::vector< SymmetricTensor<2,dim> >(this->n_q_points_bdry));
    pressure_bdry_old.resize      (n_species, std::vector<double>                  (this->n_q_points_bdry));
    shear_stress_bdry_old.resize  (n_species, std::vector< SymmetricTensor<2,dim> >(this->n_q_points_bdry));

    phi_density_bdry.resize            (n_species, std::vector< std::vector<double> >                  (this->n_q_points_bdry, std::vector<double>                  (this->dofs_per_cell)));
    grad_phi_density_bdry.resize       (n_species, std::vector< std::vector< Tensor<1,dim> > >         (this->n_q_points_bdry, std::vector< Tensor<1,dim> >         (this->dofs_per_cell)));
    phi_velocity_bdry.resize           (n_species, std::vector< std::vector< Tensor<1,dim> > >         (this->n_q_points_bdry, std::vector< Tensor<1,dim> >         (this->dofs_per_cell)));
    div_phi_velocity_bdry.resize       (n_species, std::vector< std::vector<double> >                  (this->n_q_points_bdry, std::vector<double>                  (this->dofs_per_cell)));
    grads_phi_velocity_bdry.resize     (n_species, std::vector< std::vector< SymmetricTensor<2,dim> > >(this->n_q_points_bdry, std::vector< SymmetricTensor<2,dim> >(this->dofs_per_cell)));
    delta_mass_flux_bdry.resize        (n_species, std::vector< std::vector< Tensor<1,dim> > >         (this->n_q_points_bdry, std::vector< Tensor<1,dim> >         (this->dofs_per_cell)));
    delta_momentum_flux_bdry.resize    (n_species, std::vector< std::vector< SymmetricTensor<2,dim> > >(this->n_q_points_bdry, std::vector< SymmetricTensor<2,dim> >(this->dofs_per_cell)));
    delta_pressure_bdry.resize         (n_species, std::vector< std::vector<double> >                  (this->n_q_points_bdry, std::vector<double>                  (this->dofs_per_cell)));
    delta_shear_stress_bdry.resize     (n_species, std::vector< std::vector< SymmetricTensor<2,dim> > >(this->n_q_points_bdry, std::vector< SymmetricTensor<2,dim> >(this->dofs_per_cell)));
    delta_partial_viscosity_bdry.resize(n_species, std::vector< std::vector<double> >                  (this->n_q_points_bdry, std::vector<double>                  (this->dofs_per_cell)));
    delta_bulk_viscosity_bdry.resize   (n_species, std::vector< std::vector<double> >                  (this->n_q_points_bdry, std::vector<double>                  (this->dofs_per_cell)));

    this->JxW_bdry.resize(this->n_q_points_bdry);
    this->normal_vectors.resize(this->n_q_points_bdry);
    this->tangential_vectors.resize(dim-1, std::vector< Point<dim> >(this->n_q_points_bdry));
    
    n_BY_phi_velocity_bdry_S.resize(n_species, std::vector< std::vector< SymmetricTensor<2,dim> > >(this->n_q_points_bdry,
                                                                                                    std::vector< SymmetricTensor<2,dim> >(this->dofs_per_cell)));
}

// ---                                    ---
// --- make_assemblers_cell_variable_data ---
// ---                                    ---

template<int dim>
void
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                                            FuelCellShop::Layer::BaseLayer<dim>* const                               layer)
{
    //Clear and resize
    drag_cell_old.clear();
    drag_cell_old.resize(n_species, std::vector< Tensor<1,dim> >(this->n_q_points_cell));
    
    delta_drag_cell.clear();
    delta_drag_cell.resize(n_species, std::vector< std::vector< Tensor<1,dim> > >(this->n_q_points_cell, std::vector< Tensor<1,dim> >(this->dofs_per_cell)));
    
    //Set type info
    const std::type_info& Channel                 = typeid(FuelCellShop::Layer::Channel<dim>);
    const std::type_info& GasDiffusionLayer       = typeid(FuelCellShop::Layer::GasDiffusionLayer<dim>);
    const std::type_info& MicroPorousLayer        = typeid(FuelCellShop::Layer::MicroPorousLayer<dim>);
    const std::type_info& CatalystLayer           = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    const std::type_info& info                    = layer->get_base_type();

    //Create vectors for porous media properties
    std::vector<double>                   porosity                (this->n_q_points_cell);
    std::vector< SymmetricTensor<2,dim> > permeability_INV        (this->n_q_points_cell);
    std::vector< SymmetricTensor<2,dim> > Forchheimer_permeability(this->n_q_points_cell);
    std::vector< SymmetricTensor<2,dim> > tortuosity              (this->n_q_points_cell);
    std::vector< Table< 2, SymmetricTensor<2,dim> > > maxwell_stefan_EFFECTIVE_isobaric_diffusion_coefficient_INV( this->n_q_points_cell,
                                                                                                                   Table< 2, SymmetricTensor<2,dim> >(n_species, n_species) );
    
    //-----------------------------------------
    //
    //Set Values of previous Newton iteration
    //
    //-----------------------------------------
    for(unsigned int s = 0; s < n_species; ++s)
    {
        cell_info.get_fe_val_unsplit()[ density[s]  ].get_function_values             (cell_info.global_data->vector( cell_info.global_data->n_vectors()-1 ), density_cell_old[s]);
        cell_info.get_fe_val_unsplit()[ density[s]  ].get_function_gradients          (cell_info.global_data->vector( cell_info.global_data->n_vectors()-1 ), grad_density_cell_old[s]);
        cell_info.get_fe_val_unsplit()[ velocity[s] ].get_function_values             (cell_info.global_data->vector( cell_info.global_data->n_vectors()-1 ), velocity_cell_old[s]);
        cell_info.get_fe_val_unsplit()[ velocity[s] ].get_function_divergences        (cell_info.global_data->vector( cell_info.global_data->n_vectors()-1 ), div_velocity_cell_old[s]);
        cell_info.get_fe_val_unsplit()[ velocity[s] ].get_function_symmetric_gradients(cell_info.global_data->vector( cell_info.global_data->n_vectors()-1 ), grads_velocity_cell_old[s]);
    }
    
    //-----------------------------------------
    //
    //Set Values of current Newton iteration
    //
    //-----------------------------------------
    for(unsigned int s = 0; s < n_species; ++s)
    {
        for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
        {
            for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
            {
                phi_density_cell[s][q][k]        = cell_info.get_fe_val_unsplit()[ density[s] ].value(k,q);
                grad_phi_density_cell[s][q][k]   = cell_info.get_fe_val_unsplit()[ density[s] ].gradient(k,q);

                phi_velocity_cell[s][q][k]       = cell_info.get_fe_val_unsplit()[ velocity[s] ].value(k,q);
                div_phi_velocity_cell[s][q][k]   = cell_info.get_fe_val_unsplit()[ velocity[s] ].divergence(k,q);
                grads_phi_velocity_cell[s][q][k] = cell_info.get_fe_val_unsplit()[ velocity[s] ].symmetric_gradient(k,q);
            }
        }
    }
    
    // -- Get cell data from layer:
    try
    {
        //Set properties based on layer
        if( info == Channel )
        {
            FuelCellShop::Layer::Channel<dim>* ptr = dynamic_cast< FuelCellShop::Layer::Channel<dim>* >(layer);
            
            for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
            {
                porosity[q] = 1.0; // enforce that the channel has no porosity
                tortuosity[q] = unit;
                for(unsigned int s = 0; s < n_species; ++s)
                    for(unsigned int s1 = 0; s1 < n_species; ++s1)
                        maxwell_stefan_EFFECTIVE_isobaric_diffusion_coefficient_INV[q](s,s1) = ( 1.0 / (porosity[q]*maxwell_stefan_isobaric_diffusion_coefficient(s,s1)) ) * tortuosity[q];
            }
            
            update_partial_viscosities(ptr, this->n_q_points_cell, porosity, density_cell_old, paramMatrix_old, PInv_old, partial_viscosity_old, bulk_viscosity_old);
            update_delta_partial_viscosities(ptr, partial_viscosity_old, paramMatrix_old, PInv_old, porosity, phi_density_cell, density_cell_old, delta_partial_viscosity_cell);
            update_delta_bulk_viscosities(ptr, delta_partial_viscosity_cell, delta_bulk_viscosity_cell);
        }
        else if( info == GasDiffusionLayer || info == MicroPorousLayer ||  info == CatalystLayer )
        {
            FuelCellShop::Layer::PorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::PorousLayer<dim>* >(layer);
            
            update_porosity(ptr, cell_info, porosity);
            update_permeability(ptr, cell_info, permeability_INV, Forchheimer_permeability);
            update_inv_diffusion_coefficients(ptr, maxwell_stefan_EFFECTIVE_isobaric_diffusion_coefficient_INV);
            update_partial_viscosities(ptr, this->n_q_points_cell, porosity, density_cell_old, paramMatrix_old, PInv_old, partial_viscosity_old, bulk_viscosity_old);
            update_delta_partial_viscosities(ptr, partial_viscosity_old, paramMatrix_old, PInv_old, porosity, phi_density_cell, density_cell_old, delta_partial_viscosity_cell);
            update_delta_bulk_viscosities(ptr, delta_partial_viscosity_cell, delta_bulk_viscosity_cell);
        }
        else //Unknown, throw error
        {
            std::cerr << "Layer you specified is not implemented" << std::endl;
            AssertThrow( false , ExcNotImplemented() );
        }
    }
    catch(const std::bad_cast& e)
    {
        std::cerr << "Layer object is not of appropriate type in CompressibleMultiComponentKGEquationsCoupled<dim>::make_assemblers_cell_variable_data" << std::endl;
        std::cerr << e.what() << std::endl;
    }
    
    //Simple Functions
    for(unsigned int s = 0; s < n_species; ++s)
    {
        for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
        {
            // --- mass flux ---
            mass_flux_cell_old[s][q] = density_cell_old[s][q] * velocity_cell_old[s][q];

            // --- partial pressure ---
            pressure_cell_old[s][q]  = ( 1.0 / molar_mass[s] ) * density_cell_old[s][q] * R_universal * T_mixture;
        }

        //Set values dependent on layer
        if( info == Channel )
        {
            // --- momentum flux ---
            if( inertia_in_channels[s] )
            {
                for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                {
                    Tensor<2,dim> tmp;
                    outer_product(tmp, velocity_cell_old[s][q], velocity_cell_old[s][q]);
                    momentum_flux_cell_old[s][q] = density_cell_old[s][q] * tmp;
                }
            }

            // --- shear stress ---
            if( shear_stress_in_channels[s] )
            {
                for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                {
                    shear_stress_cell_old[s][q] = 2.0 * partial_viscosity_old[q][s] * grads_velocity_cell_old[s][q]
                                                  +
                                                  bulk_viscosity_old[q][s] * div_velocity_cell_old[s][q] * unit;
                }
            }

            // --- gravity force ---
            if( gravity_in_channels[s] )
            {
                for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                {
                    gravity_cell_old[s][q] = density_cell_old[s][q] * gravity_acceleration;
                }
            }

            // --- drag force ---
            // do nothing
        }
        else if( info == GasDiffusionLayer || info == MicroPorousLayer || info == CatalystLayer )
        {
                // --- momentum flux ---
                if( inertia_in_porous_media[s] )
                {
                    for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                    {
                        Tensor<2,dim> tmp;
                        outer_product(tmp, velocity_cell_old[s][q], velocity_cell_old[s][q]);
                        momentum_flux_cell_old[s][q] = density_cell_old[s][q] * tmp;
                    }
                }

                // --- shear stress ---
                if( shear_stress_in_porous_media[s] )
                {
                    for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                    {
                        shear_stress_cell_old[s][q] = 2.0 * partial_viscosity_old[q][s] * grads_velocity_cell_old[s][q]
                                                      +
                                                      bulk_viscosity_old[q][s] * div_velocity_cell_old[s][q] * unit;
                    }
                }

                // --- gravity force ---
                if( gravity_in_porous_media[s] )
                {
                    for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                    {
                        gravity_cell_old[s][q] = density_cell_old[s][q] * gravity_acceleration;
                    }
                }

                // --- drag force ---
                if(      drag_in_porous_media[s].compare("none") == 0 )
                {
                    // do nothing
                }
                else if( drag_in_porous_media[s].compare("Darcy") == 0 )
                {
                    for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                    {
                        drag_cell_old[s][q] = - partial_viscosity_old[q][s] * permeability_INV[q] * porosity[q] * velocity_cell_old[s][q];
                    }
                }
                else if( drag_in_porous_media[s].compare("Forchheimer") == 0 )
                {
                    for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                    {
                        drag_cell_old[s][q] = - partial_viscosity_old[q][s] * permeability_INV[q] * porosity[q] * velocity_cell_old[s][q]
                                              - density_cell_old[s][q] * porosity[q] * velocity_cell_old[s][q].norm() * Forchheimer_permeability[q] * porosity[q] * velocity_cell_old[s][q];
                    }
                }
                else if( drag_in_porous_media[s].compare("Forchheimer modified") == 0 )
                {
                    for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                    {

                    }
                }
                else
                {
                    std::cerr << "Form of the drag force you specified in porous media is not implemented" << std::endl;
                    AssertThrow( false , ExcNotImplemented() );
                }
        }
        else //Unknown, throw error
        {
            std::cerr << "Layer you specified is not implemented" << std::endl;
            AssertThrow( false , ExcNotImplemented() );
        }
    }
    
    //Complex Functions
    // --- diffusion force ---
    for(unsigned int s = 0; s < n_species; ++s)
        for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
        {
            Tensor<1,dim> sum;

            for(unsigned int s1 = 0; s1 < n_species; ++s1)
            {
                if( s1 != s )
                {
                    sum += pressure_cell_old[s][q] * pressure_cell_old[s1][q] * maxwell_stefan_EFFECTIVE_isobaric_diffusion_coefficient_INV[q](s,s1) * ( velocity_cell_old[s1][q] - velocity_cell_old[s][q] );
                }
            }
            diffusion_cell_old[s][q] = sum;
        }
    
    //-----------------------------------------
    //
    //Set delta values
    //
    //-----------------------------------------
    //Simple Functions
    for(unsigned int s = 0; s < n_species; ++s)
    {
        for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
        {
            for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
            {
                // --- DELTA mass flux ---
                delta_mass_flux_cell[s][q][k] = phi_density_cell[s][q][k] * velocity_cell_old[s][q]
                                                +
                                                density_cell_old[s][q]    * phi_velocity_cell[s][q][k];

                // --- DELTA partial pressure ---
                delta_pressure_cell[s][q][k]  = ( 1.0 / molar_mass[s] ) *   phi_density_cell[s][q][k] * R_universal * T_mixture;
            }
        }

        //Set values dependent on layer
        if( info == Channel )
        {
            // --- DELTA momentum flux ---
            if( inertia_in_channels[s] )
            {
                for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                {
                    Tensor<2,dim> tmp1;
                    outer_product(tmp1, velocity_cell_old[s][q], velocity_cell_old[s][q]);

                    for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
                    {
                        Tensor<2,dim> tmp2;
                        outer_product(tmp2, phi_velocity_cell[s][q][k], velocity_cell_old[s][q]);

                        Tensor<2,dim> tmp3;
                        outer_product(tmp3, velocity_cell_old[s][q], phi_velocity_cell[s][q][k]);

                        delta_momentum_flux_cell[s][q][k] = phi_density_cell[s][q][k] * tmp1
                                                            +
                                                            density_cell_old[s][q]    * tmp2
                                                            +
                                                            density_cell_old[s][q]    * tmp3;
                    }
                }
            }

            // --- DELTA shear stress ---
            if( shear_stress_in_channels[s] )
            {
                for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                {
                    for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
                    {
                        delta_shear_stress_cell[s][q][k] = 2.0 * delta_partial_viscosity_cell[s][q][k] * grads_velocity_cell_old[s][q]
                                                           +
                                                           delta_bulk_viscosity_cell[s][q][k] * div_velocity_cell_old[s][q] * unit
                                                           +
                                                           2.0 * partial_viscosity_old[q][s] * grads_phi_velocity_cell[s][q][k]
                                                           +
                                                           bulk_viscosity_old[q][s] * div_phi_velocity_cell[s][q][k] * unit;
                    }
                }
            }

            // --- DELTA gravity force ---
            if( gravity_in_channels[s] )
            {
                for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                {
                    for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
                    {
                            delta_gravity_cell[s][q][k] = phi_density_cell[s][q][k] * gravity_acceleration;
                    }
                }
            }

            // --- DELTA drag force ---
            // do nothing
        }
        else if( info == GasDiffusionLayer || info == MicroPorousLayer || info == CatalystLayer )
        {
            // --- DELTA momentum flux ---
            if( inertia_in_porous_media[s] )
            {
                for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                {
                    Tensor<2,dim> tmp1;
                    outer_product(tmp1, velocity_cell_old[s][q], velocity_cell_old[s][q]);

                        for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
                        {
                            Tensor<2,dim> tmp2;
                            outer_product(tmp2, phi_velocity_cell[s][q][k], velocity_cell_old[s][q]);

                            Tensor<2,dim> tmp3;
                            outer_product(tmp3, velocity_cell_old[s][q], phi_velocity_cell[s][q][k]);

                            delta_momentum_flux_cell[s][q][k] = phi_density_cell[s][q][k] * tmp1
                                                                +
                                                                density_cell_old[s][q]    * tmp2
                                                                +
                                                                density_cell_old[s][q]    * tmp3;
                        }
                }
            }

            // --- DELTA shear stress ---
            if( shear_stress_in_porous_media[s] )
            {
                for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                {
                    for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
                    {
                        delta_shear_stress_cell[s][q][k] = 2.0 * delta_partial_viscosity_cell[s][q][k] * grads_velocity_cell_old[s][q]
                                                           +
                                                           delta_bulk_viscosity_cell[s][q][k] * div_velocity_cell_old[s][q] * unit
                                                           +
                                                           2.0 * partial_viscosity_old[q][s] * grads_phi_velocity_cell[s][q][k]
                                                           +
                                                           bulk_viscosity_old[q][s] * div_phi_velocity_cell[s][q][k] * unit;
                    }
                }
            }

            // --- DELTA gravity force ---
            if( gravity_in_porous_media[s] )
            {
                for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                    for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
                    {
                        delta_gravity_cell[s][q][k] = phi_density_cell[s][q][k] * gravity_acceleration;
                    }
            }

            // --- DELTA drag force ---
            if(      drag_in_porous_media[s].compare("none") == 0 )
            {
                // do nothing
            }
            else if( drag_in_porous_media[s].compare("Darcy") == 0 )
            {
                for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                    for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
                    {
                        delta_drag_cell[s][q][k] = - ( delta_partial_viscosity_cell[s][q][k] * permeability_INV[q] * porosity[q] * velocity_cell_old[s][q] 
                                                       +
                                                       partial_viscosity_old[q][s] * permeability_INV[q] * porosity[q] * phi_velocity_cell[s][q][k]);
                    }
            }
            else if( drag_in_porous_media[s].compare("Forchheimer") == 0 )
            {
                for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                    for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
                    {
                        Tensor<1,dim> tmp1;
                        tmp1  = - ( delta_partial_viscosity_cell[s][q][k] * permeability_INV[q] * porosity[q] * velocity_cell_old[s][q] 
                                    +
                                    partial_viscosity_old[q][s] * permeability_INV[q] * porosity[q] * phi_velocity_cell[s][q][k]);

                        Tensor<1,dim> tmp2;
                        tmp2  = - phi_density_cell[s][q][k] * porosity[q] * velocity_cell_old[s][q].norm() * porosity[q] * velocity_cell_old[s][q];

                        double        ratio;
                        ratio =   ( porosity[q] * velocity_cell_old[s][q] * porosity[q] * phi_velocity_cell[s][q][k] ) / ( porosity[q] * velocity_cell_old[s][q].norm() );

                        Tensor<1,dim> tmp3;
                        tmp3  = - density_cell_old[s][q] * ratio * porosity[q] * velocity_cell_old[s][q];

                        Tensor<1,dim> tmp4;
                        tmp4  = - density_cell_old[s][q] * porosity[q] * velocity_cell_old[s][q].norm() * porosity[q] * phi_velocity_cell[s][q][k];

                        delta_drag_cell[s][q][k] = tmp1 + Forchheimer_permeability[q] * ( tmp2 + tmp3 + tmp4 );
                    }
            }
            else if( drag_in_porous_media[s].compare("Forchheimer modified") == 0 )
            {
                for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                    for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
                    {  }
            }
            else
                AssertThrow( false , ExcMessage("Form of the drag force you specified in porous media is not implemented") );
        }
        else //Unknown, throw error
            AssertThrow( false , ExcMessage("Layer you specified is not implemented") );
    }

    //Complex Functions
    // --- DELTA diffusion force ---
    for(unsigned int s = 0; s < n_species; ++s)
    {
        for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
        {
            for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
            {
                Tensor<1,dim> sum1;
                Tensor<1,dim> sum2;

                for(unsigned int s1 = 0; s1 < n_species; ++s1)
                {
                    if( s1 != s )
                    {
                        sum1 += ( delta_pressure_cell[s][q][k] * pressure_cell_old[s1][q] + pressure_cell_old[s][q] * delta_pressure_cell[s1][q][k] )
                                *
                                maxwell_stefan_EFFECTIVE_isobaric_diffusion_coefficient_INV[q](s,s1) * ( velocity_cell_old[s1][q] - velocity_cell_old[s][q] );

                        sum2 += pressure_cell_old[s][q] * pressure_cell_old[s1][q]
                                *
                                maxwell_stefan_EFFECTIVE_isobaric_diffusion_coefficient_INV[q](s,s1) * ( phi_velocity_cell[s1][q][k] - phi_velocity_cell[s][q][k] );
                    }
                }

                delta_diffusion_cell[s][q][k] = sum1 + sum2;
            }
        }
    }

    //Other
    // --- JxW ---
    for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
        this->JxW_cell[q] = cell_info.get_fe_val_unsplit().JxW(q);
}

// ---                                    ---
// --- make_assemblers_bdry_variable_data ---
// ---                                    ---

template<int dim>
void
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::make_assemblers_bdry_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                                                            FuelCellShop::Layer::BaseLayer<dim>* const                               layer)
{
    //Set type info
    const std::type_info& Channel                 = typeid(FuelCellShop::Layer::Channel<dim>);
    const std::type_info& GasDiffusionLayer       = typeid(FuelCellShop::Layer::GasDiffusionLayer<dim>);
    const std::type_info& MicroPorousLayer        = typeid(FuelCellShop::Layer::MicroPorousLayer<dim>);
    const std::type_info& CatalystLayer           = typeid(FuelCellShop::Layer::CatalystLayer<dim>);

    const std::type_info& info                    = layer->get_base_type();
    
    //Create vector for porous media properties
    std::vector<double> porosity(this->n_q_points_bdry);

    //-----------------------------------------
    //
    //Set Values of previous Newton iteration
    //
    //-----------------------------------------
    for(unsigned int s = 0; s < n_species; ++s)
    {
        bdry_info.get_fe_val_unsplit()[ density[s]  ].get_function_values             (bdry_info.global_data->vector(bdry_info.global_data->n_vectors()-1 ), density_bdry_old[s]       );
        bdry_info.get_fe_val_unsplit()[ density[s]  ].get_function_gradients          (bdry_info.global_data->vector(bdry_info.global_data->n_vectors()-1 ), grad_density_bdry_old[s]  );
        bdry_info.get_fe_val_unsplit()[ velocity[s] ].get_function_values             (bdry_info.global_data->vector(bdry_info.global_data->n_vectors()-1 ), velocity_bdry_old[s]      );
        bdry_info.get_fe_val_unsplit()[ velocity[s] ].get_function_divergences        (bdry_info.global_data->vector(bdry_info.global_data->n_vectors()-1 ), div_velocity_bdry_old[s]  );
        bdry_info.get_fe_val_unsplit()[ velocity[s] ].get_function_symmetric_gradients(bdry_info.global_data->vector(bdry_info.global_data->n_vectors()-1 ), grads_velocity_bdry_old[s]);
    }

    //-----------------------------------------
    //
    //Set Values of current Newton iteration
    //
    //-----------------------------------------
    for(unsigned int s = 0; s < n_species; ++s)
    {
        for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
            for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
            {
                phi_density_bdry[s][q][k]        = bdry_info.get_fe_val_unsplit()[ density[s] ].value(k,q);
                grad_phi_density_bdry[s][q][k]   = bdry_info.get_fe_val_unsplit()[ density[s] ].gradient(k,q);

                phi_velocity_bdry[s][q][k]       = bdry_info.get_fe_val_unsplit()[ velocity[s] ].value(k,q);
                div_phi_velocity_bdry[s][q][k]   = bdry_info.get_fe_val_unsplit()[ velocity[s] ].divergence(k,q);
                grads_phi_velocity_bdry[s][q][k] = bdry_info.get_fe_val_unsplit()[ velocity[s] ].symmetric_gradient(k,q);
            }
    }
    
    try
    {
        //Set properties based on layer
        if( info == Channel )
        {
            FuelCellShop::Layer::Channel<dim>* ptr = dynamic_cast< FuelCellShop::Layer::Channel<dim>* >(layer);
            
            for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
                porosity[q] = 1.0; // enforce that the channel has no porosity

            update_partial_viscosities(ptr, this->n_q_points_bdry, porosity, density_bdry_old, paramMatrix_bdry_old, PInv_bdry_old, partial_viscosity_bdry_old, bulk_viscosity_bdry_old);
            update_delta_partial_viscosities(ptr, partial_viscosity_bdry_old, paramMatrix_bdry_old, PInv_bdry_old, porosity, phi_density_bdry, density_bdry_old, delta_partial_viscosity_bdry);
            update_delta_bulk_viscosities(ptr, delta_partial_viscosity_bdry, delta_bulk_viscosity_bdry);
        }
        else if( info == GasDiffusionLayer || info == MicroPorousLayer || info == CatalystLayer )
        {
            FuelCellShop::Layer::PorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::PorousLayer<dim>* >(layer);
            
            update_porosity(ptr, bdry_info, porosity);
            update_partial_viscosities(ptr, this->n_q_points_bdry, porosity, density_bdry_old, paramMatrix_bdry_old, PInv_bdry_old, partial_viscosity_bdry_old, bulk_viscosity_bdry_old);
            update_delta_partial_viscosities(ptr, partial_viscosity_bdry_old, paramMatrix_bdry_old, PInv_bdry_old, porosity, phi_density_bdry, density_bdry_old, delta_partial_viscosity_bdry);
            update_delta_bulk_viscosities(ptr, delta_partial_viscosity_bdry, delta_bulk_viscosity_bdry);
        }
        else //Unknown, throw error
        {
            std::cerr << "Layer you specified is not implemented" << std::endl;
            AssertThrow( false , ExcNotImplemented() );
        }
    }
    catch(const std::bad_cast& e)
    {
        std::cerr << "Layer object is not of appropriate type in CompressibleMultiComponentKGEquationsCoupled<dim>::make_assemblers_bdry_variable_data" << std::endl;
        std::cerr << e.what() << std::endl;
    }
    
    //Simple functions
    for(unsigned int s = 0; s < n_species; ++s)
    {
        for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
        {
            // --- mass flux ---
            mass_flux_bdry_old[s][q] = density_bdry_old[s][q] * velocity_bdry_old[s][q];

            // --- partial pressure ---
            pressure_bdry_old[s][q]  = ( 1.0 / molar_mass[s] ) * density_bdry_old[s][q] * R_universal * T_mixture;
        }

        //Set values dependent on layer
        if( ( (info == Channel) && inertia_in_channels[s] ) || 
            ( ( info == GasDiffusionLayer || info == MicroPorousLayer || info == CatalystLayer ) &&   inertia_in_porous_media[s] ) )
        {
            // --- momentum flux ---
            for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
                {
                    Tensor<2,dim> tmp;
                    outer_product(tmp, velocity_bdry_old[s][q], velocity_bdry_old[s][q]);
                    momentum_flux_bdry_old[s][q] = density_bdry_old[s][q] * tmp;
                }
        }
        else if ( ( (info == Channel) && shear_stress_in_channels[s] ) ||
            ( ( info == GasDiffusionLayer || info == MicroPorousLayer || info == CatalystLayer ) &&   shear_stress_in_porous_media[s] ) )
        {
            for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
            {
                shear_stress_bdry_old[s][q] = 2.0 * partial_viscosity_bdry_old[q][s] * grads_velocity_bdry_old[s][q]
                                              +
                                              bulk_viscosity_bdry_old[q][s] * div_velocity_bdry_old[s][q] * unit;
            }
        }
        else if ( info != Channel && info != GasDiffusionLayer && info != MicroPorousLayer && info != CatalystLayer )
            AssertThrow( false , ExcMessage("Layer you specified is not implemented in CompressibleMultiComponentKGEquationsCoupled<dim>::make_assemblers_bdry_variable_data") );
    }
    
    //-----------------------------------------
    //
    //Set delta values
    //
    //-----------------------------------------
    //Simple functions
    for(unsigned int s = 0; s < n_species; ++s)
    {
        for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
        {
            for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
            {
                // --- DELTA mass flux ---
                delta_mass_flux_bdry[s][q][k] = phi_density_bdry[s][q][k] * velocity_bdry_old[s][q]
                                                +
                                                density_bdry_old[s][q]    * phi_velocity_bdry[s][q][k];

                // --- DELTA partial pressure ---
                delta_pressure_bdry[s][q][k]  = ( 1.0 / molar_mass[s] ) * phi_density_bdry[s][q][k] * R_universal * T_mixture;
            }
        }

        // --- DELTA momentum flux
        if ( ( info == Channel &&  inertia_in_channels[s] ) ||
            ( ( info == GasDiffusionLayer || info == MicroPorousLayer || info == CatalystLayer ) && inertia_in_porous_media[s] ) )
        {
            for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
            {
                Tensor<2,dim> tmp1;
                outer_product(tmp1, velocity_bdry_old[s][q], velocity_bdry_old[s][q]);
                
                for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
                {
                    Tensor<2,dim> tmp2;
                    outer_product(tmp2, phi_velocity_bdry[s][q][k], velocity_bdry_old[s][q]);
                    
                    Tensor<2,dim> tmp3;
                    outer_product(tmp3, velocity_bdry_old[s][q], phi_velocity_bdry[s][q][k]);
                    
                    delta_momentum_flux_bdry[s][q][k] = phi_density_bdry[s][q][k] * tmp1
                                                        +
                                                        density_bdry_old[s][q]    * tmp2
                                                        +
                                                        density_bdry_old[s][q]    * tmp3;
                }
            }
        }
        // --- DELTA shear stress ---
        else if ( ( info == Channel && shear_stress_in_channels[s] ) ||
            ( ( info == GasDiffusionLayer || info == MicroPorousLayer || info == CatalystLayer ) && shear_stress_in_porous_media[s] ) )
        {
            for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
            {
                for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
                {
                    delta_shear_stress_bdry[s][q][k] = 2.0 * delta_partial_viscosity_bdry[s][q][k] * grads_velocity_bdry_old[s][q]
                                                       +
                                                       delta_bulk_viscosity_bdry[s][q][k] * div_velocity_bdry_old[s][q] * unit
                                                       +
                                                       2.0 * partial_viscosity_bdry_old[q][s] * grads_phi_velocity_bdry[s][q][k]
                                                       +
                                                       bulk_viscosity_bdry_old[q][s] * div_phi_velocity_bdry[s][q][k] * unit;
                }
            }
        }
        /*
        else //Unknown, throw error
        {
            std::cerr << "Layer you specified is not implemented" << std::endl;
            AssertThrow( false , ExcNotImplemented() );
        }
        */
    }
    
    // --- JxW ---
    for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
        this->JxW_bdry[q] = bdry_info.get_fe_val_unsplit().JxW(q);

    // --- normal vectors ---
    this->normal_vectors = bdry_info.get_fe_val_unsplit().get_normal_vectors();

    // --- tangential vectors ---
    FemExtras::get_tangential_vectors(this->tangential_vectors, this->normal_vectors);

    // --- symmetrized tensor product ---
    for(unsigned int s = 0; s < n_species; ++s)
    {
        for(unsigned int q = 0; q < this->n_q_points_bdry; ++q)
        {
            for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
            {
                Tensor<2,dim> n_BY_phi_velocity_bdry;
                outer_product(n_BY_phi_velocity_bdry, this->normal_vectors[q], phi_velocity_bdry[s][q][k]);

                n_BY_phi_velocity_bdry_S[s][q][k] = 0.5 * ( n_BY_phi_velocity_bdry + transpose(n_BY_phi_velocity_bdry) );
            }
        }
    }
}

//-----------------------------------------
//
//Other Functions
//
//-----------------------------------------

// ---                              ---
// --- make_internal_cell_couplings ---
// ---                              ---

template<int dim>
void
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::make_internal_cell_couplings()
{
    //Create string of equation name
    std::string eq_generic_prefix = "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - ";

    //Setup equation and variable string names
    std::vector<std::string> eq_postfixes;
    std::vector<std::string> var_postfixes;
    for(unsigned int index = 1; index <= n_species; ++index)
    {
        std::ostringstream streamOut;
        streamOut << index;
        std::string name1 = "species " + streamOut.str();
        std::string name2 = "species_" + streamOut.str();

        eq_postfixes.push_back ( name1.c_str() );
        var_postfixes.push_back( name2.c_str() );
    }

    //Initiate map between variable name and DOF coupling
    std::map< std::string, DoFTools::Coupling > tmp;
    
    //NOTE: The 6 lines below are necessary for coupling with fuel cell physics HOWEVER, cannot be contained 
    //      in if statement else errors in next if statement due to variables not existing
    //Initiate vector of string solution names
    const std::vector<std::string> solution_names = this->system_management->get_solution_names();

    bool Fs = false;
    bool Fm = false;
    bool Lm = false;

    std::vector<std::string>::const_iterator iter;

    if(coupled_with_fuel_cell_physics)
    {
        // --- Fs ---
        iter = std::find(solution_names.begin(), solution_names.end(), "electronic_electrical_potential");

        if( iter != solution_names.end() )
            Fs = true;

        // --- Fm ---
        iter = std::find(solution_names.begin(), solution_names.end(), "protonic_electrical_potential");

        if( iter != solution_names.end() )
            Fm = true;

        // --- Lm ---
        iter = std::find(solution_names.begin(), solution_names.end(), "membrane_water_content");

        if( iter != solution_names.end() )
            Lm = true;
    }

    //Create internal cell couplings
    for(unsigned int s = 1; s <= n_species; ++s)
    {
        //Create string to hold full equation name to push back into internal_cell_couplings
        std::string eq_name;
        
        //Mass Conservation
        eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[s-1];
        for(unsigned int s1 = 1; s1 <= n_species; ++s1)
        {
            std::string var_name_dens = "density_" + var_postfixes[s1-1];
            std::string var_name_velX = "velocity_" + var_postfixes[s1-1] + "_X";
            std::string var_name_velY = "velocity_" + var_postfixes[s1-1] + "_Y";
            #if deal_II_dimension == 3
                std::string var_name_velZ = "velocity_" + var_postfixes[s1-1] + "_Z";
            #endif
            
            //Set internal map for each variable of this equation
            if( s1 == s )
            {
                tmp[var_name_dens] = DoFTools::always;
                tmp[var_name_velX] = DoFTools::always;
                tmp[var_name_velY] = DoFTools::always;
                #if deal_II_dimension == 3
                    tmp[var_name_velZ] = DoFTools::always;
                #endif
            }
            else
            {
                tmp[var_name_dens] = DoFTools::none;
                tmp[var_name_velX] = DoFTools::none;
                tmp[var_name_velY] = DoFTools::none;
                #if deal_II_dimension == 3
                    tmp[var_name_velZ] = DoFTools::none;
                #endif
            }
        }
        this->internal_cell_couplings[eq_name] = tmp;
        tmp.clear();
        
        //Momentum Conservation X
        eq_name = eq_generic_prefix + "momentum conservation X - " + eq_postfixes[s-1];
        for(unsigned int s1 = 1; s1 <= n_species; ++s1)
        {
            std::string var_name_dens = "density_" + var_postfixes[s1-1];
            std::string var_name_velX = "velocity_" + var_postfixes[s1-1] + "_X";
            std::string var_name_velY = "velocity_" + var_postfixes[s1-1] + "_Y";
            #if deal_II_dimension == 3
                std::string var_name_velZ = "velocity_" + var_postfixes[s1-1] + "_Z";
            #endif

            //Set internal map for each variable of this equation
            tmp[var_name_dens] = DoFTools::always;
            tmp[var_name_velX] = DoFTools::always;
            if( s1 == s )
            {
                tmp[var_name_velY] = DoFTools::always;
                #if deal_II_dimension == 3
                    tmp[var_name_velZ] = DoFTools::always;
                #endif
            }
            else
            {
                tmp[var_name_velY] = DoFTools::none;
                #if deal_II_dimension == 3
                    tmp[var_name_velZ] = DoFTools::none;
                #endif
            }
        }
        this->internal_cell_couplings[eq_name] = tmp;
        tmp.clear();
        
        //Momentum Conservation Y
        eq_name = eq_generic_prefix + "momentum conservation Y - " + eq_postfixes[s-1];
        for(unsigned int s1 = 1; s1 <= n_species; ++s1)
        {
            std::string var_name_dens = "density_" + var_postfixes[s1-1];
            std::string var_name_velX = "velocity_" + var_postfixes[s1-1] + "_X";
            std::string var_name_velY = "velocity_" + var_postfixes[s1-1] + "_Y";
            #if deal_II_dimension == 3
                std::string var_name_velZ = "velocity_" + var_postfixes[s1-1] + "_Z";
            #endif
            
            //Set internal map for each variable of this equation
            tmp[var_name_dens] = DoFTools::always;
            tmp[var_name_velY] = DoFTools::always;
            if( s1 == s )
            {
                tmp[var_name_velX] = DoFTools::always;
                #if deal_II_dimension == 3
                    tmp[var_name_velZ] = DoFTools::always;
                #endif
            }
            else
            {
                tmp[var_name_velX] = DoFTools::none;
                #if deal_II_dimension == 3
                    tmp[var_name_velZ] = DoFTools::none;
                #endif
            }
        }
        this->internal_cell_couplings[eq_name] = tmp;
        tmp.clear();
        
        #if deal_II_dimension == 3
            //Momentum Conservation Z
            eq_name = eq_generic_prefix + "momentum conservation Z - " + eq_postfixes[s-1];
            for(unsigned int s1 = 1; s1 <= n_species; ++s1)
            {
                std::string var_name_dens = "density_" + var_postfixes[s1-1];
                std::string var_name_velX = "velocity_" + var_postfixes[s1-1] + "_X";
                std::string var_name_velY = "velocity_" + var_postfixes[s1-1] + "_Y";
                std::string var_name_velZ = "velocity_" + var_postfixes[s1-1] + "_Z";
                
                //Set internal map for each variable of this equation
                tmp[var_name_dens] = DoFTools::always;
                tmp[var_name_velZ] = DoFTools::always;
                if( s1 == s )
                {
                    tmp[var_name_velX] = DoFTools::always;
                    tmp[var_name_velY] = DoFTools::always;
                }
                else
                {
                    tmp[var_name_velX] = DoFTools::none;
                    tmp[var_name_velY] = DoFTools::none;
                }
            }
            this->internal_cell_couplings[eq_name] = tmp;
            tmp.clear();
        #endif
    }

    if(coupled_with_fuel_cell_physics)
    {
        //////////////////////////////////////////////////////
        // POST-PROCESSING ON this->internal_cell_couplings //
        //////////////////////////////////////////////////////

//         tmp.clear();
        for(couplings_map::iterator iter  = this->internal_cell_couplings.begin();
                                    iter != this->internal_cell_couplings.end();
                                    ++iter)
        {
            tmp = iter->second;

            if( Fs )
                tmp["electronic_electrical_potential"] = DoFTools::none;

            if( Fm )
                tmp["protonic_electrical_potential"]   = DoFTools::none;

            if( Lm )
                tmp["membrane_water_content"]          = DoFTools::none;

            this->internal_cell_couplings[iter->first] = tmp;
            tmp.clear();
        }
    }
}

// ---                           ---
// --- make_matrix_block_indices ---
// ---                           ---

template<int dim>
void
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::make_matrix_block_indices()
{
    //Create string of equation name
    std::string eq_generic_prefix = "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - ";

    //Setup equation and variable string names
    std::vector<std::string> eq_postfixes;
    std::vector<std::string> var_postfixes;
    for(unsigned int index = 1; index <= n_species; ++index)
    {
        std::ostringstream streamOut;
        streamOut << index;
        std::string name1 = "species " + streamOut.str();
        std::string name2 = "species_" + streamOut.str();

        eq_postfixes.push_back ( name1.c_str() );
        var_postfixes.push_back( name2.c_str() );
    }

    //Initiate integer for indices of block matrix
    unsigned int index;
    
    //Create matrix block indices
    for(unsigned int s = 1; s <= n_species; ++s)
    {
        //Create string to hold full equation name to push back into internal_cell_couplings
        std::string eq_name;
        
        //Mass Conservation
        eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[s-1];
        for(unsigned int s1 = 1; s1 <= n_species; ++s1)
        {
            std::string var_name_dens = "density_" + var_postfixes[s1-1];
            std::string var_name_velX = "velocity_" + var_postfixes[s1-1] + "_X";
            std::string var_name_velY = "velocity_" + var_postfixes[s1-1] + "_Y";
            #if deal_II_dimension == 3
                std::string var_name_velZ = "velocity_" + var_postfixes[s1-1] + "_Z";
            #endif
                
            //Set TODO for each variable of this equation
            if( s1 == s )
            {
                //For density
                index = this->system_management->matrix_block_index(eq_name, var_name_dens);
                this->matrix_block_indices.push_back(index); 
                
                //For velocity X
                index = this->system_management->matrix_block_index(eq_name, var_name_velX);
                this->matrix_block_indices.push_back(index);
                
                //For velocity Y
                index = this->system_management->matrix_block_index(eq_name, var_name_velY);
                this->matrix_block_indices.push_back(index);
                
                //For velcoity Z
                #if deal_II_dimension == 3
                    index = this->system_management->matrix_block_index(eq_name, var_name_velZ);
                    this->matrix_block_indices.push_back(index);
                #endif
            }
        }
            
        //Momentum Conservation X
        eq_name = eq_generic_prefix + "momentum conservation X - " + eq_postfixes[s-1];
        for(unsigned int s1 = 1; s1 <= n_species; ++s1)
        {
            std::string var_name_dens = "density_" + var_postfixes[s1-1];
            std::string var_name_velX = "velocity_" + var_postfixes[s1-1] + "_X";
            std::string var_name_velY = "velocity_" + var_postfixes[s1-1] + "_Y";
            #if deal_II_dimension == 3
                std::string var_name_velZ = "velocity_" + var_postfixes[s1-1] + "_Z";
            #endif
            
            //Set TODO for each variable of this equation
            
            //For density
            index = this->system_management->matrix_block_index(eq_name, var_name_dens);
            this->matrix_block_indices.push_back(index); 
            
            //For velocity X
            index = this->system_management->matrix_block_index(eq_name, var_name_velX);
            this->matrix_block_indices.push_back(index);
            
            if( s1 == s )
            {
                //For velocity Y
                index = this->system_management->matrix_block_index(eq_name, var_name_velY);
                this->matrix_block_indices.push_back(index);
                
                //For velocity Z
                #if deal_II_dimension == 3
                    index = this->system_management->matrix_block_index(eq_name, var_name_velZ);
                    this->matrix_block_indices.push_back(index);
                #endif
            }
        }
        
        //Momentum Conservation Y
        eq_name = eq_generic_prefix + "momentum conservation Y - " + eq_postfixes[s-1];
        for(unsigned int s1 = 1; s1 <= n_species; ++s1)
        {
            std::string var_name_dens = "density_" + var_postfixes[s1-1];
            std::string var_name_velX = "velocity_" + var_postfixes[s1-1] + "_X";
            std::string var_name_velY = "velocity_" + var_postfixes[s1-1] + "_Y";
            #if deal_II_dimension == 3
                std::string var_name_velZ = "velocity_" + var_postfixes[s1-1] + "_Z";
            #endif
            
            //Set TODO for each variable of this equation
            
            //For density
            index = this->system_management->matrix_block_index(eq_name, var_name_dens);
            this->matrix_block_indices.push_back(index); 
            
            //For velocity Y
            index = this->system_management->matrix_block_index(eq_name, var_name_velY);
            this->matrix_block_indices.push_back(index);
            
            if( s1 == s )
            {
                //For velocity X
                index = this->system_management->matrix_block_index(eq_name, var_name_velX);
                this->matrix_block_indices.push_back(index);
                
                //For velocity Z
                #if deal_II_dimension == 3
                    index = this->system_management->matrix_block_index(eq_name, var_name_velZ);
                    this->matrix_block_indices.push_back(index);
                #endif
            }
        }
        
        //Momentum Conservation Z
        #if deal_II_dimension == 3
            eq_name = eq_generic_prefix + "momentum conservation Z - " + eq_postfixes[s-1];
            for(unsigned int s1 = 1; s1 <= n_species; ++s1)
            {
                std::string var_name_dens = "density_" + var_postfixes[s1-1];
                std::string var_name_velX = "velocity_" + var_postfixes[s1-1] + "_X";
                std::string var_name_velY = "velocity_" + var_postfixes[s1-1] + "_Y";
                std::string var_name_velZ = "velocity_" + var_postfixes[s1-1] + "_Z";
                
                //Set TODO for each variable of this equation
                
                //For density
                index = this->system_management->matrix_block_index(eq_name, var_name_dens);
                this->matrix_block_indices.push_back(index); 
                
                if( s1 == s )
                {
                    //For velocity X
                    index = this->system_management->matrix_block_index(eq_name, var_name_velX);
                    this->matrix_block_indices.push_back(index);
                    
                    //For velocity Y
                    index = this->system_management->matrix_block_index(eq_name, var_name_velY);
                    this->matrix_block_indices.push_back(index);
                }
                //For velocity Z
                index = this->system_management->matrix_block_index(eq_name, var_name_velZ);
                this->matrix_block_indices.push_back(index);
            }
        #endif
    }
}

// ---                       ---
// --- make_residual_indices ---
// ---                       ---

template<int dim>
void
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::make_residual_indices()
{
    //Create string of equation name
    std::string eq_generic_prefix = "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - ";

    //Setup equation and variable string names
    std::vector<std::string> eq_postfixes;
    for(unsigned int index = 1; index <= n_species; ++index)
    {
        std::ostringstream streamOut;
        streamOut << index;
        std::string name = "species " + streamOut.str();

        eq_postfixes.push_back( name.c_str() );
    }

    //Initiate integer for indices of residual
    unsigned int index;

    //Create residual indices
    for(unsigned int s = 1; s <= n_species; ++s)
    {
        //Mass Conservation
        std::string eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[s-1];
        index = this->system_management->equation_name_to_index(eq_name);
        this->residual_indices.push_back(index);
        
        //Momentum Conservation X
        eq_name = eq_generic_prefix + "momentum conservation X - " + eq_postfixes[s-1];
        index = this->system_management->equation_name_to_index(eq_name);
        this->residual_indices.push_back(index);
        
        //Momentum Conservation Y
        eq_name = eq_generic_prefix + "momentum conservation Y - " + eq_postfixes[s-1];
        index = this->system_management->equation_name_to_index(eq_name);
        this->residual_indices.push_back(index);
        
        //Momentum Conservation Z
        #if deal_II_dimension == 3
            eq_name = eq_generic_prefix + "momentum conservation Z - " + eq_postfixes[s-1];
            index = this->system_management->equation_name_to_index(eq_name);
            this->residual_indices.push_back(index);
        #endif
    }
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
NAME::CompressibleMultiComponentKGEquationsCoupled<dim>::print_equation_info() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "------------------------------";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "Parameters for Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component:";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;
    
    if( !this->counter[3] ) //Channel
    {
        FcstUtilities::log << "Boolean constants:";
        FcstUtilities::log << std::endl;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Inertia in channels:";
        FcstUtilities::log << std::endl;
        for(unsigned int index = 0; index < inertia_in_channels.size(); ++index)
            FcstUtilities::log << inertia_in_channels[index] << std::endl;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Shear stress in channels:";
        FcstUtilities::log << std::endl;
        for(unsigned int index = 0; index < shear_stress_in_channels.size(); ++index)
            FcstUtilities::log << shear_stress_in_channels[index] << std::endl;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Gravity in channels:";
        FcstUtilities::log << std::endl;
        for(unsigned int index = 0; index < gravity_in_channels.size(); ++index)
            FcstUtilities::log << gravity_in_channels[index] << std::endl;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Normal velocity is suppressed weakly:";
        FcstUtilities::log << std::endl;
        for(unsigned int index = 0; index < normal_velocity_is_suppressed_weakly.size(); ++index)
            FcstUtilities::log << normal_velocity_is_suppressed_weakly[index] << std::endl;
        FcstUtilities::log << std::endl;
    }
    else //Channel and Porous Layer
    {
        FcstUtilities::log << "Boolean constants and form of the drag force in porous media:";
        FcstUtilities::log << std::endl;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Inertia in channels:";
        FcstUtilities::log << std::endl;
        for(unsigned int index = 0; index < inertia_in_channels.size(); ++index)
            FcstUtilities::log << inertia_in_channels[index] << std::endl;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Shear stress in channels:";
        FcstUtilities::log << std::endl;
        for(unsigned int index = 0; index < shear_stress_in_channels.size(); ++index)
            FcstUtilities::log << shear_stress_in_channels[index] << std::endl;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Gravity in channels:";
        FcstUtilities::log << std::endl;
        for(unsigned int index = 0; index < gravity_in_channels.size(); ++index)
            FcstUtilities::log << gravity_in_channels[index] << std::endl;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Inertia in porous media:";
        FcstUtilities::log << std::endl;
        for(unsigned int index = 0; index < inertia_in_porous_media.size(); ++index)
            FcstUtilities::log << inertia_in_porous_media[index] << std::endl;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Shear stress in porous media:";
        FcstUtilities::log << std::endl;
        for(unsigned int index = 0; index < shear_stress_in_porous_media.size(); ++index)
            FcstUtilities::log << shear_stress_in_porous_media[index] << std::endl;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Gravity in porous media:";
        FcstUtilities::log << std::endl;
        for(unsigned int index = 0; index < gravity_in_porous_media.size(); ++index)
            FcstUtilities::log << gravity_in_porous_media[index] << std::endl;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Form of the drag force in porous media:";
        FcstUtilities::log << std::endl;
        for(unsigned int index = 0; index < drag_in_porous_media.size(); ++index)
            FcstUtilities::log << drag_in_porous_media[index] << std::endl;
        FcstUtilities::log << std::endl;

        FcstUtilities::log << "Normal velocity is suppressed weakly:";
        FcstUtilities::log << std::endl;
        for(unsigned int index = 0; index < normal_velocity_is_suppressed_weakly.size(); ++index)
            FcstUtilities::log << normal_velocity_is_suppressed_weakly[index] << std::endl;
        FcstUtilities::log << std::endl;
    }

    FcstUtilities::log << "Computational constants:";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;

    for(unsigned int index = 0; index < eta.size(); ++index)
        if( normal_velocity_is_suppressed_weakly[index] )
            FcstUtilities::log << "Normal velocity suppression coefficient " << index+1 << " = " << eta[index] << std::endl;
    FcstUtilities::log << std::endl;

    for(unsigned int index = 0; index < inlet_outlet_velocity_max.size(); ++index)
        if( inlet_outlet_velocity_max[index] != 0.0 )
            FcstUtilities::log << "Maximum inlet-outlet velocity [cm/s] " << index+1 << " = " << inlet_outlet_velocity_max[index] << std::endl;
    FcstUtilities::log << std::endl;

    if( inlet_outlet_velocity_mixture_max != 0 )
        FcstUtilities::log << "Maximum inlet-outlet velocity of the mixture [cm/s] = " << inlet_outlet_velocity_mixture_max << std::endl;
    FcstUtilities::log << std::endl;

    FcstUtilities::log << "Initial data:";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;

    for(component_materialID_value_map::const_iterator iter  = this->component_materialID_value.begin();
                                                       iter != this->component_materialID_value.end();
                                                       ++iter)
    {
        std::map<types::material_id, double> tmp = iter->second;

        for(std::map<types::material_id, double>::const_iterator iter2  = tmp.begin();
                                                                 iter2 != tmp.end();
                                                                 ++iter2)
        {
            FcstUtilities::log << "Name of the solution component  = " << iter->first   << std::endl;
            FcstUtilities::log << "Material id                     = " << iter2->first  << std::endl;

            const std::string str  = iter->first;
            const std::string str2 = str.substr(0,3);

            if( str2 == "den" )
                FcstUtilities::log << "Value of the solution component = " << iter2->second << std::endl;
            if( str2 == "vel" )
                FcstUtilities::log << "Value of the solution component = " << iter2->second << std::endl;

            FcstUtilities::log << std::endl;
        }
    }

    FcstUtilities::log << "Boundary conditions:";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;
    
    //Loop Over Species
    for(unsigned int s = 0; s < n_species; ++s)
    {
        FcstUtilities::log << "Species # " << s+1;
        FcstUtilities::log << std::endl;
        FcstUtilities::log << std::endl;
        
        //Boundaries
        const std::vector< BoundaryType > tmp = this->multi_boundary_types[s];

        std::vector< BoundaryType >::const_iterator iter;
        
        //Loop Over Boundaries
        for( iter = tmp.begin(); iter != tmp.end(); ++iter )
        {
            FcstUtilities::log << iter->boundary_name      << std::endl;
            FcstUtilities::log << iter->boundary_id        << std::endl;
            FcstUtilities::log << iter->boundary_condition << std::endl;
            if( iter->boundary_condition == "Navier slip"  )
                FcstUtilities::log << "Navier slip coefficient = " << theta[s] << std::endl;
            FcstUtilities::log << std::endl;
        }
    }

    FcstUtilities::log << "Boundary data:";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;

    for(component_boundaryID_value_map::const_iterator iter  = this->component_boundaryID_value.begin();
                                                       iter != this->component_boundaryID_value.end();
                                                       ++iter)
    {
        std::map<types::boundary_id, double> tmp = iter->second;

        for(std::map<types::boundary_id, double>::const_iterator iter2  = tmp.begin();
                                                                 iter2 != tmp.end();
                                                                 ++iter2)
        {
            FcstUtilities::log << "Name of the solution component  = " << iter->first  << std::endl;
            FcstUtilities::log << "Boundary id                     = " << iter2->first << std::endl;

            const std::string str  = iter->first;
            const std::string str2 = str.substr(0,3);

            if( str2 == "den" )
                FcstUtilities::log << "Value of the solution component = " << iter2->second  << std::endl;
            if( str2 == "vel" )
                FcstUtilities::log << "Value of the solution component = " << iter2->second << std::endl;

            FcstUtilities::log << std::endl;
        }
    }

    FcstUtilities::log << "------------------------------";
    FcstUtilities::log << std::endl;
}

// ---                                              ---
// --- CompressibleMultiComponentKGEquationsCoupled ---
// ---                                              ---
//Explicit Instantiation
template class NAME::CompressibleMultiComponentKGEquationsCoupled<deal_II_dimension>;