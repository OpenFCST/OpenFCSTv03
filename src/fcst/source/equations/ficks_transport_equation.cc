//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: ficks_transport_equation.cc
//    - Description: Implementation of Ficks diffusion model.
//    - Developers: M. Secanell, M. Bhaiya, V. N. Zingan, A. Kosakian, M. Sabharwal, J. Zhou
//
//---------------------------------------------------------------------------

#include "equations/ficks_transport_equation.h"

namespace NAME = FuelCellShop::Equation;

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::FicksTransportEquation<dim>::FicksTransportEquation(FuelCell::SystemManagement& system_management,
                                                                FuelCellShop::Material::PureGas* solute,
                                                                FuelCellShop::Material::PureGas* solvent,
                                                                boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data )
:
EquationBase<dim>(system_management,data),
gas(solute),
solvent(solvent),
concentration(0.0),
electronic_pot(0.0),
protonic_pot(0.0),
temperature(353.)
{
    std::stringstream s;
    s <<"Ficks Transport Equation - "<<gas->name_material();
    this->equation_name = s.str();

    std::stringstream ss;
    ss <<gas->name_material()<<"_molar_fraction";
    this->name_base_variable = ss.str();

    FcstUtilities::log << "->FuelCellShop::Equation::FicksTransportEquation" << std::endl;

    xi.indices_exist = false;
    t_rev.indices_exist = false;
    s_liquid_water.indices_exist = false;
    p_liquid_water.indices_exist = false;

    this->counter.resize(3, true);
}

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::FicksTransportEquation<dim>::FicksTransportEquation(FuelCell::SystemManagement& system_management, std::string& name_section, boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data )
:
EquationBase<dim>(system_management,data),
concentration(0.0),
electronic_pot(0.0),
protonic_pot(0.0),
temperature(353.)
{
    
    this->equation_name = "Ficks Transport Equation - " + name_section;
    
    
    this->name_base_variable = name_section + "_molar_fraction";
    
    FcstUtilities::log << "->FuelCellShop::Equation::FicksTransportEquation" << std::endl;
    
    xi.indices_exist = false;
    t_rev.indices_exist = false;
    s_liquid_water.indices_exist = false;
    p_liquid_water.indices_exist = false;
    
    this->counter.resize(3, true);
}
// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::FicksTransportEquation<dim>::~FicksTransportEquation()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::FicksTransportEquation<dim>::declare_parameters(ParameterHandler& param)
{
    NAME::EquationBase<dim>::declare_parameters(param);
    
    param.enter_subsection("Equations");
    {
        param.enter_subsection(this->equation_name);
        {

            param.enter_subsection("Boundary data");
            {
                   param.declare_entry("species_flux",
                                       """",
                                        Patterns::Map(   Patterns::Integer(0) , Patterns::Double()   ),
                                       "Enter the molar flux in units of moles/(s*cm^2) that you would like to use as the Neumann boundary condition for each"
                                       "boundary_id. The format should be as follows: boundary_id1:value1, boundary_id2:value2. For example, using"
                                       "3:0.1, 41:0.25 will set the molar flux in boundary 3 to 0.1 and in boundary 41 to 0.25");
            }
            param.leave_subsection();

            param.enter_subsection("Boundary conditions");
            {
                param.declare_entry("Dirichlet Boundary Indicators",
                                    "",
                                    Patterns::List( Patterns::Integer(0) ),
                                    "Boundary_id(s) corresponding to Dirichlet boundaries (where molar fraction is to be specified), "
                                    "provided by comma-separated list of unsigned integers.");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    if (this->data->flag_exists("reaction")&& this->data->flag("reaction"))
    {
        FuelCellShop::Material::CatalystBase::declare_Catalyst_parameters(param);
        FuelCellShop::Material::PolymerElectrolyteBase::declare_PolymerElectrolyte_parameters(param);
        FuelCellShop::Kinetics::BaseKinetics::declare_Kinetics_parameters(param);
        param.enter_subsection("Equations");
        {
            param.enter_subsection(this->equation_name);
            {
                param.enter_subsection("Reaction data");
                {
                    param.declare_entry("Electronic potential,[V]",
                                        "0.4",
                                        Patterns::Double(),
                                        "Electronic potential to be used for the reaction boundary");
                    param.declare_entry("Protonic potential,[V]",
                                        "0.0",
                                        Patterns::Double(),
                                        "Protonic potential to be used for the reaction boundary");
                    param.declare_entry("Type of kinetics",
                                        "DoubleTrapKinetics",
                                        Patterns::Selection("DoubleTrapKinetics|ButlerVolmerKinetics|TafelKinetics|DualPathKinetics"),
                                        "Kinetics model to be used for the Reaction Boundary");
                    param.declare_entry("Reaction name",
                                        "ORR",
                                        Patterns::Selection("ORR|HOR"),
                                        "Type of reaction to be simulated at the Reaction boundary");
                    param.declare_entry("Ionomer dissolution constant,[cm/s]",
                                        "1e-1",                                               //based on Phil Wardlaw's thesis work
                                        Patterns::Double(),
                                        "Dissolution rate constant for the oxygen into the ionomer film");
                    param.declare_entry("Ionomer film thickness,[cm]",
                                        "2e-7",
                                        Patterns::Double(),
                                        "Ionomer film thickness around the catalyst");
                    param.declare_entry("A_Pt,s|g",
                                        "1.0",
                                        Patterns::Double(),
                                        "Ratio of the active platinum area to the solid-gas interface in the reconstruction");
                    param.enter_subsection("Materials");
                    {
                        param.declare_entry("Catalyst type",
                                            "Platinum",
                                            Patterns::Selection("Platinum"),
                                            "Catalyst type on the reaction boundary");
                        param.declare_entry("Electrolyte type",
                                            "Nafion",
                                            Patterns::Selection("Nafion"),
                                            "Electrolyte type on the reaction boundary");
                    }
                    param.leave_subsection();
                }
                param.leave_subsection();
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
}

// ---            ---
// --- initialize ---
// ---            ---

template<int dim>
void
NAME::FicksTransportEquation<dim>::initialize(ParameterHandler& param)
{
    NAME::EquationBase<dim>::initialize(param);
    
    param.enter_subsection("Equations");
    {
        param.enter_subsection(this->equation_name);
        {
            param.enter_subsection("Boundary data");
            {
                if( !param.get("species_flux").empty() )
                {
                    species_flux = FcstUtilities::string_to_map<unsigned int, double>( param.get("species_flux") );
                }
            }
            param.leave_subsection();
            
            param.enter_subsection("Boundary conditions");
            {
                dirichlet_bdry_ids = FcstUtilities::string_to_number<unsigned int>( Utilities::split_string_list( param.get("Dirichlet Boundary Indicators") ) );
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    if (this->data->flag_exists("reaction") && this->data->flag("reaction"))
    {
        std::string catalyst_type,electrolyte_type, reaction_kinetics, reaction_name;
        param.enter_subsection("Equations");
        {
            param.enter_subsection(this->equation_name);
            {
                param.enter_subsection("Reaction data");
                {
                    this->electronic_pot = param.get_double("Electronic potential,[V]");
                    this->protonic_pot = param.get_double("Protonic potential,[V]");
                    reaction_kinetics = param.get("Type of kinetics");
                    reaction_name = param.get("Reaction name");
                    this->k_ionomer = param.get_double("Ionomer dissolution constant,[cm/s]");
                    this->delta = param.get_double("Ionomer film thickness,[cm]");
                    this->prefactor = param.get_double("A_Pt,s|g");
                    param.enter_subsection("Materials");
                    {
                        catalyst_type = param.get("Catalyst type");
                        electrolyte_type = param.get("Electrolyte type");
                    }
                    param.leave_subsection();                   
                }
                param.leave_subsection();
            }
            param.leave_subsection();
        }
        param.leave_subsection();
        param.enter_subsection("Materials");
        param.enter_subsection("Platinum");
        param.set("Method for kinetics parameters (ORR)","Double_trap");
        param.leave_subsection();
        param.leave_subsection();
        catalyst = FuelCellShop::Material::CatalystBase::create_Catalyst(param, catalyst_type);
        electrolyte = FuelCellShop::Material::PolymerElectrolyteBase::create_PolymerElectrolyte(param, electrolyte_type);
        kinetics = FuelCellShop::Kinetics::BaseKinetics::create_Kinetics(param, reaction_kinetics);
        
        FcstUtilities::log<<"Catalyst parameter method: "<<catalyst->get_kinetic_parameter_method()<<std::endl;
        
        
        kinetics->set_catalyst(catalyst.get());
        kinetics->set_electrolyte(electrolyte.get());
        ReactionNames rxn_name;
        if (reaction_name.compare("ORR") == 0)
            rxn_name=ORR;
        else if (reaction_name.compare("HOR") == 0)
            rxn_name=HOR;
        else
            rxn_name=noReaction;
        
        kinetics->set_reaction_kinetics(rxn_name);
        catalyst->set_reaction_kinetics(rxn_name);
    }
    
    
    this->make_internal_cell_couplings();
    this->make_boundary_types();
}

// ---                             ---
// --- assemble_cell_linear_matrix ---
// ---                             ---

template<int dim>
void
NAME::FicksTransportEquation<dim>::assemble_cell_linear_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                                               const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                               FuelCellShop::Layer::BaseLayer<dim>* const                               layer
                                                              )
{
    if ( this->counter[0] )
    {
        this->make_assemblers_generic_constant_data();
        this->counter[0] = false;
    }

    if ( this->counter[1] )
    {
        this->make_assemblers_cell_constant_data(cell_info);
        this->counter[1] = false;
    }

    this->make_assemblers_cell_variable_data(cell_info, layer);
    
    Tensor<2,dim,double> conc_Deff;
    
    //-------- Looping over Quadrature points ----------------------------
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
        if (this->data->flag_exists("knudsen") && this->data->flag("knudsen"))
        {
            double radius,conc_DKnud =0;
            std::string text = boost::lexical_cast<std::string>(cell_info.cell->id());
            int cell_id =  FcstUtilities::cellId_to_index(text);
            conc_DKnud = this->concentration*this->compute_Knudsen_diffusivity(cell_id);
            conc_Deff = this->effective_diffusion_coefficient(conc_Deff_cell[q],conc_DKnud);
        }
        else
            conc_Deff=conc_Deff_cell[q];
        
        //---------------LOOP over i -----------------------------------------------------------------
        for (unsigned int i=0; i < (cell_info.fe(xi.fetype_index)).dofs_per_cell; ++i)
        {
            //--------------LOOP(s) over j-------------------------------------------------------------

            //-----------Assembling Matrix for terms corresponding to "xi" BLOCK------------------------
            for (unsigned int j=0; j < (cell_info.fe(xi.fetype_index)).dofs_per_cell; ++j)
            {
                cell_matrices[xi.block_index].matrix(i,j) += grad_phi_xi_cell[q][i] * conc_Deff * grad_phi_xi_cell[q][j] * this->JxW_cell[q];
            }
        }
    }
}

// ---                             ---
// --- compute_Knudsen_diffusivity ---
// ---                             ---

template<int dim>
double
NAME::FicksTransportEquation<dim>::compute_Knudsen_diffusivity(int cell_index)
{
    double D_Knud,radius;
    std::map< std::string, std::map<int,double>>::const_iterator iter  = this->data->field_data.find("KnudsenRadius");
    if (iter !=this->data->field_data.end())
    {
        std::map<int,double>::const_iterator iter2 = iter->second.find(cell_index);
        if (iter2 != iter->second.end())
            radius = iter2->second;
        else
        {
            FcstUtilities::log<<"Index not found!"<<std::endl;
            AssertThrow(false, ExcInternalError());
        }
    }
    else
    {
        FcstUtilities::log<<"Knudsen Radius section with name KnudsenRadius not found in the field_data"<<std::endl;
        AssertThrow(false, ExcInternalError());
    }
    D_Knud =  2. / 3 * radius * sqrt(8*Constants::R()*this->temperature/(Constants::Pi()*this->gas->get_molar_mass()))*Units::convert(1.,Units::C_UNIT, Units::UNIT);
    return D_Knud;       //Units are cm2/s
}

// ---                             ---
// --- effective_diffusion_coefficient ---
// ---                             ---

template<int dim>
Tensor<2,dim,double>
NAME::FicksTransportEquation<dim>::effective_diffusion_coefficient(Tensor<2,dim,double>& D_bulk,
                                                                   double                D_Knud)
{
    Tensor<2,dim,double> D_eff;
    for (unsigned int a=0; a<dim;++a)
        for (unsigned int b=0; b<dim;++b)
            D_eff[a][b] = D_bulk[a][b]*D_Knud/(D_Knud + D_bulk[a][b]);
    return D_eff;
}


// ---                              ---
// --- assemble_cell_Jacobian_matrix ---
// ---                              ---

template<int dim>
void
NAME::FicksTransportEquation<dim>::assemble_cell_Jacobian_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                                                   const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                   FuelCellShop::Layer::BaseLayer<dim>* const                               layer
                                                                  )
{
    this->assemble_cell_linear_matrix(cell_matrices, cell_info, layer);
    
    if (t_rev.indices_exist || s_liquid_water.indices_exist || p_liquid_water.indices_exist )
    {
        
        //-------- Looping over Quadrature points ----------------------------
        for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
        {
            //---------------LOOP over i -----------------------------------------------------------------
            for (unsigned int i=0; i < (cell_info.fe(xi.fetype_index)).dofs_per_cell; ++i)
            {
                //--------------LOOP(s) over j-------------------------------------------------------------
                
                //----------TERM CORRESPONDING TO "T" BLOCK-----------------------------------------
                if (t_rev.indices_exist)
                {
                    //-----Assembling Matrix------------------------------------------------------------
                    for (unsigned int j=0; j < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++j)
                    {
                        cell_matrices[t_rev.block_index].matrix(i,j) += grad_phi_xi_cell[q][i] * dconc_Deff_dT_cell[q] * cell_info.gradients[last_iter_cell][xi.solution_index][q] * phi_T_cell[q][j] *
                        this->JxW_cell[q];
                    }
                }
                
                //----------TERM CORRESPONDING TO "s" BLOCK-----------------------------------------
                if (s_liquid_water.indices_exist)
                {
                    //-----Assembling Matrix------------------------------------------------------------
                    for (unsigned int j=0; j < (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell; ++j)
                    {
                        cell_matrices[s_liquid_water.block_index].matrix(i,j) += grad_phi_xi_cell[q][i] * dconc_Deff_ds_cell[q] * cell_info.gradients[last_iter_cell][xi.solution_index][q] * phi_s_cell[q][j] *
                        this->JxW_cell[q];
                    }
                }
                
                //----------TERM CORRESPONDING TO "p" BLOCK-----------------------------------------
                if (p_liquid_water.indices_exist)
                {
                    //-----Assembling Matrix------------------------------------------------------------
                    for (unsigned int j=0; j < (cell_info.fe(p_liquid_water.fetype_index)).dofs_per_cell; ++j)
                    {
                        cell_matrices[p_liquid_water.block_index].matrix(i,j) += grad_phi_xi_cell[q][i] * dconc_Deff_dp_cell[q] * cell_info.gradients[last_iter_cell][xi.solution_index][q] * phi_p_cell[q][j] *
                        this->JxW_cell[q];
                    }
                }
            }
        }
    }
}


// ---                        ---
// --- assemble_cell_residual ---
// ---                        ---

template<int dim>
void
NAME::FicksTransportEquation<dim>::assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_rhs,
                                                          const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                          FuelCellShop::Layer::BaseLayer<dim>* const                               layer   )
{
    if ( this->counter[0] )
    {
        this->make_assemblers_generic_constant_data();
        this->counter[0] = false;
    }
    
    if ( this->counter[1] )
    {
        this->make_assemblers_cell_constant_data(cell_info);
        this->counter[1] = false;
    }
    
    this->make_assemblers_cell_variable_data(cell_info, layer);
    
    for (unsigned int q=0; q < this->n_q_points_cell; ++q)
    {
        for (unsigned int i=0; i < (cell_info.fe(xi.fetype_index)).dofs_per_cell; ++i)
        {
            cell_rhs.block(xi.solution_index)(i) += grad_phi_xi_cell[q][i] * conc_Deff_cell[q] * cell_info.gradients[last_iter_cell][xi.solution_index][q] * this->JxW_cell[q];
        }
    }
}

// ---                      ---
// --- assemble_bdry_matrix ---
// ---                      ---

template<int dim>
void
NAME::FicksTransportEquation<dim>::assemble_bdry_matrix(FuelCell::ApplicationCore::MatrixVector&                                 bdry_matrices,
                                                           const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                           FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
       // The boundary integral " q * [ ... ] " is always zero because of the boundary conditions we use:
       // - if xi = known then q = 0
       // - if [-C_tot * Di * grad_xi * n] = known then [ ... ] = 0
}

// ---                        ---
// --- assemble_bdry_residual ---
// ---                        ---

template<int dim>
void
NAME::FicksTransportEquation<dim>::assemble_bdry_residual(FuelCell::ApplicationCore::FEVector&                                     bdry_residual,
                                                             const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                             FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
    if ( this->counter[0] )
    {
        this->make_assemblers_generic_constant_data();
        this->counter[0] = false;
    }

    if ( this->counter[2] )
    {
        this->make_assemblers_bdry_constant_data(bdry_info);
        this->counter[2] = false;
    }
    this->make_assemblers_bdry_variable_data(bdry_info, layer);
    
    if (this->data->flag_exists("reaction") && this->data->flag("reaction"))
    {
        int coefficient;
        if (kinetics->get_reaction_name() == ORR)
            coefficient = 4;
        else if (kinetics->get_reaction_name() == HOR)
            coefficient = 2;

        std::map<unsigned int, double>::const_iterator iter = species_flux.find( bdry_info.dof_face->boundary_indicator() );
        if (iter != species_flux.end() )
        {
            
            //-------- Looping over Quadrature points ----------------------------
            for (unsigned int q = 0; q < this->n_q_points_bdry; ++q)
                for (unsigned int i = 0; i < (bdry_info.fe(xi.fetype_index)).dofs_per_cell; ++i)
                    bdry_residual.block(xi.solution_index)(i) += -phi_xi_bdry[q][i] * this->compute_current(bdry_info.values[last_iter_cell][xi.solution_index][q])/(coefficient*Constants::F()) * this->prefactor  * this->JxW_bdry[q];
        }
    }
    else 
    {
        std::map<unsigned int, double>::const_iterator iter = species_flux.find( bdry_info.dof_face->boundary_indicator() );
        if (iter != species_flux.end() )
        {
            //-------- Looping over Quadrature points ----------------------------
            for (unsigned int q = 0; q < this->n_q_points_bdry; ++q)
                for (unsigned int i = 0; i < (bdry_info.fe(xi.fetype_index)).dofs_per_cell; ++i)
                    bdry_residual.block(xi.solution_index)(i) += phi_xi_bdry[q][i] * iter->second * this->JxW_bdry[q];
            
        }
    }
}

// ---                     ---
// --- print_equation_info ---
// ---                     ---

template<int dim>
void
NAME::FicksTransportEquation<dim>::print_equation_info() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "-------------------------------------------------------------------------------" << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "INTERNAL CELL COUPLINGS FOR \"" << this->equation_name << "\":" << std::endl;
    FcstUtilities::log << std::endl;

    couplings_map::const_iterator iter;

    for( iter = this->internal_cell_couplings.begin(); iter != this->internal_cell_couplings.end(); ++iter )
    {
        FcstUtilities::log << "\"" << iter->first << "\"" << ":";

        std::map<std::string, DoFTools::Coupling> int_map = iter->second;
        std::map<std::string, DoFTools::Coupling>::const_iterator int_iter;

        for( int_iter = int_map.begin(); int_iter != int_map.end(); ++int_iter )
            FcstUtilities::log << "\"" << int_iter->first << "\"" << " is coupled as " << int_iter->second << std::endl;

        FcstUtilities::log << std::endl;
    }

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
                FcstUtilities::log << "Value of the solution component = " << iter2->second << std::endl;
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
                FcstUtilities::log << "Name of the solution component  = " << iter->first   << std::endl;
                FcstUtilities::log << "Boundary id                     = " << iter2->first  << std::endl;
                FcstUtilities::log << "Value of the solution component = " << iter2->second << std::endl;
                FcstUtilities::log << std::endl;
         }
  }

    FcstUtilities::log << std::endl;
    FcstUtilities::log << "-------------------------------------------------------------------------------" << std::endl;
    FcstUtilities::log << std::endl;
}

// ---                              ---
// --- make_internal_cell_couplings ---
// ---                              ---

template<int dim>
void
NAME::FicksTransportEquation<dim>::make_internal_cell_couplings()
{
    AssertThrow(this->system_management->solution_in_userlist(this->name_base_variable), VariableShouldExistForEquation(this->name_base_variable, this->equation_name) );
    AssertThrow(this->system_management->equation_name_to_index(this->equation_name) == this->system_management->solution_name_to_index(this->name_base_variable),
                IndexDoNotMatch(this->name_base_variable, this->equation_name) );

    std::map< std::string, DoFTools::Coupling > tmp;

    std::vector< std::string> sol_names = this->system_management->get_solution_names();

    for (unsigned int i = 0; i < sol_names.size(); ++i)
    {
        if (sol_names[i] == this->name_base_variable)
            tmp[this->name_base_variable] = DoFTools::always;

        else if (sol_names[i] == "temperature_of_REV")
            tmp["temperature_of_REV"] = DoFTools::always;
        
        else if (sol_names[i] == "liquid_water_saturation")
            tmp["liquid_water_saturation"] = DoFTools::always;
        
        else if (sol_names[i] == "capillary_pressure")
            tmp["capillary_pressure"] = DoFTools::always;
        
        else
            tmp[sol_names[i]] = DoFTools::none;
    }

    this->internal_cell_couplings[this->equation_name] = tmp;
}

// ---                     ---
// --- make_boundary_types ---
// ---                     ---

template<int dim>
void
NAME::FicksTransportEquation<dim>::make_boundary_types()
{
    for (unsigned int index = 1; index <= dirichlet_bdry_ids.size(); ++index)
    {
        BoundaryType temp_dirich;

        std::ostringstream streamOut;
        streamOut << index;
        std::string temp_name = "Dirichlet_" + streamOut.str();

        temp_dirich.boundary_name = temp_name;
        temp_dirich.boundary_id = dirichlet_bdry_ids[index-1];
        temp_dirich.boundary_condition = "Dirichlet";

        this->boundary_types.push_back(temp_dirich);
    }
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
void
NAME::FicksTransportEquation<dim>::make_assemblers_generic_constant_data()
{
    //-----------Filling VariableInfo structures------------------------------------------
    xi.solution_index = this->system_management->solution_name_to_index(this->name_base_variable);
    xi.block_index = this->system_management->matrix_block_index(this->equation_name, this->name_base_variable);
    xi.fetype_index = this->system_management->block_info->base_element[xi.solution_index];
    xi.indices_exist = true;

    //----------temperature_of_solid_phase----------------------------------------------------
    if ( this->system_management->solution_in_userlist("temperature_of_REV") )
    {
        t_rev.solution_index = this->system_management->solution_name_to_index("temperature_of_REV");
        t_rev.block_index = this->system_management->matrix_block_index(this->equation_name, "temperature_of_REV");
        t_rev.fetype_index = this->system_management->block_info->base_element[t_rev.solution_index];
        t_rev.indices_exist = true;
    }

    //----------liquid_water_saturation----------------------------------------------------
    if ( this->system_management->solution_in_userlist("liquid_water_saturation") )
    {
        s_liquid_water.solution_index = this->system_management->solution_name_to_index("liquid_water_saturation");
        s_liquid_water.block_index = this->system_management->matrix_block_index(this->equation_name, "liquid_water_saturation");
        s_liquid_water.fetype_index = this->system_management->block_info->base_element[s_liquid_water.solution_index];
        s_liquid_water.indices_exist = true;
    }
    
    //----------capillary_pressure----------------------------------------------------
    if ( this->system_management->solution_in_userlist("capillary_pressure") )
    {
        p_liquid_water.solution_index = this->system_management->solution_name_to_index("capillary_pressure");
        p_liquid_water.block_index = this->system_management->matrix_block_index(this->equation_name, "capillary_pressure");
        p_liquid_water.fetype_index = this->system_management->block_info->base_element[p_liquid_water.solution_index];
        p_liquid_water.indices_exist = true;
    }
}

// ---                                    ---
// --- make_assemblers_cell_constant_data ---
// ---                                    ---

template<int dim>
void
NAME::FicksTransportEquation<dim>::make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
    Assert( xi.indices_exist, ExcMessage("make_assemblers_generic_constant_data function not called before.") );

    this->n_q_points_cell = (cell_info.fe(xi.fetype_index)).n_quadrature_points;
    this->last_iter_cell = cell_info.global_data->find_vector(this->solution_vector_name);

    conc_Deff_cell.resize(this->n_q_points_cell);
    dconc_Deff_dT_cell.resize(this->n_q_points_cell);
    dconc_Deff_ds_cell.resize(this->n_q_points_cell);
    dconc_Deff_dp_cell.resize(this->n_q_points_cell);

    //-------------Allocation------------------------------------------
    phi_xi_cell.resize( this->n_q_points_cell, std::vector< double >( (cell_info.fe(xi.fetype_index)).dofs_per_cell ) );
    grad_phi_xi_cell.resize( this->n_q_points_cell, std::vector< Tensor<1,dim> >( (cell_info.fe(xi.fetype_index)).dofs_per_cell ) );

    if (t_rev.indices_exist)
        phi_T_cell.resize( this->n_q_points_cell, std::vector< double >( (cell_info.fe(t_rev.fetype_index)).dofs_per_cell ) );

    if (s_liquid_water.indices_exist)
        phi_s_cell.resize( this->n_q_points_cell, std::vector< double >( (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell ) );
    
    if (p_liquid_water.indices_exist)
        phi_p_cell.resize( this->n_q_points_cell, std::vector< double >( (cell_info.fe(p_liquid_water.fetype_index)).dofs_per_cell ) );

    this->JxW_cell.resize(this->n_q_points_cell);

    //-----------------------------------------------------------------
}

// ---                                    ---
// --- make_assemblers_bdry_constant_data ---
// ---                                    ---

template<int dim>
void
NAME::FicksTransportEquation<dim>::make_assemblers_bdry_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info)
{
    Assert( xi.indices_exist, ExcMessage("make_assemblers_generic_constant_data function not called before.") );

    this->n_q_points_bdry = (bdry_info.fe(xi.fetype_index)).n_quadrature_points;
    last_iter_bdry = bdry_info.global_data->find_vector(this->solution_vector_name);

    //-------------Allocation------------------------------------------
    phi_xi_bdry.resize( this->n_q_points_bdry, std::vector<double>( (bdry_info.fe(xi.fetype_index)).dofs_per_cell ) );

    //-----------------------------------------------------------------
    this->JxW_bdry.resize(this->n_q_points_bdry);
}

// ---                                    ---
// --- make_assemblers_cell_variable_data ---
// ---                                    ---

template<int dim>
void
NAME::FicksTransportEquation<dim>::make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                         FuelCellShop::Layer::BaseLayer<dim>* const layer)
{
    Assert( this->n_q_points_cell != 0, ExcMessage("make_assemblers_cell_constant_data function not called before.") );

    //---------------Effective Transport Properties---------------------------------------------------------------
    // ----- type infos -------------
    const std::type_info& GasDiffusionLayer = typeid(FuelCellShop::Layer::GasDiffusionLayer<dim>);
    const std::type_info& MicroPorousLayer  = typeid(FuelCellShop::Layer::MicroPorousLayer<dim>);
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);

    const std::type_info& base_layer = layer->get_base_type();

    // Creating some internal variables -- to be used commonly by all the dynamically casted layers
    Table< 2, Tensor<2,dim> > Deff_iso;
    double p, T;
    int index_gas, index_solvent;
    //------
    std::vector< Tensor<2,dim> > Deff_noniso;
    std::map< VariableNames, std::vector< Tensor<2,dim> > > dDeff_du;
    std::vector<VariableNames> deriv_flags;

    // ----- dynamic cast and filling the containers -----------------
    try
    {
        if (base_layer == GasDiffusionLayer)
        {
            FuelCellShop::Layer::GasDiffusionLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::GasDiffusionLayer<dim>* >(layer);

            if (!t_rev.indices_exist)         // t_rev.indices_exist = false
            {
                if (p_liquid_water.indices_exist)
                {
                    ptr->get_T_and_p(T,p);
                    SolutionVariable temperature = SolutionVariable(T, this->n_q_points_cell , temperature_of_REV);
                    ptr->set_temperature( temperature ); 
                    ptr->set_capillary_pressure( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][p_liquid_water.solution_index], capillary_pressure) );
                    deriv_flags.push_back(capillary_pressure);
                    
                    this->concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T))*Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3);   
                    std::vector< Tensor<2,dim> > Deff_capillary;
                    ptr->compute_gas_diffusion(gas, solvent);
                    ptr->effective_gas_diffusivity(Deff_capillary);
                    ptr->get_gas_index(this->gas, index_gas);
                    ptr->get_gas_index(this->solvent, index_solvent);
                    
                    ptr->set_derivative_flags(deriv_flags);
                    ptr->derivative_effective_gas_diffusivity(dDeff_du);
                    
                    for (unsigned int q=0; q<this->n_q_points_cell; ++q)
                    {
                        conc_Deff_cell[q] = this->concentration * Deff_capillary[q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);
                        dconc_Deff_dp_cell[q] = this->concentration * dDeff_du[capillary_pressure][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);
                    }   
                    
                    
                }
                else
                {
                    ptr->get_T_and_p(T,p);
                    this->concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T))*Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3);       // mol/cm^3
                    ptr->effective_gas_diffusivity(Deff_iso);               // m^2/s
                    
                    ptr->get_gas_index(this->gas, index_gas);
                    ptr->get_gas_index(this->solvent, index_solvent);
                    for (unsigned int q=0; q<this->n_q_points_cell; ++q)
                    conc_Deff_cell[q] = this->concentration * Deff_iso(index_gas,index_solvent)*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                         // mol/(cm-s)
                }
            }
            else                                // t_rev.indices_exist = true
            {
                ptr->get_p(p);
                ptr->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index], temperature_of_REV) );
                deriv_flags.push_back(temperature_of_REV);
                if (s_liquid_water.indices_exist)
                {
                    ptr->set_saturation( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][s_liquid_water.solution_index], liquid_water_saturation) );
                    deriv_flags.push_back(liquid_water_saturation);
                }
                
                else if (p_liquid_water.indices_exist)
                {
                    ptr->set_capillary_pressure( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][p_liquid_water.solution_index], capillary_pressure) );
                    deriv_flags.push_back(capillary_pressure);
                }
                
                ptr->compute_gas_diffusion(gas, solvent);
                ptr->effective_gas_diffusivity(Deff_noniso);            // m^2/s

                ptr->set_derivative_flags(deriv_flags);
                ptr->derivative_effective_gas_diffusivity(dDeff_du);
                for (unsigned int q=0; q<this->n_q_points_cell; ++q)
                {
                    T = cell_info.values[last_iter_cell][t_rev.solution_index][q];
                    this->concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T)) * Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3); // mol/cm^3
                    double dC_dT = (-1.)*((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T*T)) * Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3); // mol/(cm^3-K)
                    
                    conc_Deff_cell[q] = this->concentration * Deff_noniso[q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                                            // mol/(cm-s)
                    //std::cout<<"Effective diffusivity in the GDL for "<<this->gas->name_material()<<" in "<<this->solvent->name_material()<<" : "<<Deff_noniso[0]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2)<<std::endl;
                    dconc_Deff_dT_cell[q] = dC_dT * Deff_noniso[q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2) +
                                            this->concentration * dDeff_du[temperature_of_REV][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                       // mol/(cm-s-K)
                    
                    if (s_liquid_water.indices_exist)
                        dconc_Deff_ds_cell[q] = this->concentration * dDeff_du[liquid_water_saturation][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);              // mol/(cm-s)
                        
                        else if (p_liquid_water.indices_exist)
                            dconc_Deff_dp_cell[q] = this->concentration * dDeff_du[capillary_pressure][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);              // mol/(cm-s)
                            
                }
            }
        }

        else if (base_layer == MicroPorousLayer)
        {
            FuelCellShop::Layer::MicroPorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::MicroPorousLayer<dim>* >(layer);

            if (!t_rev.indices_exist)         // t_rev.indices_exist = false
            {
                if (p_liquid_water.indices_exist)
                {
                    ptr->get_T_and_p(T,p);
                    SolutionVariable temperature = SolutionVariable(T, this->n_q_points_cell , temperature_of_REV);
                    ptr->set_temperature( temperature );                    
                    ptr->set_capillary_pressure( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][p_liquid_water.solution_index], capillary_pressure) );                   
                    deriv_flags.push_back(capillary_pressure);

                    this->concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T))*Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3);   
                    std::vector< Tensor<2,dim> > Deff_capillary;
                    ptr->compute_gas_diffusion(gas, solvent);
                    ptr->effective_gas_diffusivity(Deff_capillary);
                    ptr->get_gas_index(this->gas, index_gas);
                    ptr->get_gas_index(this->solvent, index_solvent);
                    
                    ptr->set_derivative_flags(deriv_flags);
                    ptr->derivative_effective_gas_diffusivity(dDeff_du);
                    
                    for (unsigned int q=0; q<this->n_q_points_cell; ++q)
                    {
                        conc_Deff_cell[q] = this->concentration * Deff_capillary[q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);
                        dconc_Deff_dp_cell[q] = this->concentration * dDeff_du[capillary_pressure][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);
                    }   
                }
                else
                {
                    ptr->get_T_and_p(T,p);
                    this->concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T))*Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3);       // mol/cm^3
                    ptr->effective_gas_diffusivity(Deff_iso);               // m^2/s
                    
                    ptr->get_gas_index(this->gas, index_gas);
                    ptr->get_gas_index(this->solvent, index_solvent);
                    for (unsigned int q=0; q<this->n_q_points_cell; ++q)
                    conc_Deff_cell[q] = this->concentration * Deff_iso(index_gas,index_solvent)*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                         // mol/(cm-s)
                }
            }
            else                                // t_rev.indices_exist = true
            {
                ptr->get_p(p);
                ptr->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index], temperature_of_REV) );
                deriv_flags.push_back(temperature_of_REV);
                if (s_liquid_water.indices_exist)
                {
                    ptr->set_saturation( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][s_liquid_water.solution_index], liquid_water_saturation) );
                    deriv_flags.push_back(liquid_water_saturation);
                }
                
                else if (p_liquid_water.indices_exist)
                {
                    ptr->set_capillary_pressure( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][p_liquid_water.solution_index], capillary_pressure) );
                    deriv_flags.push_back(capillary_pressure);
                }
                
                ptr->compute_gas_diffusion(gas, solvent);
                ptr->effective_gas_diffusivity(Deff_noniso);            // m^2/s
                
                ptr->set_derivative_flags(deriv_flags);
                ptr->derivative_effective_gas_diffusivity(dDeff_du);
                for (unsigned int q=0; q<this->n_q_points_cell; ++q)
                {
                    T = cell_info.values[last_iter_cell][t_rev.solution_index][q];
                    this->concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T)) * Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3); // mol/cm^3
                    double dC_dT = (-1.)*((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T*T)) * Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3); // mol/(cm^3-K)
                    
                    conc_Deff_cell[q] = this->concentration * Deff_noniso[q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                                            // mol/(cm-s)
                    //std::cout<<"Effective diffusivity in the MPL for "<<this->gas->name_material()<<" in "<<this->solvent->name_material()<<" : "<<Deff_noniso[0]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2)<<std::endl;
                    dconc_Deff_dT_cell[q] = dC_dT * Deff_noniso[q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2) +
                                            this->concentration * dDeff_du[temperature_of_REV][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                       // mol/(cm-s-K)
                    
                    if (s_liquid_water.indices_exist)
                        dconc_Deff_ds_cell[q] = this->concentration * dDeff_du[liquid_water_saturation][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);              // mol/(cm-s)
                        
                        else if (p_liquid_water.indices_exist)
                            dconc_Deff_dp_cell[q] = this->concentration * dDeff_du[capillary_pressure][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);              // mol/(cm-s)                        
                            
                }
            }
        }
        
        else if (base_layer == CatalystLayer)
        {
            FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);
            
            if (!t_rev.indices_exist)         // t_rev.indices_exist = false
            {
                
                if (p_liquid_water.indices_exist)
                {
                    ptr->get_T_and_p(T,p);
                    SolutionVariable temperature = SolutionVariable(T, this->n_q_points_cell , temperature_of_REV);
                    ptr->set_temperature( temperature );
                    ptr->set_capillary_pressure( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][p_liquid_water.solution_index], capillary_pressure) );
                    deriv_flags.push_back(capillary_pressure);
                    
                    this->concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T))*Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3);   
                    std::vector< Tensor<2,dim> > Deff_capillary;
                    ptr->compute_gas_diffusion(gas, solvent);
                    ptr->effective_gas_diffusivity(Deff_capillary);
                    ptr->get_gas_index(this->gas, index_gas);
                    ptr->get_gas_index(this->solvent, index_solvent);
                    
                    ptr->set_derivative_flags(deriv_flags);
                    ptr->derivative_effective_gas_diffusivity(dDeff_du);
                    
                    for (unsigned int q=0; q<this->n_q_points_cell; ++q)
                    {
                        conc_Deff_cell[q] = this->concentration * Deff_capillary[q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);
                        dconc_Deff_dp_cell[q] = this->concentration * dDeff_du[capillary_pressure][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);
                    }   
                }
                else
                {
                    ptr->get_T_and_p(T,p);
                    this->concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T))*Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3);       // mol/cm^3
                    ptr->effective_gas_diffusivity(Deff_iso);               // m^2/s
                    
                    ptr->get_gas_index(this->gas, index_gas);
                    ptr->get_gas_index(this->solvent, index_solvent);
                    for (unsigned int q=0; q<this->n_q_points_cell; ++q)
                    conc_Deff_cell[q] = this->concentration * Deff_iso(index_gas,index_solvent)*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                         // mol/(cm-s)
                }
            }
            else                                // t_rev.indices_exist = true
            {
                ptr->get_p(p);
                ptr->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index], temperature_of_REV) );
                deriv_flags.push_back(temperature_of_REV);
                if (s_liquid_water.indices_exist)
                {
                    ptr->set_saturation( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][s_liquid_water.solution_index], liquid_water_saturation) );
                    deriv_flags.push_back(liquid_water_saturation);
                }
                
                else if (p_liquid_water.indices_exist)
                {
                    ptr->set_capillary_pressure( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][p_liquid_water.solution_index], capillary_pressure) );
                    deriv_flags.push_back(capillary_pressure);
                }
                
                ptr->compute_gas_diffusion(gas, solvent);
                ptr->effective_gas_diffusivity(Deff_noniso);            // m^2/s
                //std::cout<<"Effective diffusivity in the CL for "<<this->gas->name_material()<<" in "<<this->solvent->name_material()<<" : "<<Deff_noniso[0]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2)<<std::endl;
                ptr->set_derivative_flags(deriv_flags);
                ptr->derivative_effective_gas_diffusivity(dDeff_du);
                for (unsigned int q=0; q<this->n_q_points_cell; ++q)
                {
                    T = cell_info.values[last_iter_cell][t_rev.solution_index][q];
                    this->concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T)) * Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3); // mol/cm^3
                    double dC_dT = (-1.)*((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T*T)) * Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3); // mol/(cm^3-K)
                    
                    conc_Deff_cell[q] = this->concentration * Deff_noniso[q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                                            // mol/(cm-s)
                    dconc_Deff_dT_cell[q] = dC_dT * Deff_noniso[q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2) +
                                            this->concentration * dDeff_du[temperature_of_REV][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                       // mol/(cm-s-K)
                    
                    if (s_liquid_water.indices_exist)
                        dconc_Deff_ds_cell[q] = this->concentration * dDeff_du[liquid_water_saturation][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);              // mol/(cm-s)
                        
                        else if (p_liquid_water.indices_exist)
                            dconc_Deff_dp_cell[q] = this->concentration * dDeff_du[capillary_pressure][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);              // mol/(cm-s)                        
                            
                }
            }
        }
        else
            AssertThrow( false, ExcNotImplemented() );
    }
    catch(const std::bad_cast& e)
    {
        const std::type_info& info = typeid(*layer);
        FcstUtilities::log << "Object of type "<<info.name()<<" not implemented"<< std::endl;
        FcstUtilities::log << e.what() << std::endl;
    }
    this->temperature = T;

    //---------------------------------------------------------------------------------------------------------------
    //------------Looping over quadrature points in the cell --------------------------------------------------------
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
        //-------JxW----------
        this->JxW_cell[q] = (cell_info.fe(xi.fetype_index)).JxW(q);

        //------ Filling shape functions etc ----------------------------------------------------------------------
        //------ This avoids recalculating shape functions etc for efficiency -------------------------------------
        for (unsigned int k=0; k < (cell_info.fe(xi.fetype_index)).dofs_per_cell; ++k)
        {
            phi_xi_cell[q][k] = (cell_info.fe(xi.fetype_index)).shape_value(k,q);
            grad_phi_xi_cell[q][k] = (cell_info.fe(xi.fetype_index)).shape_grad(k,q);
        }

        if (t_rev.indices_exist)
            for (unsigned int k=0; k < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++k)
                phi_T_cell[q][k] = (cell_info.fe(t_rev.fetype_index)).shape_value(k,q);

        if (s_liquid_water.indices_exist)
            for (unsigned int k=0; k < (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell; ++k)
                phi_s_cell[q][k] = (cell_info.fe(s_liquid_water.fetype_index)).shape_value(k,q);
            
        if (p_liquid_water.indices_exist)
            for (unsigned int k=0; k < (cell_info.fe(p_liquid_water.fetype_index)).dofs_per_cell; ++k)
                phi_p_cell[q][k] = (cell_info.fe(p_liquid_water.fetype_index)).shape_value(k,q);
    }
}

// ---                                    ---
// --- make_assemblers_bdry_variable_data ---
// ---                                    ---

template<int dim>
void
NAME::FicksTransportEquation<dim>::make_assemblers_bdry_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                                         FuelCellShop::Layer::BaseLayer<dim>* const layer)
{
    Assert( this->n_q_points_bdry != 0, ExcMessage("make_assemblers_bdry_constant_data function not called before.") );

    //---------------------------------------------------------------------------------------------------------------
    //------------Looping over quadrature points in the cell --------------------------------------------------------
    for (unsigned int q=0; q < this->n_q_points_bdry; ++q)
    {
        this->JxW_bdry[q] = (bdry_info.fe(xi.fetype_index)).JxW(q);

        for (unsigned int k=0; k < (bdry_info.fe(xi.fetype_index)).dofs_per_cell; ++k)
            phi_xi_bdry[q][k] = (bdry_info.fe(xi.fetype_index)).shape_value(k,q);
    }
    
    /*
     * If boundary assemblers are called before cell assemblers, it is possible that this->concentration will be zero.
     * In that case, we compute the concentration here.
     */
    if (this->concentration<1e-12)
        {
        const std::type_info& GasDiffusionLayer = typeid(FuelCellShop::Layer::GasDiffusionLayer<dim>);
        const std::type_info& MicroPorousLayer  = typeid(FuelCellShop::Layer::MicroPorousLayer<dim>);
        const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);

        const std::type_info& base_layer = layer->get_base_type(); 
    
        double p, T;
    
        try
        {
            if (base_layer == GasDiffusionLayer)
            {
                FuelCellShop::Layer::GasDiffusionLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::GasDiffusionLayer<dim>* >(layer); 
                ptr->get_T_and_p(T,p);
            }

            else if (base_layer == MicroPorousLayer)
            {
                FuelCellShop::Layer::MicroPorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::MicroPorousLayer<dim>* >(layer);
                ptr->get_T_and_p(T,p);
            }
            
            else if (base_layer == CatalystLayer)
            {
                FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);
                ptr->get_T_and_p(T,p);
            }
            else
                AssertThrow( false, ExcNotImplemented() );
            
            this->concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T))*Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3);             
        }
        catch(const std::bad_cast& e)
        {
            const std::type_info& info = typeid(*layer);
            FcstUtilities::log << "Object of type "<<info.name()<<" not implemented"<< std::endl;
            FcstUtilities::log << e.what() << std::endl;
        }        
    }
}

template<int dim>
double
NAME::FicksTransportEquation<dim>::compute_current(double x_gas)
{
    double c_eq = this->concentration* x_gas * Constants::R() * this->temperature/(this->electrolyte->get_H_O2()*1.0e-6);
    double c_react = 0.5 * c_eq;
    double f_x = get_ICCP_residual(c_react, c_eq);
    double f_x_; //First derivative
    unsigned int loops = 0;
    
    while(std::abs(f_x) > 1e-12)
    {
        double perc_change = 1.0e-5;
        f_x_ = (get_ICCP_residual(c_react*(1.0+perc_change), c_eq) - get_ICCP_residual(c_react*(1-perc_change), c_eq))/(2.0*perc_change*c_react);
        
        double step_factor = 1.0;
        
        // Apply over-relaxation for first 20 steps
        if (loops < 20)
            step_factor *= 0.05*loops;
        
        while(c_react - step_factor*(f_x/f_x_) < 0.0)
        {
            //if the step will take c_inner below 0 then reduce it
            step_factor /=2.0;
            
            if(step_factor < 1e-6) //Something has gone really wrong... (used to be 1e-5)
            {
                f_x = get_ICCP_residual(c_react, c_eq);
                FcstUtilities::log<<"Inner C_O2_REACT value failed to converge (Step reduction loop). Residual is "<<std::abs(f_x)<<" C_react is "<<c_react<<std::endl;
                AssertThrow (false, ExcMessage ("BC inner C_O2_REACT value failed to converge (Step reduction loop)."));
            }
        }
        
        c_react -= step_factor*(f_x/f_x_);
        f_x = get_ICCP_residual(c_react, c_eq);
        
        if(loops++ > 2000) {
            FcstUtilities::log<<"Inner C_O2_REACT value failed to converge (Step reduction loop). Residual is "<<std::abs(f_x)<<" C_react is "<<c_react<<std::endl;
            AssertThrow (false, ExcMessage ("BC inner C_O2_REACT value failed to converge (Step reduction loop)."));
        }
    }
    std::vector<double> J(1,0.0);
    kinetics->current_density(J);
    return J[0];
}

template<int dim>
double
NAME::FicksTransportEquation<dim>::get_ICCP_residual(double c_react, double c_eq)
{
    std::vector<double> J(1,0.0);
    
    std::vector<SolutionVariable> temp_react;
    
    SolutionVariable proton_pot = SolutionVariable(this->protonic_pot, 1, protonic_electrical_potential);
    SolutionVariable electron_pot = SolutionVariable(this->electronic_pot, 1, electronic_electrical_potential);
    SolutionVariable temp = SolutionVariable(this->temperature, 1, temperature_of_REV);
    
    
    kinetics->set_electrolyte_potential(proton_pot);
    kinetics->set_solid_potential(electron_pot);
    kinetics->set_temperature (temp);
    int coefficient;
    if (kinetics->get_reaction_name() == ORR)
    {
        coefficient = 4;
        temp_react.push_back(SolutionVariable(c_react, 1, oxygen_concentration));
    }
    
    else if (kinetics->get_reaction_name() == HOR)
    {
        coefficient = 2;
        temp_react.push_back(SolutionVariable(c_react, 1, hydrogen_concentration));
    }
    kinetics->set_reactant_concentrations(temp_react);
    kinetics->current_density(J);
    double D_O2;
    this->electrolyte->set_T(this->temperature);
    this->electrolyte->oxygen_diffusivity(D_O2);
    double D_O2_delta= D_O2/this->delta;
    float term2 = this->k_ionomer * (D_O2_delta * (c_eq - c_react))/( this->k_ionomer + D_O2_delta);
    float term1 = J[0]/(coefficient * Constants::F());
    double residual = term1 - term2;
    return residual;
}
// ---            ---
// --- class_test ---
// ---            ---

template<int dim>
void
NAME::FicksTransportEquation<dim>::class_test()
{
  FuelCellShop::Material::Hydrogen H2;
  FuelCellShop::Material::Oxygen O2;
  FuelCell::ApplicationCore::BlockInfo block_info;
  Table< 2, DoFTools::Coupling > cell_couplings;
  Table< 2, DoFTools::Coupling > flux_couplings;
  FuelCell::SystemManagement sys(block_info, cell_couplings, flux_couplings);
  FcstUtilities::log<<"Create object FicksTransportEquation:"<<std::endl;
  FuelCellShop::Equation::FicksTransportEquation<dim> test(sys,&O2,&H2);
  test.print_equation_info();
}

// ---                         ---
// --- explicit instantiations ---
// ---                         ---

template class NAME::FicksTransportEquation<deal_II_dimension>;
