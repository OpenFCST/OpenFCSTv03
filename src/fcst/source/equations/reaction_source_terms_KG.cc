// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: reaction_source_terms.cc
// - Description: This class is used to assemble both cell matrix and cell residual
//                for reaction source terms in the catalyst layers for various equation classes
// - Developers: Madhur Bhaiya, Marc Secanell, Valentin N. Zingan
//
// ----------------------------------------------------------------------------

#include "equations/reaction_source_terms_KG.h"

namespace NAME = FuelCellShop::Equation;

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::ReactionSourceTermsKG<dim>::ReactionSourceTermsKG(FuelCell::SystemManagement& system_management,
                                                        boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data)
:
NAME::ReactionSourceTermsBase<dim>(system_management,data),
electronic_electrical_potential_extractor(0),
protonic_electrical_potential_extractor(0)
{
    FcstUtilities::log << "->FuelCellShop::Equation::ReactionSourceTermsKG" << std::endl;
    
    eq_generic_prefix = "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - ";    

    this->counter.resize(2, true);
    
    anode_kinetics = NULL;
    cathode_kinetics = NULL;

    indexO2  = -1; // does not exist in the problem
    indexH2  = -1; // does not exist in the problem
    indexH2O = -1; // does not exist in the problem
}

// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::ReactionSourceTermsKG<dim>::~ReactionSourceTermsKG()
{}

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::ReactionSourceTermsKG<dim>::declare_parameters(ParameterHandler& param) const
{

}

// ---            ---
// --- initialize ---
// ---            ---

template<int dim>
void
NAME::ReactionSourceTermsKG<dim>::initialize(ParameterHandler& param)
{
    Assert( cathode_kinetics!=NULL || anode_kinetics!=NULL, ExcMessage("At least one of either cathode/anode kinetics should be set before ReactionSourceTermsKG::initialize.") );

    //-- Assertion checks for "Required" solution variables for Cathode/Anode kinetics
    AssertThrow( this->system_management->solution_in_userlist("electronic_electrical_potential"), VariableNotFoundForKinetics("Cathode/Anode", "electronic_electrical_potential") );
    AssertThrow( this->system_management->solution_in_userlist("protonic_electrical_potential"), VariableNotFoundForKinetics("Cathode/Anode", "protonic_electrical_potential") );
}

// ---                                ---
// --- adjust_internal_cell_couplings ---
// ---                                ---
// TODO This routine is not general. It will not work with appPEMFC or with anything else other than cathodeKG... Needs to be redone!!!
template<int dim>
void
NAME::ReactionSourceTermsKG<dim>::adjust_internal_cell_couplings(std::vector< couplings_map >& equation_map,
                                                                 const std::vector< FuelCellShop::Material::PureGas* >& gases)
{
    AssertThrow( equation_map.size() != 0, ExcMessage("Vector size should be greater than zero in ReactionSourceTermsKG::adjust_internal_cell_couplings") );
    AssertThrow( gases.size() <= 5,
                 ExcMessage("The total number of gases used in this problem can not exceed 5") );
    
    // set number of species and gas indices:
    n_species = gases.size();    
    for(unsigned int g = 1; g <= n_species; ++g)
    {
        if( gases[g-1]->name_material() == "oxygen"   )
            indexO2  = g;
        else if( gases[g-1]->name_material() == "hydrogen" )
            indexH2  = g;
        else if( gases[g-1]->name_material() == "water"    )
            indexH2O = g;
    }
    
    //-- Set the multiplier to be used in for the reactions 
    // TODO: (THIS IS NOT THE PLACE FOR IT) !!!!!!!!!!!!!!
    // Move to make_assemblers_generic_constant_data
    molar_mass.resize(n_species);
    for(unsigned int g = 0; g < n_species; ++g) 
        molar_mass[g] = Units::convert(gases[g]->get_molar_mass(), Units::UNIT, Units::K_UNIT); // Note: 1.0e3 is used to go from kg/mol to g/mol
    if( indexO2  != -1 )
        multiplierO2  = - molar_mass[indexO2-1]  / ( 4.0 * Constants::F() );
    if( indexH2  != -1 )
        multiplierH2  = - molar_mass[indexH2-1]  / ( 2.0 * Constants::F() );
    if( indexH2O != -1 )
        multiplierH2O =   molar_mass[indexH2O-1] / ( 2.0 * Constants::F() );
    
    //-- TODO:  why can I not use blocks here?
    std::vector<unsigned int> density_indices(n_species);   
    for(unsigned int g = 0; g < n_species; ++g)
    {
        density_indices[g] = g*(dim+1) + this->system_management->solution_name_to_index("density_species_1");
        density_extractors.push_back( FEValuesExtractors::Scalar(density_indices[g]) );
    }
    electronic_electrical_potential_extractor.component = this->system_management->solution_name_to_index("electronic_electrical_potential");
    protonic_electrical_potential_extractor.component   = this->system_management->solution_name_to_index("protonic_electrical_potential");
    
    //--
    for(unsigned int g = 1; g <= n_species; ++g)
    {
        std::ostringstream streamOut;
        streamOut << g;
        std::string name = "species " + streamOut.str();
        eq_postfixes.push_back( name.c_str() );
        
        name = "species_" + streamOut.str();
        var_postfixes.push_back( name.c_str() );
    }
    
    //-- Find coupling_map that corresponds to the mass conservation equation for the species we would like to react
    // - index stores the vector component with the coupling_map corresponding to that equation.
    // TODO Now it is assumed species 1 is the only reacting species
    unsigned int index;
    bool         flag = false;
    for(unsigned int i = 0; i < equation_map.size(); ++i)
    {
        const couplings_map tmp = equation_map[i];      
        std::string reacting_species_equation = "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - mass conservation - species 1";
        couplings_map::const_iterator iter = tmp.find(reacting_species_equation.c_str());        
        if( iter != tmp.end() )
        {
            index = i;
            flag  = true;
            break;
        }
    }
    AssertThrow( flag, ExcMessage("Kerkhof equations are not found in this problem") );
    
    //-- 
    for(unsigned int g = 1; g <= n_species; ++g)
    {
        if( g == indexO2 )
        {
            if( (anode_kinetics!=NULL   && anode_kinetics->get_reaction_name()   == ORR) ||
                (cathode_kinetics!=NULL && cathode_kinetics->get_reaction_name() == ORR) )
                set_species_couplings(index,g,equation_map);
        }
        
        if( g == indexH2 )
        {
            if( (anode_kinetics!=NULL   && anode_kinetics->get_reaction_name()   == HOR) ||
                (cathode_kinetics!=NULL && cathode_kinetics->get_reaction_name() == HOR) )
                set_species_couplings(index,g,equation_map);
        }
        if( g == indexH2O )
        {
            if( (anode_kinetics!=NULL   && anode_kinetics->get_reaction_name()   == ORR) ||
                (cathode_kinetics!=NULL && cathode_kinetics->get_reaction_name() == ORR) ) 
            {
                set_species_couplings(index,g,equation_map);
                
                // also add coupling wiht oxygen:
                eq_name = "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - mass conservation - " + eq_postfixes[g-1];                          
                couplings_map::iterator iter = equation_map[index].find(eq_name);
                set_density_couplings (indexO2, iter);
            }
        }
    }
    
    //-- OTHER EQUATIONS  
    for(unsigned int i = 0; i < equation_map.size(); ++i)
        if( i != index )
            for( couplings_map::iterator iter  = equation_map[i].begin(); iter != equation_map[i].end(); ++iter )
            {
                if( iter->first == "Electron Transport Equation" )
                {                    
                    iter->second["protonic_electrical_potential"] = DoFTools::always;
                    
                    if( (anode_kinetics!=NULL   && anode_kinetics->get_reaction_name()   == ORR ) ||
                        (cathode_kinetics!=NULL && cathode_kinetics->get_reaction_name() == ORR ) )
                        set_density_couplings (indexO2, iter);
                    
                    if( (anode_kinetics!=NULL   && anode_kinetics->get_reaction_name()   == HOR) ||
                        (cathode_kinetics!=NULL && cathode_kinetics->get_reaction_name() == HOR) )
                        set_density_couplings (indexH2, iter);
                }
                
                if( iter->first == "Proton Transport Equation" )
                {
                    iter->second["electronic_electrical_potential"] = DoFTools::always;
                    
                    if( (anode_kinetics!=NULL   && anode_kinetics->get_reaction_name()   == ORR ) ||
                        (cathode_kinetics!=NULL && cathode_kinetics->get_reaction_name() == ORR ) )
                        set_density_couplings (indexO2, iter);
                    
                    if( (anode_kinetics!=NULL   && anode_kinetics->get_reaction_name()   == HOR ) ||
                        (cathode_kinetics!=NULL && cathode_kinetics->get_reaction_name() == HOR ) )
                        set_density_couplings (indexH2, iter);
                }
            }
}

// ---                      ---
// --- assemble_cell_matrix ---
// ---                      ---

template<int dim>
void
NAME::ReactionSourceTermsKG<dim>::assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                                     const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                     FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
    
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    const std::type_info& info          = layer->get_base_type();
    
    if( info == CatalystLayer )
    {
        FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);
        const FuelCellShop::Material::GasMixture* gas_mixture = ptr->get_gas_mixture();
        
        if( this->counter[0] )
        {
            this->make_assemblers_generic_constant_data();
            this->counter[0] = false;
        }
        if( this->counter[1] )
        {
            this->make_assemblers_cell_constant_data(cell_info);
            this->counter[1] = false;
        }
        
        //-- CELL VARIABLE DATA 
        this->make_assemblers_cell_variable_data(cell_info, layer);
        
        //-- LOCAL CELL MATRIX 
        FullMatrix<double> local_matrix(this->dofs_per_cell, this->dofs_per_cell);
        
        /////////// ///////////// /////////
        // ANODE // // CATHODE // // HOR //
        /////////// ///////////// /////////        
        if( ptr->get_kinetics() == anode_kinetics   && ptr->get_kinetics()->get_reaction_name() == HOR ||
            ptr->get_kinetics() == cathode_kinetics && ptr->get_kinetics()->get_reaction_name() == HOR )
        {            
            for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                    for(unsigned int j = 0; j < this->dofs_per_cell; ++j)
                    {
                        for(unsigned int s = 0; s < n_species; ++s)
                        {
                            if( s+1 == indexH2 )
                            {
                                local_matrix(i,j) += this->JxW_cell[q]*phi_density[s][q][i]*multiplierH2*( DHOR_current_density_Dhydrogen_concentration[q]*(1.0/molar_mass[s])*phi_density[s][q][j]
                                + DHOR_current_density_Delectronic_electrical_potential[q]*phi_electronic_electrical_potential[q][j]
                                + DHOR_current_density_Dprotonic_electrical_potential[q]*phi_protonic_electrical_potential[q][j] );
                            }
                        }
                        
                        //-- Electronic poitential:
                        local_matrix(i,j) += this->JxW_cell[q]*phi_electronic_electrical_potential[q][i]*(1.0)*( DHOR_current_density_Dhydrogen_concentration[q]*(1.0/molar_mass[indexH2-1])*phi_density[indexH2-1][q][j]
                        + DHOR_current_density_Delectronic_electrical_potential[q]*phi_electronic_electrical_potential[q][j]
                        + DHOR_current_density_Dprotonic_electrical_potential[q]*phi_protonic_electrical_potential[q][j] );
                        //-- Protonic potential:
                        local_matrix(i,j) += this->JxW_cell[q]*phi_protonic_electrical_potential[q][i]*(-1.0)*( DHOR_current_density_Dhydrogen_concentration[q]*(1.0/molar_mass[indexH2-1])*phi_density[indexH2-1][q][j]
                        + DHOR_current_density_Delectronic_electrical_potential[q]*phi_electronic_electrical_potential[q][j]
                        + DHOR_current_density_Dprotonic_electrical_potential[q]*phi_protonic_electrical_potential[q][j] );
                    }
        }
        /////////// ///////////// /////////
        // ANODE // // CATHODE // // ORR //
        /////////// ///////////// ///////// 
        if( ptr->get_kinetics() == anode_kinetics   && ptr->get_kinetics()->get_reaction_name() == ORR ||
            ptr->get_kinetics() == cathode_kinetics && ptr->get_kinetics()->get_reaction_name() == ORR )
        {
            for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                    for(unsigned int j = 0; j < this->dofs_per_cell; ++j)
                    {
                        for(unsigned int s = 0; s < n_species; ++s)
                        {
                            if( s+1 == indexO2 )
                            {
                                local_matrix(i,j) += this->JxW_cell[q]*phi_density[s][q][i]*multiplierO2*( DORR_current_density_Doxygen_concentration[q]*(1.0/molar_mass[s])*phi_density[s][q][j]
                                + DORR_current_density_Delectronic_electrical_potential[q]*phi_electronic_electrical_potential[q][j]
                                + DORR_current_density_Dprotonic_electrical_potential[q]*phi_protonic_electrical_potential[q][j] );
                            }
                            else if ( s+1 == indexH2O )
                            {
                                local_matrix(i,j) += this->JxW_cell[q]*phi_density[s][q][i]*multiplierH2O*( DORR_current_density_Doxygen_concentration[q]*(1.0/molar_mass[indexO2-1])*phi_density[indexO2-1][q][j]
                                + DORR_current_density_Delectronic_electrical_potential[q]*phi_electronic_electrical_potential[q][j]
                                + DORR_current_density_Dprotonic_electrical_potential[q]*phi_protonic_electrical_potential[q][j] );
                            }
                        }
                        // -- ELECTRON TRANSPORT EQUATION
                        local_matrix(i,j) += this->JxW_cell[q]*phi_electronic_electrical_potential[q][i]*(-1.0)*( DORR_current_density_Doxygen_concentration[q]*(1.0/molar_mass[indexO2-1])*phi_density[indexO2-1][q][j]
                        + DORR_current_density_Delectronic_electrical_potential[q]*phi_electronic_electrical_potential[q][j]
                        + DORR_current_density_Dprotonic_electrical_potential[q]*phi_protonic_electrical_potential[q][j] );
                        
                        // -- PROTON TRANSPORT EQUATION                        
                        local_matrix(i,j) += this->JxW_cell[q]*phi_protonic_electrical_potential[q][i]*(1.0)*( DORR_current_density_Doxygen_concentration[q]*(1.0/molar_mass[indexO2-1])*phi_density[indexO2-1][q][j]
                        + DORR_current_density_Delectronic_electrical_potential[q]*phi_electronic_electrical_potential[q][j]
                        + DORR_current_density_Dprotonic_electrical_potential[q]*phi_protonic_electrical_potential[q][j] );
                        
                    }
        }
        //-- DEALII -> APPFRAME
         this->dealII_to_appframe(cell_matrices,
                                  local_matrix,
                                  this->matrix_block_indices);
    }
    else
    {
        FcstUtilities::log << "Layer you specified is not FuelCellShop::Layer::CatalystLayer<dim>" << std::endl;
        AssertThrow( false , ExcInternalError() );
    }

}

// ---                        ---
// --- assemble_cell_residual ---
// ---                        ---

template<int dim>
void
NAME::ReactionSourceTermsKG<dim>::assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
                                                       const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                       FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
    
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    const std::type_info& info          = layer->get_base_type();
    
    ////////////////////
    // CATALYST LAYER //
    ////////////////////
    
    if( info == CatalystLayer )
    {
        FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);
        const FuelCellShop::Material::GasMixture* gas_mixture = ptr->get_gas_mixture();
        
        // -- GENERIC CONSTANT DATA 
        if( this->counter[0] )
        {
            this->make_assemblers_generic_constant_data();
            this->counter[0] = false;
        }
        //-- CELL CONSTANT DATA
        if( this->counter[1] )
        {
            this->make_assemblers_cell_constant_data(cell_info);
            this->counter[1] = false;
        }
        
        //-- CELL VARIABLE DATA 
        this->make_assemblers_cell_variable_data(cell_info, layer);
        
        //-- LOCAL CELL RESIDUAL
        Vector<double> local_residual(this->dofs_per_cell);                   
        
        /////////// ///////////// /////////
        // ANODE // // CATHODE // // HOR //
        /////////// ///////////// /////////
        if( ptr->get_kinetics() == anode_kinetics   && ptr->get_kinetics()->get_reaction_name() == HOR
            ||
            ptr->get_kinetics() == cathode_kinetics && ptr->get_kinetics()->get_reaction_name() == HOR )
        {                       
            
            for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                {
                    //-- LOOP OVER SPECIES
                    for(unsigned int s = 0; s < n_species; ++s)
                    {
                        if( s+1 == indexH2 )        {                            
                            local_residual(i) += phi_density[s][q][i]*multiplierH2*HOR_current_density[q]*this->JxW_cell[q];
                        }
                    }
                    //-- Electronic potential
                    local_residual(i) += phi_electronic_electrical_potential[q][i]*(1.0)*HOR_current_density[q]*this->JxW_cell[q];
                    //-- Protonic potential
                    local_residual(i) += phi_protonic_electrical_potential[q][i]*(-1.0)*HOR_current_density[q]*this->JxW_cell[q];
                }
        }
        /////////// ///////////// /////////
        // ANODE // // CATHODE // // ORR //
        /////////// ///////////// /////////                       
        if( ptr->get_kinetics() == anode_kinetics   && ptr->get_kinetics()->get_reaction_name() == ORR
            ||
            ptr->get_kinetics() == cathode_kinetics && ptr->get_kinetics()->get_reaction_name() == ORR )
        {
            
            for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                {
                    //-- Mass transport eq.                          
                    for(unsigned int s = 0; s < n_species; ++s)
                    {
                        if( s+1 == indexO2 )
                            local_residual(i) += phi_density[s][q][i]*multiplierO2*ORR_current_density[q]*this->JxW_cell[q];                              
                        if( s+1 == indexH2O )
                            local_residual(i) += phi_density[s][q][i]*multiplierH2O*ORR_current_density[q]*this->JxW_cell[q];
                    }
                    //-- 
                    local_residual(i) += phi_electronic_electrical_potential[q][i]*(-1.0)*ORR_current_density[q]*this->JxW_cell[q];
                    //--
                    local_residual(i) += phi_protonic_electrical_potential[q][i]*(1.0)*ORR_current_density[q]*this->JxW_cell[q];
                }
        }

        //-- DEALII -> APPFRAME                    
         this->dealII_to_appframe(cell_residual,
                                  local_residual,
                                  this->residual_indices);
    }
    //-- OTHER LAYERS
    else
    {
        FcstUtilities::log << "Layer you specified is not FuelCellShop::Layer::CatalystLayer<dim>" << std::endl;
        AssertThrow( false , ExcInternalError() );
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
NAME::ReactionSourceTermsKG<dim>::make_assemblers_generic_constant_data()
{
       this->make_matrix_block_indices();
       this->make_residual_indices();
}


// ---                                     ---
// --- make_assemblers_cell_constant_data ---
// ---                                     ---

template<int dim>
void
NAME::ReactionSourceTermsKG<dim>::make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
       this->dofs_per_cell   = cell_info.get_fe_val_unsplit().dofs_per_cell;
       this->n_q_points_cell = cell_info.get_fe_val_unsplit().n_quadrature_points;
       
       this->JxW_cell.resize(this->n_q_points_cell);

       T_mixture.resize(this->n_q_points_cell);
               
       ORR_current_density.resize(this->n_q_points_cell);
       DORR_current_density_Doxygen_concentration.resize(this->n_q_points_cell);
       DORR_current_density_Delectronic_electrical_potential.resize(this->n_q_points_cell);
       DORR_current_density_Dprotonic_electrical_potential.resize(this->n_q_points_cell);

       HOR_current_density.resize(this->n_q_points_cell);
       DHOR_current_density_Dhydrogen_concentration.resize(this->n_q_points_cell);
       DHOR_current_density_Delectronic_electrical_potential.resize(this->n_q_points_cell);
       DHOR_current_density_Dprotonic_electrical_potential.resize(this->n_q_points_cell);

       density_old.resize(this->n_q_points_cell);
       electronic_electrical_potential_old.resize(this->n_q_points_cell);
       protonic_electrical_potential_old.resize(this->n_q_points_cell);

       phi_density.resize( n_species, std::vector< std::vector<double> >( this->n_q_points_cell, std::vector<double>( this->dofs_per_cell ) ) );
       phi_electronic_electrical_potential.resize( this->n_q_points_cell, std::vector<double>( this->dofs_per_cell ) );
       phi_protonic_electrical_potential.resize( this->n_q_points_cell, std::vector<double>( this->dofs_per_cell ) );

       last_iter_cell = cell_info.global_data->find_vector(this->solution_vector_name);
}

// ---                                     ---
// --- make_assemblers_cell_variable_data2 ---
// ---                                     ---

template<int dim>
void
NAME::ReactionSourceTermsKG<dim>::make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                    FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
    for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
        this->JxW_cell[q] = cell_info.get_fe_val_unsplit().JxW(q);
    }
    
    for(unsigned int g = 1; g <= n_species; ++g)
    {
        var_name = "density_" + var_postfixes[g-1];
        density_old[g-1] = cell_info.values[last_iter_cell][this->system_management->solution_name_to_index(var_name)];
    }
    
    electronic_electrical_potential_old = cell_info.values[last_iter_cell][this->system_management->solution_name_to_index("electronic_electrical_potential")];
    protonic_electrical_potential_old   = cell_info.values[last_iter_cell][this->system_management->solution_name_to_index("protonic_electrical_potential")];
    
    for(unsigned int q = 0; q < this->n_q_points_cell; ++q) 
    {
        for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
        {
            for(unsigned int g = 0; g < n_species; ++g)
            {
                phi_density[g][q][k] = cell_info.get_fe_val_unsplit()[ density_extractors[g] ].value(k,q);
            }
            
            phi_electronic_electrical_potential[q][k] = cell_info.get_fe_val_unsplit()[ electronic_electrical_potential_extractor ].value(k,q);
            phi_protonic_electrical_potential[q][k] = cell_info.get_fe_val_unsplit()[ protonic_electrical_potential_extractor ].value(k,q);
        }
    }
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    const std::type_info& info          = layer->get_base_type();
       
    if( info == CatalystLayer )
    {
        FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);
        
        //-- For the case of constant temperature, get the temperature from gas mixture:
        for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
            T_mixture[q] = ptr->get_gas_mixture()->get_temperature();
                     
        //-- HOR
        if( ptr->get_kinetics() == anode_kinetics   && ptr->get_kinetics()->get_reaction_name() == HOR  ||
            ptr->get_kinetics() == cathode_kinetics && ptr->get_kinetics()->get_reaction_name() == HOR )
        {
            std::vector<double> hydrogen_concentration_old(this->n_q_points_cell);
            
            if( indexH2 != -1 )
                for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                    hydrogen_concentration_old[q] = density_old[indexH2-1][q]/molar_mass[indexH2-1];
                else
                    AssertThrow( false, ExcInternalError() );
                
                std::vector<FuelCellShop::SolutionVariable> solution_variables;
            solution_variables.push_back( FuelCellShop::SolutionVariable( hydrogen_concentration_old,
                                                                          VariableNames::hydrogen_concentration) );
            
            solution_variables.push_back( FuelCellShop::SolutionVariable( electronic_electrical_potential_old,
                                                                          VariableNames::electronic_electrical_potential) );
            
            solution_variables.push_back( FuelCellShop::SolutionVariable( protonic_electrical_potential_old,
                                                                          VariableNames::protonic_electrical_potential) );
            
            std::vector<VariableNames> derivative_flags;
            derivative_flags.push_back(VariableNames::hydrogen_concentration);
            derivative_flags.push_back(VariableNames::electronic_electrical_potential);
            derivative_flags.push_back(VariableNames::protonic_electrical_potential);
            
            ptr->set_solution(solution_variables);
            ptr->set_cell_id(cell_info.dof_active_cell->index());
            ptr->current_density( HOR_current_density );
            
            std::map< VariableNames, std::vector<double> > Dcurrent;
            ptr->set_derivative_flags(derivative_flags);
            ptr->derivative_current_density(Dcurrent);
            
            DHOR_current_density_Dhydrogen_concentration = Dcurrent.at(VariableNames::hydrogen_concentration);
            
            for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                DHOR_current_density_Dhydrogen_concentration[q] *= (Constants::R()*T_mixture[q])/(ptr->get_electrolyte()->get_H_H2()*1.0e-6);
            
            DHOR_current_density_Delectronic_electrical_potential = Dcurrent.at(VariableNames::electronic_electrical_potential);
            DHOR_current_density_Dprotonic_electrical_potential   = Dcurrent.at(VariableNames::protonic_electrical_potential);
        }
        
        //-- ORR        
        if( ptr->get_kinetics() == anode_kinetics   && ptr->get_kinetics()->get_reaction_name() == ORR ||
            ptr->get_kinetics() == cathode_kinetics && ptr->get_kinetics()->get_reaction_name() == ORR )
        {
            std::vector<double> oxygen_concentration_old(this->n_q_points_cell);
            
            if( indexO2 != -1 )
                for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                    oxygen_concentration_old[q] = density_old[indexO2-1][q]/molar_mass[indexO2-1];
                else
                    AssertThrow( false, ExcInternalError() );
                
                std::vector<FuelCellShop::SolutionVariable> solution_variables;
            
            solution_variables.push_back( FuelCellShop::SolutionVariable( oxygen_concentration_old,
                                                                          VariableNames::oxygen_concentration) );
            
            solution_variables.push_back( FuelCellShop::SolutionVariable( electronic_electrical_potential_old,
                                                                          VariableNames::electronic_electrical_potential) );
            
            solution_variables.push_back( FuelCellShop::SolutionVariable( protonic_electrical_potential_old,
                                                                          VariableNames::protonic_electrical_potential) );
            
            std::vector<VariableNames> derivative_flags;
            derivative_flags.push_back(VariableNames::oxygen_concentration);
            derivative_flags.push_back(VariableNames::electronic_electrical_potential);
            derivative_flags.push_back(VariableNames::protonic_electrical_potential);
            
            ptr->set_solution(solution_variables);
            ptr->set_cell_id(cell_info.dof_active_cell->index());
            ptr->current_density( ORR_current_density );
            
            std::map< VariableNames, std::vector<double> > Dcurrent;
            ptr->set_derivative_flags(derivative_flags);
            ptr->derivative_current_density(Dcurrent);
            
            DORR_current_density_Doxygen_concentration = Dcurrent.at(VariableNames::oxygen_concentration);
            
            for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                DORR_current_density_Doxygen_concentration[q] *= (Constants::R()*T_mixture[q])/(ptr->get_electrolyte()->get_H_O2()*1.0e-6);
            
            DORR_current_density_Delectronic_electrical_potential = Dcurrent.at(VariableNames::electronic_electrical_potential);
            DORR_current_density_Dprotonic_electrical_potential   = Dcurrent.at(VariableNames::protonic_electrical_potential);
        }
    }
    
    //////////////////
    // OTHER LAYERS //
    //////////////////    
    else
    {
        FcstUtilities::log << "Layer you specified is not FuelCellShop::Layer::CatalystLayer<dim>" << std::endl;
        AssertThrow( false , ExcInternalError() );
    }
}

// ---                     ---
// --- print_equation_info ---
// ---                     ---

template<int dim>
void
NAME::ReactionSourceTermsKG<dim>::print_equation_info() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "-------------------------------------------------------------------------------" << std::endl;
    FcstUtilities::log << std::endl;



    FcstUtilities::log << std::endl;
    FcstUtilities::log << "-------------------------------------------------------------------------------" << std::endl;
    FcstUtilities::log << std::endl;
}

/////////////////////////////////////////////////////////////
//
// Functions necessary for dealII_to_appframe routines to 
// split local matrix and residual into blocks.
//
/////////////////////////////////////////////////////////////
// ---                           ---
// --- make_matrix_block_indices ---
// ---                           ---
template<int dim>
void
NAME::ReactionSourceTermsKG<dim>::make_matrix_block_indices()
{

       unsigned int index;

       //-- KERKHOF EQUATIONS
       for(unsigned int g = 1; g <= n_species; ++g)
       {
           eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[g-1];           
           
           if( g == indexO2 )
           {
               if( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR
                   ||
                   cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )
               {
                   //-- DENSITY OF OXYGEN
                   var_name = "density_" + var_postfixes[g-1];
                   index = this->system_management->matrix_block_index(eq_name, var_name);
                   this->matrix_block_indices.push_back(index);
               }
           }           
           if( g == indexH2 )
           {
               if( anode_kinetics   && anode_kinetics->get_reaction_name()   == HOR
                   ||
                   cathode_kinetics && cathode_kinetics->get_reaction_name() == HOR )
               {
                   eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[g-1];
                   
                   // DENSITY OF HYDROGEN 
                   var_name = "density_" + var_postfixes[g-1];
                   index = this->system_management->matrix_block_index(eq_name, var_name);
                   this->matrix_block_indices.push_back(index);
                   
               }
           }
           if( g == indexH2O )
           {
               if( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR
                   ||
                   cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )
               {                  
                   //-- DENSITY OF OXYGEN 
                   var_name = "density_" + var_postfixes[indexO2-1];
                   index = this->system_management->matrix_block_index(eq_name, var_name);
                   this->matrix_block_indices.push_back(index);
                   //-- Density of Water vapour
                   var_name = "density_" + var_postfixes[indexH2O-1];
                   index = this->system_management->matrix_block_index(eq_name, var_name);
                   this->matrix_block_indices.push_back(index);   
               }
           }
           //-- All equatiions that have reactions are also coupled with potentials:
           if ( (g == indexH2 && ( anode_kinetics   && anode_kinetics->get_reaction_name()   == HOR || cathode_kinetics && cathode_kinetics->get_reaction_name() == HOR )) ||
                ( (g == indexO2 || g == indexH2O) && ( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR || cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )) )
           {
               //-- Fs 
               var_name = "electronic_electrical_potential";
               index = this->system_management->matrix_block_index(eq_name, var_name);
               this->matrix_block_indices.push_back(index);
               //-- Fm 
               var_name = "protonic_electrical_potential";
               index = this->system_management->matrix_block_index(eq_name, var_name);
               this->matrix_block_indices.push_back(index);
           }
       }

       /////////////////////////////////
       // ELECTRON TRANSPORT EQUATION //
       /////////////////////////////////
       eq_name = "Electron Transport Equation";

       if( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR
           ||
           cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )
       {
              var_name = "density_" + var_postfixes[indexO2-1];
              index = this->system_management->matrix_block_index(eq_name, var_name);
              this->matrix_block_indices.push_back(index);
       }

       if( anode_kinetics   && anode_kinetics->get_reaction_name()   == HOR
           ||
           cathode_kinetics && cathode_kinetics->get_reaction_name() == HOR )
       {
              var_name = "density_" + var_postfixes[indexH2-1];
              index = this->system_management->matrix_block_index(eq_name, var_name);
              this->matrix_block_indices.push_back(index);
       }

       var_name = "electronic_electrical_potential";
       index = this->system_management->matrix_block_index(eq_name, var_name);
       this->matrix_block_indices.push_back(index);

       var_name = "protonic_electrical_potential";
       index = this->system_management->matrix_block_index(eq_name, var_name);
       this->matrix_block_indices.push_back(index);

       ///////////////////////////////
       // PROTON TRANSPORT EQUATION //
       ///////////////////////////////
       eq_name = "Proton Transport Equation";

       if( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR
           ||
           cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )
       {
              var_name = "density_" + var_postfixes[indexO2-1];
              index = this->system_management->matrix_block_index(eq_name, var_name);
              this->matrix_block_indices.push_back(index);
       }

       if( anode_kinetics   && anode_kinetics->get_reaction_name()   == HOR
           ||
           cathode_kinetics && cathode_kinetics->get_reaction_name() == HOR )
       {
              var_name = "density_" + var_postfixes[indexH2-1];
              index = this->system_management->matrix_block_index(eq_name, var_name);
              this->matrix_block_indices.push_back(index);
       }

       var_name = "electronic_electrical_potential";
       index = this->system_management->matrix_block_index(eq_name, var_name);
       this->matrix_block_indices.push_back(index);

       var_name = "protonic_electrical_potential";
       index = this->system_management->matrix_block_index(eq_name, var_name);
       this->matrix_block_indices.push_back(index);
}

// ---                       ---
// --- make_residual_indices ---
// ---                       ---

template<int dim>
void
NAME::ReactionSourceTermsKG<dim>::make_residual_indices()
{
    unsigned int index;
    
    //-- KERKHOF EQUATIONS 
    for(unsigned int g = 1; g <= n_species; ++g)
    {
        
        if( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR
            ||
            cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )
        {
            if( g == indexO2 )
            {
                eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[g-1];
                index = this->system_management->equation_name_to_index(eq_name);
                this->residual_indices.push_back(index);
            }
            if( g == indexH2O )
            {
                eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[g-1];
                index = this->system_management->equation_name_to_index(eq_name);
                this->residual_indices.push_back(index);
            }
        }
        
        
        if( anode_kinetics   && anode_kinetics->get_reaction_name()   == HOR
            ||
            cathode_kinetics && cathode_kinetics->get_reaction_name() == HOR )
        {
            if( g == indexH2 )
            {
                eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[g-1];
                index = this->system_management->equation_name_to_index(eq_name);
                this->residual_indices.push_back(index);
            }
        }
    }

    // ELECTRON TRANSPORT EQUATION 
    eq_name = "Electron Transport Equation";
    index = this->system_management->equation_name_to_index(eq_name);
    this->residual_indices.push_back(index);
    
    // PROTON TRANSPORT EQUATION 
    eq_name = "Proton Transport Equation";
    index = this->system_management->equation_name_to_index(eq_name);
    this->residual_indices.push_back(index);
}

// ---                           ---
// ---  EXPLICIT INSTANTIATIONS  ---
// ---                           ---

template class NAME::ReactionSourceTermsKG<deal_II_dimension>;