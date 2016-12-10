// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: response_reaction_heat.cc
// - Description: It contains definition of reaction heat response evaluator classes.
// - Developers: Madhur Bhaiya
// - $Id:
//
// ----------------------------------------------------------------------------

#include <postprocessing/response_reaction_heat.h>

namespace NAME = FuelCellShop::PostProcessing;

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::ORRReactionHeatResponse<dim>::initialize(ParameterHandler&param )
{    
    if ( this->system_management->solution_in_userlist("oxygen_molar_fraction") )
    {
        xOxygen.solution_index = this->system_management->solution_name_to_index("oxygen_molar_fraction"); 
        xOxygen.fetype_index = this->system_management->block_info->base_element[xOxygen.solution_index];
        xOxygen.indices_exist = true;
    }
    else
        throw std::runtime_error("oxygen_molar_fraction variable is required for ORRReactionHeatResponse");
        
    if ( this->system_management->solution_in_userlist("electronic_electrical_potential") )
    {
        phiS.solution_index = this->system_management->solution_name_to_index("electronic_electrical_potential"); 
        phiS.fetype_index = this->system_management->block_info->base_element[phiS.solution_index];
        phiS.indices_exist = true;
    }
    else
        throw std::runtime_error("electronic_electrical_potential variable is required for ORRReactionHeatResponse");
        
    if ( this->system_management->solution_in_userlist("protonic_electrical_potential") )
    {
        phiM.solution_index = this->system_management->solution_name_to_index("protonic_electrical_potential"); 
        phiM.fetype_index = this->system_management->block_info->base_element[phiM.solution_index];
        phiM.indices_exist = true;
    }
    else
        throw std::runtime_error("protonic_electrical_potential variable is required for ORRReactionHeatResponse");
        
    if ( this->system_management->solution_in_userlist("temperature_of_REV") )
    {
        tRev.solution_index = this->system_management->solution_name_to_index("temperature_of_REV"); 
        tRev.fetype_index = this->system_management->block_info->base_element[tRev.solution_index];
        tRev.indices_exist = true;
    }
    else
        throw std::runtime_error("temperature_of_REV variable is required for ORRReactionHeatResponse");
        
    reaction_heat = new FuelCellShop::Equation::ReactionHeat;
    reaction_heat->set_kinetics( reaction_source->get_cathode_kinetics() );
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::ORRReactionHeatResponse<dim>::compute_responses(const typename DoFApplication<dim>::CellInfo& info,
                                                      FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                                      std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& respMap) const
{
    respMap.clear();
    
    // Get the total number of quadrature points in the cell
    unsigned int n_quad_cell = (info.fe(phiS.fetype_index)).n_quadrature_points;
    
    unsigned int solIndex = info.global_data->find_vector("Solution");
    
    // Determining responses based on which layer we are in.      
    const std::type_info& base_layer = layer->get_base_type();

    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    
    try
    {
        if (base_layer == CatalystLayer)
        {
            FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);
            ptr->set_cell_id(info.dof_active_cell->index());
            
            // Creating a vector to store current and a map to store heat values.
            std::vector<double> current_vec;
            std::map< FuelCellShop::PostProcessing::ResponsesNames, std::vector<double> > heat_map;
            
            // Creating solution variable vector to passed to catalyst layer classes
            std::vector< FuelCellShop::SolutionVariable > solution_variables;
            solution_variables.push_back( FuelCellShop::SolutionVariable(&info.values[solIndex][xOxygen.solution_index], oxygen_molar_fraction) );
            solution_variables.push_back( FuelCellShop::SolutionVariable(&info.values[solIndex][phiM.solution_index], protonic_electrical_potential) );
            solution_variables.push_back( FuelCellShop::SolutionVariable(&info.values[solIndex][phiS.solution_index], electronic_electrical_potential) );
            solution_variables.push_back( FuelCellShop::SolutionVariable(&info.values[solIndex][tRev.solution_index], temperature_of_REV) );
            if ( this->system_management->solution_in_userlist("membrane_water_content") )
                solution_variables.push_back(FuelCellShop::SolutionVariable(&info.values[solIndex][this->system_management->solution_name_to_index("membrane_water_content")], membrane_water_content));
                
            // Setting solution in kinetics and reaction heat objects.
            ptr->set_solution(solution_variables);
            reaction_heat->set_electrolyte_potential( FuelCellShop::SolutionVariable(&info.values[solIndex][phiM.solution_index], protonic_electrical_potential) );
            reaction_heat->set_solid_potential( FuelCellShop::SolutionVariable(&info.values[solIndex][phiS.solution_index], electronic_electrical_potential) );
            reaction_heat->set_temperature( FuelCellShop::SolutionVariable(&info.values[solIndex][tRev.solution_index], temperature_of_REV) );
            
            // Computing current density
            ptr->current_density(current_vec);
            
            // Computing various reaction heat responses by initializing factors and then filling the heat_map
            std::vector<double> response_vec; // A temporary vector to store computed heat values before transferring them to heat_map
            
            // 1. ORR_reaction_heat
            reaction_heat->initialize_factors(reaction_source->get_irrev_heat_ccl(),
                                              reaction_source->get_irrev_heat_acl(),
                                              reaction_source->get_rev_heat(),
                                              reaction_source->get_factor_rev_heat_ccl(),
                                              reaction_source->get_water_vap_heat_ccl());
            reaction_heat->heat_source(response_vec, current_vec);
            heat_map[FuelCellShop::PostProcessing::ResponsesNames::ORR_reaction_heat] = response_vec;
            
            // 2. ORR_irrev_heat
            reaction_heat->initialize_factors(reaction_source->get_irrev_heat_ccl(),
                                              false,
                                              false,
                                              0.0,
                                              false);
            reaction_heat->heat_source(response_vec, current_vec);
            heat_map[FuelCellShop::PostProcessing::ResponsesNames::ORR_irrev_heat] = response_vec;
            
            // 3. ORR_rev_heat
            reaction_heat->initialize_factors(false,
                                              false,
                                              reaction_source->get_rev_heat(),
                                              reaction_source->get_factor_rev_heat_ccl(),
                                              false);
            reaction_heat->heat_source(response_vec, current_vec);
            heat_map[FuelCellShop::PostProcessing::ResponsesNames::ORR_rev_heat] = response_vec;
            
            // 4. ORR_watervap_heat
            reaction_heat->initialize_factors(false,
                                              false,
                                              false,
                                              0.0,
                                              reaction_source->get_water_vap_heat_ccl());
            reaction_heat->heat_source(response_vec, current_vec);
            heat_map[FuelCellShop::PostProcessing::ResponsesNames::ORR_watervap_heat] = response_vec;

            // Integrating various responses
            for (unsigned int q = 0; q < n_quad_cell; ++q)
            {
                respMap[FuelCellShop::PostProcessing::ResponsesNames::ORR_reaction_heat] += (heat_map[FuelCellShop::PostProcessing::ResponsesNames::ORR_reaction_heat][q]*info.fe(phiS.fetype_index).JxW(q));
                respMap[FuelCellShop::PostProcessing::ResponsesNames::ORR_irrev_heat] += (heat_map[FuelCellShop::PostProcessing::ResponsesNames::ORR_irrev_heat][q]*info.fe(phiS.fetype_index).JxW(q));
                respMap[FuelCellShop::PostProcessing::ResponsesNames::ORR_rev_heat] += (heat_map[FuelCellShop::PostProcessing::ResponsesNames::ORR_rev_heat][q]*info.fe(phiS.fetype_index).JxW(q));
                respMap[FuelCellShop::PostProcessing::ResponsesNames::ORR_watervap_heat] += (heat_map[FuelCellShop::PostProcessing::ResponsesNames::ORR_watervap_heat][q]*info.fe(phiS.fetype_index).JxW(q));
            }
        }
        else
            AssertThrow(false, ExcMessage("ORRReactionHeatResponse can only be used with the CatalystLayer object"));
    }
    catch(const std::bad_cast& e)
    {
        const std::type_info& info = typeid(*layer);  
        FcstUtilities::log << "Object of type "<<info.name()<<" not implemented"<< std::endl;
        FcstUtilities::log << e.what() << std::endl;
    }    
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


template <int dim>
void 
NAME::HORReactionHeatResponse<dim>::initialize(ParameterHandler&param )
{
    if ( this->system_management->solution_in_userlist("electronic_electrical_potential") )
    {
        phiS.solution_index = this->system_management->solution_name_to_index("electronic_electrical_potential"); 
        phiS.fetype_index = this->system_management->block_info->base_element[phiS.solution_index];
        phiS.indices_exist = true;
    }
    else
        throw std::runtime_error("electronic_electrical_potential variable is required for HORReactionHeatResponse");
        
    if ( this->system_management->solution_in_userlist("protonic_electrical_potential") )
    {
        phiM.solution_index = this->system_management->solution_name_to_index("protonic_electrical_potential"); 
        phiM.fetype_index = this->system_management->block_info->base_element[phiM.solution_index];
        phiM.indices_exist = true;
    }
    else
        throw std::runtime_error("protonic_electrical_potential variable is required for HORReactionHeatResponse");
        
    if ( this->system_management->solution_in_userlist("temperature_of_REV") )
    {
        tRev.solution_index = this->system_management->solution_name_to_index("temperature_of_REV"); 
        tRev.fetype_index = this->system_management->block_info->base_element[tRev.solution_index];
        tRev.indices_exist = true;
    }
    else
        throw std::runtime_error("temperature_of_REV variable is required for HORReactionHeatResponse");
        
    if ( this->system_management->solution_in_userlist("hydrogen_molar_fraction") )
    {
        xHydrogen.solution_index = this->system_management->solution_name_to_index("hydrogen_molar_fraction"); 
        xHydrogen.fetype_index = this->system_management->block_info->base_element[xHydrogen.solution_index];
        xHydrogen.indices_exist = true;
    }
    else if ( this->system_management->solution_in_userlist("water_molar_fraction") )
    {
        xWater.solution_index = this->system_management->solution_name_to_index("water_molar_fraction"); 
        xWater.fetype_index = this->system_management->block_info->base_element[xWater.solution_index];
        xWater.indices_exist = true;
    }
    else
        throw std::runtime_error("Either hydrogen_molar_fraction or water_molar_fraction, as an variable is required for HORReactionHeatResponse");
        
    reaction_heat = new FuelCellShop::Equation::ReactionHeat;
    reaction_heat->set_kinetics( reaction_source->get_anode_kinetics() );
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::HORReactionHeatResponse<dim>::compute_responses(const typename DoFApplication<dim>::CellInfo& info,
                                                      FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                                      std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& respMap) const
{
    respMap.clear();
    
    // Get the total number of quadrature points in the cell
    unsigned int n_quad_cell = (info.fe(phiS.fetype_index)).n_quadrature_points;
    
    unsigned int solIndex = info.global_data->find_vector("Solution");
    
    // Determining responses based on which layer we are in.      
    const std::type_info& base_layer = layer->get_base_type();

    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    
    try
    {
        if (base_layer == CatalystLayer)
        {
            FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);
            ptr->set_cell_id(info.dof_active_cell->index());
            
            // Creating a vector to store current and a map to store heat values.
            std::vector<double> current_vec;
            std::map< FuelCellShop::PostProcessing::ResponsesNames, std::vector<double> > heat_map;
            
            // Creating solution variable vector to passed to catalyst layer classes
            std::vector<double> x_H2(n_quad_cell, 0.0);
            if (xHydrogen.indices_exist)
            {
                x_H2 = info.values[solIndex][xHydrogen.solution_index];
            }
            else if (xWater.indices_exist)
            {
                for(unsigned int q = 0; q < n_quad_cell; ++q)
                    x_H2[q] = 1.0 - info.values[solIndex][xWater.solution_index][q];
            }
            else // Some bug otherwise it should not enter here
                AssertThrow( false, ExcNotImplemented() );
            
            std::vector< FuelCellShop::SolutionVariable > solution_variables;
            solution_variables.push_back( FuelCellShop::SolutionVariable(&x_H2, hydrogen_molar_fraction) );
            solution_variables.push_back( FuelCellShop::SolutionVariable(&info.values[solIndex][phiM.solution_index], protonic_electrical_potential) );
            solution_variables.push_back( FuelCellShop::SolutionVariable(&info.values[solIndex][phiS.solution_index], electronic_electrical_potential) );
            solution_variables.push_back( FuelCellShop::SolutionVariable(&info.values[solIndex][tRev.solution_index], temperature_of_REV) );
            if ( this->system_management->solution_in_userlist("membrane_water_content") )
                solution_variables.push_back(FuelCellShop::SolutionVariable(&info.values[solIndex][this->system_management->solution_name_to_index("membrane_water_content")], membrane_water_content));
                
            // Setting solution in kinetics and reaction heat objects.
            ptr->set_solution(solution_variables);
            reaction_heat->set_electrolyte_potential( FuelCellShop::SolutionVariable(&info.values[solIndex][phiM.solution_index], protonic_electrical_potential) );
            reaction_heat->set_solid_potential( FuelCellShop::SolutionVariable(&info.values[solIndex][phiS.solution_index], electronic_electrical_potential) );
            reaction_heat->set_temperature( FuelCellShop::SolutionVariable(&info.values[solIndex][tRev.solution_index], temperature_of_REV) );
            
            // Computing current density
            ptr->current_density(current_vec);
            
            // Computing various reaction heat responses by initializing factors and then filling the heat_map
            std::vector<double> response_vec; // A temporary vector to store computed heat values before transferring them to heat_map
            
            // 1. HOR_reaction_heat
            reaction_heat->initialize_factors(reaction_source->get_irrev_heat_ccl(),
                                              reaction_source->get_irrev_heat_acl(),
                                              reaction_source->get_rev_heat(),
                                              reaction_source->get_factor_rev_heat_ccl(),
                                              reaction_source->get_water_vap_heat_ccl());
            reaction_heat->heat_source(response_vec, current_vec);
            heat_map[FuelCellShop::PostProcessing::ResponsesNames::HOR_reaction_heat] = response_vec;
            
            // 2. HOR_irrev_heat
            reaction_heat->initialize_factors(false,
                                              reaction_source->get_irrev_heat_acl(),
                                              false,
                                              0.0,
                                              false);
            reaction_heat->heat_source(response_vec, current_vec);
            heat_map[FuelCellShop::PostProcessing::ResponsesNames::HOR_irrev_heat] = response_vec;
            
            // 3. HOR_rev_heat
            reaction_heat->initialize_factors(false,
                                              false,
                                              reaction_source->get_rev_heat(),
                                              reaction_source->get_factor_rev_heat_ccl(),
                                              false);
            reaction_heat->heat_source(response_vec, current_vec);
            heat_map[FuelCellShop::PostProcessing::ResponsesNames::HOR_rev_heat] = response_vec;

            // Integrating various responses
            for (unsigned int q = 0; q < n_quad_cell; ++q)
            {
                respMap[FuelCellShop::PostProcessing::ResponsesNames::HOR_reaction_heat] += (heat_map[FuelCellShop::PostProcessing::ResponsesNames::HOR_reaction_heat][q]*info.fe(phiS.fetype_index).JxW(q));
                respMap[FuelCellShop::PostProcessing::ResponsesNames::HOR_irrev_heat] += (heat_map[FuelCellShop::PostProcessing::ResponsesNames::HOR_irrev_heat][q]*info.fe(phiS.fetype_index).JxW(q));
                respMap[FuelCellShop::PostProcessing::ResponsesNames::HOR_rev_heat] += (heat_map[FuelCellShop::PostProcessing::ResponsesNames::HOR_rev_heat][q]*info.fe(phiS.fetype_index).JxW(q));
            }
        }
        else
            AssertThrow(false, ExcMessage("HORReactionHeatResponse can only be used with the CatalystLayer object"));
    }
    catch(const std::bad_cast& e)
    {
        const std::type_info& info = typeid(*layer);  
        FcstUtilities::log << "Object of type "<<info.name()<<" not implemented"<< std::endl;
        FcstUtilities::log << e.what() << std::endl;
    }    
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::ORRReactionHeatResponse<deal_II_dimension>;
template class NAME::HORReactionHeatResponse<deal_II_dimension>;