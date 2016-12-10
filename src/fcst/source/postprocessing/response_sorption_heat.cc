// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: response_sorption_heat.cc
// - Description: It contains definition sorption heat response evaluator class.
// - Developers: Madhur Bhaiya
// - $Id:
//
// ----------------------------------------------------------------------------

#include <postprocessing/response_sorption_heat.h>

namespace NAME = FuelCellShop::PostProcessing;

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::SorptionHeatResponse<dim>::initialize(ParameterHandler&param )
{    
    if ( this->system_management->solution_in_userlist("water_molar_fraction") )
    {
        xWater.solution_index = this->system_management->solution_name_to_index("water_molar_fraction"); 
        xWater.fetype_index = this->system_management->block_info->base_element[xWater.solution_index];
        xWater.indices_exist = true;
    }
    else
        throw std::runtime_error("water_molar_fraction variable is required for SorptionHeatResponse");
        
    if ( this->system_management->solution_in_userlist("membrane_water_content") )
    {
        lambda.solution_index = this->system_management->solution_name_to_index("membrane_water_content"); 
        lambda.fetype_index = this->system_management->block_info->base_element[lambda.solution_index];
        lambda.indices_exist = true;
    }
    else
        throw std::runtime_error("membrane_water_content variable is required for SorptionHeatResponse");
        
    if ( this->system_management->solution_in_userlist("temperature_of_REV") )
    {
        tRev.solution_index = this->system_management->solution_name_to_index("temperature_of_REV"); 
        tRev.fetype_index = this->system_management->block_info->base_element[tRev.solution_index];
        tRev.indices_exist = true;
    }
    else
        throw std::runtime_error("temperature_of_REV variable is required for SorptionHeatResponse");
        
    time_constant = sorption_source->get_time_constant();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::SorptionHeatResponse<dim>::compute_responses(const typename DoFApplication<dim>::CellInfo& info,
                                                        FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                                        std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& respMap) const
{
    respMap.clear();
    
    // Get the total number of quadrature points in the cell
    unsigned int n_quad_cell = (info.fe(xWater.fetype_index)).n_quadrature_points;
    
    unsigned int solIndex = info.global_data->find_vector("Solution");
    
    // Determining factor and computing effective properties based on which layer we are in.      
    const std::type_info& base_layer = layer->get_base_type();

    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    
    try
    {
        if (base_layer == CatalystLayer)
        {
            if ( sorption_source->get_flag_sorp_heat_cl() == false )
            {
                respMap[FuelCellShop::PostProcessing::ResponsesNames::sorption_heat] = 0.0;
                return;
            }
            else
            {
                FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);
                
                double rhoDry = ptr->get_electrolyte()->get_density();
                double EW = ptr->get_electrolyte()->get_EW();
                ptr->get_electrolyte()->set_water_molar_fraction( FuelCellShop::SolutionVariable(&info.values[solIndex][xWater.solution_index], water_molar_fraction) );
                ptr->get_electrolyte()->set_temperature( FuelCellShop::SolutionVariable(&info.values[solIndex][tRev.solution_index], temperature_of_REV) );
                std::vector<double> lambdaVec( info.values[solIndex][lambda.solution_index] );
                std::vector<double> lambdaEq( n_quad_cell, 0.0 );
                ptr->get_electrolyte()->sorption_isotherm(lambdaEq);
                std::vector<double> hSorp(n_quad_cell, 0.0);
                ptr->get_electrolyte()->sorption_enthalpy(hSorp);
                
                for (unsigned int q = 0; q < n_quad_cell; ++q)
                    respMap[FuelCellShop::PostProcessing::ResponsesNames::sorption_heat] += ( ((time_constant*rhoDry)/EW)*(lambdaEq[q]-lambdaVec[q])*hSorp[q]*info.fe(xWater.fetype_index).JxW(q) ) ;
            }
        }
        else
            AssertThrow(false, ExcMessage("SorptionHeatResponse can only be used with the CatalystLayer object"));
    
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
template class NAME::SorptionHeatResponse<deal_II_dimension>;
