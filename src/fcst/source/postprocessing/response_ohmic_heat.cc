// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: response_ohmic_heat.cc
// - Description: It contains definitions of electronic and protonic ohmic heat response evaluator classes.
// - Developers: Madhur Bhaiya
// - $Id:
//
// ----------------------------------------------------------------------------

#include <postprocessing/response_ohmic_heat.h>

namespace NAME = FuelCellShop::PostProcessing;

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::ElectronOhmicHeatResponse<dim>::initialize(ParameterHandler&param )
{    
    if ( this->system_management->solution_in_userlist("electronic_electrical_potential") )
    {
        phiS.solution_index = this->system_management->solution_name_to_index("electronic_electrical_potential"); 
        phiS.fetype_index = this->system_management->block_info->base_element[phiS.solution_index];
        phiS.indices_exist = true;
    }
    else
        throw std::runtime_error("electronic_electrical_potential variable is required for ElectronOhmicHeatResponse");
    
    if ( thermal_equation->get_electron_ohmic_heat_gdl() )
        factor_GDL = 1;
    
    if ( thermal_equation->get_electron_ohmic_heat_mpl() )
        factor_MPL = 1;
        
    if ( thermal_equation->get_electron_ohmic_heat_cl() )
        factor_CL = 1;
    
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::ElectronOhmicHeatResponse<dim>::compute_responses(const typename DoFApplication<dim>::CellInfo& info,
                                                        FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                                        std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& respMap) const
{
    respMap.clear();
    Tensor<2,dim> sigmaSeff;
    
    // Get the total number of quadrature points in the cell
    unsigned int n_quad_cell = (info.fe(phiS.fetype_index)).n_quadrature_points;
    
    unsigned int solIndex = info.global_data->find_vector("Solution");
    
    // A temporary variable which will store the factor values based on the layer we are in.
    unsigned int factor_cell = 0;
    
    // Determining factor and computing effective properties based on which layer we are in.      
    const std::type_info& base_layer = layer->get_base_type();
    
    const std::type_info& GasDiffusionLayer = typeid(FuelCellShop::Layer::GasDiffusionLayer<dim>);
    const std::type_info& MicroPorousLayer  = typeid(FuelCellShop::Layer::MicroPorousLayer<dim>);
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    
    try
    {
        if (base_layer == GasDiffusionLayer)
        {
            FuelCellShop::Layer::GasDiffusionLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::GasDiffusionLayer<dim>* >(layer);
            ptr->effective_electron_conductivity(sigmaSeff);
            
            factor_cell = factor_GDL;
        }
        else if (base_layer == MicroPorousLayer)
        {
            FuelCellShop::Layer::MicroPorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::MicroPorousLayer<dim>* >(layer);
            ptr->effective_electron_conductivity(sigmaSeff);
            
            factor_cell = factor_MPL;
        }
        else if (base_layer == CatalystLayer)
        {
            FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);
            double sigmaSeff_iso;
            ptr->effective_electron_conductivity(sigmaSeff_iso);
            for (unsigned int i = 0; i < dim; ++i)
                sigmaSeff[i][i] = sigmaSeff_iso;
            
            factor_cell = factor_CL;
        }
        else
            AssertThrow(false, ExcMessage("ElectronOhmicHeatResponse can only be used with either GasDiffusionLayer, MicroPorousLayer or CatalystLayer object"));
    
    }
    catch(const std::bad_cast& e)
    {
        const std::type_info& info = typeid(*layer);  
        FcstUtilities::log << "Object of type "<<info.name()<<" not implemented"<< std::endl;
        FcstUtilities::log << e.what() << std::endl;
    }
    
    for (unsigned int q=0; q<n_quad_cell; ++q)
    {
        respMap[FuelCellShop::PostProcessing::ResponsesNames::electron_ohmic_heat] += ( (info.gradients[solIndex][phiS.solution_index][q] * (sigmaSeff * info.gradients[solIndex][phiS.solution_index][q])) 
                                                                                     * info.fe(phiS.fetype_index).JxW(q) * factor_cell );
    }
    
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

template <int dim>
void 
NAME::ProtonOhmicHeatResponse<dim>::initialize(ParameterHandler&param )
{    
    if ( this->system_management->solution_in_userlist("protonic_electrical_potential") )
    {
        phiM.solution_index = this->system_management->solution_name_to_index("protonic_electrical_potential"); 
        phiM.fetype_index = this->system_management->block_info->base_element[phiM.solution_index];
        phiM.indices_exist = true;
    }
    else
        throw std::runtime_error("protonic_electrical_potential variable is required for ProtonOhmicHeatResponse");
        
    if ( this->system_management->solution_in_userlist("temperature_of_REV") )
    {
        tRev.solution_index = this->system_management->solution_name_to_index("temperature_of_REV"); 
        tRev.fetype_index = this->system_management->block_info->base_element[tRev.solution_index];
        tRev.indices_exist = true;
    }
    else
        throw std::runtime_error("temperature_of_REV variable is required for ProtonOhmicHeatResponse");
        
    if ( this->system_management->solution_in_userlist("membrane_water_content") )
    {
        lambda.solution_index = this->system_management->solution_name_to_index("membrane_water_content"); 
        lambda.fetype_index = this->system_management->block_info->base_element[lambda.solution_index];
        lambda.indices_exist = true;
    } 
    
    if ( thermal_equation->get_proton_ohmic_heat_cl() )
        factor_CL = 1;
    
    if ( thermal_equation->get_proton_ohmic_heat_ml() )
        factor_ML = 1;    
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::ProtonOhmicHeatResponse<dim>::compute_responses(const typename DoFApplication<dim>::CellInfo& info,
                                                        FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                                        std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& respMap) const
{
    respMap.clear();
    
    // Get the total number of quadrature points in the cell
    unsigned int n_quad_cell = (info.fe(phiM.fetype_index)).n_quadrature_points;
    
    std::vector<double> sigmaMeff(n_quad_cell, 0.0);
    
    unsigned int solIndex = info.global_data->find_vector("Solution");
    
    // A temporary variable which will store the factor values based on the layer we are in.
    unsigned int factor_cell = 0;
    
    // Determining factor and computing effective properties based on which layer we are in.      
    const std::type_info& base_layer = layer->get_base_type();
    
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    const std::type_info& MembraneLayer     = typeid(FuelCellShop::Layer::MembraneLayer<dim>);
    
    try
    {
        if (base_layer == CatalystLayer)
        {
            FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);
            
            ptr->get_electrolyte()->set_temperature( FuelCellShop::SolutionVariable(&info.values[solIndex][tRev.solution_index], temperature_of_REV) );
            if (lambda.indices_exist)
                ptr->get_electrolyte()->set_membrane_water_content( FuelCellShop::SolutionVariable(&info.values[solIndex][lambda.solution_index], membrane_water_content) );
                
            ptr->effective_proton_conductivity(sigmaMeff);
            factor_cell = factor_CL;

        }
        else if (base_layer == MembraneLayer)
        {
            FuelCellShop::Layer::MembraneLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::MembraneLayer<dim>* >(layer);
            
            ptr->get_electrolyte()->set_temperature( FuelCellShop::SolutionVariable(&info.values[solIndex][tRev.solution_index], temperature_of_REV) );
            if (lambda.indices_exist)
                ptr->get_electrolyte()->set_membrane_water_content( FuelCellShop::SolutionVariable(&info.values[solIndex][lambda.solution_index], membrane_water_content) );
                
            ptr->effective_proton_conductivity(sigmaMeff);            
            factor_cell = factor_ML;
        }
        else
            AssertThrow(false, ExcMessage("ProtonOhmicHeatResponse can only be used with either Membrane or CatalystLayer object"));
    
    }
    catch(const std::bad_cast& e)
    {
        const std::type_info& info = typeid(*layer);  
        FcstUtilities::log << "Object of type "<<info.name()<<" not implemented"<< std::endl;
        FcstUtilities::log << e.what() << std::endl;
    }
    
    for (unsigned int q=0; q<n_quad_cell; ++q)
    {
        respMap[FuelCellShop::PostProcessing::ResponsesNames::proton_ohmic_heat] += ( (info.gradients[solIndex][phiM.solution_index][q] * (sigmaMeff[q] * info.gradients[solIndex][phiM.solution_index][q])) 
                                                                                     * info.fe(phiM.fetype_index).JxW(q) * factor_cell );
    }
    
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::ElectronOhmicHeatResponse<deal_II_dimension>;
template class NAME::ProtonOhmicHeatResponse<deal_II_dimension>;
