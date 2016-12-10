//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: solid_layer.cc
//    - Description: Solid layer class
//    - Developers: M. Secanell and Jie Zhou
//    - $Id: solid_layer.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <solid_layer.h>

namespace NAME = FuelCellShop::Layer;

//---------------------------------------------------------------------------

template <int dim>
NAME::SolidLayer<dim>::SolidLayer(const std::string &section_layer_name, ParameterHandler& param)
: NAME::BaseLayer<dim> (section_layer_name)
{
    temperature_enum = VariableNames::temperature_of_REV;
    this->  Solutions[temperature_enum] = T_vector;
    set_solution(Solutions);
    solution_iterator = Solutions.find(temperature_enum);
    solid = FuelCellShop::Material::PureSolid::create_PureSolid (section_layer_name, param);
}

//---------------------------------------------------------------------------
template <int dim>
NAME::SolidLayer<dim>::~SolidLayer()
{}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::SolidLayer<dim>::declare_parameters (const std::string& name,ParameterHandler &param) const
{
    NAME::BaseLayer<dim>::declare_parameters(param);
    
    std::string list_of_materials;
    
    for (typename FuelCellShop::Material::PureSolid::_mapFactory::iterator iterator = FuelCellShop::Material::PureSolid::get_mapFactory()->begin(); 
         iterator != FuelCellShop::Material::PureSolid::get_mapFactory()->end(); 
         iterator++)
         {
            list_of_materials += iterator->first;
            if (iterator != FuelCellShop::Material::PureSolid::get_mapFactory()->begin())
                list_of_materials += "|";
         }
         
         
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(this->name);
        {
            
            param.declare_entry ("Solid material type",
                                 "Graphite",
                                 Patterns::Selection(list_of_materials),
                                 "Solid Material");
                                 
                                 //-- Geometry (here or in the mesh)
                                 param.declare_entry ("Boundary id of the solid-ionomer interface",
                                 "0",
                                 Patterns::Integer(),
                                 "ID used to read mesh");
                                 param.declare_entry ("Boundary id of the solid-pore interface",
                                 "1",
                                 Patterns::Integer(),
                                 "ID used to read mesh");

                                 
        }
        param.leave_subsection();
    }
    param.leave_subsection();  
    
    FuelCellShop::Material::PureSolid::declare_PureSolid_parameters(param);
   
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::SolidLayer<dim>::initialize (ParameterHandler &param)
{
    
    NAME::BaseLayer<dim>::initialize(param);
       
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(this->name);
        {
            concrete_solid_name = param.get("Solid material type");
            boundary_ids.push_back(param.get_integer("Boundary id of the solid-ionomer interface"));
            boundary_ids.push_back(param.get_integer("Boundary id of the solid-pore interface"));
            
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    solid = FuelCellShop::Material::PureSolid::create_PureSolid(concrete_solid_name, param);
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::SolidLayer<dim>::effective_electron_conductivity(Tensor<2,dim>& sigma_eff ) const {
    
    for (int j = 0; j < dim; j++)            
        for (int k = 0; k < dim; k++)                
                sigma_eff[j][k] = solid->get_electrical_conductivity(solution_iterator->second.get_data()->at(0));
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::SolidLayer<dim>::effective_electron_conductivity(std::vector<Tensor<2,dim>>& sigma_eff ) const {
    
    for (int i = 0; i <= solution_iterator->second.get_data()->size() ; i++ )        
        for (int j = 0; j < dim; j++)            
            for (int k = 0; k < dim; k++)                
                sigma_eff[i][j][k] = solid->get_electrical_conductivity(solution_iterator->second.get_data()->at(i));            
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::SolidLayer<dim>::derivative_effective_electron_conductivity(std::map<VariableNames, Tensor<2,dim>>& Dsigma_eff ) const {
    
    typename std::map<VariableNames, Tensor<2,dim>>::iterator iterator_tensor = Dsigma_eff.find(temperature_enum);
    for (int j = 0; j < dim; j++)
        for (int k = 0; k < dim; k++)
                iterator_tensor->second[j][k] = solid->get_Delectrical_conductivity_Dtemperature(solution_iterator->second.get_data()->at(0)); 
}
//---------------------------------------------------------------------------
template <int dim>
void 
NAME::SolidLayer<dim>::derivative_effective_electron_conductivity(std::vector< std::map< VariableNames, Tensor<2,dim> > >& Dsigma_eff ) const
{
    for (int i = 0; i <= solution_iterator->second.get_data()->size() ; i++ ){
        typename std::map<VariableNames, Tensor<2,dim>>::iterator iterator_tensor = Dsigma_eff.at(i).find(temperature_enum);
        for (int j = 0; j <=dim; j++)
            for (int k = 0; k <=dim; k++)
                
                iterator_tensor->second[j][k] = solid->get_Delectrical_conductivity_Dtemperature(solution_iterator->second.get_data()->at(i));
    }
}
//---------------------------------------------------------------------------
template <int dim>
void 
NAME::SolidLayer<dim>::effective_thermal_conductivity(Tensor<2,dim>& Omega_eff ) const
{
    for (int j = 0; j < dim; j++)
        for (int k = 0; k < dim; k++)
            Omega_eff[j][k] = solid->get_thermal_conductivity(solution_iterator->second.get_data()->at(0));
}
//---------------------------------------------------------------------------
template <int dim>
void 
NAME::SolidLayer<dim>::effective_thermal_conductivity(std::vector<Tensor<2,dim>>& Omega_eff ) const
{
    for (int i = 0; i <= solution_iterator->second.get_data()->size() ; i++ )        
        for (int j = 0; j < dim; j++)            
            for (int k = 0; k < dim; k++)                
                Omega_eff[i][j][k] = solid->get_thermal_conductivity(solution_iterator->second.get_data()->at(i));
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::SolidLayer<dim>::derivative_effective_thermal_conductivity(std::map<VariableNames, Tensor<2,dim> >& DOmega_eff ) const {
    
    typename std::map<VariableNames, Tensor<2,dim>>::iterator iterator_tensor = DOmega_eff.find(temperature_enum);
    
    for (int j = 0; j < dim; j++)
        for (int k = 0; k < dim; k++)
            
            iterator_tensor->second[j][k] = solid->get_Dthermal_conductivity_Dtemperature(solution_iterator->second.get_data()->at(0));
}
//---------------------------------------------------------------------------
template <int dim>
void 
NAME::SolidLayer<dim>::derivative_effective_thermal_conductivity(std::vector< std::map<VariableNames, Tensor<2,dim> > >& DOmega_eff ) const {
    
    for (int i = 0; i <= solution_iterator->second.get_data()->size() ; i++ ){
        typename std::map<VariableNames, Tensor<2,dim>>::iterator iterator_tensor = DOmega_eff.at(i).find(temperature_enum);
        
        for (int j = 0; j <=dim; j++)
            for (int k = 0; k <=dim; k++)
                iterator_tensor->second[j][k] = solid->get_Dthermal_conductivity_Dtemperature(solution_iterator->second.get_data()->at(i));
    }
}
// //---------------------------------------------------------------------------
// template <int dim>
// void 
// NAME::SolidLayer<dim>::test_class()
// {
//     NAME:SolidLayer<dim> tina("name of subsection in parameter file");
//     FuelCellShop::Material::CarbonFiber cf;
//     tina.set_solid_and_compute(&cf,298);
//     double val;
//     tina.effective_electron_conductivity(val);
//     FcstUtilities::log<<"The computed conductivity is: "<<val<<std::endl;
// } 
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::SolidLayer<deal_II_dimension>;
