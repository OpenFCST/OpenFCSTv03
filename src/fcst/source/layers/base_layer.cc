// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: base_layer.cc
// - Description: This is a base class for all available FCST layers
// - Developers: Marc Secanell Gallart,    University of Alberta
//               Madhur Bhaiya, University of Alberta
//               Valentin N. Zingan, U of A
// - $Id: base_layer.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <layers/base_layer.h>

namespace NAME = FuelCellShop::Layer;

//---------------------------------------------------------------------------
template <int dim>
NAME::BaseLayer<dim>::BaseLayer(const std::string& name)
  : name(name)
{
  //FcstUtilities::log<<"Creating layer " << name<<std::endl;
}

//---------------------------------------------------------------------------
template <int dim>
NAME::BaseLayer<dim>::~BaseLayer()
{}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::BaseLayer<dim>::initialize (ParameterHandler &param)
{
  param.enter_subsection("Fuel cell data");
  {
    param.enter_subsection(this->name);
    {
      material_ids = FcstUtilities::string_to_number<unsigned int>( Utilities::split_string_list( param.get("Material id") ) );
           
      this->set_local_material_id(material_ids.at(0));

    }
    param.leave_subsection();
  }
  param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
bool
NAME::BaseLayer<dim>::belongs_to_material(const unsigned int material)
{
    
    //bool result = material_ids.end() != std::find(material_ids.begin(), material_ids.end(), material);
    std::vector<unsigned int>::iterator ind = std::find(material_ids.begin(), material_ids.end(), material);
    bool result = false;
    if (ind != material_ids.end())
        result = true;
    
    
    if (result)
        this->set_local_material_id(material);
    
    return result;
}

//---------------------------------------------------------------------------
//template <int dim>
//void
//NAME::BaseLayer<dim>::set_local_material_id(const unsigned int material)
//{
//    bool result = material_ids.end() != std::find(material_ids.begin(), material_ids.end(), material);
//
//    if (result)
//        this->local_material_id = material;
//    else
//        AssertThrow(false, ExcMessage("set_local_material_id: The local material_ID you want to assign does not exist in material_ids for the layer"));
//}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::BaseLayer<dim>::print_layer_properties() const
{
	  const std::type_info& info = typeid(*this);
	  FcstUtilities::log << "Pure function " << __FUNCTION__
		  << " called in Class "
		  << info.name() << std::endl;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::BaseLayer<deal_II_dimension>;
