//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: micro_porous_layer.h
//    - Description: Base Micro-Porous Layer Class. It implements the interface for other micro-porous layer classes
//        and some common methods.
//    - Developers: M. Secanell
//
//---------------------------------------------------------------------------

#include <layers/micro_porous_layer.h>

namespace NAME = FuelCellShop::Layer;

template <int dim>
const std::string NAME::MicroPorousLayer<dim>::concrete_name ("MicroPorousLayer");

//---------------------------------------------------------------------------
template <int dim>
NAME::MicroPorousLayer<dim>::MicroPorousLayer()
  : NAME::PorousLayer<dim>()
{}

//---------------------------------------------------------------------------
template <int dim>
NAME::MicroPorousLayer<dim>::MicroPorousLayer(const std::string& name)
  : NAME::PorousLayer<dim>(name)
{}

//---------------------------------------------------------------------------
template <int dim>
NAME::MicroPorousLayer<dim>::~MicroPorousLayer()
{}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::MicroPorousLayer<dim>::declare_parameters (const std::string& name, 
                                                 ParameterHandler &param) const
{
    FuelCellShop::Layer::PorousLayer<dim>::declare_parameters(name,param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(name);
        {
            param.declare_entry("Micro porous layer type",
                                "SGL24BC",
                                Patterns::Selection("SGL24BC | DesignMPL"),
                                "The type of the Micro Porous Layer ");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}   

//---------------------------------------------------------------------------
template <int dim>
void
NAME::MicroPorousLayer<dim>::initialize (ParameterHandler &param)
{
  NAME::PorousLayer<dim>::initialize(param);
  
  param.enter_subsection("Fuel cell data"); 
  {
    param.enter_subsection(this->name); 
    { 

    }
    param.leave_subsection();
  }
  param.leave_subsection();
}

//---------------------------------------------------------------------------
//-------------------- OPTIMIZATION ROUTINES --------------------------------
//---------------------------------------------------------------------------
/*
template <int dim>
Tensor<2,dim>
NAME::MicroPorousLayer<dim>::Deffective_transport_property_pores_Dporosity(const double prop) const
{

}
*/

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations. 
template class NAME::MicroPorousLayer<deal_II_dimension>;
