//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: PSD_base.h
//    - Description: Base class for pore size distribution model.
//    - Developers: 2009-13 by Marc Secanell, University of Alberta
//                  2013-14 by Jie Zhou, University of Alberta
//    - $ $
//
//---------------------------------------------------------------------------

#include <microscale/PSD_none.h>

namespace NAME = FuelCellShop::MicroScale;

template<int dim>
const std::string NAME::NonePSD<dim>::concrete_name("NonePSD");

template<int dim>
NAME::NonePSD<dim> const* NAME::NonePSD<dim>::PROTOTYPE = new NAME::NonePSD<dim>();

//---------------------------------------------------------------------------
template<int dim>
NAME::NonePSD<dim>::NonePSD() :
NAME::BasePSD<dim>() 
{
    this->get_mapFactory()->insert(
            std::pair<std::string, FuelCellShop::MicroScale::BasePSD<dim>*>(
                    concrete_name, this));
}

//---------------------------------------------------------------------------
template<int dim>
NAME::NonePSD<dim>::NonePSD(std::string name) :
NAME::BasePSD<dim>(name) 
{
}

template class NAME::NonePSD<deal_II_dimension>;