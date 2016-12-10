//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: carbon.cc
//    - Description: Class characterizing a carbon black support
//    - Developers: M. Secanell and Madhur Bhaiya
//    - $Id: carbon.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <materials/carbon.h>

namespace NAME = FuelCellShop::Material;

const std::string NAME::CarbonBlack::concrete_name ("CarbonBlack");

NAME::CarbonBlack const* NAME::CarbonBlack::PROTOTYPE = new NAME::CarbonBlack(true);

//---------------------------------------------------------------------------
NAME::CarbonBlack::CarbonBlack()
:
CatalystSupportBase ()
{
    this->electrical_conductivity = 88.84;
    this->thermal_conductivity = 1.0;
    this->density = 2.0;
}

//---------------------------------------------------------------------------
NAME::CarbonBlack::CarbonBlack(const bool create_replica)
:
CatalystSupportBase ()
{
    if (create_replica)
        this->get_mapFactory()->insert(std::pair<std::string, CatalystSupportBase* > (this->concrete_name, this) );           
}

//---------------------------------------------------------------------------
NAME::CarbonBlack::~CarbonBlack()
{}

//---------------------------------------------------------------------------
void
NAME::CarbonBlack::declare_parameters(ParameterHandler &param) const
{
    param.enter_subsection("Materials");
    {
        param.enter_subsection(NAME::CarbonBlack::concrete_name);
        {
            
            param.declare_entry ("Density [g/cm^3]",
                                 "2.0", 
                                 Patterns::Double());
            param.declare_entry ("Electrical conductivity [S/cm]",
                                 "88.84", 
                                 Patterns::Double(),
                                 "Electrical conductivity of a group of particles packed to have 0% porosity. In S/cm");
            param.declare_entry ("Thermal conductivity [W/(cm-K)]",
                                 "1.0", 
                                 Patterns::Double());
            
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
void
NAME::CarbonBlack::initialize (ParameterHandler &param)
{
    
    param.enter_subsection("Materials");
    {
        param.enter_subsection(NAME::CarbonBlack::concrete_name);
        {
            this->electrical_conductivity = param.get_double("Electrical conductivity [S/cm]");
            this->thermal_conductivity = param.get_double("Thermal conductivity [W/(cm-K)]");
            this->density = param.get_double("Density [g/cm^3]");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

