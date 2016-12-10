//---------------------------------------------------------------------------
//
//    Copyright (C) 2011 by Marc Secanell, University of Alberta
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <materials/carbon_fiber.h>

namespace NAME=FuelCellShop::Material;

//---------------------------------------------------------------------------
NAME::CarbonFiber::CarbonFiber(std::string name)
:
FiberBase (name)
{
	this->electrical_conductivity=88.84;
	this->thermal_conductivity=1.0;
	this->density=2.0;
}

//---------------------------------------------------------------------------
NAME::CarbonFiber::~CarbonFiber()
{}

//---------------------------------------------------------------------------
void
NAME::CarbonFiber::declare_parameters(ParameterHandler &param) const
{
  param.enter_subsection("Fuel cell data"); 
  {
    param.enter_subsection("Materials");
    {
      param.enter_subsection(name);
      {
			param.declare_entry ("Electrical conductivity",
						"88.84", 
						Patterns::Double(),
						"Electrical conductivity of a group of particles packed to have 0% porosity. In S/cm");
			param.declare_entry ("Thermal conductivity",
						"1.0", 
						Patterns::Double());
			param.declare_entry ("Density [g/cm^3]",
						"2.0", 
						Patterns::Double());
      }
      param.leave_subsection();
    }
    param.leave_subsection();
  }
  param.leave_subsection();
}

//---------------------------------------------------------------------------
void
NAME::CarbonFiber::initialize (ParameterHandler &param)
{
  param.enter_subsection("Fuel cell data"); 
  {
    param.enter_subsection("Materials");
    {
      param.enter_subsection(name);
      {
		this->electrical_conductivity = param.get_double("Electrical conductivity");
		this->thermal_conductivity = param.get_double("Thermal conductivity");
		this->density = param.get_double("Density [g/cm^3]");
      }
      param.leave_subsection();
    }
    param.leave_subsection();
  }
  param.leave_subsection();
}

double 
NAME::CarbonFiber::get_electrical_conductivity(double temperature) const
{
	double elec_conductivity(0.);
	elec_conductivity = 100.+100000.*(1./353.0-1./temperature);
	return elec_conductivity;
}




