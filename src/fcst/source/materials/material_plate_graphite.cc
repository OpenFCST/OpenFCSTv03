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

// Include FCST classes
#include <materials/material_plate_graphite.h>

namespace NAME = FuelCellShop::Material;

//---------------------------------------------------------------------------
NAME::MaterialPlateGraphite::MaterialPlateGraphite(std::string name)
	: FuelCellShop::Material::MaterialPlateBase(name)
{
  //implement routine  
}

//---------------------------------------------------------------------------
NAME::MaterialPlateGraphite::~MaterialPlateGraphite()
{}

//---------------------------------------------------------------------------
void 
NAME::MaterialPlateGraphite::declare_parameters (ParameterHandler &param) const
{
  FuelCellShop::Material::MaterialPlateBase::declare_parameters(param);
  
  param.enter_subsection("Fuel cell data"); 
  {
    param.enter_subsection("Materials");
    {
      param.enter_subsection(name);
      {
	param.declare_entry ("Electron conductivity (S/cm)",
			     "1e200",
			     Patterns::Double(),
			     "Electrical conductivity of material. In S/cm");
	param.declare_entry ("Thermal conductivity (W/(cm K))",
			     "1.0", 
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
NAME::MaterialPlateGraphite::initialize (ParameterHandler &param)
{
  FuelCellShop::Material::BaseMaterial::initialize(param);

  param.enter_subsection("Fuel cell data"); 
  {
    param.enter_subsection("Materials");
    {
      param.enter_subsection(name);
      {
	this->electron_conductivity = param.get_double("Electron conductivity (S/cm)");
	this->thermal_conductivity = param.get_double("Thermal conductivity (W/(cm K))");
      }
      param.leave_subsection();
    }
    param.leave_subsection();
  }
  param.leave_subsection();
}

//---------------------------------------------------------------------------
double
NAME::MaterialPlateGraphite::get_electron_conductivity() const
{
  return this->electron_conductivity;
}
	
//---------------------------------------------------------------------------
void
NAME::MaterialPlateGraphite::get_electron_conductivity_derivative(double &E, std::vector<double>& dE) const
{
  E = this->electron_conductivity;
  for (unsigned int i=0; i<derivative_flags.size(); i++)
  {
    dE[i] = 0.0;
  }
}


/*
//---------------------------------------------------------------------------
double 
NAME::MaterialPlateGraphite::get_thermal_conductivity() const
{
};

//---------------------------------------------------------------------------
void 
NAME::MaterialPlateGraphite::get_thermal_conductivity_derivative(double &, std::vector<double>&) const
{
};

//---------------------------------------------------------------------------
double 
NAME::MaterialPlateGraphite::get_youngs_modulus() const
{

  
};

//---------------------------------------------------------------------------
void 
NAME::MaterialPlateGraphite::get_youngs_modulus_derivative(double &, std::vector<double>&) const
{
};

//---------------------------------------------------------------------------
double 
NAME::MaterialPlateGraphite::get_poissons_ratio() const
{
};

//---------------------------------------------------------------------------
void 
NAME::MaterialPlateGraphite::get_poissons_modulus_derivative(double &, std::vector<double>&) const
{
};

//---------------------------------------------------------------------------
double 
NAME::MaterialPlateGraphite::get_expansion_coefficient() const
{
};


//---------------------------------------------------------------------------
void 
NAME::MaterialPlateGraphite::get_expansion_coefficient_derivative(double &E, std::vector<double>& dE) const
{
};
*/