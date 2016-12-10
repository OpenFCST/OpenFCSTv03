//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: PureLiquid.cc
//    - Description: Material class for pure liquids.
//    - Developers: Madhur Bhaiya and Phil Wardlaw
//    - Id: $Id: PureLiquid.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <materials/PureLiquid.h>

namespace NAME = FuelCellShop::Material;

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//LIQUID WATER-------------------------------------------------------------------
//---------------------------------------------------------------------------
NAME::LiquidWater::LiquidWater()
{
  M   = 18.015; // g/mol
  molecular_width = 0.28e-9;
  HenryO2 = 5.146e+5*exp(-500/353); //TODO: REPLACE TEMP Hard COde
}
//---------------------------------------------------------------------------
NAME::LiquidWater::~LiquidWater()
{
}
//---------------------------------------------------------------------------
char* 
NAME::LiquidWater::get_name()
{
  return "Liquid Water";
}
//---------------------------------------------------------------------------
char* 
NAME::LiquidWater::get_formula()
{
  return "H2O (liquid)";
}

//---------------------------------------------------------------------------
double
NAME::LiquidWater::latentVap_heat(const double& Temp)
{
    double T_c = Temp - 273.15;
    double h_fg = (2500.304 - (2.2521025*T_c) - (0.021465847*(std::pow(T_c,1.5))) + (3.1750136e-4*(std::pow(T_c,2.5))) - (2.8607959e-5*(std::pow(T_c,3.0)))) * 18.015; // J/mol
    return h_fg;
}

//---------------------------------------------------------------------------
double
NAME::LiquidWater::deriv_latentVap_heat(const double& Temp)
{
    double T_c = Temp - 273.15;
    double dh_fg_dT = ((3.1750136e-4*2.5*(std::pow(T_c,1.5))) - 2.2521025 - (0.021465847*1.5*(std::pow(T_c,0.5))) - (2.8607959e-5*3.0*(std::pow(T_c,2.0)))) * 18.015; //J/mol-(K)
    return dh_fg_dT;
}

//---------------------------------------------------------------------------
double
NAME::LiquidWater::viscosity(const double& Temp)
{
    double exponent = 247.8/(Temp - 140.0);
    double mu = 2.414e-4 * std::pow(10.0, exponent); // g/(cm-s)
    return mu;
}
//---------------------------------------------------------------------------
double
NAME::LiquidWater::surface_tension(const double& Temp)
{
    double B = 235.8 * 1.0e-3;  // N/m
    double mu = 1.256;
    double b_small = -0.625;
    double T_c = (647.15 - Temp) / 647.15;
    double theta =  B * std::pow(T_c,mu) * (1.0 + b_small * T_c);   // N/m
    return theta;
}
//---------------------------------------------------------------------------
double
NAME::LiquidWater::deriv_viscosity(const double& Temp)
{
    double exponent = 247.8/(Temp - 140.0);
    double dmu_dT = 2.414e-4 * std::pow(10.0, exponent) * std::log(10.0) * ( ((-1.0)*247.8)/((Temp-140.0)*(Temp-140.0)) );
    return dmu_dT;
}

//---------------------------------------------------------------------------
void
NAME::LiquidWater::declare_parameters (ParameterHandler &param) const
{


    param.enter_subsection("Materials");
    {
        param.enter_subsection("Water");
        {
            param.declare_entry("Oxygen diffusion coefficient [cm^2/s]",
                    "9.19e-5",
                    Patterns::Double());
            param.declare_entry("Proton diffusion coefficient [cm^2/s]",
                    "9.2e-5",
                    Patterns::Double());
            param.declare_entry("Relative permittivity",
                    "60",
                    Patterns::Double());
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}


//---------------------------------------------------------------------------
void
NAME::LiquidWater::initialize(ParameterHandler &param)
{

    param.enter_subsection("Materials");
    {
        param.enter_subsection("Water");
        {
            oxygen_diffusion_coeff = param.get_double("Oxygen diffusion coefficient [cm^2/s]");
            proton_diffusion_coeff = param.get_double("Proton diffusion coefficient [cm^2/s]");
            relative_permittivity = param.get_double("Relative permittivity");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}


//---------------------------------------------------------------------------
