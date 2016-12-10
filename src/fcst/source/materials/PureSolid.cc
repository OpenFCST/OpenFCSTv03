//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: PureSolid.cc
//    - Description: Material class for pure solids.
//    - Developers: Marc Secanell and Jie Zhou
//    - $Id: PureSolid.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <materials/PureSolid.h>

namespace NAME = FuelCellShop::Material;

//---------------------------------------------------------------------------
const std::string NAME::Graphite::concrete_name = "Graphite";
NAME::Graphite const* NAME::Graphite::PROTOTYPE = new NAME::Graphite();

//---------------------------------------------------------------------------
const std::string NAME::DummySolid::concrete_name = "DummySolid";
NAME::DummySolid const* NAME::DummySolid::PROTOTYPE = new NAME::DummySolid();


// ---             ---
// --- Constructor ---
// ---             ---

NAME::PureSolid::PureSolid(const std::string& name)
 :
NAME::BaseMaterial(name)
{}
 
// ---            ---
// --- Destructor ---
// ---            ---
 
NAME::PureSolid::~PureSolid()
{ }

//---------------------------------------------------------------------------
// ---             ---
// --- Constructor ---
// ---             ---
NAME::Graphite::Graphite()
:
NAME::PureSolid()
{
    
}

//---------------------------------------------------------------------------
NAME::Graphite::Graphite(const std::string& name)
:
NAME::PureSolid(name)
{
    this->get_mapFactory()->insert(std::pair<std::string, PureSolid* >(concrete_name, this));
}

// ---            ---
// --- Destructor ---
// ---            ---

NAME::Graphite::~Graphite()
{ }

//---------------------------------------------------------------------------
void 
NAME::Graphite::declare_parameters (ParameterHandler &param) const {
    
    param.enter_subsection(concrete_name);
    {
        /*
         * all the values are cited from 
         * http://www.poco.com/MaterialsandServices/Graphite/SemiconductorGrades/HPD.aspx   
         */   
        param.declare_entry("Density [g/cm^3]",
                            "1.77",
                            Patterns::Double(),
                            "This is the density of the graphite material");
        
        param.declare_entry("electrical_conductivity [S/M]",
                            "5.0e5",
                            Patterns::Double(),
                            "This is the electrical conductivity of the graphite material");
        
        param.declare_entry("thermal_conductivity [watts/m K]",
                            "85.0",
                            Patterns::Double(),
                            "This is the thermal conductivity of the graphite material");
        
        param.declare_entry("coefficient_thermal_expansion [microns/m °C]",
                            "8.1",
                            Patterns::Double(),
                            "This is the coefficient_thermal_expansion of the graphite material");
        
        param.declare_entry("compressive_strength [N/mm2]",
                            "150.0",
                            Patterns::Double(),
                            "This is the compressive_strength of the graphite material");
        
        param.declare_entry("H2_permeability [cm3*cm-2*s-1]",
                            "2.0e-6",
                            Patterns::Double(),
                            "This is the H2_permeability of the graphite material");
        
        param.declare_entry("Poissons_ratio",
                            "0.17",
                            Patterns::Double(),
                            "This is the Poissons_ratio of the graphite material");
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
void 
NAME::Graphite::initialize (ParameterHandler &param){
    
    param.enter_subsection(concrete_name);
    {
        this->density = param.get_double("Density [g/cm^3]");
        this->electrical_conductivity = param.get_double("electrical_conductivity [S/M]");
        this->thermal_conductivity = param.get_double("thermal_conductivity [watts/m K]");
        this->coefficient_thermal_expansion = param.get_double("coefficient_thermal_expansion [microns/m °C]");
        this->compressive_strength = param.get_double("compressive_strength [N/mm2]");
        this->H2_permeability = param.get_double("H2_permeability [cm3*cm-2*s-1]");
        this->Poissons_ratio = param.get_double("Poissons_ratio");
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
double
NAME::Graphite::get_density() const
{
    return density;
}

//---------------------------------------------------------------------------
double
NAME::Graphite::get_electrical_conductivity(double temperature) const
{
    return electrical_conductivity;
}

//---------------------------------------------------------------------------
void
NAME::Graphite::get_electrical_conductivity(std::vector<double> temperature,std::vector<double>& dst) const
{
    for (int i = 0; i <= temperature.size();i++)
    dst[i] = electrical_conductivity;
}
//---------------------------------------------------------------------------
double
NAME::Graphite::get_thermal_conductivity(double temperature) const
{
    return thermal_conductivity;
}

//---------------------------------------------------------------------------
void
NAME::Graphite::get_thermal_conductivity(std::vector<double> temperature,std::vector<double>& dst) const
{
    for (int i = 0; i <= temperature.size();i++)
    dst[i] = thermal_conductivity;
}

//---------------------------------------------------------------------------
double
NAME::Graphite::get_coefficient_thermal_expansion(double) const
{
    return coefficient_thermal_expansion;
}

//---------------------------------------------------------------------------
double
NAME::Graphite::get_compressive_strength(double) const
{
    return compressive_strength;
}

//---------------------------------------------------------------------------
double
NAME::Graphite::get_H2_permeability(double) const
{
    return H2_permeability;
}

//---------------------------------------------------------------------------
double
NAME::Graphite::get_Poissons_ratio(double) const
{
    return Poissons_ratio;
}

//---------------------------------------------------------------------------
// ---             ---
// --- Constructor ---
// ---             ---
NAME::DummySolid::DummySolid()
:
NAME::PureSolid()
{}

NAME::DummySolid::DummySolid(const std::string& name)
:
NAME::PureSolid(name)
{
    this->get_mapFactory()->insert(std::pair<std::string, NAME::PureSolid* >(concrete_name, this));
}
// ---            ---
// --- Destructor ---
// ---            ---

NAME::DummySolid::~DummySolid()
{ }

//---------------------------------------------------------------------------
void 
NAME::DummySolid::declare_parameters (ParameterHandler &param) const {
    
    param.enter_subsection(concrete_name);
    {
        param.declare_entry("Density [g/cm^3]",
                            "1.70",
                            Patterns::Double(),
                            "This is the density of the dummy material");
        
        param.declare_entry("coefficient_thermal_expansion [microns/m °C]",
                            "8.1",
                            Patterns::Double(),
                            "This is the coefficient_thermal_expansion of the graphite material");
        
        param.declare_entry("compressive_strength [N/mm2]",
                            "150.0",
                            Patterns::Double(),
                            "This is the compressive_strength of the graphite material");
        
        param.declare_entry("H2_permeability [cm3*cm-2*s-1]",
                            "2.0e-6",
                            Patterns::Double(),
                            "This is the H2_permeability of the graphite material");
        
        param.declare_entry("Poissons_ratio",
                            "0.17",
                            Patterns::Double(),
                            "This is the Poissons_ratio of the graphite material");
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
void 
NAME::DummySolid::initialize (ParameterHandler &param){
    
    param.enter_subsection(concrete_name);
    {
        this->density = param.get_double("Density [g/cm^3]");
        this->coefficient_thermal_expansion = param.get_double("coefficient_thermal_expansion [microns/m °C]");
        this->compressive_strength = param.get_double("compressive_strength [N/mm2]");
        this->H2_permeability = param.get_double("H2_permeability [cm3*cm-2*s-1]");
        this->Poissons_ratio = param.get_double("Poissons_ratio");
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
double
NAME::DummySolid::get_density() const
{
    return density;
}
//---------------------------------------------------------------------------
double
NAME::DummySolid::get_electrical_conductivity(double temperature) const
{
    double electrical_conductivity = temperature * 2.0 / 3.0;
    return electrical_conductivity;
}
//---------------------------------------------------------------------------
void
NAME::DummySolid::get_electrical_conductivity(std::vector<double> temperature,std::vector<double>& dst) const
{
    for (int i=0;i<=temperature.size();i++)
        dst [i] = temperature [i] * 2.0/3.0; 
}
//---------------------------------------------------------------------------
double
NAME::DummySolid::get_Delectrical_conductivity_Dtemperature(double temperature) const
{
    double Delectrical_conductivity_Dtemperature = 2.0/3.0;
    return Delectrical_conductivity_Dtemperature;
}
//---------------------------------------------------------------------------
void
NAME::DummySolid::get_Delectrical_conductivity_Dtemperature(std::vector<double> temperature,std::vector<double>& dst) const
{
    for (int i=0;i<=temperature.size();i++)
        dst [i] = 2.0/3.0; 
}
//---------------------------------------------------------------------------
double
NAME::DummySolid::get_thermal_conductivity(double temperature) const
{
    double thermal_conductivity = temperature*temperature * 2.0 / 3.0;
    return thermal_conductivity;
}
//---------------------------------------------------------------------------

void
NAME::DummySolid::get_thermal_conductivity(std::vector<double> temperature,std::vector<double>& dst) const
{
    for (int i=0;i<=temperature.size();i++)
        dst [i] = temperature [i] * temperature [i] * 2.0/3.0; 
}
//---------------------------------------------------------------------------

double
NAME::DummySolid::get_Dthermal_conductivity_Dtemperature(double temperature) const
{
    double Dthermal_conductivity_Dtemperature = temperature * 2.0/3.0;
    return Dthermal_conductivity_Dtemperature;
}
//---------------------------------------------------------------------------
void
NAME::DummySolid::get_Dthermal_conductivity_Dtemperature(std::vector<double> temperature,std::vector<double>& dst) const
{
    for (int i=0;i<=temperature.size();i++)
        dst [i] = temperature [i] * 2.0/3.0; 
}

//---------------------------------------------------------------------------
double
NAME::DummySolid::get_coefficient_thermal_expansion(double) const
{
    return coefficient_thermal_expansion;
}

//---------------------------------------------------------------------------
double
NAME::DummySolid::get_compressive_strength(double) const
{
    return compressive_strength;
}

//---------------------------------------------------------------------------
double
NAME::DummySolid::get_H2_permeability(double) const
{
    return H2_permeability;
}

//---------------------------------------------------------------------------
double
NAME::DummySolid::get_Poissons_ratio(double) const
{
    return Poissons_ratio;
}

