// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: experimental_fluid.cc
// - Description: This class describes a fluid
// - Developers: Valentin N. Zingan, University of Alberta
// - Id: $Id: experimental_fluid.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <materials/experimental_fluid.h>

namespace NAME = FuelCellShop::Material;

       //////////////////////////////////////////////////
       //////////////////////////////////////////////////
       // CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
       //////////////////////////////////////////////////
       //////////////////////////////////////////////////

// ---             ---
// --- Constructor ---
// ---             ---

NAME::ExperimentalFluid::ExperimentalFluid(const std::string& name)
:
NAME::BaseMaterial(name)
{ }

// ---            ---
// --- Destructor ---
// ---            ---

NAME::ExperimentalFluid::~ExperimentalFluid()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

void
NAME::ExperimentalFluid::declare_parameters(ParameterHandler& param) const
{
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Materials");
        {
            param.enter_subsection(this->name);
            {
                param.declare_entry("Density",
                                    "1.0",
                                    Patterns::Double(),
                                    " ");

                param.declare_entry("Dynamic viscosity",
                                    "1.0e-2",
                                    Patterns::Double(),
                                    " ");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

// ---            ---
// --- initialize ---
// ---            ---

void
NAME::ExperimentalFluid::initialize(ParameterHandler& param)
{
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Materials");
        {
            param.enter_subsection(this->name);
            {
                density           = param.get_double("Density");
                dynamic_viscosity = param.get_double("Dynamic viscosity");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

       ////////////////////////
       ////////////////////////
       // ACCESSORS AND INFO //
       ////////////////////////
       ////////////////////////

// ---                           ---
// --- print_material_properties ---
// ---                           ---

void
NAME::ExperimentalFluid::print_material_properties() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "------------------------------";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "Parameters for " << this->name << ":";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "Density [kg/m^3]:           " << density;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "Dynamic viscosity [Pa sec]: " << dynamic_viscosity;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "------------------------------";
    FcstUtilities::log << std::endl;
}