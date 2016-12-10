// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: experimental_solid.cc
// - Description: This class describes a solid
// - Developers: Valentin N. Zingan, University of Alberta
// - Id: $Id: experimental_solid.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <materials/experimental_solid.h>

namespace NAME = FuelCellShop::Material;

       //////////////////////////////////////////////////
       //////////////////////////////////////////////////
       // CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
       //////////////////////////////////////////////////
       //////////////////////////////////////////////////

// ---             ---
// --- Constructor ---
// ---             ---

NAME::ExperimentalSolid::ExperimentalSolid(const std::string& name)
:
NAME::BaseMaterial(name)
{ }

// ---            ---
// --- Destructor ---
// ---            ---

NAME::ExperimentalSolid::~ExperimentalSolid()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

void
NAME::ExperimentalSolid::declare_parameters(ParameterHandler& param) const
{
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Materials");
        {
            param.enter_subsection(this->name);
            {

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
NAME::ExperimentalSolid::initialize(ParameterHandler& param)
{
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Materials");
        {
            param.enter_subsection(this->name);
            {

            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}