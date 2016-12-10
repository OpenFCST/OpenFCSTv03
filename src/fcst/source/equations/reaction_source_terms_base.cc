// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: reaction_source_terms_base.cc
// - Description: This class is used to assemble both cell matrix and cell residual
//                for reaction source terms in the catalyst layers for various equation classes
// - Developers: Madhur Bhaiya, Marc Secanell
//
// ----------------------------------------------------------------------------

#include "equations/reaction_source_terms_base.h"

namespace NAME = FuelCellShop::Equation;

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::ReactionSourceTermsBase<dim>::ReactionSourceTermsBase(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data)
:
NAME::EquationBase<dim>(system_management,data)
{

}

// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::ReactionSourceTermsBase<dim>::~ReactionSourceTermsBase()
{
}

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::ReactionSourceTermsBase<dim>::declare_parameters(ParameterHandler& param) const
{
}

// ---                           ---
// ---  EXPLICIT INSTANTIATIONS  ---
// ---                           ---

template class NAME::ReactionSourceTermsBase<deal_II_dimension>;