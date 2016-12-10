// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: initial_and_boundary_data.cc
// - Description: This namespace contains data and methods
//                that handle initial and boundary data of
//                a problem at hand
// - Developers: Valentin N. Zingan, University of Alberta
// - Id: $Id: initial_and_boundary_data.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <application_core/initial_and_boundary_data.h>

namespace NAME = FuelCell::InitialAndBoundaryData;

       //////////////////////////////////////////////////
       //////////////////////////////////////////////////
       // CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
       //////////////////////////////////////////////////
       //////////////////////////////////////////////////

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::InitialOrBoundaryDataBase<dim>::InitialOrBoundaryDataBase(boost::shared_ptr<ApplicationData> data,
                                                                const unsigned int                 n_components)
:
Function<dim>(n_components),
data(data)
{ }

// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::InitialOrBoundaryDataBase<dim>::~InitialOrBoundaryDataBase()
{ }

       ///////////////////////
       ///////////////////////
       // SERVICE FUNCTIONS //
       ///////////////////////
       ///////////////////////

// ---       ---
// --- value ---
// ---       ---

template<int dim>
double
NAME::InitialOrBoundaryDataBase<dim>::value(const Point<dim>&  point,
                                            const unsigned int no_component) const
{
  AssertThrow( no_component < this->n_components , ExcIndexRange( no_component , 0 , this->n_components ) );

  return this->math_expression(point,
                               no_component);
}

// ---                 ---
// --- math_expression ---
// ---                 ---

template<int dim>
double
NAME::InitialOrBoundaryDataBase<dim>::math_expression(const Point<dim>&  point,
                                                      const unsigned int no_component) const
{
  print_caller_name(__FUNCTION__);
}

       /////////////////////
       /////////////////////
       // MINOR FUNCTIONS //
       /////////////////////
       /////////////////////

// ---                   ---
// --- print_caller_name ---
// ---                   ---

template<int dim>
void
NAME::InitialOrBoundaryDataBase<dim>::print_caller_name(const std::string& caller_name) const
{
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << caller_name << " called in Class " << info.name() << std::endl;
}

       /////////////////////////////
       /////////////////////////////
       // EXPLICIT INSTANTIATIONS //
       /////////////////////////////
       /////////////////////////////

// ---                           ---
// --- InitialOrBoundaryDataBase ---
// ---                           ---

template class NAME::InitialOrBoundaryDataBase<deal_II_dimension>;