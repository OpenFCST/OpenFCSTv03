// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: channel.cc
// - Description: This class describes a channel
// - Developers: Valentin N. Zingan, University of Alberta
// - Id: $Id$
//
// ----------------------------------------------------------------------------

#include <layers/channel.h>

namespace NAME = FuelCellShop::Layer;

//////////////////////////////////////////////////
//////////////////////////////////////////////////
// CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
//////////////////////////////////////////////////
//////////////////////////////////////////////////

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::Channel<dim>::Channel(const std::string& name)
:
NAME::BaseLayer<dim>(name),
fluid(nullptr),
gas_mixture(nullptr)
{ }

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::Channel<dim>::Channel(const std::string&                         name,
                            FuelCellShop::Material::ExperimentalFluid& fluid)
:
NAME::BaseLayer<dim>(name),
fluid(&fluid),
gas_mixture(nullptr)
{ }

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::Channel<dim>::Channel(const std::string&                  name,
                            FuelCellShop::Material::GasMixture& gas_mixture)
:
NAME::BaseLayer<dim>(name),
fluid(nullptr),
gas_mixture(&gas_mixture)
{ }

// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::Channel<dim>::~Channel()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::Channel<dim>::declare_parameters(ParameterHandler& param) const
{
  NAME::BaseLayer<dim>::declare_parameters(this->name,
                                           param);

  param.enter_subsection("Fuel cell data");
  {
    param.enter_subsection(this->name);
    {
           param.declare_entry("Roughness [cm]",
                               "1.0e-4",
                                Patterns::Double(),
                               " ");

           param.declare_entry("Effective electronic conductivity [S/cm]",
                               "100.0",
                                Patterns::Double(),
                               " ");
    }
    param.leave_subsection();
  }
  param.leave_subsection();
}

// ---            ---
// --- initialize ---
// ---            ---

template<int dim>
void
NAME::Channel<dim>::initialize(ParameterHandler& param)
{
  NAME::BaseLayer<dim>::initialize(param);

  param.enter_subsection("Fuel cell data");
  {
    param.enter_subsection(this->name);
    {
           roughness = param.get_double("Roughness [cm]");

           for(unsigned int i = 0; i < dim; ++i)
             for(unsigned int j = 0; j < dim; ++j)
               if( i == j )
                 effective_electronic_conductivity[i][j] = param.get_double("Effective electronic conductivity [S/cm]");
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

// ---                        ---
// --- print_layer_properties ---
// ---                        ---

template<int dim>
void
NAME::Channel<dim>::print_layer_properties() const
{
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;
  FcstUtilities::log << "------------------------------";
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;
  FcstUtilities::log << "Parameters for "  << this->name << ":";
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;
  for(unsigned int i = 0; i < this->material_ids.size(); ++i)
    FcstUtilities::log << "Material ids: " << this->material_ids.at(i) << std::endl;
  FcstUtilities::log << std::endl;
  FcstUtilities::log << "Roughness [cm]: "  << roughness;
  FcstUtilities::log << std::endl;
  FcstUtilities::log << "Effective electronic conductivity [S/cm]: " << effective_electronic_conductivity;
  FcstUtilities::log << std::endl;

  if( fluid )
    fluid->print_material_properties();

  if( gas_mixture )
    gas_mixture->print_material_properties();
}

/////////////////////////////
/////////////////////////////
// EXPLICIT INSTANTIATIONS //
/////////////////////////////
/////////////////////////////

// ---         ---
// --- Channel ---
// ---         ---

template class NAME::Channel<deal_II_dimension>;