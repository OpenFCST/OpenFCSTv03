//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: membrane_layer.h
//    - Description: (base) class for membrane layers
//    - Developers: Madhur Bhaiya (2012-13)
//    - Id: $Id: membrane_layer.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <layers/membrane_layer.h>

namespace NAME = FuelCellShop::Layer;

//---------------------------------------------------------------------------
template <int dim>
NAME::MembraneLayer<dim>::MembraneLayer()
  : NAME::BaseLayer<dim> ()
{}

//---------------------------------------------------------------------------
template <int dim>
NAME::MembraneLayer<dim>::MembraneLayer(std::string name)
  : NAME::BaseLayer<dim>(name)
{
    electrolyte = boost::shared_ptr< FuelCellShop::Material::PolymerElectrolyteBase > ();
}

//---------------------------------------------------------------------------
template <int dim>
NAME::MembraneLayer<dim>::~MembraneLayer()
{}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::declare_parameters(const std::string& pem_section_name, 
                                             ParameterHandler &param) const
{
    
    NAME::BaseLayer<dim>::declare_parameters(pem_section_name,param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(pem_section_name);
        {
            param.declare_entry("Membrane layer type",
                                "NafionMembrane",
                                Patterns::Selection("NafionMembrane"),
                                " ");
            param.declare_entry("Electrolyte type",
                                "Nafion",
                                Patterns::Selection("Nafion"),
                                " ");
            // Declare electrolyte material parameters:
            FuelCellShop::Material::PolymerElectrolyteBase::declare_PolymerElectrolyte_parameters(param);
        }     
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::initialize (ParameterHandler &param)
{
    FuelCellShop::Layer::BaseLayer<dim>::initialize(param);
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection(this->name); 
        { 
            electrolyte_type = param.get("Electrolyte type");
            electrolyte = FuelCellShop::Material::PolymerElectrolyteBase::create_PolymerElectrolyte(param, electrolyte_type);
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::effective_proton_conductivity(double&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
	    << " called in Class "
	    << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::effective_proton_conductivity(std::vector<double>&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
	    << " called in Class "
	    << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::derivative_effective_proton_conductivity(std::map< VariableNames, std::vector<double> >&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
	    << " called in Class "
	    << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::effective_water_diffusivity(double&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
	    << " called in Class "
	    << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::effective_water_diffusivity(std::vector<double>&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
	    << " called in Class "
	    << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::derivative_effective_water_diffusivity(std::map< VariableNames, std::vector<double> >&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
	    << " called in Class "
	    << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::effective_oxygen_diffusivity(double&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
	    << " called in Class "
	    << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::effective_oxygen_diffusivity(std::vector<double>&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
	    << " called in Class "
	    << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::derivative_effective_oxygen_diffusivity(std::map< VariableNames, std::vector<double> >&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
	    << " called in Class "
	    << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::effective_thermal_conductivity(double&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
	    << " called in Class "
	    << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::effective_thermal_conductivity(std::vector<double>&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
	    << " called in Class "
	    << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::effective_thermal_conductivity(std::vector< Tensor<2,dim> >&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
	    << " called in Class "
	    << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::derivative_effective_thermal_conductivity(std::vector< std::vector<double> >&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
	    << " called in Class "
	    << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::derivative_effective_thermal_conductivity(std::vector< std::vector< Tensor<2,dim> > >&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
	    << " called in Class "
	    << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::effective_thermoosmotic_diffusivity(std::vector<double>&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
            << " called in Class "
            << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::MembraneLayer<dim>::derivative_effective_thermoosmotic_diffusivity(std::map< VariableNames, std::vector<double> >&) const
{
  // Print info:
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << __FUNCTION__
            << " called in Class "
            << info.name() << std::endl;
  
  // Throw an exception:
  ExcNotImplemented ();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations. 
template class NAME::MembraneLayer<deal_II_dimension>;
