//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: nafion_membrane.cc
//    - Description: Class representing Nafion membrane layer class - returning effective transport properties
//    - Developers: Madhur Bhaiya (2012-13)
//    - Id: $Id: nafion_membrane.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include <layers/nafion_membrane.h>

namespace NAME = FuelCellShop::Layer;

template <int dim>
const std::string NAME::NafionMembrane<dim>::concrete_name ("NafionMembrane");

template <int dim>
NAME::NafionMembrane<dim> const* NAME::NafionMembrane<dim>::PROTOTYPE = new NAME::NafionMembrane<dim>();

//---------------------------------------------------------------------------
template <int dim>
NAME::NafionMembrane<dim>::NafionMembrane()
  : NAME::MembraneLayer<dim>()
{
    //FcstUtilities::log<<" Register NafionMembrane to FactoryMap"<<std::endl;
    this->get_mapFactory()->insert(std::pair<std::string, NafionMembrane<dim>* > (this->concrete_name, this) );     
}

//---------------------------------------------------------------------------
template <int dim>
NAME::NafionMembrane<dim>::NafionMembrane(std::string name)
  : NAME::MembraneLayer<dim>(name)
{
    FcstUtilities::log<<" Created a NafionMembrane"<<std::endl;
}

//---------------------------------------------------------------------------
template <int dim>
NAME::NafionMembrane<dim>::~NafionMembrane()
{}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::declare_parameters(const std::string& pem_section_name,
                                              ParameterHandler &param) const
{
    
    NAME::MembraneLayer<dim>::declare_parameters(pem_section_name, param);
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection(pem_section_name); 
        {
            param.enter_subsection(this->concrete_name); 
            {
                param.declare_entry ("Method effective thermal conductivity",
                                     "Given",
                                     Patterns::Selection("Given"),
                                     "Method used to compute effective thermal conductivity");
                param.declare_entry ("Thermal conductivity, [W/(cm K)]",
                                     "0.015",
                                     Patterns::Double());
                
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::initialize (ParameterHandler &param)
{
    FuelCellShop::Layer::MembraneLayer<dim>::initialize(param);
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection(this->name); 
        {
            param.enter_subsection(concrete_name);
            {
                method_thermal_conductivity = param.get("Method effective thermal conductivity");
                thermal_conductivity = param.get_double("Thermal conductivity, [W/(cm K)]");
            }
            param.leave_subsection();
            
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::effective_proton_conductivity(double& sigma_eff) const
{
    this->electrolyte->proton_conductivity(sigma_eff);
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::effective_proton_conductivity(std::vector<double>& sigma_eff) const
{
    this->electrolyte->proton_conductivity(sigma_eff);
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::derivative_effective_proton_conductivity(std::map< VariableNames, std::vector<double> >& dSigma_eff) const
{
    Assert(this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called before NafionMembrane::derivative_effective_proton_conductivity."));
    this->electrolyte->proton_conductivity_derivative(dSigma_eff);
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::effective_water_diffusivity(double& D_w_eff) const
{
    this->electrolyte->water_diffusivity(D_w_eff);
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::effective_water_diffusivity(std::vector<double>& D_w_eff) const
{
    this->electrolyte->water_diffusivity(D_w_eff);
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::derivative_effective_water_diffusivity(std::map< VariableNames, std::vector<double> >& dD_w_eff) const
{
    Assert(this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called before NafionMembrane::derivative_effective_water_diffusivity."));
    this->electrolyte->water_diffusivity_derivative(dD_w_eff);
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::effective_oxygen_diffusivity(double& D_O_eff) const
{
    this->electrolyte->oxygen_diffusivity(D_O_eff);
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::effective_oxygen_diffusivity(std::vector<double>& D_O_eff) const
{
    this->electrolyte->oxygen_diffusivity(D_O_eff);
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::derivative_effective_oxygen_diffusivity(std::map< VariableNames, std::vector<double> >& dD_O_eff) const
{
    Assert(this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called before NafionMembrane::derivative_effective_oxygen_diffusivity."));
    this->electrolyte->oxygen_diffusivity_derivative(dD_O_eff);
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::effective_thermal_conductivity(double& prop_eff) const
{
    if (method_thermal_conductivity == "Given")
        prop_eff = thermal_conductivity;
    
    else
    {
        FcstUtilities::log << "Unknown method to compute effective thermal conductivity in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::effective_thermal_conductivity(std::vector<double>& prop) const
{
    Assert (prop.size() != 0, ExcMessage("Vector not initialized before passing as argument in NafionMembrane::effective_thermal_conductivity."));

    if (method_thermal_conductivity == "Given")
    {
        for (unsigned int i = 0; i<prop.size(); ++i)
        prop[i] = thermal_conductivity;
    }
    else
    {
        FcstUtilities::log << "Unknown method to compute effective thermal conductivity in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::NafionMembrane<dim>::effective_thermal_conductivity(std::vector< Tensor<2,dim> >& prop) const
{
    Assert (prop.size() != 0, ExcMessage("Vector not initialized before passing as argument in NafionMembrane::effective_thermal_conductivity."));
    
    if (method_thermal_conductivity == "Given")
    {
        for (unsigned int i = 0; i<prop.size(); ++i)
            for (unsigned int j = 0; j<dim; ++j)
                prop[i][j][j] = thermal_conductivity;
    }
    else
    {
        FcstUtilities::log << "Unknown method to compute effective thermal conductivity in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
    
}
//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::derivative_effective_thermal_conductivity(std::vector< std::vector<double> >& prop) const
{
    Assert (prop.size() != 0, ExcMessage("Vector not initialized at outer level before passing as argument in NafionMembrane::derivative_effective_thermal_conductivity."));
    Assert (prop[0].size() != 0, ExcMessage("Vector not initialized at inner level before passing as argument in NafionMembrane::derivative_effective_thermal_conductivity."));

    if (method_thermal_conductivity == "Given")
    {
        for (unsigned int i = 0; i<prop.size(); ++i)
            for (unsigned int j = 0; j<prop[i].size(); ++j)
                prop[i][j] = 0.0;
    }

    else
    {
        FcstUtilities::log << "Unknown method to compute derivative of effective thermal conductivity in "<<__FILE__ <<" line "<<__LINE__<<std::endl;
        abort();
    }
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::effective_thermoosmotic_diffusivity(std::vector<double>& D_T_eff) const
{
    this->electrolyte->thermoosmotic_coeff(D_T_eff);
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::NafionMembrane<dim>::derivative_effective_thermoosmotic_diffusivity(std::map< VariableNames, std::vector<double> >& dD_T_eff) const
{
    Assert(this->derivative_flags.size()!=0, ExcMessage("set_derivative_flags has not been probably called before NafionMembrane::derivative_effective_thermoosmotic_diffusivity."));
    this->electrolyte->thermoosmotic_coeff_derivative(dD_T_eff);
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations. 
template class NAME::NafionMembrane<deal_II_dimension>;
