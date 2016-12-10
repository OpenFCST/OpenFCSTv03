//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: PSD_base.h
//    - Description: Base class for pore size distribution model.
//    - Developers: 2009-13 by Marc Secanell, University of Alberta
//                  2013-14 by Jie Zhou, University of Alberta
//    - $ $
//
//---------------------------------------------------------------------------

#include <microscale/PSD_dual.h>

namespace NAME = FuelCellShop::MicroScale;

template<int dim>
const std::string NAME::DualPSD<dim>::concrete_name("DualPSD");

template<int dim>
NAME::DualPSD<dim> const* NAME::DualPSD<dim>::PROTOTYPE = new NAME::DualPSD<dim>();

//---------------------------------------------------------------------------
template<int dim>
NAME::DualPSD<dim>::DualPSD() :
NAME::BasePSD<dim>(),
psd_hi("HIPSD"),
psd_ho("HOPSD")
{
    this->get_mapFactory()->insert(
            std::pair<std::string, FuelCellShop::MicroScale::BasePSD<dim>*>(
                    concrete_name, this));
}

//---------------------------------------------------------------------------
template<int dim>
NAME::DualPSD<dim>::DualPSD(std::string name) :
NAME::BasePSD<dim>(name) 
{
}

//---------------------------------------------------------------------------
template<int dim>
NAME::DualPSD<dim>::~DualPSD() 
{
}

//---------------------------------------------------------------------------

template <int dim>
void
NAME::DualPSD<dim>::declare_parameters (ParameterHandler &param) const
{
  psd_hi.declare_parameters(param);
  
  psd_ho.declare_parameters(param);
  
  FuelCellShop::MicroScale::BasePSD<dim>::declare_parameters(param);
  
  param.enter_subsection("PSD parameters");
  {
      param.enter_subsection("BasePSD");
      {
          param.enter_subsection("DualPSD"); 
          {
              param.declare_entry("capillay pressure",
                                  "0.0",
                                  Patterns::Double(0.0),
                                  "capillay pressure in Pascal");
              
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
NAME::DualPSD<dim>::initialize (ParameterHandler &param)
{
    
    psd_hi.initialize(param);
    
    psd_ho.initialize(param);
    
    FuelCellShop::MicroScale::BasePSD<dim>::_initialize(param);
    
    param.enter_subsection("PSD parameters");
    {
        param.enter_subsection("BasePSD");
        {
            param.enter_subsection("DualPSD"); 
            {
                pressure_c = param.get_double("capillay pressure");
                
            }
            param.leave_subsection(); 
        }
        param.leave_subsection();
    }
  param.leave_subsection();
  
}

//---------------------------------------------------------------------------

template<int dim>
void
NAME::DualPSD<dim>::get_maximum_cross_sectional_areas(double& a_max) const
{
    double a_max_HI,a_max_HO;
    
    psd_hi.get_maximum_cross_sectional_areas(a_max_HI);
    
    psd_ho.get_maximum_cross_sectional_areas(a_max_HO);
    
    a_max =  a_max_HI + a_max_HO;
}

// ---              ---
// --- get_saturation ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::DualPSD<dim>::get_saturation(std::vector<double>& S) const
{
    std::vector<double> S_HI;
    psd_hi.get_saturation(S_HI); 
    std::vector<double> S_HO;
    psd_ho.get_saturation(S_HO);
    S.clear();
    S.resize(S_HI.size());
    
    for(unsigned int i = 0; i < S_HO.size(); ++i)
        
        S[i] = S_HI[i]+S_HO[i];
}
//---------------------------------------------------------------------------
template<int dim>
void
NAME::DualPSD<dim>::get_derivative_saturation(std::vector<double>& derivate_S) const
{
    std::vector<double> S_HI;
    psd_hi.get_derivative_saturation(S_HI); 
    std::vector<double> S_HO;
    psd_ho.get_derivative_saturation(S_HO);
    derivate_S.clear();
    derivate_S.resize(S_HI.size());
    
    for(unsigned int i = 0; i < S_HO.size(); ++i)
        
        derivate_S[i] = S_HI[i]+S_HO[i];
}
// ---              ---
// --- get_pore_saturated_permeability ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::DualPSD<dim>::get_global_saturated_permeability(double& saturated_permeability) const
{
    saturated_permeability = 0.0;
    double p_hi = 0.0;
    double p_ho = 0.0;
    psd_hi.get_global_saturated_permeability(this->get_porosity(),p_hi);
    psd_ho.get_global_saturated_permeability(this->get_porosity(),p_ho);
    saturated_permeability = p_hi + p_ho;
    
}
// ---              ---
// --- get_pore_liquid_permeability ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::DualPSD<dim>::get_relative_liquid_permeability(std::vector<double>& liquid_permeability) const
{
    std::vector<double> S;
    S.clear();
    this->get_saturation(S);
    
    std::vector<double> liquid_HI_permeability;
    psd_hi.get_pore_HI_liquid_saturated_permeability(this->get_porosity(),S,liquid_HI_permeability);
    
    std::vector<double> liquid_HO_permeability;
    psd_ho.get_pore_HO_liquid_saturated_permeability(this->get_porosity(),S,liquid_HO_permeability);
    
    double saturated_permeability;
    this->get_global_saturated_permeability(saturated_permeability);
    
    liquid_permeability.clear();
    liquid_permeability.resize(liquid_HI_permeability.size());
    
    for(unsigned int q = 0; q < liquid_HI_permeability.size(); ++q)
        
        liquid_permeability[q] = (liquid_HI_permeability[q] + liquid_HO_permeability[q]) / saturated_permeability;
    
}

//---------------------------------------------------------------------------
template<int dim>
void
NAME::DualPSD<dim>::get_derivative_relative_liquid_permeability(std::vector<double>& derivative_liquid_permeability) const
{
    std::vector<double> S;
    this->get_saturation(S);
    std::vector<double> ds_dp;
    this->get_derivative_saturation(ds_dp);
    
    std::vector<double> derivative_liquid_HI_permeability;
    psd_hi.get_derivative_pore_HI_liquid_saturated_permeability(this->get_porosity(),S,ds_dp,derivative_liquid_HI_permeability);
    
    std::vector<double> derivative_sliquid_HO_permeability;
    psd_ho.get_derivative_pore_HO_liquid_saturated_permeability(this->get_porosity(),S,ds_dp,derivative_sliquid_HO_permeability);
    
    double saturated_permeability;
    this->get_global_saturated_permeability(saturated_permeability);
    
    derivative_liquid_permeability.clear();
    derivative_liquid_permeability.resize(derivative_liquid_HI_permeability.size());
    
    for(unsigned int q = 0; q < derivative_liquid_HI_permeability.size(); ++q)
        
        derivative_liquid_permeability[q] = (derivative_liquid_HI_permeability[q] + derivative_sliquid_HO_permeability[q]) / saturated_permeability;
    
}
// ---              ---
// --- get_pore_gas_permeability ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::DualPSD<dim>::get_relative_gas_permeability(std::vector<double>& gas_permeability) const
{
    
    std::vector<double> S;
    S.clear();
    get_saturation(S);
    
    std::vector<double> HI_gas_permeability;
    psd_hi.get_pore_HI_liquid_saturated_permeability(this->get_porosity(),S ,HI_gas_permeability);
    
    std::vector<double> HO_gas_permeability;
    psd_ho.get_pore_HO_gas_saturated_permeability(this->get_porosity(),S,HO_gas_permeability);
    
    double saturated_permeability;
    get_global_saturated_permeability(saturated_permeability);
    
    gas_permeability.clear();
    gas_permeability.resize(HI_gas_permeability.size());
    
    for(unsigned int q = 0; q < HI_gas_permeability.size(); ++q)
        
        gas_permeability[q] = (HI_gas_permeability[q] + HO_gas_permeability[q]) / saturated_permeability;
}
// ---              ---
// --- get_pore_liquid_gas_interfacial_surface ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::DualPSD<dim>::get_liquid_gas_interfacial_surface(std::vector<double>& alv) const
{
    
    alv.clear();
    std::vector<double> HI_liquid_gas_interfacial_surface;
    std::vector<double> HO_liquid_gas_interfacial_surface;
    std::vector<double> liquid_gas_interfacial_surface;
    
    double a_max;
    this->get_maximum_cross_sectional_areas(a_max);
    HI_liquid_gas_interfacial_surface.clear();
    HO_liquid_gas_interfacial_surface.clear();
    liquid_gas_interfacial_surface.clear();
    psd_hi.get_liquid_gas_interfacial_surface_withoutPb (HI_liquid_gas_interfacial_surface);
    psd_ho.get_liquid_gas_interfacial_surface_withoutPb (HO_liquid_gas_interfacial_surface);
    liquid_gas_interfacial_surface.resize(HI_liquid_gas_interfacial_surface.size());

    alv.resize(HI_liquid_gas_interfacial_surface.size());
    
    for(unsigned int i = 0; i < HI_liquid_gas_interfacial_surface.size(); ++i)
    {
        liquid_gas_interfacial_surface[i] = HI_liquid_gas_interfacial_surface[i] + HO_liquid_gas_interfacial_surface[i];
    }
    
    for(unsigned int i = 0; i < liquid_gas_interfacial_surface.size(); ++i)
        
        alv[i] = (liquid_gas_interfacial_surface[i]/a_max)
                 *(1- liquid_gas_interfacial_surface[i]/a_max)
                 *liquid_gas_interfacial_surface[i];
}
//---------------------------------------------------------------------------
template<int dim>
void
NAME::DualPSD<dim>::get_derivative_liquid_gas_interfacial_surface(std::vector<double>& derivative_alv) const
{
    derivative_alv.clear();
    std::vector<double> HI_alv;
    std::vector<double> HO_alv;
    psd_hi.get_derivative_liquid_gas_interfacial_surface (HI_alv);
    psd_ho.get_derivative_liquid_gas_interfacial_surface (HO_alv);
    derivative_alv.resize(HI_alv.size());
    
    for (unsigned i=0; i<HI_alv.size(); ++i)
        derivative_alv[i] = HI_alv[i] + HO_alv[i];
    
    
}
// ---              ---
// --- get_pore_wetted_wall ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::DualPSD<dim>::get_wetted_wall_surface_area(std::vector<double>& wetted_wall_surface_area) const
{
    std::vector<double> HI_wetted_wall_surface_area;
    std::vector<double> HO_wetted_wall_surface_area;
    
    psd_hi.get_pore_HI_wetted_wall_surface_area(HI_wetted_wall_surface_area);
    psd_ho.get_pore_HO_wetted_wall_surface_area(HO_wetted_wall_surface_area);
    
    wetted_wall_surface_area.clear();
    wetted_wall_surface_area.resize(HI_wetted_wall_surface_area.size());
    
    for(unsigned int i = 0; i < HI_wetted_wall_surface_area.size(); ++i)
        
        wetted_wall_surface_area[i] = HI_wetted_wall_surface_area[i] + HO_wetted_wall_surface_area[i];
    
}
// ---              ---
// --- get_pore_knudsen_radius ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::DualPSD<dim>::get_knudsen_radius(std::vector<double>& knudsen_radius) const
{
    std::vector<double> knudsen_radius_C1;
    std::vector<double> knudsen_radius_C2;
    std::vector<double> knudsen_radius_C3;
    std::vector<double> knudsen_radius_C4;
    
    psd_hi.get_pore_knudsen_radius_C1(knudsen_radius_C1);
    psd_ho.get_pore_knudsen_radius_C2(knudsen_radius_C2);
    psd_hi.get_pore_knudsen_radius_C3(knudsen_radius_C3);
    psd_ho.get_pore_knudsen_radius_C4(knudsen_radius_C4);
    
    knudsen_radius.clear();
    knudsen_radius.resize(knudsen_radius_C1.size());
    
    for(unsigned int i = 0; i < knudsen_radius_C1.size(); ++i)
        
        knudsen_radius [i] = (knudsen_radius_C1[i] + knudsen_radius_C2[i]) / (knudsen_radius_C3[i] + knudsen_radius_C4[i]);
    
    
}
// ---              ---
// --- get_pore_diffusivity ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::DualPSD<dim>::get_diffusivity() const
{
    
    
}
//---------------------------------------------------------------------------
template <int dim>
void
NAME::DualPSD<dim>::get_PSD_plot (const std::vector<double> critical_radius,std::vector<double>& dx_dr_HI,std::vector<double>& dx_dr_HO, std::vector<double>& dx_dr) const
{
    psd_hi.get_PSD_plot(critical_radius,dx_dr_HI);
    psd_ho.get_PSD_plot(critical_radius,dx_dr_HO);
    for (unsigned int i = 0; i < dx_dr_HI.size(); ++i)
    dx_dr[i] = dx_dr_HI[i] + dx_dr_HO[i];
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::DualPSD<deal_II_dimension>;