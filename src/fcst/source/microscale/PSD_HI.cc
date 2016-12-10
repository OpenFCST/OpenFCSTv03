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

#include <microscale/PSD_HI.h>

namespace NAME = FuelCellShop::MicroScale;

template<int dim>
const std::string NAME::HIPSD<dim>::concrete_name("HIPSD");

template<int dim>
NAME::HIPSD<dim> const* NAME::HIPSD<dim>::PROTOTYPE = new NAME::HIPSD<dim>();

//---------------------------------------------------------------------------
template<int dim>
NAME::HIPSD<dim>::HIPSD() 
:
FuelCellShop::MicroScale::BasePSD<dim>() 
{
    this->get_mapFactory()->insert(std::pair<std::string, FuelCellShop::MicroScale::BasePSD<dim>* >(concrete_name, this));
} 

//---------------------------------------------------------------------------
template<int dim>
NAME::HIPSD<dim>::HIPSD(std::string name) 
:
NAME::BasePSD<dim>(name) 
{}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::HIPSD<dim>::declare_parameters (ParameterHandler &param) const
{
  FuelCellShop::MicroScale::BasePSD<dim>::declare_parameters(param);
  
  param.enter_subsection("PSD parameters");
  {
      param.enter_subsection("BasePSD");
      {
          param.enter_subsection("HIPSD"); 
          {
            
            param.declare_entry("Hydrophilic Mode probability global",
                                "0.0",
                                Patterns::List(Patterns::Double(0.0)),
                                "Contribution of the distribution mode into the PSD ");
            
            param.declare_entry("Hydrophilic Mode characteristic radius global",
                                "0.0",
                                Patterns::List(Patterns::Double(0.0)),
                                "Characteristic pore size of the distribution mode into the PSD ");
            
            param.declare_entry("Hydrophilic Mode width global",
                                "0.0",
                                Patterns::List(Patterns::Double(0.0)),
                                "Characteristic pore size of the distribution mode into the PSD ");
            
            param.declare_entry("capillay pressure",
                                "1.0",
                                Patterns::Double(0.0),
                                "capillay pressure in Pascal");
            
            param.declare_entry("Static Contact Angle HI",
                                "80.0",
                                Patterns::Double(0.0),
                                "Static Contact Angle for hydrophilic");

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
NAME::HIPSD<dim>::initialize (ParameterHandler &param)
{
    
    FuelCellShop::MicroScale::BasePSD<dim>::_initialize(param);
    
    param.enter_subsection("PSD parameters");
    {
        param.enter_subsection("BasePSD");
        {
            param.enter_subsection("HIPSD"); 
            {
                contact_angle_HI = param.get_double("Static Contact Angle HI");
                
                pressure_c = param.get_double("capillay pressure");
                
                fHI_k = FcstUtilities::string_to_number<double>( Utilities::split_string_list( param.get("Hydrophilic Mode probability global") ) );
                
                rHI_k = FcstUtilities::string_to_number<double>( Utilities::split_string_list( param.get("Hydrophilic Mode characteristic radius global") ) );
                
                sHI_k = FcstUtilities::string_to_number<double>( Utilities::split_string_list( param.get("Hydrophilic Mode width global") ) );
            }
            param.leave_subsection();  
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}


// ---              ---
// --- get_critical_radius ---
// ---              ---
//---------------------------------------------------------------------------

template<int dim>
void
NAME::HIPSD<dim>::get_critical_radius(std::vector<double>& critical_radius) const
{
    std::vector<double> p_c;
    p_c.clear();
    critical_radius.clear();
    
    if ( critical_radius_is_initialized )
    {
        critical_radius = critical_radius_computed;
    }
    else 
    {
        if (Capillary_pressure_vector.is_initialized ())
        {
            p_c.resize(this->Capillary_pressure_vector.size());
            
            for (unsigned int q = 0; q < p_c.size(); ++q)
                p_c[q] = this->Capillary_pressure_vector[q];
        }
        else
        {
            p_c.push_back (pressure_c);
        }
        
        critical_radius.resize(p_c.size());
        
        for(unsigned int q = 0; q < p_c.size() ; ++q)
        {
            
            if (p_c[q] >= 0)
            {
                critical_radius[q] = std::numeric_limits<double>::infinity();
            }
            else 
            {
                critical_radius[q] =  -1.0 * 2.0 * this->gamma 
                                      * (   std::cos(contact_angle_HI * Constants::Pi()/180.0) /p_c[q]);  
            }   
        }
    }
}

//---------------------------------------------------------------------------

template<int dim>
void
NAME::HIPSD<dim>::get_derivate_critical_radius(std::vector<double>& derivate_critical_radius) const
{
     Assert( Capillary_pressure_vector.is_initialized (), ExcMessage("Capillary_pressure is not initialized") );
    
    derivate_critical_radius.clear();
    derivate_critical_radius.resize(this->Capillary_pressure_vector.size());

    std::vector<double> p_c,critical_radius;
    
    p_c.resize(this->Capillary_pressure_vector.size());
    critical_radius.resize(this->Capillary_pressure_vector.size());
    
    for (unsigned int q = 0; q < p_c.size(); ++q)
    {
        p_c[q] = this->Capillary_pressure_vector[q];
    }
    
    for(unsigned int q = 0; q < p_c.size() ; ++q)
    {
        if (p_c[q] >= 0) {
            critical_radius[q] = std::numeric_limits<double>::infinity();            
            derivate_critical_radius[q] = 0.0;
        }
        else {
            derivate_critical_radius[q] = 2.0 * this->gamma * (std::cos(contact_angle_HI * Constants::Pi()/180.0) /std::pow(p_c[q],2.0));        
        }
    } 
}

//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_maximum_cross_sectional_areas(double &a_max_HI) const
{
    a_max_HI=0.0;
    
    for(unsigned int i = 0; i < this->sHI_k.size(); ++i){
        a_max_HI +=  this->F_HI * this->fHI_k[i] * std::exp (this->sHI_k[i] * this->sHI_k[i] /2.0 )/(4.0*this->rHI_k[i]);
    }
}
// ---              ---
// --- get_saturation ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_saturation(std::vector<double>& S) const
{
    std::vector<double> critical_radius;
    S.clear();
    critical_radius.clear();
    this->get_critical_radius(critical_radius);
    S.resize(critical_radius.size());
    
    if (saturation_is_initialized)  {
        S = saturation_computed;
    }
    else {
        for(unsigned int q = 0; q < critical_radius.size(); ++q)
            for(unsigned int i = 0; i < rHI_k.size(); ++i) 
                    S [q] += this->F_HI * 0.5 * fHI_k[i] * ( 1.0 + std::erf (   (std::log(critical_radius[q]) - std::log(rHI_k[i]) )/(  sHI_k[i] * std::pow(2.0,0.5)  )));        
    }
}
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_derivative_saturation(std::vector<double>& derivate_S) const
{
    derivate_S.clear();
    derivate_S.resize(this->Capillary_pressure_vector.size());
    
    std::vector<double> critical_radius;
    std::vector<double> derivative_critical_radius;
    this->get_critical_radius(critical_radius);
    this->get_derivate_critical_radius(derivative_critical_radius);
    
    for(unsigned int i = 0; i < critical_radius.size(); ++i)
        for(unsigned int j = 0; j < fHI_k.size(); ++j)
            
                derivate_S[i] +=  this->F_HI 
                                * 0.5 
                                * fHI_k[j] 
                                * (1.0 / critical_radius[i])
                                * (1.0 / (sHI_k[j]*std::sqrt(2.0)))
                                * derivative_critical_radius[i]
                                * 2.0
                                * (1.0 / Constants::Pi())
                                * std::exp (- std::pow( (std::log(critical_radius[i]) - std::log(rHI_k[j])) / (sHI_k[j]*std::sqrt(2.0)),2.0 )); 
      
}
// ---              ---
// --- get_pore_saturated_permeability ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_global_saturated_permeability(double& saturated_permeability) const
{
    saturated_permeability = 0.0;
    
    for(unsigned int i = 0; i < rHI_k.size(); ++i)

        saturated_permeability +=  this->F_HI
                                   *std::pow(this->get_porosity() / this->lamda , 2.0) 
                                   *( std::exp( (-2.0) * sHI_k[i] * sHI_k[i])      
                                   *rHI_k[i] * rHI_k[i] * fHI_k[i])   
                                   *(1.0/8.0); 
}

//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_global_saturated_permeability(const double porosity,
                                                    double& saturated_permeability) const
{
    for(unsigned int i = 0; i < rHI_k.size(); ++i)
        
        saturated_permeability +=  this->F_HI
                                   *std::pow( porosity / this->lamda , 2.0) 
                                   *( std::exp( (-2.0) * sHI_k[i] * sHI_k[i])      
                                   *rHI_k[i] * rHI_k[i] * fHI_k[i])   
                                   *(1.0/8.0);
}
// ---              ---
// --- get_pore_liquid_permeability ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_pore_HI_liquid_saturated_permeability(std::vector<double>& saturated_HI_permeability) const
{
    std::vector<double> S;
    S.clear();
    get_saturation(S);
    
    std::vector<double> critical_radius;
    critical_radius.clear();
    this->get_critical_radius(critical_radius);
    
    saturated_HI_permeability.clear();
    
    saturated_HI_permeability.resize(S.size());
    
    for(unsigned int q = 0; q < S.size(); ++q)
    {
        for(unsigned int i = 0; i < sHI_k.size(); ++i)
        {
            saturated_HI_permeability [q] +=  this->F_HI 
                                              *std::pow(this->get_porosity()* S[q] /this->lamda,2.0) 
                                              *std::exp ( (-2.0) * sHI_k[i] * sHI_k[i] ) 
                                              * rHI_k[i] * rHI_k[i] * fHI_k[i] 
                                              *( std::erf ( (std::log(critical_radius[q]) - std::log(rHI_k[i])) / (sHI_k[i] * std::pow(2.0,0.5)) - sHI_k[i] * std::pow(2.0,0.5) ) + 1)  
                                              * (1.0/16.0) ;
        }
    }   

     
        
}

//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_pore_HI_liquid_saturated_permeability(const double porosity,
                                                            const std::vector<double> S,
                                                            std::vector<double>& saturated_HI_permeability) const
{
    std::vector<double> critical_radius;
    this->get_critical_radius(critical_radius);
    
    saturated_HI_permeability.clear();
    saturated_HI_permeability.resize(S.size());
    
    for(unsigned int q = 0; q < S.size(); ++q)        
        for(unsigned int i = 0; i < sHI_k.size(); ++i) {
                saturated_HI_permeability [q] +=  this->F_HI * (1.0/16.0)
                                              * std::pow(porosity * S[q] /this->lamda,2.0) 
                                              * std::exp ( (-2.0) * sHI_k[i] * sHI_k[i] ) 
                                              * rHI_k[i] * rHI_k[i] * fHI_k[i]
                                              * ( std::erf ( (std::log(critical_radius[q]) - std::log(rHI_k[i])) / (sHI_k[i] * std::pow(2.0,0.5)) - sHI_k[i] * std::pow(2.0,0.5) ) + 1.0);                                    
            }
            
}
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_derivative_pore_HI_liquid_saturated_permeability(std::vector<double>& derivative_kr_l) const
{
}

//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_derivative_pore_HI_liquid_saturated_permeability(const double porosity, const std::vector<double> S,
                                                                       const std::vector<double> ds_dp,
                                                                        std::vector<double>& derivative_kr_l) const
{
    std::vector<double> critical_radius;
    std::vector<double> derivative_critical_radius;
    this->get_critical_radius(critical_radius);
    this->get_derivate_critical_radius(derivative_critical_radius);
    
    derivative_kr_l.clear();
    derivative_kr_l.resize(S.size());
    
    for(unsigned int q = 0; q < S.size(); ++q)        
        for(unsigned int i = 0; i < sHI_k.size(); ++i)            
            derivative_kr_l [q] += this->F_HI * (1.0/16.0)*std::pow(porosity* S[q] /this->lamda,2.0) 
                                   *std::exp ( (-2.0) * sHI_k[i] * sHI_k[i] ) * rHI_k[i] * rHI_k[i] * fHI_k[i] 
                                   * (1.0/ (sHI_k[i] * std::sqrt(2.0))) * (2.0/std::sqrt(Constants::Pi())) * (1.0/critical_radius[q])
                                   * derivative_critical_radius[q] 
                                   * std::exp ( - (  std::pow( (sHI_k[i] * std::sqrt(2.0)) - (std::log(critical_radius[q]) -std::log(rHI_k[i]))/ (sHI_k[i]*std::sqrt(2.0)) ,2.0) ) )
                                   +
                                   this->F_HI * (1.0/16.0) * 2.0 * (porosity * porosity * S[q] / (this->lamda * this->lamda) )  
                                   * std::exp ( (-2.0) * sHI_k[i] * sHI_k[i] ) * rHI_k[i] * rHI_k[i] * fHI_k[i] 
                                   * ( 1 - std::erf ( sHI_k[i] * std::sqrt(2.0) - (std::log(critical_radius[q]) - std::log(rHI_k[i])) / (sHI_k[i] * std::sqrt(2.0))  ) ) * ds_dp[q];
}



//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_relative_liquid_permeability(std::vector<double>& liquid_permeability) const
{
    std::vector<double> liquid_HI_permeability;
    liquid_HI_permeability.clear();
    this->get_pore_HI_liquid_saturated_permeability(liquid_HI_permeability);
    
    double saturated_permeability;
    this->get_global_saturated_permeability(saturated_permeability);
    
    liquid_permeability.clear();    
    liquid_permeability.resize(liquid_HI_permeability.size());
    
    for(unsigned int q = 0; q < liquid_HI_permeability.size(); ++q)
    {
        liquid_permeability[q] = liquid_HI_permeability[q] / saturated_permeability;
    }
}


//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_derivative_relative_liquid_permeability(std::vector<double>& derivative_liquid_permeability) const
{
    std::vector<double> derivative_liquid_HI_permeability;
    derivative_liquid_HI_permeability.clear();
    this->get_derivative_pore_HI_liquid_saturated_permeability(derivative_liquid_HI_permeability);
    
    double saturated_permeability;
    this->get_global_saturated_permeability(saturated_permeability);
    
    derivative_liquid_permeability.clear();
    derivative_liquid_permeability.resize(derivative_liquid_HI_permeability.size());
    
    for(unsigned int q = 0; q < derivative_liquid_HI_permeability.size(); ++q)
    {
        derivative_liquid_permeability[q] = derivative_liquid_HI_permeability[q] / saturated_permeability;
    }
}
// ---              ---
// --- get_pore_gas_permeability ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_pore_HI_gas_saturated_permeability(std::vector<double>& saturated_HI_permeability) const
{
    std::vector<double> S;
    S.clear();
    get_saturation(S);
    
    std::vector<double> critical_radius;
    critical_radius.clear();
    get_critical_radius(critical_radius);
    
    saturated_HI_permeability.clear();
    saturated_HI_permeability.resize(S.size());
    
    for(unsigned int q = 0; q < S.size(); ++q)        
        for(unsigned int i = 0; i < sHI_k.size(); ++i)            
            saturated_HI_permeability[q] += this->F_HI
                                            * std::pow(this->get_porosity() * (1.0 - S[q])/this->lamda,2.0) 
                                            * (std::exp ((-2.0) * sHI_k[i] * sHI_k[i]) * rHI_k[i] * rHI_k[i] * fHI_k[i]) 
                                            * (1.0/16.0) 
                                            * ( -std::erf ( (std::log(critical_radius[q]) - std::log(rHI_k[i]))/(sHI_k[i] * std::pow(2.0,0.5)) - sHI_k[i] * std::pow(2.0,0.5) ) + 1.0 );
        
        
}
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_pore_HI_gas_saturated_permeability(const double porosity,
                                                         const std::vector<double> S,
                                                         std::vector<double>& saturated_HI_permeability) const
{
    std::vector<double> critical_radius;
    critical_radius.clear();
    get_critical_radius(critical_radius);
    
    saturated_HI_permeability.clear();
    saturated_HI_permeability.resize(S.size());
    
    for(unsigned int q = 0; q < S.size(); ++q)
        
        for(unsigned int i = 0; i < sHI_k.size(); ++i)
            
            saturated_HI_permeability[q] += this->F_HI
                                            * std::pow( porosity * (1.0 - S[q])/this->lamda,2.0) 
                                            * (std::exp ((-2.0) * sHI_k[i] * sHI_k[i]) * rHI_k[i] * rHI_k[i] * fHI_k[i]) 
                                            * (1.0/16.0) 
                                            * ( -std::erf ( (std::log(critical_radius[q]) - std::log(rHI_k[i]))/(sHI_k[i] * std::pow(2.0,0.5)) - sHI_k[i] * std::pow(2.0,0.5) ) + 1.0 );
          
}
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_relative_gas_permeability(std::vector<double>& gas_permeability) const
{
    std::vector<double> HI_gas_permeability;
    this->get_pore_HI_gas_saturated_permeability(HI_gas_permeability);
    
    double saturated_permeability;
    this->get_global_saturated_permeability(saturated_permeability);
    
    gas_permeability.clear();
    gas_permeability.resize(HI_gas_permeability.size());
    
    for(unsigned int q = 0; q < HI_gas_permeability.size(); ++q)
        
        gas_permeability[q] = HI_gas_permeability[q]  / saturated_permeability;
}

// ---              ---
// --- get_pore_liquid_gas_interfacial_surface ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_liquid_gas_interfacial_surface_withoutPb(std::vector<double>& HI_liquid_gas_interfacial_surface_a) const
{
    std::vector<double> critical_radius;
    this->get_critical_radius(critical_radius);
    HI_liquid_gas_interfacial_surface_a.clear();
    HI_liquid_gas_interfacial_surface_a.resize(critical_radius.size());
    
    for(unsigned int q = 0; q <critical_radius.size(); ++q)
    
        for(unsigned int i = 0; i < sHI_k.size(); ++i)
        
            HI_liquid_gas_interfacial_surface_a[q] +=   this->F_HI 
                                                        * fHI_k[i] 
                                                        * std::exp(sHI_k[i] * sHI_k[i] /2.0 )  
                                                        * (1.0/ rHI_k[i])  
                                                        * ( 1.0 + std::erf ( (std::log(critical_radius[q]) - std::log(rHI_k[i])) / (sHI_k[i] * std::sqrt(2.0)) + sHI_k[i] * std::sqrt(2.0) / 2.0) )
                                                        / 8.0;
}
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_liquid_gas_interfacial_surface(std::vector<double>& HI_liquid_gas_interfacial_surface) const
{
    HI_liquid_gas_interfacial_surface.clear();
    double a_max ;
    this->get_maximum_cross_sectional_areas(a_max);
    
    std::vector<double> HI_liquid_gas_interfacial_surface_a;
    this->get_liquid_gas_interfacial_surface_withoutPb(HI_liquid_gas_interfacial_surface_a);
    HI_liquid_gas_interfacial_surface.resize(HI_liquid_gas_interfacial_surface_a.size());
    
    for(unsigned int i = 0; i <HI_liquid_gas_interfacial_surface_a.size(); ++i)
    
        HI_liquid_gas_interfacial_surface [i] =  HI_liquid_gas_interfacial_surface_a[i]/a_max 
                                                * (1 - HI_liquid_gas_interfacial_surface_a[i]/a_max ) 
                                                * HI_liquid_gas_interfacial_surface_a[i];
}
//---------------------------------------------------------------------------

template<int dim>
void
NAME::HIPSD<dim>::get_derivative_liquid_gas_interfacial_surface_increment(std::vector<double>& derivate_alv) const
{
    derivate_alv.clear();
    std::vector<double> critical_radius;
    std::vector<double> derivative_critical_radius;
    this->get_critical_radius(critical_radius);
    this->get_derivate_critical_radius(derivative_critical_radius);
    derivate_alv.resize(critical_radius.size());
    
    for (unsigned i = 0 ; i < critical_radius.size(); ++i)
    {
        for(unsigned int k = 0; k < sHI_k.size(); ++k)
        {
                derivate_alv[i] += this->F_HI 
                                   * fHI_k[k] 
                                   * std::exp(sHI_k[k] * sHI_k[k] /2.0 )  
                                   * (1.0/(rHI_k[k]*8.0))  
                                   
                                   * (1.0/(sHI_k[k] * std::sqrt(2.0)))
                                   
                                   * (2.0/std::sqrt(Constants::Pi()))
                                   
                                   * std::exp(- std::pow( (std::log(critical_radius[i]) - std::log(rHI_k[k])) / (std::log(sHI_k[k]) * std::sqrt(2.0)) + (std::log(sHI_k[k])/2.0) ,2.0) )
                                   
                                   * (1.0/critical_radius[i])
                                   * derivative_critical_radius[i];
        }
   }
}
//---------------------------------------------------------------------------

template<int dim>
void
NAME::HIPSD<dim>::get_derivative_liquid_gas_interfacial_surface(std::vector<double>& derivative_alv) const
{
    derivative_alv.clear();
    std::vector<double> alv;
    std::vector<double> derivative;
    this->get_liquid_gas_interfacial_surface_withoutPb(alv);
    this->get_derivative_liquid_gas_interfacial_surface_increment(derivative);
    derivative_alv.resize(derivative.size());
    
    double a_max;
    this->get_maximum_cross_sectional_areas(a_max);
    
    for (unsigned i = 0; i<alv.size() ; ++i)
    {                                            
        derivative_alv[i] = 2.0 * alv[i] * derivative[i] / a_max
                            -
                            3.0 * alv[i] * alv[i] * derivative[i] / (a_max*a_max);
    }
}

// ---              ---
// --- get_pore_wetted_wall ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_pore_HI_wetted_wall_surface_area(std::vector<double>& HI_wetted_wall_surface_area) const
{
    std::vector<double> critical_radius;
    get_critical_radius(critical_radius);
    
    HI_wetted_wall_surface_area.clear();
    HI_wetted_wall_surface_area.resize(critical_radius.size());
    
    for(unsigned int i = 0; i < critical_radius.size(); ++i)
        
        for(unsigned int q = 0; q < rHI_k.size(); ++q)
            
            HI_wetted_wall_surface_area[i] = this->F_HI 
                                             * (fHI_k[q] * std::exp (sHI_k[q]*sHI_k[q]/2.0) / rHI_k[q]) 
                                             * (1.0 + std::erf ( (std::log(critical_radius[i]) - std::log (rHI_k[q])) / (sHI_k[q] * std::sqrt(2.0)) + sHI_k[q] / std::sqrt(2.0)));
        
}
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_wetted_wall_surface_area(std::vector<double>& wetted_wall_surface_area) const
{
    std::vector<double> HI_wetted_wall_surface_area;
    get_pore_HI_wetted_wall_surface_area(HI_wetted_wall_surface_area);
    
    wetted_wall_surface_area.clear();
    wetted_wall_surface_area.resize(HI_wetted_wall_surface_area.size());
    
    for(unsigned int i = 0; i < HI_wetted_wall_surface_area.size(); ++i)
        
        wetted_wall_surface_area[i] = HI_wetted_wall_surface_area[i];
    
}
// ---              ---
// --- get_pore_knudsen_radius ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_pore_knudsen_radius_C1(std::vector<double>& knudsen_radius_C1) const
{
    
    std::vector<double> critical_radius;
    this->get_critical_radius(critical_radius);
    
    knudsen_radius_C1.clear();
    knudsen_radius_C1.resize(critical_radius.size());
    
    for(unsigned int i = 0; i < critical_radius.size(); ++i)
        
        for(unsigned int q = 0; q < rHI_k.size(); ++q)
            {
                knudsen_radius_C1 [i] += this->F_HI 
                                        * 0.5 
                                        * fHI_k[q] 
                                        *( 1.0 - std::erf ((std::log(critical_radius[i]) - std::log(rHI_k[q]))/(sHI_k[q] * std::pow(2.0,0.5)))); 
            }  
        
}
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_pore_knudsen_radius_C3(std::vector<double>& knudsen_radius_C3) const
{
    std::vector<double> critical_radius;
    this->get_critical_radius(critical_radius);
    
    knudsen_radius_C3.clear();
    knudsen_radius_C3.resize(critical_radius.size());
    
    for(unsigned int i = 0; i < critical_radius.size(); ++i)
        
        for(unsigned int q = 0; q < rHI_k.size(); ++q)
            
            knudsen_radius_C3[i] += this->F_HI 
                                   * (fHI_k[q] * std::exp (sHI_k[q]*sHI_k[q]/2.0) / (rHI_k[q])) 
                                   * ( 1.0 - std::erf ( (std::log(critical_radius[i]) - std::log (rHI_k[q]) ) / (sHI_k[q] * std::sqrt(2.0)) + sHI_k[q] / std::sqrt(2.0)));
        
        
}
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_knudsen_radius(std::vector<double>& knudsen_radius) const
{
    std::vector<double> knudsen_radius_C1;
    std::vector<double> knudsen_radius_C3;
    
    this->get_pore_knudsen_radius_C1(knudsen_radius_C1);
    this->get_pore_knudsen_radius_C3(knudsen_radius_C3);
    
    knudsen_radius.clear();
    knudsen_radius.resize(knudsen_radius_C1.size());
    
    for(unsigned int i = 0; i < knudsen_radius_C1.size(); ++i)
        
        knudsen_radius [i] =  knudsen_radius_C1[i] / knudsen_radius_C3[i];
    
    
}
// ---              ---
// --- get_pore_diffusivity ---
// ---              ---
//---------------------------------------------------------------------------
template<int dim>
void
NAME::HIPSD<dim>::get_diffusivity() const
{
    
    
}
//---------------------------------------------------------------------------
template <int dim>
void
NAME::HIPSD<dim>::get_PSD_plot (const std::vector<double> critical_radius, std::vector<double>& dx_dr_HI) const
{
    std::vector<double> E_HI;
    E_HI.resize(critical_radius.size());
    dx_dr_HI.resize(critical_radius.size());
    for(unsigned int i = 0; i < critical_radius.size(); ++i)
    {
        for (unsigned int k = 0; k < this->rHI_k.size(); ++k)
        {
            E_HI[i] += std::exp(-1.0 * std::pow( (std::log(critical_radius[i]) - std::log(this->rHI_k[k]) )/ (sHI_k[k]*std::sqrt(2.0)),2.0)   )  ;
        }
        for (unsigned int k = 0; k < rHI_k.size(); ++k)
        {
            dx_dr_HI[i] =  this->F_HI * E_HI[i] * fHI_k[k] / (critical_radius[i] * sHI_k[k] * std::sqrt(2.0*Constants::Pi()));
        }
    }
    
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::HIPSD<deal_II_dimension>;