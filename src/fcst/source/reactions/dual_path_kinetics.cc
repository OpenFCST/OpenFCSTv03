//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: dual_path_kinetics.cc
//    - Description: Dual Path Kinetics model for hydrogen oxidation reaction
//    - Developers: M. Secanell, M. Moore and Madhur Bhaiya
//    - $Id: dual_path_kinetics.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include "reactions/dual_path_kinetics.h"

namespace NAME = FuelCellShop::Kinetics;

const std::string NAME::DualPathKinetics::concrete_name ("DualPathKinetics");

NAME::DualPathKinetics const* NAME::DualPathKinetics::PROTOTYPE = new NAME::DualPathKinetics(true);

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
NAME::DualPathKinetics::DualPathKinetics(const bool construct_replica)
:
BaseKinetics ()
{
    if (construct_replica)
        this->get_mapFactory()->insert(std::pair<std::string, BaseKinetics* > (this->concrete_name, this) );           
}

//---------------------------------------------------------------------------
NAME::DualPathKinetics::DualPathKinetics()
:
BaseKinetics ()
{
    //FcstUtilities::log << "Dual Path kinetics." << std::endl;
    
    // Default values for kinetic parameters.
    j_0T = 0.47;        // A/cm^2
    j_0H = 0.01;        // A/cm^2
    potential_constant = 1.2;
    ref_potential = 0.;
}

//---------------------------------------------------------------------------
NAME::DualPathKinetics::~DualPathKinetics()
{}

//---------------------------------------------------------------------------
void
NAME::DualPathKinetics::declare_parameters(ParameterHandler &param) const
{
    param.enter_subsection("Kinetics");
    {
        param.enter_subsection(NAME::DualPathKinetics::concrete_name);
        {
            param.declare_entry ("TV exchange current density, [A/cm^2]",
                                 "0.47",
                                 Patterns::Double(),
                                 "Exchange current density for TV pathway");
            param.declare_entry ("HV exchange current density, [A/cm^2]",
                                 "0.01",
                                 Patterns::Double(),
                                 "Exchange current density for HV pathway");
            param.declare_entry ("Potential range constant",
                                 "1.2",
                                 Patterns::Double(),
                                 "Potential range constant");
            param.declare_entry ("Reference potential",
                                 "0.0",
                                 Patterns::Double(),
                                 "Reference potential");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
void
NAME::DualPathKinetics::initialize(ParameterHandler &param)
{
    param.enter_subsection("Kinetics"); 
    {
        param.enter_subsection(NAME::DualPathKinetics::concrete_name); 
        {
            j_0T = param.get_double("TV exchange current density, [A/cm^2]");
            j_0H = param.get_double("HV exchange current density, [A/cm^2]");
            potential_constant = param.get_double("Potential range constant");
            ref_potential = param.get_double("Reference potential");
        }
        param.leave_subsection(); 
    }
    param.leave_subsection();
}  

//---------------------------------------------------------------------------

void
NAME::DualPathKinetics::current_density (std::vector<double>& coef)
{
    if (!kin_param_initialized)
        init_kin_param();
    
    n_quad = phi_m.size();
    coef.resize(n_quad);
    
    // Compute the coefficient, first loop over the quadrature points
    for (unsigned int i=0; i<n_quad; ++i)
    {        
        if (reactants_map.at(hydrogen_concentration)[i] <= 0.0)
            coef[i] = 0.0;
        else
        {
            double eta = phi_s[i]-phi_m[i]-catalyst->voltage_cell_th(T[i]);
            
            // Now compute the current density, assuming that the backward reaction is negligible.
            // See equation [25]
            coef[i] = (reactants_map.at(hydrogen_concentration)[i]/ref_conc_H2)*
            ( j_0T*( 1.0 - exp( ((-2.0)*F*eta)/(potential_constant*R*T[i]) ) ) + j_0H*( exp( (F*eta)/(2.0*R*T[i]) ) - ( exp( (-F*eta)/(potential_constant*R*T[i]) ) )*(exp( (-F*eta)/(2.0*R*T[i]) ))  ));
        }
        
        if (std::isinf(coef[i]))
            coef[i] = std::numeric_limits<double>::max()*1e-150;
        else if (std::isnan(coef[i]))
            coef[i] = 0.0;
    }
}

//---------------------------------------------------------------------------

void
NAME::DualPathKinetics::derivative_current (std::map< VariableNames, std::vector<double> >& dcoef_du)
{
    Assert( derivative_flags.size() != 0, ExcMessage("Derivative flags are not set using set_derivative_flags method before the DualPathKinetics::derivative_current.") );
    Assert( electrolyte->get_H_H2() != 0.0, ExcMessage("Derivatives at the moment are defined only for Thin film (Henry's law) case."));
    
    if (!kin_param_initialized)
        init_kin_param();
    
    n_quad = phi_m.size();
    // Loop over the flags
    for (unsigned int i=0; i<derivative_flags.size(); ++i)
    {
        std::vector<double> dcurrent(n_quad, 0.);         // Defaulting all derivatives to 0.0
        
        if ( (derivative_flags[i] == hydrogen_molar_fraction) || (derivative_flags[i] == hydrogen_concentration))
        {
            for (unsigned int j=0; j<n_quad; ++j)
            {
                if (reactants_map.at(hydrogen_concentration)[j] <= 0.0)
                    dcurrent[j] = 0.0;
                
                else
                {
                    double eta = phi_s[j]-phi_m[j]-catalyst->voltage_cell_th(T[j]);
                    dcurrent[j] =  (1.0/ref_conc_H2)*(j_0T*( 1.0 - exp( ((-2.0)*F*eta)/(potential_constant*R*T[j]) ) ) +
                    j_0H*( exp( (F*eta)/(2.0*R*T[j]) ) - ( exp( (-F*eta)/(potential_constant*R*T[j]) ) )*(exp( (-F*eta)/(2.0*R*T[j]) ))  ));
                    
                    if (derivative_flags[i] == hydrogen_molar_fraction)
                        dcurrent[j] *= (p_total/(electrolyte->get_H_H2()));
                }
                
                if (std::isinf(dcurrent[j]))
                    dcurrent[j] = std::numeric_limits<double>::max()*1e-150;
                else if (std::isnan(dcurrent[j]))
                    dcurrent[j] = 0.0;
            }
             if (derivative_flags[i] == hydrogen_molar_fraction)
                 dcoef_du[hydrogen_molar_fraction] = dcurrent;
             else
                 dcoef_du[hydrogen_concentration] = dcurrent;
        }//x_H2
        
        else if (derivative_flags[i] == electronic_electrical_potential)
        { 
            for (unsigned int j=0; j<n_quad; ++j)
            {
                if (reactants_map.at(hydrogen_concentration)[j] <= 0.0)
                    dcurrent[j] = 0.0;
                
                else
                {
                    double eta = phi_s[j]-phi_m[j]-catalyst->voltage_cell_th(T[j]);
                    dcurrent[j] = (reactants_map.at(hydrogen_concentration)[j]/(ref_conc_H2))*
                    ( j_0T*((2.0*F)/(potential_constant*R*T[j]))*exp((-2.0*F*(eta))/(potential_constant*R*T[j])) + 
                    j_0H*( (F/(2.0*R*T[j]))*exp((F*(eta))/(2.0*R*T[j])) + 
                    (F/(R*T[j]))*(0.5+(1.0/potential_constant))*exp((-F*(eta))/(potential_constant*R*T[j]))*exp((-F*(eta))/(2.0*R*T[j]))) );
                }
                
                if (std::isinf(dcurrent[j]))
                    dcurrent[j] = std::numeric_limits<double>::max()*1e-150;
                else if (std::isnan(dcurrent[j]))
                    dcurrent[j] = 0.0;
            }
            
            dcoef_du[electronic_electrical_potential] = dcurrent;
        }//phi_s
        
        else if (derivative_flags[i] == protonic_electrical_potential)
        {
            for (unsigned int j=0; j<n_quad; ++j)
            {
                if (reactants_map.at(hydrogen_concentration)[j] <= 0.0)
                    dcurrent[j] = 0.0;
                
                else
                {
                    double eta = phi_s[j]-phi_m[j]-catalyst->voltage_cell_th(T[j]);
                    dcurrent[j] = (reactants_map.at(hydrogen_concentration)[j]/(ref_conc_H2))*
                    ( j_0T*((-2.0*F)/(potential_constant*R*T[j]))*exp((-2.0*F*(eta))/(potential_constant*R*T[j])) + 
                    j_0H*( (-F/(2.0*R*T[j]))*exp((F*(eta))/(2.0*R*T[j])) + 
                    (F/(R*T[j]))*(-0.5-(1.0/potential_constant))*exp((-F*(eta))/(potential_constant*R*T[j]))*exp((-F*(eta))/(2.0*R*T[j]))));
                }
                
                if (std::isinf(dcurrent[j]))
                    dcurrent[j] = std::numeric_limits<double>::max()*1e-150;
                else if (std::isnan(dcurrent[j]))
                    dcurrent[j] = 0.0;    
            }
            
            dcoef_du[protonic_electrical_potential] = dcurrent;
        }//phi_m
        
        else if (derivative_flags[i] == temperature_of_REV)
        {
            for (unsigned int j=0; j<n_quad; ++j)
            {
                if (reactants_map.at(hydrogen_concentration)[j] <= 0.0)
                    dcurrent[j] = 0.0;
                
                else
                {
                    double eta = phi_s[j]-phi_m[j]-catalyst->voltage_cell_th(T[j]);
                    
                    dcurrent[j] = (reactants_map.at(hydrogen_concentration)[j]/ref_conc_H2)*(F/R)*
                    ( ((2.0*j_0T)/potential_constant)*exp((-2.0)*F*eta/(potential_constant*R*T[j])) +
                    j_0H*( 0.5*exp(F*eta/(2.0*R*T[j])) + 
                    exp(-F*eta/(potential_constant*R*T[j]))*exp(-F*eta/(2.0*R*T[j]))*(0.5 + (1.0/potential_constant)) ) )*
                    (-1.0)*((T[j]*catalyst->dvoltage_cell_th_dT(T[j])  + eta)/(T[j]*T[j]));
                }
                
                if (std::isinf(dcurrent[j]))
                    dcurrent[j] = std::numeric_limits<double>::max()*1e-150;
                else if (std::isnan(dcurrent[j]))
                    dcurrent[j] = 0.0;   
            }
            
            dcoef_du[temperature_of_REV] = dcurrent;
        }//T
        
        else
        {
            dcoef_du[ derivative_flags[i] ] = dcurrent;
        }//everything else
    }//end of derivative_flags loop  
}

// ---------------This is to test derivative_current method using numerical differentiation -----------------------------
//----------------Based on old format --needs to be changed..Madhur -----------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------
// 
//      //size checking for vectors
//      dcoef_du.clear();
//      dcoef_du.resize(this->derivative_flags.size(), std::vector<double>(n_quad, 0.0));
//      
//      unsigned int location = 9;
//      
//      for (unsigned int k=0; k<species_names.size(); ++k)
//      {
//              //Find the location of Hydrogen molar fraction in the species_names list
//              if (species_names[k] == "Hydrogen molar fraction")
//              {
//                      location = k;
//              }
//      }
//      
//      //Loop over the flags
//      for (unsigned int i=0; i<this->derivative_flags.size(); ++i)
//      {
//              if (this->derivative_flags[i] == "Hydrogen molar fraction")
//              {
//                      //Compute the current using the current solution.
//                      std::vector<double> coef1(std::vector<double>(n_quad, 0.0));
// 
//                      current_density(coef1);
//                      double delta = pow(10, -7);
//                      //Modify solution
//                      for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                      {
//                              species_solution[location][j] += delta;
//                      }
//                      
//                      //Compute the new current
//                      std::vector<double> coef2(std::vector<double>(n_quad, 0.0));
//                      current_density(coef2);
//                      
//                      //Compute the coefficient
//                      for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                      {
//                              dcoef_du[i][j] = ((coef2[j] - coef1[j])*(this->p_t/(this->electrolyte->get_H_H2())))/delta;
//                              //FcstUtilities::log << "hydrogen numerical " << dcoef_du[i][j] << std::endl;
//                      }
//                      for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                      {
//                              species_solution[location][j] -= delta;
//                      }      
//              }//x_H2
//              
//              
//              else if (this->derivative_flags[i] == "Electronic electrical potential")
//              { 
//                      //Compute the current using the current solution.
//                      std::vector<double> coef1(std::vector<double>(n_quad, 0.0));
//                      current_density(coef1);
//                      double delta = pow(10, -7);
//                      
//                      //Modify solution
//                      for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                      {
//                              phi_s[j] += delta;
//                      }
//                      
//                      //Compute the new current
//                      std::vector<double> coef2(std::vector<double>(n_quad, 0.0));
//                      current_density(coef2);
//                      
//                      //Compute the coefficient
//                      for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                      {
//                              dcoef_du[i][j] = (coef2[j] - coef1[j])/delta;
//                              //FcstUtilities::log << "electron numerical " << dcoef_du[i][j] << std::endl;
//                      }
//                      //Restore the solution to the correct value.
//                      for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                      {
//                              phi_s[j] -= delta;
//                      }        
//                      
//              }//phi_s
//              
//              else if (this->derivative_flags[i] == "Protonic electrical potential")
//              {
//                      //Compute the current using the current solution.
//                      std::vector<double> coef1(std::vector<double>(n_quad, 0.0));
//                      current_density(coef1);
//                      
//                      double delta = pow(10, -7);
//                      
//                      //Modify solution
//                      for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                      {
//                              phi_m[j] += delta;
//                      }
//                      
//                      //Compute the new current
//                      std::vector<double> coef2(std::vector<double>(n_quad, 0.0));
//                      current_density(coef2);
//                      
//                      //Compute the coefficient
//                      for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                      {
//                              dcoef_du[i][j] = (coef2[j] - coef1[j])/delta;
//                              //FcstUtilities::log << "proton numerical " << dcoef_du[i][j] << std::endl;
//                      }    
//                      for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                      {
//                              phi_m[j] -= delta;
//                      }       
//              }//phi_m
//              
//              else if (this->derivative_flags[i] == "Temperature")
//              {
//                      //Compute the current using the current solution.
//                      std::vector<double> coef1(std::vector<double>(n_quad, 0.0));
//                      current_density(coef1);
//                      
//                      double delta = pow(10.0, -10);
//                      
//                      //Modify solution
//                      for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                      {
//                              T[j] += delta;
//                      }
//                      
//                      //Compute the new current
//                      std::vector<double> coef2(std::vector<double>(n_quad, 0.0));
//                      current_density(coef2);
//                      
//                      //Compute the coefficient
//                      for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                      {
//                              dcoef_du[i][j] = (coef2[j] - coef1[j])/delta;
//                              FcstUtilities::log << "temperature numerical " << dcoef_du[i][j] << std::endl;
//                      }    
//                      for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                      {
//                              T[j] -= delta;
//                      }       
//              }//T
//              
//              else
//              {
//                      //FcstUtilities::log << "Derivative flag " << this->derivative_flags[i] << " not implemented." << std::endl;
//              }//everything else
//              
//      }//end of derivative_flags loop
