//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: bv_kinetics.cc
//    - Description: Butler-Volmer Kinetics model class
//    - Developers: M. Secanell, M. Moore and M. Bhaiya
//    - $Id: bv_kinetics.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include "reactions/bv_kinetics.h"

namespace NAME = FuelCellShop::Kinetics;

const std::string NAME::ButlerVolmerKinetics::concrete_name ("ButlerVolmerKinetics");

NAME::ButlerVolmerKinetics const* NAME::ButlerVolmerKinetics::PROTOTYPE = new NAME::ButlerVolmerKinetics(true);

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
NAME::ButlerVolmerKinetics::ButlerVolmerKinetics(const bool construct_replica)
:
BaseKinetics ()
{
    if (construct_replica)
        this->get_mapFactory()->insert(std::pair<std::string, BaseKinetics* > (this->concrete_name, this) );           
}

//---------------------------------------------------------------------------

NAME::ButlerVolmerKinetics::ButlerVolmerKinetics()
:
BaseKinetics ()
{
    FcstUtilities::log << "Butler Volmer kinetics." << std::endl;
}

//---------------------------------------------------------------------------
NAME::ButlerVolmerKinetics::~ButlerVolmerKinetics()
{}

//---------------------------------------------------------------------------
void
NAME::ButlerVolmerKinetics::current_density (std::vector<double>& coef) 
{
    if (!kin_param_initialized)
        init_kin_param();
    
    n_quad = phi_m.size();
    coef.resize(n_quad);
    double species_comp;
    
    // Compute the coefficient, first loop over the quadrature points
    for (unsigned int i=0; i<n_quad; ++i)
    {
        species_comp = 1.;
        //loop over each species in the reaction and compute its contribution to the current density
        for (std::map< VariableNames, SolutionVariable >::const_iterator iter=reactants_map.begin(); iter != reactants_map.end(); ++iter)
        {
            if (iter->second[i] < 0.0)
                species_comp = 0.0;

            else
                species_comp = species_comp*( pow(iter->second[i]/ref_conc.at(iter->first), gamma.at(iter->first)));
        }
        
        coef[i] = (catalyst->exchange_current_density(T[i]))*species_comp*(  exp( alpha_a*F*(phi_s[i]-phi_m[i]-(catalyst->voltage_cell_th(T[i])))/(R*T[i]) ) 
                                                                                            - 
                                                                             exp((-alpha_c*F*(phi_s[i]-phi_m[i]-(catalyst->voltage_cell_th(T[i]))))/(R*T[i])) );
        
        if (std::isinf(coef[i]))
            coef[i] = std::numeric_limits<double>::max()*1e-150;
        else if (std::isnan(coef[i]))
            coef[i] = 0.0;
    }
}

//---------------------------------------------------------------------------
void
NAME::ButlerVolmerKinetics::derivative_current (std::map< VariableNames, std::vector<double> >& dcoef_du) 
{
    Assert( derivative_flags.size() != 0, ExcMessage("Derivative flags are not set using set_derivative_flags method before the ButlerVolmerKinetics::derivative_current.") );
    
    if (!kin_param_initialized)
        init_kin_param();
    
    n_quad = phi_m.size();
    double species_comp;

    // Loop over the flags
    for (unsigned int i=0; i<derivative_flags.size(); ++i)
    {
        std::vector<double> dcurrent(n_quad, 0.);         // Defaulting all derivatives to 0.0
        
        if ( (derivative_flags[i] == oxygen_molar_fraction)  || 
             (derivative_flags[i] == hydrogen_molar_fraction)||
             (reactants_map.find(derivative_flags[i])!=reactants_map.end()) )
        {
            for (unsigned int j=0; j<n_quad; ++j)
            {
                species_comp = 1.;
                // Loop over each species in the reaction and compute its contribution
                for (std::map< VariableNames, SolutionVariable >::const_iterator iter=reactants_map.begin(); iter != reactants_map.end(); ++iter)
                {
                    if (iter->second[j] < 0.0)
                    {
                        species_comp *= 0.0;
                    }
                    
                    else if ( derivative_flags[i] == oxygen_molar_fraction && iter->first == oxygen_concentration )
                    {
                        Assert( electrolyte->get_H_O2() != 0.0, ExcMessage("Derivatives at the moment are defined only using Thin fim (Henry's law) case for Oxygen molar fraction."));
                        
                        if ( gamma.at(iter->first)-1. >= 0. )
                            species_comp *= ( (gamma.at(iter->first)*pow(1./(ref_conc.at(iter->first)),gamma.at(iter->first))*(p_total/electrolyte->get_H_O2())) *pow(iter->second[j],gamma.at(iter->first)-1.) );
                        else
                            species_comp *= ( (gamma.at(iter->first)*pow(1./(ref_conc.at(iter->first)),gamma.at(iter->first))*(p_total/electrolyte->get_H_O2())) /pow(iter->second[j],1.-gamma.at(iter->first)) );
                    }
                    
                    else if ( derivative_flags[i] == hydrogen_molar_fraction && iter->first == hydrogen_concentration )
                    {
                        Assert( electrolyte->get_H_H2() != 0.0, ExcMessage("Derivatives at the moment are defined only using Thin fim (Henry's law) case for Hydrogen molar fraction."));
                        
                        if ( gamma.at(iter->first)-1. >= 0. )
                            species_comp *= ( (gamma.at(iter->first)*pow(1./(ref_conc.at(iter->first)),gamma.at(iter->first))*(p_total/electrolyte->get_H_H2())) *pow(iter->second[j],gamma.at(iter->first)-1.) );
                        else
                            species_comp *= ( (gamma.at(iter->first)*pow(1./(ref_conc.at(iter->first)),gamma.at(iter->first))*(p_total/electrolyte->get_H_H2())) /pow(iter->second[j],1.-gamma.at(iter->first)) );
                    }
                    
                    else if ( iter->first == derivative_flags[i] )
                    {                        
                        if ( gamma.at(iter->first)-1. >= 0. )
                            species_comp *= ( (gamma.at(iter->first)*pow(1./(ref_conc.at(iter->first)),gamma.at(iter->first))) *pow(iter->second[j],gamma.at(iter->first)-1.) );
                        else
                            species_comp *= ( (gamma.at(iter->first)*pow(1./(ref_conc.at(iter->first)),gamma.at(iter->first))) /pow(iter->second[j],1.-gamma.at(iter->first)) );
                    }
                    
                    else
                        species_comp *= ( pow(iter->second[j]/ref_conc.at(iter->first), gamma.at(iter->first)) );
                }

                dcurrent[j] = (catalyst->exchange_current_density(T[j]))*species_comp* ( exp( alpha_a*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*T[j]) ) 
                                                                                                            -
                                                                                         exp(-alpha_c*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*T[j]) ) );

                if (std::isinf(dcurrent[j]))
                    dcurrent[j] = std::numeric_limits<double>::max()*1e-150;
                else if (std::isnan(dcurrent[j]))
                    dcurrent[j] = 0.0;
            }
            
            dcoef_du[ derivative_flags[i] ] = dcurrent;
        }// x_O2 / x_H2 / Reactant concentrations

        else if (derivative_flags[i] == protonic_electrical_potential)
        {
            for (unsigned int j=0; j<n_quad; ++j)
            {
                species_comp = 1.;
                // Loop over each species in the reaction and compute its contribution
                for (std::map< VariableNames, SolutionVariable >::const_iterator iter=reactants_map.begin(); iter != reactants_map.end(); ++iter)
                {
                    if (iter->second[j] < 0.0)
                        species_comp *= 0.0;
                    else
                        species_comp *= ( pow(iter->second[j]/ref_conc.at(iter->first), gamma.at(iter->first)) );
                }
                              
                dcurrent[j] = (catalyst->exchange_current_density(T[j]))*species_comp*( exp( alpha_a*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*T[j]) )*( -alpha_a*F/(R*T[j]) ) 
                                                                                                        - 
                                                                                        exp(-alpha_c*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*T[j]) )*(  alpha_c*F/(R*T[j]) ) );

                if (std::isinf(dcurrent[j]))
                    dcurrent[j] = std::numeric_limits<double>::max()*1e-150;
                else if (std::isnan(dcurrent[j]))
                    dcurrent[j] = 0.0;
            }
            
            dcoef_du[protonic_electrical_potential] = dcurrent;
        }//phi_m
        
        else if (derivative_flags[i] == electronic_electrical_potential)
        {
            for (unsigned int j=0; j<n_quad; ++j)
            {
                species_comp = 1.;
                // Loop over each species in the reaction and compute its contribution
                for (std::map< VariableNames, SolutionVariable >::const_iterator iter=reactants_map.begin(); iter != reactants_map.end(); ++iter)
                {
                    if (iter->second[j] < 0.0)
                        species_comp *= 0.0;
                    else
                        species_comp *= ( pow(iter->second[j]/ref_conc.at(iter->first), gamma.at(iter->first)) );
                }

                dcurrent[j] = (catalyst->exchange_current_density(T[j]))*species_comp*( exp( alpha_a*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*T[j]) )*( alpha_a*F/(R*T[j]) ) 
                                                                                                        - 
                                                                                        exp(-alpha_c*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*T[j]) )*(-alpha_c*F/(R*T[j]) ) ); 
                
                if (std::isinf(dcurrent[j]))
                    dcurrent[j] = std::numeric_limits<double>::max()*1e-150;
                else if (std::isnan(dcurrent[j]))
                    dcurrent[j] = 0.0;	
            }
            
            dcoef_du[electronic_electrical_potential] = dcurrent;
        }//phi_s
        
        else if (derivative_flags[i] == temperature_of_REV)
        {
            for (unsigned int j=0; j<n_quad; ++j)
            {
                species_comp = 1.;
                // Loop over each species in the reaction and compute its contribution
                for (std::map< VariableNames, SolutionVariable >::const_iterator iter=reactants_map.begin(); iter != reactants_map.end(); ++iter)
                {
                    if (iter->second[j] < 0.0)
                        species_comp *= 0.0;
                    else
                        species_comp *= ( pow(iter->second[j]/ref_conc.at(iter->first), gamma.at(iter->first)) );
                }

                double first_term = (catalyst->derivative_exchange_current_density(T[j]))*(  exp( alpha_a*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*T[j]) ) 
                                                                                                                - 
                                                                                             exp((-alpha_c*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j]))))/(R*T[j])) );
                double intermed_c = ( alpha_c*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*pow(T[j],2.))) + (alpha_c*F*(catalyst->dvoltage_cell_th_dT(T[j]))/(R*T[j]));
                double intermed_a = (-alpha_a*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*pow(T[j],2.))) - (alpha_a*F*(catalyst->dvoltage_cell_th_dT(T[j]))/(R*T[j]));
                double second_term = (catalyst->exchange_current_density(T[j])) * ( exp( alpha_a*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*T[j]) )*intermed_a
                                                                                                        -
                                                                                    exp((-alpha_c*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j]))))/(R*T[j]))*intermed_c );
                
                dcurrent[j] = species_comp*(first_term + second_term);
                
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