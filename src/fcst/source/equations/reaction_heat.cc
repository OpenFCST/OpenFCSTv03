//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: reaction_heat.cc
//    - Description: Class for computing heat source terms due to electrochemical reaction
//    - Developers: Madhur Bhaiya
//
//---------------------------------------------------------------------------

#include "equations/reaction_heat.h"

namespace NAME = FuelCellShop::Equation;

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
NAME::ReactionHeat::ReactionHeat()
{
    //FcstUtilities::log << "->FuelCellShop::Equation::ReactionHeat" << std::endl;

    kinetics = NULL;
    factors_initialized = false;
}

//---------------------------------------------------------------------------
NAME::ReactionHeat::~ReactionHeat()
{}

//---------------------------------------------------------------------------
void
NAME::ReactionHeat::heat_source (std::vector<double>& heat,
                                 const std::vector<double>& current) const
{
    Assert( kinetics != NULL, ExcMessage("Kinetics object not initialized using set_kinetics method in ReactionHeat object"));
    Assert( phi_m.is_initialized() && phi_s.is_initialized() && T.is_initialized(), ExcMessage("Either phi_m/phi_s/T is not set before the ReactionHeat::heat_source.") );
    Assert( kinetics->get_reaction_name()==ORR || kinetics->get_reaction_name()==HOR, ExcMessage("Reaction should be either ORR or HOR only.") );
    Assert( factors_initialized, ExcMessage("initialize_factors method not called before the ReactionHeat::heat_source.") );
    
    heat.resize(phi_m.size());
    
    if (kinetics->get_reaction_name() == ORR)
    {
        for (unsigned int i=0; i<heat.size(); ++i)
        {
            heat[i] = current.at(i) * ( (phi_s[i] - phi_m[i] - kinetics->get_cat()->voltage_cell_th(T[i]))*factor_irrev_ORR 
                                            + 
                                        T[i]*entropy_rxn(T[i])*factor_rev_ORR 
                                            + 
                                        FuelCellShop::Material::LiquidWater::latentVap_heat(T[i])*factor_vap_ORR );
            
            if (std::isinf(heat[i]))
                heat[i] = std::numeric_limits<double>::max()*1e-150;
            else if (std::isnan(heat[i]))
                heat[i] = 0.0;
        }
    }
    
    else if (kinetics->get_reaction_name() == HOR)
    {        
        for (unsigned int i=0; i<heat.size(); ++i)
        {
            heat[i] = current.at(i) * ( (phi_s[i] - phi_m[i] - kinetics->get_cat()->voltage_cell_th(T[i]))*factor_irrev_HOR 
                                            + 
                                        T[i]*entropy_rxn(T[i])*factor_rev_HOR );
            
            if (std::isinf(heat[i]))
                heat[i] = std::numeric_limits<double>::max()*1e-150;
            else if (std::isnan(heat[i]))
                heat[i] = 0.0;
        }
    }
}

//---------------------------------------------------------------------------
void
NAME::ReactionHeat::derivative_heat_source (std::map< VariableNames, std::vector<double> >& heat_derived,
                                            const std::map< VariableNames, std::vector<double> >& current_derived,
                                            const std::vector<double>& current) const
{
    Assert( kinetics != NULL, ExcMessage("Kinetics object not initialized using set_kinetics method in ReactionHeat object"));
    Assert( phi_m.is_initialized() && phi_s.is_initialized() && T.is_initialized(), ExcMessage("Either phi_m/phi_s/T is not set before the ReactionHeat::derivative_heat_source.") );
    Assert( kinetics->get_reaction_name()==ORR || kinetics->get_reaction_name()==HOR, ExcMessage("Reaction should be either ORR or HOR only.") );
    Assert( derivative_flags.size() != 0, ExcMessage("Derivative flags are not set using set_derivative_flags method before the ReactionHeat::derivative_heat_source.") );
    Assert( factors_initialized, ExcMessage("initialize_factors method not called before the ReactionHeat::derivative_heat_source.") );
    
    if (kinetics->get_reaction_name() == ORR)
    {   
        // Loop over the flags
        for (unsigned int i=0; i<derivative_flags.size(); ++i)
        {
            Assert( current_derived.find(derivative_flags[i])!=current_derived.end(), ExcMessage("Corresponding derivative does not exist in the current derivative map from kinetics.") );
            
            std::vector<double> dheat(phi_m.size(), 0.);         // Defaulting all derivatives to 0.0
            
            if (derivative_flags[i] == electronic_electrical_potential)
            {
                for (unsigned int j=0; j<dheat.size(); ++j)
                {
                    dheat[j] = ( current_derived.at(electronic_electrical_potential)[j] * ( (phi_s[j] - phi_m[j] - kinetics->get_cat()->voltage_cell_th(T[j]))*factor_irrev_ORR 
                                        +
                                T[j]*entropy_rxn(T[j])*factor_rev_ORR 
                                        + 
                                FuelCellShop::Material::LiquidWater::latentVap_heat(T[j])*factor_vap_ORR ) )    +  ( current[j] * factor_irrev_ORR );
                    
                    if (std::isinf(dheat[j]))
                        dheat[j] = std::numeric_limits<double>::max()*1e-150;
                    else if (std::isnan(dheat[j]))
                        dheat[j] = 0.0;
                }
                
                heat_derived[electronic_electrical_potential] = dheat;
            }//phi_s
            
            else if (derivative_flags[i] == protonic_electrical_potential)
            {
                for (unsigned int j=0; j<dheat.size(); ++j)
                {
                    dheat[j] = ( current_derived.at(protonic_electrical_potential)[j] * ( (phi_s[j] - phi_m[j] - kinetics->get_cat()->voltage_cell_th(T[j]))*factor_irrev_ORR 
                                    + 
                                T[j]*entropy_rxn(T[j])*factor_rev_ORR 
                                    + 
                                FuelCellShop::Material::LiquidWater::latentVap_heat(T[j])*factor_vap_ORR ) )    -  ( current[j] * factor_irrev_ORR );
                    
                    if (std::isinf(dheat[j]))
                        dheat[j] = std::numeric_limits<double>::max()*1e-150;
                    else if (std::isnan(dheat[j]))
                        dheat[j] = 0.0;
                }
                
                heat_derived[protonic_electrical_potential] = dheat;
            }//phi_m
            
            else if (derivative_flags[i] == temperature_of_REV)
            {
                for (unsigned int j=0; j<dheat.size(); ++j)
                {
                    dheat[j] = ( current_derived.at(temperature_of_REV)[j] * ( (phi_s[j] - phi_m[j] - kinetics->get_cat()->voltage_cell_th(T[j]))*factor_irrev_ORR 
                                    + 
                                T[j]*entropy_rxn(T[j])*factor_rev_ORR 
                                    + 
                                FuelCellShop::Material::LiquidWater::latentVap_heat(T[j])*factor_vap_ORR ) ) 
                                                + 
                               ( current[j] * ( (-1.)*kinetics->get_cat()->dvoltage_cell_th_dT(T[j])*factor_irrev_ORR 
                                                        +
                                                 entropy_rxn(T[j])*factor_rev_ORR 
                                                        +
                                                 T[j]*deriv_entropy_rxn(T[j])*factor_rev_ORR
                                                        +
                                                 FuelCellShop::Material::LiquidWater::deriv_latentVap_heat(T[j])*factor_vap_ORR ) );
                    
                    if (std::isinf(dheat[j]))
                        dheat[j] = std::numeric_limits<double>::max()*1e-150;
                    else if (std::isnan(dheat[j]))
                        dheat[j] = 0.0;
                }
                
                heat_derived[temperature_of_REV] = dheat;
            }//T
            
            else
            {
                for (unsigned int j=0; j<dheat.size(); ++j)
                {
                    dheat[j] = ( current_derived.at(derivative_flags[i])[j] * ( (phi_s[j] - phi_m[j] - kinetics->get_cat()->voltage_cell_th(T[j]))*factor_irrev_ORR 
                                    + 
                                T[j]*entropy_rxn(T[j])*factor_rev_ORR 
                                    + 
                                FuelCellShop::Material::LiquidWater::latentVap_heat(T[j])*factor_vap_ORR ) );
                    
                    if (std::isinf(dheat[j]))
                        dheat[j] = std::numeric_limits<double>::max()*1e-150;
                    else if (std::isnan(dheat[j]))
                        dheat[j] = 0.0;
                }
                
                heat_derived[ derivative_flags[i] ] = dheat;
            }// Other variables including x_O2 etc. 
        }//end of derivative_flags loop
    }
            
    else if (kinetics->get_reaction_name() == HOR)
    {
        // Loop over the flags
        for (unsigned int i=0; i<derivative_flags.size(); ++i)
        {
            Assert( current_derived.find(derivative_flags[i])!=current_derived.end(), ExcMessage("Corresponding derivative does not exist in the current derivative map from kinetics.") );
            
            std::vector<double> dheat(phi_m.size(), 0.);         // Defaulting all derivatives to 0.0
                       
            if (derivative_flags[i] == electronic_electrical_potential)
            {
                for (unsigned int j=0; j<dheat.size(); ++j)
                {
                    dheat[j] = ( current_derived.at(electronic_electrical_potential)[j] * ( (phi_s[j] - phi_m[j] - kinetics->get_cat()->voltage_cell_th(T[j]))*factor_irrev_HOR 
                                        +
                                T[j]*entropy_rxn(T[j])*factor_rev_HOR ) )    +  ( current[j] * factor_irrev_HOR );
                    
                    if (std::isinf(dheat[j]))
                        dheat[j] = std::numeric_limits<double>::max()*1e-150;
                    else if (std::isnan(dheat[j]))
                        dheat[j] = 0.0;
                }
                
                heat_derived[electronic_electrical_potential] = dheat;
            }//phi_s
            
            else if (derivative_flags[i] == protonic_electrical_potential)
            {
                for (unsigned int j=0; j<dheat.size(); ++j)
                {
                    dheat[j] = ( current_derived.at(protonic_electrical_potential)[j] * ( (phi_s[j] - phi_m[j] - kinetics->get_cat()->voltage_cell_th(T[j]))*factor_irrev_HOR 
                                    + 
                                T[j]*entropy_rxn(T[j])*factor_rev_HOR ) )    -  ( current[j] * factor_irrev_HOR );
                    
                    if (std::isinf(dheat[j]))
                        dheat[j] = std::numeric_limits<double>::max()*1e-150;
                    else if (std::isnan(dheat[j]))
                        dheat[j] = 0.0;
                }
                
                heat_derived[protonic_electrical_potential] = dheat;
            }//phi_m
            
            else if (derivative_flags[i] == temperature_of_REV)
            {
                for (unsigned int j=0; j<dheat.size(); ++j)
                {
                    dheat[j] = ( current_derived.at(temperature_of_REV)[j] * ( (phi_s[j] - phi_m[j] - kinetics->get_cat()->voltage_cell_th(T[j]))*factor_irrev_HOR 
                                    + 
                                T[j]*entropy_rxn(T[j])*factor_rev_HOR ) ) 
                                                + 
                               ( current[j] * ( (-1.)*kinetics->get_cat()->dvoltage_cell_th_dT(T[j])*factor_irrev_HOR 
                                                        +
                                                 entropy_rxn(T[j])*factor_rev_HOR 
                                                        +
                                                 T[j]*deriv_entropy_rxn(T[j])*factor_rev_HOR ) );
                    
                    if (std::isinf(dheat[j]))
                        dheat[j] = std::numeric_limits<double>::max()*1e-150;
                    else if (std::isnan(dheat[j]))
                        dheat[j] = 0.0;	
                }
                
                heat_derived[temperature_of_REV] = dheat;
            }//T
            
            else
            {
                for (unsigned int j=0; j<dheat.size(); ++j)
                {
                    dheat[j] = ( current_derived.at(derivative_flags[i])[j] * ( (phi_s[j] - phi_m[j] - kinetics->get_cat()->voltage_cell_th(T[j]))*factor_irrev_HOR 
                                    + 
                                T[j]*entropy_rxn(T[j])*factor_rev_HOR ) );
                    
                    if (std::isinf(dheat[j]))
                        dheat[j] = std::numeric_limits<double>::max()*1e-150;
                    else if (std::isnan(dheat[j]))
                        dheat[j] = 0.0;   
                }
                
                heat_derived[ derivative_flags[i] ] = dheat;
            }// Other variables including x_H2 etc.          
        }//end of derivative_flags loop
    } // end of ORR/HOR if-else
    
}
