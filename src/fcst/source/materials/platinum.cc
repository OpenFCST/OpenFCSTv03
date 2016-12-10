//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: platinum.cc
//    - Description: Class representing Platinum material class
//    - Developers: Peter Dobson(2011), Madhur Bhaiya(2012-13) and M. Secanell(2013)
//    - $Id: platinum.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <materials/platinum.h>

namespace NAME = FuelCellShop::Material;

const std::string NAME::Platinum::concrete_name ("Platinum");

NAME::Platinum const* NAME::Platinum::PROTOTYPE = new NAME::Platinum(true);

//---------------------------------------------------------------------------
NAME::Platinum::Platinum(std::string name)
:
CatalystBase (name)
{
    // Default Values
    density = 21.5;
    alpha_a_ORR = 1.0;
    alpha_c_ORR = 1.0;
    i_0_ref_ORR = 2.47e-2;
    c_O2_ref_ORR = 0.725e-5;
    c_H_ref_ORR = 1.818e-3;
    gamma_O2_ORR = 1.0;
    gamma_H_ORR = 1.0;
    method_kinetics_ORR = "Given";
    given_OCV_ORR = 1.229;

    alpha_a_HOR = 0.5;
    alpha_c_HOR = 0.5;
    i_0_ref_HOR = 1e6;
    c_H2_ref_HOR = 5.64e-5;
    gamma_H2_HOR = 0.25;
    given_OCV_HOR = 0.;
}

//---------------------------------------------------------------------------
NAME::Platinum::Platinum(const bool create_replica)
:
CatalystBase ()
{
    if (create_replica)
        this->get_mapFactory()->insert(std::pair<std::string, CatalystBase* > (this->concrete_name, this) );           
}
//---------------------------------------------------------------------------
NAME::Platinum::~Platinum()
{
}

//---------------------------------------------------------------------------
void
NAME::Platinum::declare_parameters(ParameterHandler &param) const
{
    
    param.enter_subsection("Materials");
    {
        param.enter_subsection(NAME::Platinum::concrete_name);
        {
            param.declare_entry ("Density [g/cm^3]",
                                 "21.5", 
                                 Patterns::Double());
            
            //############ Oxygen Reduction Reaction Parameters ############
            param.declare_entry ("Method for kinetics parameters (ORR)",
                                 "Parthasarathy",
                                 Patterns::Selection("Given|Parthasarathy|Parthasarathy_hcd|Double_trap|Neyerlin"),
                                 "Method used to compute specific kinetics parameters such as OCV, reference exchange current density and reference concenetrations.");
                        
            param.declare_entry ("Anodic transfer coefficient (ORR)",
                                 "0.5", 
                                 Patterns::Double());
            
            param.declare_entry ("Cathodic transfer coefficient (ORR)",
                                 "1.0", 
                                 Patterns::Double());
                                             
            param.declare_entry ("Given Open Cell Voltage (ORR) [V]",
                                 "1.229", 
                                 Patterns::Double());
            
            param.declare_entry ("Given Open Cell Voltage (HOR) [V]",
                                 "0.0", 
                                 Patterns::Double());
            
            param.declare_entry ("Reference exchange current density (ORR) [uA/cm2]",
                                 "2.707e-2", 
                                 Patterns::Double());
            
            param.declare_entry ("Reference oxygen concentration (ORR)",
                                 "0.725e-5", 
                                 Patterns::Double());
            
            param.declare_entry ("Reference proton concentration (ORR)",
                                 "1.818e-3", 
                                 Patterns::Double());
            
            param.declare_entry ("Oxygen reaction order (ORR)",
                                 "1.0", 
                                 Patterns::Double());
            
            param.declare_entry ("Proton reaction order (ORR)",
                                 "1.0", 
                                 Patterns::Double());
            
            //############ Hydrogen Oxidation Reaction Parameters ############
            param.declare_entry ("Anodic transfer coefficient (HOR)",
                                 "0.5", 
                                 Patterns::Double());
            
            param.declare_entry ("Cathodic transfer coefficient (HOR)",
                                 "0.5", 
                                 Patterns::Double());
            
            param.declare_entry ("Reference exchange current density (HOR) [uA/cm2]",
                                 "1e6", 
                                 Patterns::Double());
            
            param.declare_entry ("Reference hydrogen concentration (HOR)",
                                 "5.64e-5", 
                                 Patterns::Double());
            
            param.declare_entry ("Hydrogen reaction order (HOR)",
                                 "0.25", 
                                 Patterns::Double());
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
void
NAME::Platinum::initialize (ParameterHandler& param)
{
    
    param.enter_subsection("Materials");
    {
        param.enter_subsection(NAME::Platinum::concrete_name);
        {
            density = param.get_double("Density [g/cm^3]");
            //############ Oxygen Reduction Reaction Parameters ############            
            method_kinetics_ORR = param.get("Method for kinetics parameters (ORR)");
            alpha_a_ORR = param.get_double("Anodic transfer coefficient (ORR)");
            alpha_c_ORR = param.get_double("Cathodic transfer coefficient (ORR)");
            // overwrite if method given:
            if ((method_kinetics_ORR == "Parthasarathy") ||  (method_kinetics_ORR == "Neyerlin"))
                alpha_c_ORR = 1.0;
            else if (method_kinetics_ORR == "Parthasarathy_hcd")
                alpha_c_ORR = 0.5;
           
            given_OCV_ORR = param.get_double("Given Open Cell Voltage (ORR) [V]");
            i_0_ref_ORR = param.get_double("Reference exchange current density (ORR) [uA/cm2]");
            // overwrite if method given, however since it is T dependent, overwrite done later.
            c_O2_ref_ORR = param.get_double("Reference oxygen concentration (ORR)");
            // overwrite if method given:
            if ( (method_kinetics_ORR == "Parthasarathy") || (method_kinetics_ORR == "Parthasarathy_hcd") ){
                const double H_O2 = 3.1664e10;      // Henry's constant
                c_O2_ref_ORR = (5.0*101325.0)/H_O2;
            }
            else if (method_kinetics_ORR == "Double_trap"){
                const double H_O2 = 3.1664e10;      // Henry's constant
                //The 1.05atm curve in Parthasarathy's pressure dependent paper was chosen as the reference condition
                c_O2_ref_ORR = (1.05*101325.0)/H_O2;
            }
            else if (method_kinetics_ORR == "Neyerlin"){
                const double H_O2 = 3.1664e10;      // Henry's constant
                //Neyerlin assumes 1atm to be the reference condition
                c_O2_ref_ORR = (1.0*101325.0)/H_O2;
            }
            c_H_ref_ORR = param.get_double("Reference proton concentration (ORR)");
            gamma_O2_ORR = param.get_double("Oxygen reaction order (ORR)");
            // overwrite if method given:
            if ((method_kinetics_ORR == "Parthasarathy") ||  (method_kinetics_ORR == "Parthasarathy_hcd"))
                gamma_O2_ORR = 1.0;
            else if (method_kinetics_ORR == "Neyerlin")
                gamma_O2_ORR = 0.54;
            
            gamma_H_ORR = param.get_double("Proton reaction order (ORR)");
            
            
            
            //############ Hydrogen Oxidation Reaction Parameters ############
            given_OCV_HOR = param.get_double("Given Open Cell Voltage (HOR) [V]");
            alpha_a_HOR = param.get_double("Anodic transfer coefficient (HOR)");
            alpha_c_HOR = param.get_double("Cathodic transfer coefficient (HOR)");
            i_0_ref_HOR = param.get_double("Reference exchange current density (HOR) [uA/cm2]");
            c_H2_ref_HOR = param.get_double("Reference hydrogen concentration (HOR)");
            gamma_H2_HOR = param.get_double("Hydrogen reaction order (HOR)");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
void
NAME::Platinum::alpha_anodic(double& alpha_a) const
{
    if (name_reaction_kinetics == ORR)
        alpha_a = alpha_a_ORR;
    
    else if (name_reaction_kinetics == HOR)
        alpha_a = alpha_a_HOR;
    
    else
        Assert(false, ExcNotImplemented());
}

//---------------------------------------------------------------------------
void
NAME::Platinum::alpha_cathodic(double& alpha_c) const
{
    if (name_reaction_kinetics == ORR)
        alpha_c = alpha_c_ORR;
    else if (name_reaction_kinetics == HOR)
        alpha_c = alpha_c_HOR;
    
    else 
        Assert(false, ExcNotImplemented());
}

//---------------------------------------------------------------------------
double
NAME::Platinum::exchange_current_density(const double& Temp) const
{
    if (name_reaction_kinetics == ORR)
    {
        if (method_kinetics_ORR == "Given")
        {
            return (1.e-6*i_0_ref_ORR);
        }

        else if (method_kinetics_ORR == "Parthasarathy")
        {
            //double A = 1.69278e-10; //<- NOTE: <- From Parthasarathy et al. Temperature dependance (@ 30C, 5 atm). This values conflict with data from Pressure dependance paper.
            double A = 2.48e-8; //<- NOTE: <- i0 average between the two papers from Parthasarthy at 50C and 5 atm
            double E_a = 80987.61;
            return (A * exp( (-E_a/Constants::R()) * ( 1.0/Temp - 1.0/323.15 ) ));
        }
        else if (method_kinetics_ORR == "Parthasarathy_hcd")
        {
            //double A = 2.83792e-7; //<- NOTE: <- From Parthasarathy et al. Temperature dependance (@ 30C, 5 atm). This values conflict with data from Pressure dependance paper. 
            double A = 3.08e-6; //<- NOTE: <- i0 average between the two papers from Parthasarthy at 50C and 5 atm
            double E_a = 28920.95;
            return (A * exp( (-E_a/Constants::R()) * ( 1.0/Temp - 1.0/323.15 ) ));
        }
        else if (method_kinetics_ORR == "Neyerlin")
        {
            double i_ref_standard = 2.47e-8;
            double E_rev = 67000.;
            double exponential = -E_rev/(Constants::R())*( 1.0/Temp - 1.0/353.0);
            return (i_ref_standard*exp(exponential));
        }
    }
    
    else if (name_reaction_kinetics == HOR)
    {
        return (1.e-6*i_0_ref_HOR);
    }
    
    else 
        Assert(false, ExcNotImplemented());
}

//---------------------------------------------------------------------------
double
NAME::Platinum::derivative_exchange_current_density(const double& Temp) const
{
    if (name_reaction_kinetics == ORR)
    {
        if (method_kinetics_ORR == "Given")
        {
            return 0.0;
        }
        else if (method_kinetics_ORR == "Parthasarathy")
        {
            double A = 2.48e-8; //<- NOTE: <- To match Parthasarathy's data this value should be x by 5.2 since this is the roughness factor of his electrode (M. Secanell, 2015)
            double E_a = 80987.61;
            return (A * exp( (-E_a/Constants::R()) * ( 1.0/Temp - 1.0/323.15 ) ) * (E_a/(Constants::R() * Temp * Temp)));
        }
        else if (method_kinetics_ORR == "Parthasarathy_hcd")
        {
            double A = 3.08e-6; //<- NOTE: <- To match Parthasarathy's data this value should be x by 5.2 since this is the roughness factor of his electrode (M. Secanell, 2015)
            double E_a = 28920.95;
            return (A * exp( (-E_a/Constants::R()) * ( 1.0/Temp - 1.0/323.15 ) ) * (E_a/(Constants::R() * Temp * Temp)));
        }
        else if (method_kinetics_ORR == "Neyerlin")
        {
            double i_ref_standard = 2.47e-8;
            double E_rev = 67000.;
            double exponential = (-E_rev/Constants::R()) * ( 1.0/Temp - 1.0/353.0 );
            return (i_ref_standard*exp(exponential)*(E_rev/(Constants::R() * Temp * Temp)) );
        }
    }
    
    else if (name_reaction_kinetics == HOR)
    {
        return 0.0;
    }
    
    else 
        Assert(false, ExcNotImplemented());
}

//---------------------------------------------------------------------------
double
NAME::Platinum::voltage_cell_th(const double& Temp) const
{
    if (name_reaction_kinetics == ORR)
    {
        if (method_kinetics_ORR == "Given")
        {
            return given_OCV_ORR;
        }
        else if ( (method_kinetics_ORR == "Parthasarathy") || (method_kinetics_ORR == "Parthasarathy_hcd") )
        {
            double E_T = ((70650.0+(8.0*Temp*std::log(Temp))-(92.84*Temp))*(4.184/(2.0*Constants::F())))+(std::log(std::pow(5.0,0.5))*(Constants::R()*Temp/(2.0*Constants::F())));
            return E_T;
        }
        else if (method_kinetics_ORR == "Double_trap")
        {
            double E_T = ((70650.0+(8.0*Temp*std::log(Temp))-(92.84*Temp))*(4.184/(2.0*Constants::F())))+(std::log(std::pow(1.05,0.5))*(Constants::R()*Temp/(2.0*Constants::F())));
            return E_T;
        }
        else if (method_kinetics_ORR == "Neyerlin")
        {
            //Note that the value in the log10() function should be p_o2/p_o2^ref where p_o2^ref is the reference 
            //pressure of 1 atm and p_o2 operating pressure. As the platinum class does not have access to the 
            // operating conditions class this has been set to 1.5atm, as most simulations are run at 1-2atm, leading
            //to errors on the order of a mV.
            double E_T = 1.23 -0.9e-3*(Temp - 298) + ((2.303*Constants::R()*Temp)/(4.0*Constants::F()))*log10(1.5);
            return E_T;
        }
    }
    else if (name_reaction_kinetics == HOR)
    {
        return given_OCV_HOR;
    }
    
    else 
        Assert (false, ExcNotImplemented ());
}

//---------------------------------------------------------------------------
double
NAME::Platinum::dvoltage_cell_th_dT(const double& Temp) const
{
    if (name_reaction_kinetics == ORR)
    {
        if (method_kinetics_ORR == "Given")
        {
            return 0.0;
        }
        else if ( (method_kinetics_ORR == "Parthasarathy") || (method_kinetics_ORR == "Parthasarathy_hcd") )
        {
            double dE_T = (((8.0*std::log(Temp)) + 8.0 - 92.84)*(4.184/(2.0*Constants::F()))) + (std::log(std::pow(5.0,0.5))*(Constants::R()/(2.0*Constants::F())));			
            return dE_T;
        }
        else if (method_kinetics_ORR == "Double_trap")
        {
            double dE_T = (((8.0*std::log(Temp)) + 8.0 - 92.84)*(4.184/(2.0*Constants::F()))) + (std::log(std::pow(1.05,0.5))*(Constants::R()/(2.0*Constants::F())));
            return dE_T;
        }
        else if (method_kinetics_ORR == "Neyerlin")
        {
            double dE_T = -0.9e-3 + ((2.303*Constants::R())/(4.0*Constants::F()))*log10(1.5);
            return dE_T;
        }
    }
    
    else if (name_reaction_kinetics == HOR)
    {
        return 0.0;
    }
    
    else 
        Assert (false, ExcNotImplemented ());
}

//---------------------------------------------------------------------------
void
NAME::Platinum::reference_concentration(const std::vector<VariableNames>& names,
                                        std::map<VariableNames, double>& cref_map) const
{    
    if (name_reaction_kinetics == ORR)
    {
        for (unsigned int i=0; i<names.size();++i)
        {
            if (names[i] == oxygen_concentration)
                cref_map[oxygen_concentration] = c_O2_ref_ORR;
            
            else if (names[i] == proton_concentration)
                cref_map[proton_concentration] = c_H_ref_ORR;
            else 
               Assert (false, ExcNotImplemented ()); 
        }
    }
    
    else if (name_reaction_kinetics == HOR)
    {
        for (unsigned int i=0; i<names.size();++i)
        {
            if (names[i] == hydrogen_concentration)
            {
                cref_map[hydrogen_concentration] = c_H2_ref_HOR;
            }
            
            else 
                Assert (false, ExcNotImplemented ()); 
        }
    }
}

//---------------------------------------------------------------------------
void
NAME::Platinum::reaction_order(const std::vector<VariableNames>& names,
                               std::map<VariableNames, double>& gamma_map) const
{    
    if (name_reaction_kinetics == ORR)
    {
        for (unsigned int i=0; i<names.size();++i)
        {
            if (names[i] == oxygen_concentration)
                gamma_map[oxygen_concentration] = gamma_O2_ORR;
            else if (names[i] == proton_concentration)
                gamma_map[proton_concentration] = gamma_H_ORR;            
            else
                gamma_map[ names[i] ] = 0.0;
        }
    }
    
    else if (name_reaction_kinetics == HOR)
    {
        for (unsigned int i=0; i<names.size();++i)
        {
            if (names[i] == hydrogen_concentration)
                gamma_map[hydrogen_concentration] = gamma_H2_HOR;
            
            else
                gamma_map[ names[i] ] = 0.0;
        }
    }
}

//---------------------------------------------------------------------------
