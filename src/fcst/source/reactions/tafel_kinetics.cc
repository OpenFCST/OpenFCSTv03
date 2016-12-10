//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: tafel_kinetics.cc
//    - Description: Tafel Kinetics model class
//    - Developers: M. Moore, M. Bhaiya, M. Secanell and V. Zingan
//    - $Id: tafel_kinetics.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include "reactions/tafel_kinetics.h"

namespace NAME = FuelCellShop::Kinetics;

const std::string NAME::TafelKinetics::concrete_name ("TafelKinetics");

NAME::TafelKinetics const* NAME::TafelKinetics::PROTOTYPE = new NAME::TafelKinetics(true);

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

NAME::TafelKinetics::TafelKinetics(const bool construct_replica)
:
BaseKinetics ()
{
    if (construct_replica)
        this->get_mapFactory()->insert(std::pair<std::string, BaseKinetics* > (this->concrete_name, this) );
}

//---------------------------------------------------------------------------
NAME::TafelKinetics::TafelKinetics()
:
BaseKinetics ()
{
    //FcstUtilities::log << "Tafel kinetics." << std::endl;
}

//---------------------------------------------------------------------------
NAME::TafelKinetics::~TafelKinetics()
{}

//---------------------------------------------------------------------------
void
NAME::TafelKinetics::declare_parameters(ParameterHandler& param) const
{}

//---------------------------------------------------------------------------
void
NAME::TafelKinetics::initialize(ParameterHandler& param)
{}

//---------------------------------------------------------------------------
void
NAME::TafelKinetics::current_density (std::vector<double> &coef)
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
        for(std::map< VariableNames, SolutionVariable >::const_iterator iter=reactants_map.begin(); iter != reactants_map.end(); ++iter)
        {
            species_comp *= negative_concentration_correction(iter->second[i], gamma.at(iter->first), false)
                            * pow( std::fabs(iter->second[i]) / ref_conc.at(iter->first), gamma.at(iter->first) );
        }

        coef[i] = (catalyst->exchange_current_density(T[i])) * species_comp * exp( -alpha_c*F*(phi_s[i]-phi_m[i]-(catalyst->voltage_cell_th(T[i])))/(R*T[i]) );
    }
}

//---------------------------------------------------------------------------
void
NAME::TafelKinetics::derivative_current (std::map< VariableNames, std::vector<double> >& dcoef_du)
{
     Assert( derivative_flags.size() != 0, ExcMessage("Derivative flags are not set using set_derivative_flags method before the TafelKinetics::derivative_current.") );

    if (!kin_param_initialized)
        init_kin_param();

    n_quad = phi_m.size();
    double species_comp;

    // Loop over the flags
    for (unsigned int i=0; i<derivative_flags.size(); ++i)
    {
        std::vector<double> dcurrent(n_quad, 0.);

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
                    if ( derivative_flags[i] == oxygen_molar_fraction && iter->first == oxygen_concentration )
                    {
                        Assert( electrolyte->get_H_O2() != 0.0, ExcMessage("Derivatives at the moment are defined only using Thin film (Henry's law) case for Oxygen molar fraction."));
                        
                        species_comp *= negative_concentration_correction(iter->second[j], gamma.at(iter->first)) * 
                        gamma.at(iter->first) * pow(p_total/(ref_conc.at(iter->first)*electrolyte->get_H_O2()),gamma.at(iter->first))*pow(std::fabs(iter->second[j]),gamma.at(iter->first)-1.);
                    }

                    else if ( derivative_flags[i] == hydrogen_molar_fraction && iter->first == hydrogen_concentration )
                    {
                        Assert( electrolyte->get_H_H2() != 0.0, ExcMessage("Derivatives at the moment are defined only using Thin film (Henry's law) case for Hydrogen molar fraction."));
                        
                        species_comp *= negative_concentration_correction(iter->second[j], gamma.at(iter->first)) * 
                        gamma.at(iter->first) * pow(p_total/(ref_conc.at(iter->first)*electrolyte->get_H_H2()),gamma.at(iter->first))*pow(std::fabs(iter->second[j]),gamma.at(iter->first)-1.);
                    }
                    else if ( derivative_flags[i] == iter->first )
                    {                        
                        species_comp *= negative_concentration_correction(iter->second[j], gamma.at(iter->first)) * 
                                        gamma.at(iter->first) * pow( std::fabs(iter->second[j])/(ref_conc.at(iter->first)) , gamma.at(iter->first) - 1.0 )  * (1.0/ref_conc.at(iter->first));                            
                    }
                    // For most cases I want the derivative with respect to the reactant concentration. Then:
                    else
                    {
                        species_comp *= negative_concentration_correction(iter->second[j], gamma.at(iter->first))
                                        * pow( std::fabs(iter->second[j]) / ref_conc.at(iter->first), gamma.at(iter->first) );
                    }
                }

                dcurrent[j] = (catalyst->exchange_current_density(T[j])) * species_comp * exp( -alpha_c*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*T[j]) );
            }

            dcoef_du[ derivative_flags[i] ] = dcurrent;

        } // w.r.t molar fractions and reactant concentrations

        else if (derivative_flags[i] == electronic_electrical_potential)
        {
            for (unsigned int j=0; j<n_quad; ++j)
            {
                species_comp = 1.;
                // Loop over each species in the reaction and compute its contribution
                for (std::map< VariableNames, SolutionVariable >::const_iterator iter=reactants_map.begin(); iter != reactants_map.end(); ++iter)
                {
                    species_comp *= negative_concentration_correction(iter->second[j], gamma.at(iter->first)) 
                            * pow( std::fabs(iter->second[j]) / ref_conc.at(iter->first), gamma.at(iter->first) );
                }

                dcurrent[j] = (catalyst->exchange_current_density(T[j])) * species_comp * exp(-alpha_c*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*T[j]) ) * (-alpha_c*F/(R*T[j]));
            }

            dcoef_du[electronic_electrical_potential] = dcurrent;

        } // phi_s

        else if (derivative_flags[i] == protonic_electrical_potential)
        {
            for (unsigned int j=0; j<n_quad; ++j)
            {
                species_comp = 1.;
                // Loop over each species in the reaction and compute its contribution
                for (std::map< VariableNames, SolutionVariable >::const_iterator iter=reactants_map.begin(); iter != reactants_map.end(); ++iter)
                {
                    species_comp *= negative_concentration_correction(iter->second[j], gamma.at(iter->first)) 
                            * pow( std::fabs(iter->second[j]) / ref_conc.at(iter->first), gamma.at(iter->first) );
                }

                dcurrent[j] = (catalyst->exchange_current_density(T[j])) * species_comp * exp(-alpha_c*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*T[j]) ) * (alpha_c*F/(R*T[j]));
            }

            dcoef_du[protonic_electrical_potential] = dcurrent;

        } // phi_m

        else if (derivative_flags[i] == temperature_of_REV)
        {
            for (unsigned int j=0; j<n_quad; ++j)
            {
                species_comp = 1.;
                // Loop over each species in the reaction and compute its contribution
                for (std::map< VariableNames, SolutionVariable >::const_iterator iter=reactants_map.begin(); iter != reactants_map.end(); ++iter)
                {
                    species_comp *= negative_concentration_correction(iter->second[j], gamma.at(iter->first)) 
                           * pow( std::fabs(iter->second[j]) / ref_conc.at(iter->first), gamma.at(iter->first) );
                }

                double first_term = (catalyst->derivative_exchange_current_density(T[j]))*exp( -alpha_c*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*T[j]) );
                double intermed = (alpha_c*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*pow(T[j],2.))) + (alpha_c*F*(catalyst->dvoltage_cell_th_dT(T[j]))/(R*T[j]));
                double second_term = (catalyst->exchange_current_density(T[j])) * (exp( -alpha_c*F*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])))/(R*T[j]) )) * intermed;
            }

            dcoef_du[temperature_of_REV] = dcurrent;

        } // temperature

        else
        {
            FcstUtilities::log << "Wrong w.r.t derivative in Tafel kinetics - flags" << std::endl;
            AssertThrow( false , ExcNotImplemented() );
        }
    }
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//----------Based on old format --needs to be changed..Madhur ---------------
//----------Using numerical differentiation to test Tafel kinetics ----------
//
//   //size checking for vectors
//   Assert (phi_m.size() == dcoef_du[0].size(), ExcDimensionMismatch(phi_m.size(), dcoef_du[0].size()));
//
//   //Loop over the flags
//   for (unsigned int i=0; i<this->derivative_flags.size(); ++i)
//   {
//     if (this->derivative_flags[i] == "Protonic electrical potential")
//     {
//       //Compute the current using the current solution.
//       std::vector<double> coef1(std::vector<double>(dcoef_du[i].size(), 0.0));
//       current_density(coef1);
//
//       double delta = pow(10, -7);
//
//       //Modify solution
//       for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//       {
//              phi_m[j] += delta;
//       }
//
//       //Compute the new current
//       std::vector<double> coef2(std::vector<double>(dcoef_du[i].size(), 0.0));
//       current_density(coef2);
//
//       //Compute the coefficient
//       for (unsigned int j=0; j<dcoef_du[0].size(); ++j)
//       {
//      dcoef_du[i][j] = (coef2[j] - coef1[j])/delta;
//      FcstUtilities::log << "Numerical phi_m derivative: " << dcoef_du[i][j] << std::endl;
//       }
//       for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//       {
//      phi_m[j] -= delta;
//       }
//     }//T
//
//     else
//     {
//
//     }//everything else
//
//    }//end of derivative_flags loop
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------