//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2010-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: double_trap_kinetics.h
//    - Description: DoubleTrap Kinetics implements the double trap model proposed by
//      Wang et al. and re-developed by Moore et al.
//    - Developers: M. Moore, M. Bhaiya and M. Secanell
//    - $Id: double_trap_kinetics.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include "reactions/double_trap_kinetics.h"

namespace NAME = FuelCellShop::Kinetics;

const std::string NAME::DoubleTrapKinetics::concrete_name ("DoubleTrapKinetics");

NAME::DoubleTrapKinetics const* NAME::DoubleTrapKinetics::PROTOTYPE = new NAME::DoubleTrapKinetics(true);

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
NAME::DoubleTrapKinetics::DoubleTrapKinetics(const bool construct_replica)
:
BaseKinetics ()
{
    if (construct_replica)
        this->get_mapFactory()->insert(std::pair<std::string, BaseKinetics* > (this->concrete_name, this) );           
}

//---------------------------------------------------------------------------
NAME::DoubleTrapKinetics::DoubleTrapKinetics()
:
BaseKinetics ()
{
    //FcstUtilities::log << "DoubleTrap kinetics." << std::endl;
    
    //Values for the free energies of activation and adsorption based on the PEMFC values from 
    //Wang's 2008 paper.
    G_DA_0 = 0.358; //overridden by initialize()
    G_RA_0 = 0.531; //overridden by initialize()
    G_RT_0 = 0.578; //overridden by initialize()
    G_RD_0 = 0.5;   //overridden by initialize()
    G_O = -0.524;   //overridden by initialize()
    G_OH = -0.12;   //overridden by initialize()
    prefactor = 1000; //overridden by initialize()  NOTE: <- To match Parthasarathy's data this value should be x by 5.2 since this is the roughness factor of his electrode (M. Secanell, 2015)
    alpha = 0.5;    //overridden by initialize()
}

//---------------------------------------------------------------------------
NAME::DoubleTrapKinetics::~DoubleTrapKinetics()
{}

//---------------------------------------------------------------------------
void
NAME::DoubleTrapKinetics::declare_parameters(ParameterHandler &param) const
{
    param.enter_subsection("Kinetics");
    {
        param.enter_subsection(NAME::DoubleTrapKinetics::concrete_name);
        {
            param.declare_entry ("DA equilibrium free energy, [eV]",    // Hard coding Michael's data saves having to pass values through text file.
                                 "0.3907",                              // DoubleTrap Results from 2008 paper "0.358", Michael Moore's fitted data "0.3907" (more accurate)
                                 Patterns::Double(),
                                 "Equilibrium free energy for dissociative adsorption");
            param.declare_entry ("RA equilibrium free energy, [eV]",
                                 "0.6094",                              // DoubleTrap Results from 2008 paper "0.531", Michael Moore's fitted data "0.6094" (more accurate)
                                 Patterns::Double(),
                                 "Equilibrium free energy for reductive adsorption");
            param.declare_entry ("RT equilibrium free energy, [eV]",
                                 "0.5904",                              // DoubleTrap Results from 2008 paper "0.578", Michael Moore's fitted data "0.5904" (more accurate)
                                 Patterns::Double(),
                                 "Equilibrium free energy for reductive desorption");
            param.declare_entry ("RD equilibrium free energy, [eV]",
                                 "0.278",                               // DoubleTrap Results from 2008 paper "0.5", Michael Moore's fitted data "0.278" (more accurate)
                                 Patterns::Double(),
                                 "Equilibrium free energy for reductive desorption");
            param.declare_entry ("O coverage equilibrium free energy, [eV]",
                                 "-0.3431",                             // DoubleTrap Results from 2008 paper "-0.524", Michael Moore's fitted data "-0.3431" (more accurate)
                                 Patterns::Double(),
                                 "Equilibrium free energy for the adsorbed O intermediate");
            param.declare_entry ("OH coverage equilibrium free energy, [eV]",
                                 "-0.376",                              // DoubleTrap Results from paper "-0.12", Michael Moore's fitted data "-0.376" (more accurate)
                                 Patterns::Double(),
                                 "Equilibrium free energy for the adsorbed OH intermediate");
            param.declare_entry ("Reference prefactor, [A/cm^2]",
                                 "1000.",
                                 Patterns::Double(),
                                 "Reference prefactor that sets the scale of the activation free energy");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
void
NAME::DoubleTrapKinetics::initialize(ParameterHandler &param)
{
    param.enter_subsection("Kinetics"); 
    {
        param.enter_subsection(NAME::DoubleTrapKinetics::concrete_name); 
        {
            G_DA_0 = param.get_double("DA equilibrium free energy, [eV]");
            G_RA_0 = param.get_double("RA equilibrium free energy, [eV]");
            G_RT_0 = param.get_double("RT equilibrium free energy, [eV]");
            G_RD_0 = param.get_double("RD equilibrium free energy, [eV]");
            G_O = param.get_double("O coverage equilibrium free energy, [eV]");
            G_OH = param.get_double("OH coverage equilibrium free energy, [eV]");
            prefactor= param.get_double("Reference prefactor, [A/cm^2]");
        }
        param.leave_subsection(); 
    }
    param.leave_subsection();
}  

//---------------------------------------------------------------------------
void
NAME::DoubleTrapKinetics::compute_energies()
{
    if (!kin_param_initialized)
        init_kin_param();
    
    n_quad = phi_m.size();
    //Initialise the vectors
    free_energies.resize(8, std::vector<double>(n_quad, 0.0));
    intrinsic_rates.resize(8, std::vector<double>(n_quad, 0.0));
    theta_OH.resize(n_quad,0.0);
    theta_O.resize(n_quad,0.0);
    
    for (unsigned int j = 0; j<n_quad; ++j)
    {
        //Compute the activation energies.          
        free_energies[0][j] = G_DA_0;
        free_energies[1][j] = G_DA_0 - G_O;
        free_energies[2][j] = G_RA_0 + alpha*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])));
        free_energies[3][j] = G_RA_0 - G_OH - (1.0-alpha)*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j]))); 
        free_energies[4][j] = G_RT_0 + alpha*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])));
        free_energies[5][j] = G_RT_0 - G_OH + G_O - (1.0-alpha)*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j]))); 
        free_energies[6][j] = G_RD_0 + alpha*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j])));
        free_energies[7][j] = G_RD_0 + G_OH - (1.0-alpha)*(phi_s[j]-phi_m[j]-(catalyst->voltage_cell_th(T[j]))); 
                
        for(unsigned int k = 0; k<8; ++k)
            intrinsic_rates[k][j] = exp(-free_energies[k][j]/(K*T[j]));
        
        //Compute the oxygen fraction. Note that 0.5 is hard coded as the reaction order, as changing 
        //it would require the entire system of equations to be changed.
        double C = oxygen_concentration_ratio(j);
        
        double Ch = proton_concentration_ratio(j);


        //Compute the coverages
        theta_OH[j] = ( C*intrinsic_rates[0][j]*(Ch*C*intrinsic_rates[2][j]+intrinsic_rates[7][j]-Ch*intrinsic_rates[4][j]) - (Ch*C*intrinsic_rates[2][j]+intrinsic_rates[7][j])*(C*intrinsic_rates[0][j]+intrinsic_rates[1][j]+Ch*intrinsic_rates[4][j]) ) /
        ( (C*intrinsic_rates[0][j]-intrinsic_rates[5][j])*(Ch*C*intrinsic_rates[2][j]+intrinsic_rates[7][j]-Ch*intrinsic_rates[4][j]) - (Ch*C*intrinsic_rates[2][j]+intrinsic_rates[3][j]+intrinsic_rates[5][j]+Ch*intrinsic_rates[6][j]+intrinsic_rates[7][j])*(C*intrinsic_rates[0][j]+intrinsic_rates[1][j]+Ch*intrinsic_rates[4][j]) );

        theta_O[j] = ( C*intrinsic_rates[0][j]*(Ch*C*intrinsic_rates[2][j]+intrinsic_rates[3][j]+intrinsic_rates[5][j]+Ch*intrinsic_rates[6][j]+intrinsic_rates[7][j]) - (Ch*C*intrinsic_rates[2][j]+intrinsic_rates[7][j])*(C*intrinsic_rates[0][j]-intrinsic_rates[5][j]) ) /
        ( (C*intrinsic_rates[0][j]+intrinsic_rates[1][j]+Ch*intrinsic_rates[4][j])*(Ch*C*intrinsic_rates[2][j]+intrinsic_rates[3][j]+intrinsic_rates[5][j]+Ch*intrinsic_rates[6][j]+intrinsic_rates[7][j]) - (Ch*C*intrinsic_rates[2][j]+intrinsic_rates[7][j]-Ch*intrinsic_rates[4][j])*(C*intrinsic_rates[0][j]-intrinsic_rates[5][j]) );
    }
}

//---------------------------------------------------------------------------
void
NAME::DoubleTrapKinetics::current_density (std::vector<double>& coef)
{  
    //Call compute_energies function to compute all the necessary parameters
    compute_energies();
    coef.resize(n_quad);
    
    // Compute the coefficient, first loop over the quadrature points
    for (unsigned int i=0; i<n_quad; ++i)
    {
        double corrector = negative_concentration_correction(i);
        
        double Ch = proton_concentration_ratio(i);
        
        coef[i] = corrector*(Ch*prefactor*intrinsic_rates[6][i]*theta_OH[i] - prefactor*intrinsic_rates[7][i]*(1.0 - theta_O[i] - theta_OH[i]));
        
    }
}

//---------------------------------------------------------------------------
void 
NAME::DoubleTrapKinetics::d_OH_coverage_du (std::map< VariableNames, std::vector<double> >& dtheta_du) const
{
    Assert(kin_param_initialized, ExcInternalError());

    //Need to use the quotient rule to compute the derivative. Start with the value of the 
    //numerator and denominator.
    std::vector<double> num(n_quad, 0.0);
    std::vector<double> denom(n_quad, 0.0);
    
    for (unsigned int j = 0; j<n_quad; ++j)
    {
        //Compute the oxygen fraction
        double C = oxygen_concentration_ratio(j);
        
        double Ch = proton_concentration_ratio(j);


        num[j] = (C*intrinsic_rates[0][j]*(Ch*C*intrinsic_rates[2][j]+intrinsic_rates[7][j]-Ch*intrinsic_rates[4][j])
        -(Ch*C*intrinsic_rates[2][j]+intrinsic_rates[7][j])*(C*intrinsic_rates[0][j]+intrinsic_rates[1][j]+Ch*intrinsic_rates[4][j]) );
        
        denom[j] = (C*intrinsic_rates[0][j]-intrinsic_rates[5][j])*(Ch*C*intrinsic_rates[2][j]+intrinsic_rates[7][j]-Ch*intrinsic_rates[4][j])
        -(Ch*C*intrinsic_rates[2][j]+intrinsic_rates[3][j]+intrinsic_rates[5][j]+Ch*intrinsic_rates[6][j]+intrinsic_rates[7][j])*(C*intrinsic_rates[0][j]+intrinsic_rates[1][j]+Ch*intrinsic_rates[4][j]);
    }
    
    std::vector<std::vector<double> > drates_du(8, std::vector<double>(n_quad, 0.0));
    
    //Loop over the flags
    for (unsigned int i=0; i<derivative_flags.size(); ++i)
    {
        std::vector<double> dtheta(n_quad, 0.);         // Defaulting all derivatives to 0.0
        
        if (derivative_flags[i] == oxygen_molar_fraction)
        {
            //First compute the derivative of the intrinsic rates wrt to oxygen
            drate_du(drates_du, derivative_flags[i]); //variable not used
            //Now compute the derivative of the numerator (of the coverage term) and denominator
            // with respect to oxygen.
            std::vector<double> dnum_co2(n_quad, 0.0);
            std::vector<double> ddenom_co2(n_quad, 0.0);
            
            for (unsigned int j = 0; j<n_quad; ++j)
            {
                double C = oxygen_concentration_ratio(j);
                double dC_dco2;
                if (reactants_map.at(oxygen_concentration)[j] < 0.0)
                {
                    dC_dco2 = 1.0e-12;
                }
                else
                {
                    dC_dco2 = pow(1/c_ref_oxygen, 0.5)*0.5*(1/ pow((reactants_map.at(oxygen_concentration)[j]),0.5))*(p_total/electrolyte->get_H_O2());
                }
                
                double Ch = proton_concentration_ratio(j);
                
                dnum_co2[j] = -intrinsic_rates[0][j]*Ch*intrinsic_rates[4][j]*dC_dco2 - intrinsic_rates[1][j]*Ch*intrinsic_rates[2][j]*dC_dco2 - Ch*intrinsic_rates[2][j]*Ch*intrinsic_rates[4][j]*dC_dco2;
                
                ddenom_co2[j] = -intrinsic_rates[0][j]*Ch*intrinsic_rates[4][j]*dC_dco2 - Ch*intrinsic_rates[2][j]*intrinsic_rates[5][j]*dC_dco2 - intrinsic_rates[0][j]*intrinsic_rates[3][j]*dC_dco2
                - intrinsic_rates[0][j]*intrinsic_rates[5][j]*dC_dco2 - intrinsic_rates[0][j]*Ch*intrinsic_rates[6][j]*dC_dco2 - intrinsic_rates[1][j]*Ch*intrinsic_rates[2][j]*dC_dco2
                - Ch*intrinsic_rates[2][j]*Ch*intrinsic_rates[4][j]*dC_dco2;
                
                dtheta[j]= (denom[j]*dnum_co2[j] - num[j]*ddenom_co2[j])/(denom[j]*denom[j]);
            }
            
            dtheta_du[oxygen_molar_fraction] = dtheta;
        }//end derivative oxygen
        

        else if (derivative_flags[i] == proton_concentration)
        {
            std::vector<double> dnum_ch(n_quad, 0.0);
            std::vector<double> ddenom_ch(n_quad, 0.0);

            for (unsigned int j = 0; j<n_quad; ++j)
            {
                double C = oxygen_concentration_ratio(j);

                double Ch = proton_concentration_ratio(j);
                double dCh = 1/c_ref_protons;

                dnum_ch[j] = C*intrinsic_rates[0][j]*(dCh*C*intrinsic_rates[2][j] -dCh*intrinsic_rates[4][j]) -
                        ((Ch*C*intrinsic_rates[2][j] +intrinsic_rates[7][j])*(dCh*intrinsic_rates[4][j]) +
                        (dCh*C*intrinsic_rates[2][j])*(C*intrinsic_rates[0][j] + intrinsic_rates[1][j] +Ch*intrinsic_rates[4][j]));

                ddenom_ch[j] = ((C*intrinsic_rates[0][j] - intrinsic_rates[5][j])*(dCh*C*intrinsic_rates[2][j] - dCh*intrinsic_rates[4][j]))-(
                        (Ch*C*intrinsic_rates[2][j] + intrinsic_rates[3][j]+ intrinsic_rates[5][j] + Ch*intrinsic_rates[6][j] + intrinsic_rates[7][j])*(dCh*intrinsic_rates[4][j]) +
                        (dCh*C*intrinsic_rates[2][j] +dCh*intrinsic_rates[6][j])*(C*intrinsic_rates[0][j] + intrinsic_rates[1][j] + Ch*intrinsic_rates[4][j]));

                dtheta[j]= (denom[j]*dnum_ch[j] - num[j]*ddenom_ch[j])/(denom[j]*denom[j]);
            }

            dtheta_du[proton_concentration] = dtheta;

        }

        else if (derivative_flags[i] == electronic_electrical_potential)
        {
            //First compute the derivative of the intrinsic rates wrt to the electronic potential
            drate_du(drates_du, derivative_flags[i]);
            //Now compute the derivative of the numerator (of the coverage term) and denominator
            // with respect to electronic potential.
            std::vector<double> dnum_phi_s(n_quad, 0.0);
            std::vector<double> ddenom_phi_s(n_quad, 0.0);
            for (unsigned int j = 0; j<n_quad; ++j)
            {
                double C = oxygen_concentration_ratio(j);
                double Ch = proton_concentration_ratio(j);
                
                dnum_phi_s[j] = -C*intrinsic_rates[0][j]*Ch*drates_du[4][j] - C*intrinsic_rates[1][j]*Ch*drates_du[2][j] - C*Ch*intrinsic_rates[2][j]*Ch*drates_du[4][j] - C*Ch*intrinsic_rates[4][j]*Ch*drates_du[2][j]
                - intrinsic_rates[1][j]*drates_du[7][j] - Ch*intrinsic_rates[4][j]*drates_du[7][j] - intrinsic_rates[7][j]*Ch*drates_du[4][j];
                
                ddenom_phi_s[j] = -C*intrinsic_rates[0][j]*Ch*drates_du[4][j] - Ch*C*intrinsic_rates[2][j]*drates_du[5][j] - C*intrinsic_rates[5][j]*Ch*drates_du[2][j] - intrinsic_rates[5][j]*drates_du[7][j]
                - intrinsic_rates[7][j]*drates_du[5][j] - C*intrinsic_rates[0][j]*drates_du[3][j] - C*intrinsic_rates[0][j]*drates_du[5][j] - C*intrinsic_rates[0][j]*Ch*drates_du[6][j]
                - C*intrinsic_rates[1][j]*Ch*drates_du[2][j] - intrinsic_rates[1][j]*drates_du[3][j] - intrinsic_rates[1][j]*drates_du[5][j]- intrinsic_rates[1][j]*Ch*drates_du[6][j]
                - intrinsic_rates[1][j]*drates_du[7][j] - Ch*C*intrinsic_rates[2][j]*Ch*drates_du[4][j] - Ch*C*intrinsic_rates[4][j]*Ch*drates_du[2][j] - intrinsic_rates[3][j]*Ch*drates_du[4][j]
                - Ch*intrinsic_rates[4][j]*drates_du[3][j] - Ch*intrinsic_rates[4][j]*Ch*drates_du[6][j] - Ch*intrinsic_rates[6][j]*Ch*drates_du[4][j] - Ch*intrinsic_rates[4][j]*drates_du[7][j]
                - intrinsic_rates[7][j]*Ch*drates_du[4][j];
                
                dtheta[j]= (denom[j]*dnum_phi_s[j] - num[j]*ddenom_phi_s[j])/(denom[j]*denom[j]);
            }
            
            dtheta_du[electronic_electrical_potential] = dtheta;
        }//end derivative electronic
        
        else if (derivative_flags[i] == protonic_electrical_potential)
        {
            //First compute the derivative of the intrinsic rates wrt the protonic potential
            drate_du(drates_du, derivative_flags[i]);
            //Now compute the derivative of the numerator (of the coverage term) and denominator
            // with respect to protonic potential.
            std::vector<double> dnum_phi_m(n_quad, 0.0);
            std::vector<double> ddenom_phi_m(n_quad, 0.0);
            for (unsigned int j = 0; j<n_quad; ++j)
            {
                double C = oxygen_concentration_ratio(j);
                double Ch = proton_concentration_ratio(j);
                
                dnum_phi_m[j] = -C*intrinsic_rates[0][j]*Ch*drates_du[4][j] - C*intrinsic_rates[1][j]*Ch*drates_du[2][j] - Ch*C*intrinsic_rates[2][j]*Ch*drates_du[4][j] - Ch*C*intrinsic_rates[4][j]*Ch*drates_du[2][j]
                - intrinsic_rates[1][j]*drates_du[7][j] - Ch*intrinsic_rates[4][j]*drates_du[7][j] - intrinsic_rates[7][j]*Ch*drates_du[4][j];

                ddenom_phi_m[j] = -C*intrinsic_rates[0][j]*Ch*drates_du[4][j] - Ch*C*intrinsic_rates[2][j]*drates_du[5][j] - C*intrinsic_rates[5][j]*Ch*drates_du[2][j] - intrinsic_rates[5][j]*drates_du[7][j]
                - intrinsic_rates[7][j]*drates_du[5][j] - C*intrinsic_rates[0][j]*drates_du[3][j] - C*intrinsic_rates[0][j]*drates_du[5][j] - C*intrinsic_rates[0][j]*Ch*drates_du[6][j]
                - C*intrinsic_rates[1][j]*Ch*drates_du[2][j] - intrinsic_rates[1][j]*drates_du[3][j] - intrinsic_rates[1][j]*drates_du[5][j]- intrinsic_rates[1][j]*Ch*drates_du[6][j]
                - intrinsic_rates[1][j]*drates_du[7][j] - Ch*C*intrinsic_rates[2][j]*Ch*drates_du[4][j] - C*Ch*intrinsic_rates[4][j]*Ch*drates_du[2][j] - intrinsic_rates[3][j]*Ch*drates_du[4][j]
                - Ch*intrinsic_rates[4][j]*drates_du[3][j] - Ch*intrinsic_rates[4][j]*Ch*drates_du[6][j] - Ch*intrinsic_rates[6][j]*Ch*drates_du[4][j] - Ch*intrinsic_rates[4][j]*drates_du[7][j]
                - intrinsic_rates[7][j]*Ch*drates_du[4][j];

                dtheta[j]= (denom[j]*dnum_phi_m[j] - num[j]*ddenom_phi_m[j])/(denom[j]*denom[j]);
            }
            
            dtheta_du[protonic_electrical_potential] = dtheta;
        }//end derivative protonic
        
        else if (derivative_flags[i] == temperature_of_REV)
        {
            //First compute the derivative of the intrinsic rates wrt the protonic potential
            drate_du(drates_du, derivative_flags[i]);
            //Now compute the derivative of the numerator (of the coverage term) and denominator
            // with respect to Temperature.
            std::vector<double> dnum_T(n_quad, 0.0);
            std::vector<double> ddenom_T(n_quad, 0.0);
            for (unsigned int j = 0; j<n_quad; ++j)
            {
                double C = oxygen_concentration_ratio(j);
                double Ch = proton_concentration_ratio(j);

                dnum_T[j] = (C*intrinsic_rates[0][j]*(Ch*C*drates_du[2][j]+drates_du[7][j]-Ch*drates_du[4][j])) +
                ((Ch*C*intrinsic_rates[2][j]+intrinsic_rates[7][j]-Ch*intrinsic_rates[4][j])*C*drates_du[0][j]) -
                ((C*Ch*intrinsic_rates[2][j]+intrinsic_rates[7][j])*(C*drates_du[0][j]+drates_du[1][j]+Ch*drates_du[4][j])) -
                ((C*intrinsic_rates[0][j]+intrinsic_rates[1][j]+Ch*intrinsic_rates[4][j])*(C*Ch*drates_du[2][j]+drates_du[7][j]));
                
                ddenom_T[j] = ((C*intrinsic_rates[0][j]-intrinsic_rates[5][j])*(C*Ch*drates_du[2][j]+drates_du[7][j]-Ch*drates_du[4][j])) +
                ((C*Ch*intrinsic_rates[2][j]+intrinsic_rates[7][j]-Ch*intrinsic_rates[4][j])*(C*drates_du[0][j]-drates_du[5][j])) -
                ((C*Ch*intrinsic_rates[2][j]+intrinsic_rates[3][j]+intrinsic_rates[5][j]+Ch*intrinsic_rates[6][j]+intrinsic_rates[7][j])*(C*drates_du[0][j]+drates_du[1][j]+Ch*drates_du[4][j])) -
                ((C*intrinsic_rates[0][j]+intrinsic_rates[1][j]+Ch*intrinsic_rates[4][j])*(C*Ch*drates_du[2][j]+drates_du[3][j]+drates_du[5][j]+Ch*drates_du[6][j]+drates_du[7][j]));
                
                dtheta[j]= (denom[j]*dnum_T[j] - num[j]*ddenom_T[j])/(denom[j]*denom[j]);
            }
            
            dtheta_du[temperature_of_REV] = dtheta;
        }//end derivative Temperature
        
        else
        {
            dtheta_du[ derivative_flags[i] ] = dtheta;
        }
    }
}

//---------------------------------------------------------------------------
void 
NAME::DoubleTrapKinetics::d_O_coverage_du (std::map< VariableNames, std::vector<double> >& dtheta_du) const
{
    Assert(kin_param_initialized, ExcInternalError());
    
    //Need to use the quotient rule to compute the derivative. Start with the value of the numerator and denominator.
    std::vector<double> num(n_quad, 0.0);
    std::vector<double> denom(n_quad, 0.0);
    
    for (unsigned int j = 0; j<n_quad; ++j)
    {
        //Compute the oxygen fraction
        double C = oxygen_concentration_ratio(j);
        double Ch = proton_concentration_ratio(j);
        
        num[j] = ( (C*intrinsic_rates[0][j])*(Ch*C*intrinsic_rates[2][j]+intrinsic_rates[3][j]+intrinsic_rates[5][j]+Ch*intrinsic_rates[6][j]+intrinsic_rates[7][j])
        - ( (C*Ch*intrinsic_rates[2][j]+intrinsic_rates[7][j])*(C*intrinsic_rates[0][j]-intrinsic_rates[5][j]) ) );
        
        denom[j] = ( (C*intrinsic_rates[0][j]+intrinsic_rates[1][j]+Ch*intrinsic_rates[4][j])*(C*Ch*intrinsic_rates[2][j]+intrinsic_rates[3][j]+intrinsic_rates[5][j]+Ch*intrinsic_rates[6][j]+intrinsic_rates[7][j]) )
        -( (C*Ch*intrinsic_rates[2][j]+intrinsic_rates[7][j]-Ch*intrinsic_rates[4][j])*(C*intrinsic_rates[0][j]-intrinsic_rates[5][j]) );
    }
    
    std::vector<std::vector<double> > drates_du(8, std::vector<double>(n_quad, 0.0));
    
    //Loop over the flags
    for (unsigned int i=0; i<derivative_flags.size(); ++i)
    {
        std::vector<double> dtheta(n_quad, 0.);         // Defaulting all derivatives to 0.0
        
        if (derivative_flags[i] == oxygen_molar_fraction)
        {
            //First compute the derivative of the intrinsic rates wrt to oxygen
            drate_du(drates_du, derivative_flags[i]);
            //Now compute the derivative of the numerator (of the coverage term) and denominator
            // with respect to oxygen.
            std::vector<double> dnum_co2(n_quad, 0.0);
            std::vector<double> ddenom_co2(n_quad, 0.0);
            
            for (unsigned int j = 0; j<n_quad; ++j)
            {
                double C = oxygen_concentration_ratio(j);
                double Ch = proton_concentration_ratio(j);
                double dC_dco2;
                if (reactants_map.at(oxygen_concentration)[j] < 0.0)
                {
                    dC_dco2 = 1.0e-12;
                }
                else
                {
                    dC_dco2 = pow(1./c_ref_oxygen, 0.5)*0.5*(1./ pow((reactants_map.at(oxygen_concentration)[j]),0.5))*(p_total/electrolyte->get_H_O2());
                }
                               
                dnum_co2[j] = intrinsic_rates[0][j]*intrinsic_rates[3][j]*dC_dco2 + intrinsic_rates[0][j]*intrinsic_rates[5][j]*dC_dco2 + intrinsic_rates[0][j]*Ch*intrinsic_rates[6][j]*dC_dco2
                        + Ch*intrinsic_rates[2][j]*intrinsic_rates[5][j]*dC_dco2;
                
                ddenom_co2[j] = +intrinsic_rates[0][j]*intrinsic_rates[3][j]*dC_dco2 + intrinsic_rates[0][j]*intrinsic_rates[5][j]*dC_dco2 + intrinsic_rates[0][j]*Ch*intrinsic_rates[6][j]*dC_dco2
                + intrinsic_rates[1][j]*Ch*intrinsic_rates[2][j]*dC_dco2 + Ch*intrinsic_rates[2][j]*Ch*intrinsic_rates[4][j]*dC_dco2 + Ch*intrinsic_rates[2][j]*intrinsic_rates[5][j]*dC_dco2
                + intrinsic_rates[0][j]*Ch*intrinsic_rates[4][j]*dC_dco2;
                
                dtheta[j]= (denom[j]*dnum_co2[j] - num[j]*ddenom_co2[j])/(denom[j]*denom[j]);
            }
            
            dtheta_du[oxygen_molar_fraction] = dtheta;
        }//end derivative oxygen
        
        else if (derivative_flags[i] == proton_concentration)
        {
            std::vector<double> dnum_ch(n_quad, 0.0);
            std::vector<double> ddenom_ch(n_quad, 0.0);

            for (unsigned int j = 0; j<n_quad; ++j)
            {
                double C = oxygen_concentration_ratio(j);
                double Ch = proton_concentration_ratio(j);
                double dCh = 1/c_ref_protons;

                dnum_ch[j] = C*intrinsic_rates[0][j]*(dCh*C*intrinsic_rates[2][j] + dCh*intrinsic_rates[6][j]) -((dCh*C*intrinsic_rates[2][j])*(C*intrinsic_rates[0][j] -intrinsic_rates[5][j]));

                ddenom_ch[j] = (((C*intrinsic_rates[0][j]+intrinsic_rates[1][j]+Ch*intrinsic_rates[4][j])*(dCh*C*intrinsic_rates[2][j] + dCh*intrinsic_rates[6][j])) +(dCh*intrinsic_rates[4][j]*
                        (Ch*C*intrinsic_rates[2][j] + intrinsic_rates[3][j] + intrinsic_rates[5][j] +Ch*intrinsic_rates[6][j] + intrinsic_rates[7][j])) -
                        ((dCh*C*intrinsic_rates[2][j] -dCh*intrinsic_rates[4][j])*(C*intrinsic_rates[0][j] - intrinsic_rates[5][j])));

                dtheta[j]= (denom[j]*dnum_ch[j] - num[j]*ddenom_ch[j])/(denom[j]*denom[j]);
            }

            dtheta_du[proton_concentration] = dtheta;

        }


        else if (derivative_flags[i] == electronic_electrical_potential)
        {
            //First compute the derivative of the intrinsic rates wrt to the electronic potential
            drate_du(drates_du, derivative_flags[i]);
            //Now compute the derivative of the numerator (of the coverage term) and denominator
            // with respect to electronic potential.
            std::vector<double> dnum_phi_s(n_quad, 0.0);
            std::vector<double> ddenom_phi_s(n_quad, 0.0);
            for (unsigned int j = 0; j<n_quad; ++j)
            {
                double C = oxygen_concentration_ratio(j);
                
                double Ch = proton_concentration_ratio(j);

                dnum_phi_s[j] = C*intrinsic_rates[0][j]*drates_du[3][j] + C*intrinsic_rates[0][j]*drates_du[5][j] + C*intrinsic_rates[0][j]*Ch*drates_du[6][j] + Ch*C*intrinsic_rates[2][j]*drates_du[5][j]
                + C*intrinsic_rates[5][j]*Ch*drates_du[2][j] + intrinsic_rates[5][j]*drates_du[7][j] + intrinsic_rates[7][j]*drates_du[5][j];
                
                ddenom_phi_s[j] = C*intrinsic_rates[0][j]*Ch*drates_du[4][j] + C*Ch*intrinsic_rates[2][j]*drates_du[5][j] + C*intrinsic_rates[5][j]*Ch*drates_du[2][j] + intrinsic_rates[5][j]*drates_du[7][j]
                + intrinsic_rates[7][j]*drates_du[5][j] + C*intrinsic_rates[0][j]*drates_du[3][j] + C*intrinsic_rates[0][j]*drates_du[5][j] + C*intrinsic_rates[0][j]*Ch*drates_du[6][j]
                + C*intrinsic_rates[1][j]*Ch*drates_du[2][j] + intrinsic_rates[1][j]*drates_du[3][j] + intrinsic_rates[1][j]*drates_du[5][j]+ intrinsic_rates[1][j]*Ch*drates_du[6][j]
                + intrinsic_rates[1][j]*drates_du[7][j] + C*Ch*intrinsic_rates[2][j]*Ch*drates_du[4][j] + C*Ch*intrinsic_rates[4][j]*Ch*drates_du[2][j] + intrinsic_rates[3][j]*Ch*drates_du[4][j]
                + Ch*intrinsic_rates[4][j]*drates_du[3][j] + Ch*intrinsic_rates[4][j]*Ch*drates_du[6][j] + Ch*intrinsic_rates[6][j]*Ch*drates_du[4][j] + Ch*intrinsic_rates[4][j]*drates_du[7][j]
                + intrinsic_rates[7][j]*Ch*drates_du[4][j];
                
                dtheta[j]= (denom[j]*dnum_phi_s[j] - num[j]*ddenom_phi_s[j])/(denom[j]*denom[j]);
            }
            
            dtheta_du[electronic_electrical_potential] = dtheta;
        }//end derivative electronic
        
        else if (derivative_flags[i] == protonic_electrical_potential)
        {
            //First compute the derivative of the intrinsic rates wrt the protonic potential
            drate_du(drates_du, derivative_flags[i]);
            //Now compute the derivative of the numerator (of the coverage term) and denominator
            // with respect to protonic potential.
            std::vector<double> dnum_phi_m(n_quad, 0.0);
            std::vector<double> ddenom_phi_m(n_quad, 0.0);
            for (unsigned int j = 0; j<n_quad; ++j)
            {
                double C = oxygen_concentration_ratio(j);
                
                double Ch = proton_concentration_ratio(j);
                
                dnum_phi_m[j] = C*intrinsic_rates[0][j]*drates_du[3][j] + C*intrinsic_rates[0][j]*drates_du[5][j] + C*intrinsic_rates[0][j]*Ch*drates_du[6][j] + C*Ch*intrinsic_rates[2][j]*drates_du[5][j]
                + C*intrinsic_rates[5][j]*Ch*drates_du[2][j] + intrinsic_rates[5][j]*drates_du[7][j] + intrinsic_rates[7][j]*drates_du[5][j];

                ddenom_phi_m[j] = C*intrinsic_rates[0][j]*Ch*drates_du[4][j] + C*Ch*intrinsic_rates[2][j]*drates_du[5][j] + C*intrinsic_rates[5][j]*Ch*drates_du[2][j] + intrinsic_rates[5][j]*drates_du[7][j]
                + intrinsic_rates[7][j]*drates_du[5][j] + C*intrinsic_rates[0][j]*drates_du[3][j] + C*intrinsic_rates[0][j]*drates_du[5][j] + C*intrinsic_rates[0][j]*Ch*drates_du[6][j]
                + C*intrinsic_rates[1][j]*Ch*drates_du[2][j] + intrinsic_rates[1][j]*drates_du[3][j] + intrinsic_rates[1][j]*drates_du[5][j]+ intrinsic_rates[1][j]*Ch*drates_du[6][j]
                + intrinsic_rates[1][j]*drates_du[7][j] + C*Ch*intrinsic_rates[2][j]*Ch*drates_du[4][j] + C*Ch*intrinsic_rates[4][j]*Ch*drates_du[2][j] + intrinsic_rates[3][j]*Ch*drates_du[4][j]
                + Ch*intrinsic_rates[4][j]*drates_du[3][j] + Ch*intrinsic_rates[4][j]*Ch*drates_du[6][j] + Ch*intrinsic_rates[6][j]*Ch*drates_du[4][j] + Ch*intrinsic_rates[4][j]*drates_du[7][j]
                + intrinsic_rates[7][j]*Ch*drates_du[4][j];

                dtheta[j]= (denom[j]*dnum_phi_m[j] - num[j]*ddenom_phi_m[j])/(denom[j]*denom[j]);
            }
            
            dtheta_du[protonic_electrical_potential] = dtheta;
        }//end derivative protonic
        
        else if (derivative_flags[i] == temperature_of_REV)
        {
            //First compute the derivative of the intrinsic rates wrt the Temperature
            drate_du(drates_du, derivative_flags[i]);
            //Now compute the derivative of the numerator (of the coverage term) and denominator
            // with respect to Temperature.
            std::vector<double> dnum_T(n_quad, 0.0);
            std::vector<double> ddenom_T(n_quad, 0.0);
            for (unsigned int j = 0; j<n_quad; ++j)
            {
                double C = oxygen_concentration_ratio(j);
                
                double Ch = proton_concentration_ratio(j);
 
                dnum_T[j] = (C*intrinsic_rates[0][j]*(C*Ch*drates_du[2][j]+drates_du[3][j]+drates_du[5][j]+Ch*drates_du[6][j]+drates_du[7][j])) +
                (C*drates_du[0][j]*(C*Ch*intrinsic_rates[2][j]+intrinsic_rates[3][j]+intrinsic_rates[5][j]+Ch*intrinsic_rates[6][j]+intrinsic_rates[7][j])) -
                ((C*Ch*intrinsic_rates[2][j]+intrinsic_rates[7][j])*(C*drates_du[0][j]-drates_du[5][j])) -
                ((C*intrinsic_rates[0][j]-intrinsic_rates[5][j])*(C*Ch*drates_du[2][j]+drates_du[7][j]));
                
                ddenom_T[j] = ((C*intrinsic_rates[0][j]+intrinsic_rates[1][j]+Ch*intrinsic_rates[4][j])*(C*Ch*drates_du[2][j]+drates_du[3][j]+drates_du[5][j]+Ch*drates_du[6][j]+drates_du[7][j])) +
                ((C*Ch*intrinsic_rates[2][j]+intrinsic_rates[3][j]+intrinsic_rates[5][j]+Ch*intrinsic_rates[6][j]+intrinsic_rates[7][j])*(C*drates_du[0][j]+drates_du[1][j]+Ch*drates_du[4][j])) -
                ((C*Ch*intrinsic_rates[2][j]+intrinsic_rates[7][j]-Ch*intrinsic_rates[4][j])*(C*drates_du[0][j]-drates_du[5][j])) -
                ((C*intrinsic_rates[0][j]-intrinsic_rates[5][j])*(C*Ch*drates_du[2][j]+drates_du[7][j]-Ch*drates_du[4][j]));

                dtheta[j]= (denom[j]*dnum_T[j] - num[j]*ddenom_T[j])/(denom[j]*denom[j]);
            }
            
            dtheta_du[temperature_of_REV] = dtheta;
        }//end derivative Temperature
        
        else
        {
            dtheta_du[ derivative_flags[i] ] = dtheta;
        }
    }
}

//---------------------------------------------------------------------------
void
NAME::DoubleTrapKinetics::drate_du(std::vector<std::vector<double> >& drates_du, const VariableNames& flag) const
{
    Assert(kin_param_initialized, ExcInternalError());

    if (flag == oxygen_molar_fraction)
    {
        //Loop over the first 4 rates only, as the others do not depend on oxygen
        for (unsigned int i = 0; i<4; ++i)
        {
            //Loop for quadrature points
            for(unsigned int j = 0; j<n_quad; ++j)
                drates_du[i][j] = 0.0;
        }
        //Set the remaining derivatives to zero.
        for(unsigned int i = 4; i<8; ++i)
        {
            for (unsigned int j = 0; j<n_quad; ++j)
                drates_du[i][j] = 0.0;
        }
    }//end derivative oxygen

    if (flag == proton_concentration)
    {
        //All derivatives are 0
        for (unsigned int i = 0; i<8; ++i)
        {
            //Loop for quadrature points
            for(unsigned int j = 0; j<n_quad; ++j)
                drates_du[i][j] = 0.0;
        }
    }//end derivative proton_concentration

    
    else if (flag == electronic_electrical_potential)
    {
        //Loop over all the rates, bar the first two as the forward and backward reactions for the 
        //disassociative adsorption does not depend on the potential. In the first loop we go over
        //the even number rates (i.e. the forward rates) as they depend on alpha...
        for (unsigned int i = 2; i<8; i=i+2)
        {
            //Loop for quadrature points
            for(unsigned int j = 0; j<n_quad; ++j)
                drates_du[i][j] = ((-1.0)/(K*T[j]))*intrinsic_rates[i][j]*alpha;
        }
        
        //whereas the odd numbered rates (i.e. the backward rates) depend on -(1-alpha)
        for (unsigned int i = 3; i<8; i=i+2)
        {
            //Loop for quadrature points
            for(unsigned int j = 0; j<n_quad; ++j)
                drates_du[i][j] = ((-1.0)/(K*T[j]))*intrinsic_rates[i][j]*(-1.0)*(1.0-alpha);
        }
        
        //Loop over the first two rates and set the values to zero. 
        for (unsigned int i = 0; i<2; ++i)
        {
            //Loop for quadrature points
            for(unsigned int j = 0; j<n_quad; ++j)
                drates_du[i][j] = 0.0;
        }     
    }//end derivative electronic
    
    //The protonic is the same as the electronic but with the signs switched.
    else if (flag == protonic_electrical_potential)
    {
        //Loop over all the rates, bar the first two as the forward and backward reactions for the 
        //disassociative adsorption does not depend on the potential. In the first loop we go over
        //the even number rates (i.e. the forward rates) as they depend on alpha...
        for (unsigned int i = 2; i<8; i=i+2)
        {
            //Loop for quadrature points
            for(unsigned int j = 0; j<n_quad; ++j)
                drates_du[i][j] = ((-1.0)/(K*T[j]))*intrinsic_rates[i][j]*(-1.0)*alpha;
            
        } 
        
        //whereas the odd numbered rates (i.e. the backward rates) depend on -(1-alpha)
        for (unsigned int i = 3; i<8; i=i+2)
        {
            //Loop for quadrature points
            for(unsigned int j = 0; j<n_quad; ++j)
                drates_du[i][j] = ((-1.0)/(K*T[j]))*intrinsic_rates[i][j]*(1.0-alpha);
        }
        
        //Loop over the first two rates and set the values to zero. 
        for (unsigned int i = 0; i<2; ++i)
        {
            //Loop for quadrature points
            for(unsigned int j = 0; j<n_quad; ++j)
                drates_du[i][j] = 0.0;
        }     
    }//end derivative protonic
    
    else if (flag == temperature_of_REV)
    {
        //Loop over all the rates, bar the first two as the forward and backward reactions for the 
        //disassociative adsorption does not depend on the potential. In the first loop we go over
        //the even number rates (i.e. the forward rates) as they depend on alpha...
        for (unsigned int i = 0; i<8; ++i)
        {
            for (unsigned int j=0; j<n_quad; ++j)
            {
                double deltagDT;
                if ( (i==0) || (i==1) )
                    deltagDT = 0.0;
                else if ( (i==2) || (i==4) || (i==6) )
                    deltagDT = (-1.0)*alpha*(catalyst->dvoltage_cell_th_dT(T[j]));
                else if ( (i==3) || (i==5) || (i==7) )
                    deltagDT = (1.0-alpha)*(catalyst->dvoltage_cell_th_dT(T[j]));
                
                drates_du[i][j] = intrinsic_rates[i][j]*((-1.0)/(K))*((deltagDT*T[j] - free_energies[i][j])/(T[j]*T[j]));
            }
        }     
    }//end derivative temperature
}

//---------------------------------------------------------------------------
void
NAME::DoubleTrapKinetics::derivative_current (std::map< VariableNames, std::vector<double> >& dcoef_du)
{
    Assert( derivative_flags.size() != 0, ExcMessage("Derivative flags are not set using set_derivative_flags method before the DoubleTrapKinetics::derivative_current.") );
    Assert( std::find(derivative_flags.begin(), derivative_flags.end(), oxygen_concentration) == derivative_flags.end(),
            ExcMessage("Derivative w.r.t. oxygen concentration is not currently implemented, instead it returns the derivative w.r.t. oxygen molar fraction.") );
    Assert( electrolyte->get_H_O2() != 0.0, ExcMessage("Derivatives at the moment are defined only for Thin fim (Henry's law) case."));

    compute_energies();  
    
    //Compute the derivatives of the OH coverage wrt the solution
    std::map< VariableNames, std::vector<double> > d_OH_theta_du;
    d_OH_coverage_du(d_OH_theta_du);
    
    std::map< VariableNames, std::vector<double> > d_O_theta_du;
    d_O_coverage_du(d_O_theta_du);
    
    // Loop over the flags
    for (unsigned int i=0; i<derivative_flags.size(); ++i)
    {
        std::vector<double> dcurrent(n_quad, 0.);         // Defaulting all derivatives to 0.0

        if (derivative_flags[i] == oxygen_molar_fraction) //#######Add Ch constant
        {
            for (unsigned int j=0; j<n_quad; ++j)
            {
                double corrector = negative_concentration_correction(j);
                double Ch = proton_concentration_ratio(j);

                dcurrent[j] = corrector*(Ch*prefactor*intrinsic_rates[6][j]*d_OH_theta_du[oxygen_molar_fraction][j] + prefactor*intrinsic_rates[7][j]*d_O_theta_du[oxygen_molar_fraction][j] + prefactor*intrinsic_rates[7][j]*d_OH_theta_du[oxygen_molar_fraction][j]);
            }
            
            dcoef_du[oxygen_molar_fraction] = dcurrent;
        }//x_o2
        
        if (derivative_flags[i] == proton_concentration)
        {
            for (unsigned int j=0; j<n_quad; ++j)
            {
                double corrector = negative_concentration_correction(j);                
                double Ch = proton_concentration_ratio(j);
                
                dcurrent[j] = corrector*((1/c_ref_protons)*prefactor*intrinsic_rates[6][j]*theta_OH[j] + Ch*prefactor*intrinsic_rates[6][j]*d_OH_theta_du[proton_concentration][j] +  prefactor*intrinsic_rates[7][j]*(d_O_theta_du[proton_concentration][j] + d_OH_theta_du[proton_concentration][j]));
            }

            dcoef_du[proton_concentration] = dcurrent;
        }//Cp

        else if (derivative_flags[i] == electronic_electrical_potential) //#######Add Ch constant
        {
            for (unsigned int j=0; j<n_quad; ++j)
            {
                double corrector = negative_concentration_correction(j);                
                double Ch = proton_concentration_ratio(j);
                
                dcurrent[j] = corrector*(Ch*prefactor*intrinsic_rates[6][j]*d_OH_theta_du[electronic_electrical_potential][j] +prefactor*theta_OH[j]*(-1./(K*T[j]))*Ch*intrinsic_rates[6][j]*alpha
                    - prefactor*(-1./(K*T[j]))*intrinsic_rates[7][j]*-(1.-alpha) 
                    + prefactor*(-1./(K*T[j]))*intrinsic_rates[7][j]*-(1.-alpha)*theta_O[j] + prefactor*intrinsic_rates[7][j]*d_O_theta_du[electronic_electrical_potential][j]
                    + prefactor*(-1./(K*T[j]))*intrinsic_rates[7][j]*-(1.-alpha)*theta_OH[j] + prefactor*intrinsic_rates[7][j]*d_OH_theta_du[electronic_electrical_potential][j]);
            }
            
            dcoef_du[electronic_electrical_potential] = dcurrent;
        }//phi_s
        
        else if (derivative_flags[i] == protonic_electrical_potential)//#######Add Ch constant
        {
            for (unsigned int j=0; j<n_quad; ++j)
            {
                double corrector = negative_concentration_correction(j);                
                double Ch = proton_concentration_ratio(j);
                
                dcurrent[j] = corrector*(Ch*prefactor*intrinsic_rates[6][j]*d_OH_theta_du[protonic_electrical_potential][j] +prefactor*theta_OH[j]*(-1./(K*T[j]))*Ch*intrinsic_rates[6][j]*-alpha
                    - prefactor*(-1./(K*T[j]))*intrinsic_rates[7][j]*(1.-alpha) 
                    + prefactor*(-1./(K*T[j]))*intrinsic_rates[7][j]*(1.-alpha)*theta_O[j] + prefactor*intrinsic_rates[7][j]*d_O_theta_du[protonic_electrical_potential][j]
                    + prefactor*(-1./(K*T[j]))*intrinsic_rates[7][j]*(1.-alpha)*theta_OH[j] + prefactor*intrinsic_rates[7][j]*d_OH_theta_du[protonic_electrical_potential][j]); 

            }
            
            dcoef_du[protonic_electrical_potential] = dcurrent;
        }//phi_m
        
        else if (derivative_flags[i] == temperature_of_REV) //#######Add Ch constant
        {
            std::vector< std::vector<double> > dratesDu(8, std::vector<double>(n_quad, 0.0));
            drate_du(dratesDu, derivative_flags[i]);
            
            for (unsigned int j=0; j<n_quad; ++j)
            {
                double corrector = negative_concentration_correction(j);                
                double Ch = proton_concentration_ratio(j);
                
                dcurrent[j] = corrector*(( Ch*dratesDu[6][j]*theta_OH[j] + Ch*intrinsic_rates[6][j]*d_OH_theta_du[temperature_of_REV][j] -
                       dratesDu[7][j]*(1.0-theta_O[j]-theta_OH[j]) + intrinsic_rates[7][j]*(d_O_theta_du[temperature_of_REV][j]+d_OH_theta_du[temperature_of_REV][j]) )*prefactor);
 
            }
            
            dcoef_du[temperature_of_REV] = dcurrent;
        }//T
        
        else
        {
            dcoef_du[ derivative_flags[i] ] = dcurrent;
        }
    }//end of derivative_flags loop  
}

// ----------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------
// ---------------This is to test derivative_current method using numerical differentiation -----------------------------
//              dcoef_du.clear();
//              dcoef_du.resize(this->derivative_flags.size(), std::vector<double>(n_quad, 0.0));
//              
//              //size checking for vectors
//              
//              
//              unsigned int location = 9;
//              
//              for (unsigned int k=0; k<species_names.size(); ++k)
//              {
//                      //Find the location of oxygen molar fraction in the species_names list
//                      if (species_names[k] == "Oxygen molar fraction")
//                      {
//                              location = k;
//                      }
//              }
//              
//              //Loop over the flags
//              for (unsigned int i=0; i<this->derivative_flags.size(); ++i)
//              {
//                      if (this->derivative_flags[i] == "Oxygen molar fraction")
//                      {
//                              //Compute the current using the current solution.
//                              std::vector<double> coef1(std::vector<double>(dcoef_du[i].size(), 0.0));
//                              //std::vector<double> theta_O_1(std::vector<double>(dcoef_du[i].size(), 0.0));
//                              //std::vector<double> theta_OH_1(std::vector<double>(dcoef_du[i].size(), 0.0));
//                              dcoef_du.resize(this->derivative_flags.size(), std::vector<double>(dcoef_du[i].size(), 0.0));
// 
//                              current_density(coef1);
//                              //theta_O_1 = theta_O;
//                              //theta_OH_1 = theta_OH;
//                              double delta = std::pow(10.0, -7);
//                              //Modify solution
//                              for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                              {
//                                      double xO2 = (species_solution[location][j]*this->electrolyte->get_H_O2())/p_total;
//                                      double xO2_new = xO2 + delta;
//                                      species_solution[location][j] = xO2_new*(p_total/this->electrolyte->get_H_O2());
//                              }
//                              
//                              //Compute the new current
//                              std::vector<double> coef2(std::vector<double>(dcoef_du[i].size(), 0.0));
//                              current_density(coef2);
//                              
//                              //Compute the coefficient
//                              for (unsigned int j=0; j<dcoef_du[0].size(); ++j)
//                              {
//                                      dcoef_du[i][j] = (coef2[j] - coef1[j])/delta;
//                                      //FcstUtilities::log << "Numerical grad O oxygen " << (theta_O[j] - theta_O_1[j])/delta << std::endl;
//                                      //FcstUtilities::log << "Numerical grad OH oxygen " << (theta_OH[j] - theta_OH_1[j])/delta << std::endl;
//                                      //FcstUtilities::log << "Numerical grad oxygen " << dcoef_du[i][j] << std::endl;
//                                      //FcstUtilities::log << "numerical oxygen: " << dcoef_du[i][j] << std::endl;
//                              }
//                              
//                              for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                              {
//                                      double xO2_new = (species_solution[location][j]*this->electrolyte->get_H_O2())/p_total;
//                                      double xO2 = xO2_new - delta;
//                                      species_solution[location][j] = xO2*(p_total/this->electrolyte->get_H_O2());
//                              }      
//                      }//x_o2
//                      
//                      
//                      else if (this->derivative_flags[i] == "Electronic electrical potential")
//                      { 
//                              //Compute the current using the current solution.
//                              std::vector<double> coef1(std::vector<double>(dcoef_du[i].size(), 0.0));
// //                           std::vector<double> theta_O_1(std::vector<double>(dcoef_du[i].size(), 0.0));
// //                           std::vector<double> theta_OH_1(std::vector<double>(dcoef_du[i].size(), 0.0));
//                              current_density(coef1);
// //                           theta_O_1 = theta_O;
// //                           theta_OH_1 = theta_OH;
//                              double delta = std::pow(10.0, -7);
//                              
//                              //Modify solution
//                              for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                              {
//                                      phi_s[j] += delta;
//                              }
//                              
//                              //Compute the new current
//                              std::vector<double> coef2(std::vector<double>(dcoef_du[i].size(), 0.0));
//                              current_density(coef2);
//                              
//                              //Compute the coefficient
//                              for (unsigned int j=0; j<dcoef_du[0].size(); ++j)
//                              {
//                                      dcoef_du[i][j] = (coef2[j] - coef1[j])/delta;
//                                      //FcstUtilities::log << "Numerical grad electron " << dcoef_du[i][j] << std::endl;
// //                                   FcstUtilities::log << "Numerical grad O electron " << (theta_O_1[j] - theta_O[j])/delta << std::endl;
// //                                   FcstUtilities::log << "Numerical grad OH electron " << (theta_OH_1[j] - theta_OH[j])/delta << std::endl;
//                                      //FcstUtilities::log << "numerical electron: " << dcoef_du[i][j] << std::endl;
//                              }
//                              
//                              //Restore the solution to the correct value.
//                              for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                              {
//                                      phi_s[j] -= delta;
//                              }        
//                              
//                      }//phi_s
//                      
//                      else if (this->derivative_flags[i] == "Protonic electrical potential")
//                      {
//                              //Compute the current using the current solution.
//                              std::vector<double> coef1(std::vector<double>(dcoef_du[i].size(), 0.0));
//                              current_density(coef1);
//                              
//                              double delta = std::pow(10.0, -7);
//                              
//                              //Modify solution
//                              for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                              {
//                                      phi_m[j] += delta;
//                              }
//                              
//                              //Compute the new current
//                              std::vector<double> coef2(std::vector<double>(dcoef_du[i].size(), 0.0));
//                              current_density(coef2);
//                              
//                              //Compute the coefficient
//                              for (unsigned int j=0; j<dcoef_du[0].size(); ++j)
//                              {
//                                      dcoef_du[i][j] = (coef2[j] - coef1[j])/delta;
//                                      //FcstUtilities::log << "Numerical grad proton " << dcoef_du[i][j] << std::endl;
//                                      //FcstUtilities::log << "numerical proton: " << dcoef_du[i][j] << std::endl;
//                              }
//                              for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                              {
//                                      phi_m[j] -= delta;
//                              }       
//                      }//phi_m
//                      
//                      else if (this->derivative_flags[i] == "Temperature")
//                      {
//                              //Compute the current using the current solution.
//                              std::vector<double> coef1(std::vector<double>(dcoef_du[i].size(), 0.0));
//                              std::vector<double> theta_O_1(std::vector<double>(dcoef_du[i].size(), 0.0));
//                              std::vector<double> theta_OH_1(std::vector<double>(dcoef_du[i].size(), 0.0));
//                              current_density(coef1);
//                              theta_O_1 = theta_O;
//                              theta_OH_1 = theta_OH;
//                              
//                              double delta = std::pow(10.0, -6);
//                              
//                              //Modify solution
//                              for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                              {
//                                      T[j] += delta;
//                              }
//                              
//                              //Compute the new current
//                              std::vector<double> coef2(std::vector<double>(dcoef_du[i].size(), 0.0));
//                              current_density(coef2);
//                              
//                              //Compute the coefficient
//                              for (unsigned int j=0; j<dcoef_du[0].size(); ++j)
//                              {
//                                      dcoef_du[i][j] = (coef2[j] - coef1[j])/delta;
//                                      //FcstUtilities::log << "Numerical grad O temperature " << (theta_O[j] - theta_O_1[j])/delta << std::endl;
//                                      //FcstUtilities::log << "Numerical grad OH electron " << (theta_OH[j] - theta_OH_1[j])/delta << std::endl;
//                                      FcstUtilities::log << "numerical temperature: " << dcoef_du[i][j] << std::endl;
//                              }
//                              
//                              for (unsigned int j=0; j<dcoef_du[i].size(); ++j)
//                              {
//                                      T[j] -= delta;
//                              }       
//                      }//T
//                      
//                      else
//                      {
//                              //FcstUtilities::log << "Derivative flag " << this->derivative_flags[i] << " not implemented." << std::endl;
//                      }//everything else
//                      
//              }//end of derivative_flags loop
//      }
