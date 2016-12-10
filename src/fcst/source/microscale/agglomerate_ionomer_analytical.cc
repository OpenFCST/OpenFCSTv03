// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: agglomerate_ionomer_analytical.cc
// - Description: Used to solve a system of equations representing a spherical ionomer-filled agglomerate.
// - Developers: Peter Dobson <pdobson@ualberta.ca>
//               Marc Secanell Gallart, University of Alberta
// - $Id: agglomerate_ionomer_analytical.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <microscale/agglomerate_ionomer_analytical.h>

namespace NAME = FuelCellShop::MicroScale;

const std::string NAME::IonomerAgglomerateAnalytical::concrete_name("IonomerAgglomerateAnalytical");
NAME::IonomerAgglomerateAnalytical const* NAME::IonomerAgglomerateAnalytical::PROTOTYPE = new NAME::IonomerAgglomerateAnalytical(concrete_name);


NAME::IonomerAgglomerateAnalytical::IonomerAgglomerateAnalytical(std::string concrete_name){
    //FcstUtilities::log <<" Register " + concrete_name +  " to FactoryMap"<<std::endl;
       this->get_mapFactory()->insert(
               std::pair<std::string, NAME::IonomerAgglomerateAnalytical*>(
                       concrete_name, this));
}



//---------------------------------------------------------------------------
NAME::IonomerAgglomerateAnalytical::IonomerAgglomerateAnalytical()
{
    this->has_derivatives_ = true;
    checked_kinetics = false;
}

//---------------------------------------------------------------------------
void
NAME::IonomerAgglomerateAnalytical::check_kinetics(){

    bool error_state = false;
    std::string msg;

    //Check kinetics type
    if(typeid(*kinetics.get()).name() != typeid(FuelCellShop::Kinetics::TafelKinetics).name()
            and typeid(*kinetics.get()).name() != typeid(FuelCellShop::Kinetics::DualPathKinetics).name())
    {
        msg = "Incorrect kinetics model chosen for IonomerAgglomerateAnalytical";
	msg = msg + "\n\tSelected kinetics: " + typeid(*kinetics.get()).name();
	msg = msg + "\n\tExpected kinetics: " + typeid(FuelCellShop::Kinetics::TafelKinetics).name();
	msg = msg + " or " + typeid(FuelCellShop::Kinetics::DualPathKinetics).name();
        error_state = true;
    }
    else
    {
        //Else check ORR order

        if(this->reactant == oxygen_molar_fraction){
            std::vector<VariableNames> names;
            names.push_back(oxygen_concentration);

            std::map<VariableNames, double> gamma_map;
            boost::shared_ptr<FuelCellShop::Material::CatalystBase> catalyst;
            catalyst = this->layer->get_resource<FuelCellShop::Material::CatalystBase>();
            catalyst->reaction_order( names,  gamma_map);
            double reaction_order = gamma_map[oxygen_concentration];

            if (reaction_order < 0.999999 || reaction_order > 1.000001) {
                msg = "Incorrect ORR order chosen for IonomerAgglomerateAnalytical";
                error_state = true;
            }
        }

    }


    if(error_state){

        throw std::runtime_error(msg);
    }

    checked_kinetics = true;
}

//---------------------------------------------------------------------------
void
NAME::IonomerAgglomerateAnalytical::set_structure()
{


    CL_Properties props = this->layer->get_properties();

    electrolyte = this->layer->get_resource<FuelCellShop::Material::PolymerElectrolyteBase>();
    kinetics = this->layer->get_resource<FuelCellShop::Kinetics::BaseKinetics>();
    AV = props[CLPropNames::active_area_scaled]; // convert to 1/m
    P= props[CLPropNames::pressure];

    r_agg*=1e-7;// convert to cm
    delta_agg *=1e-7; // convert to cm

    interface = r_agg/(r_agg+delta_agg);//set to 1.0 for no thin film


}

//---------------------------------------------------------------------------
FuelCellShop::SolutionMap
NAME::IonomerAgglomerateAnalytical::compute_current ()
{

    if(not checked_kinetics)
        check_kinetics();

    //True pressure is now available from the layer.
    //Was 0 at set structure time.
    P= this->layer->get_properties()[CLPropNames::pressure];

    double E_r = 0.0;
    double I = 0.0;
    if(this->solutions.at(this->reactant)[this->sol_index]<= 0.0){
        SolutionMap sols;
        sols.push_back(SolutionVariable(0, 1, current_density));
        sols.push_back(SolutionVariable(0,1, CL_effectiveness));

        return sols; //To stabalize FEM solution
    }

        
        //Get the electrolyte effective properties for the current temperature
        electrolyte->set_T(this->solutions.at(temperature_of_REV)[this->sol_index]);
    electrolyte->set_lambda(this->solutions.at(membrane_water_content)[this->sol_index]);
    
    double molarFactor;
    if (this->solutions.at(this->reactant).get_variablename() == oxygen_molar_fraction)
    {
        electrolyte->oxygen_diffusivity(D_R_N);
        H_R_N = electrolyte->get_H_O2();
        tempReactantName = oxygen_concentration;
        molarFactor = 4.0;
    }
    else if(this->solutions.at(this->reactant).get_variablename() == hydrogen_molar_fraction)
    {
        
        //AssertThrow(false, ExcMessage("Hydrogen reaction not yet implemented in numerical ionomer agglomerate."));
        electrolyte->hydrogen_diffusivity(D_R_N);
        H_R_N = electrolyte->get_H_H2();
        tempReactantName = hydrogen_concentration;
        molarFactor = 2.0;
        
    }
    
    double D_eff = pow(epsilon_agg,1.5) * D_R_N;
       
    // Set the oxygen concentration and membrane potential at the boundary
    // Set the solid phase potential for the domain
    c_R = (this->solutions.at(this->reactant)[this->sol_index]*P)/H_R_N;
    phi_S = this->solutions.at(electronic_electrical_potential)[this->sol_index];
    phi_M = this->solutions.at(protonic_electrical_potential)[this->sol_index];
    
    double unit_size = 1;
    double agg_size = interface/unit_size;
    double volume = (4.0/3.0)*pi*pow(unit_size,3.0);
    double agg_volume = (4.0/3.0)*pi*pow(agg_size,3.0);
    
    //get the solution at the quadrature points
    std::vector<double> r_rate(1,0.0);
    
    std::vector<SolutionVariable> c_reactants;
    c_reactants.push_back(SolutionVariable (c_R,1, tempReactantName));
    SolutionVariable v_membrane(phi_M, 1, protonic_electrical_potential);
    SolutionVariable v_solid(phi_S, 1, electronic_electrical_potential);
    SolutionVariable temp_(this->solutions.at(temperature_of_REV)[this->sol_index], 1, temperature_of_REV);
    
    this->kinetics->set_reactant_concentrations(c_reactants);
    this->kinetics->set_electrolyte_potential(v_membrane);
    this->kinetics->set_solid_potential(v_solid);
    this->kinetics->set_temperature(temp_);
    this->kinetics->set_p_t(P);
    
    //this->kinetics->molar_reaction_rate(r_rate, "Oxygen molar fraction");
    this->kinetics->current_density(r_rate);
    r_rate[0]/=c_R;
    double k_c = AV*r_rate[0] / (molarFactor*F);
    k_c *= volume/agg_volume;
    double E_agg = compute_Er(k_c, D_eff);
    
    //compute the concentration at the boundary
    double cR_boundary= c_R * pow(1.0 + ((pow(r_agg,2)*delta_agg*E_agg*k_c)/(3*(r_agg + delta_agg)*D_R_N)), -1);
    
    //Analytical expression for oxygen concentration across domain
    /*co2_final.clear();
    co2_final.resize(mesh_final.size(),0.0);
    //Compute oxygen concentration profile across the agglomerate
    double temp = sqrt(k_c/D_eff);
    
    for(int i =0; i<mesh_final.size();i++)
    {
        
        co2_final[i]=cR_boundary*(r_agg*sinh(r_agg*mesh_final[i]*temp))/(r_agg*mesh_final[i]*sinh(r_agg*temp));
    }*/
    
    I = (molarFactor*F) * c_R * pow(1/(E_agg*k_c) + ((pow(r_agg,2)*delta_agg)/(3*(r_agg + delta_agg)*D_R_N)), -1);
    double I_max = (molarFactor*F) * c_R * k_c * agg_volume/volume;
    
    double I_avg = I*agg_volume/volume;
    E_r = I_avg/I_max;
    
    if (std::isnan(E_r))
        E_r = 0.0;
    if (std::isnan(I_avg))
        I_avg = 0.0;
    




    SolutionMap sols;
    sols.push_back(SolutionVariable(I_avg, 1, current_density));
    sols.push_back(SolutionVariable(E_r,1, CL_effectiveness));

    if(this->kinetics->has_coverage(OH_coverage)){
        std::vector<double> OH_c;
        this->kinetics->OH_coverage(OH_c);
        sols.push_back(SolutionVariable(OH_c,OH_coverage));
    }

    if(this->kinetics->has_coverage(O_coverage)){
        std::vector<double> O_c;
        this->kinetics->OH_coverage(O_c);
        sols.push_back(SolutionVariable(O_c,O_coverage));
    }

    return sols;
    
}


//---------------------------------------------------------------------------
std::vector<double>
NAME::IonomerAgglomerateAnalytical::compute_derivative_current ()
{
    std::vector<double> dI(3, 0.0);
    if(this->solutions.at(this->reactant)[this->sol_index]< 0.0 < 0.0)
        return dI; // Quit on non-real values
        
        //Get the electrolyte effective properties for the current temperature
        electrolyte->set_T(this->solutions.at(temperature_of_REV)[this->sol_index]);
    electrolyte->set_lambda(this->solutions.at(membrane_water_content)[this->sol_index]);
    
    double molarFactor;
    if (this->solutions.at(this->reactant).get_variablename() == oxygen_molar_fraction)
    {
        electrolyte->oxygen_diffusivity(D_R_N);
        H_R_N = electrolyte->get_H_O2();
        tempReactantName = oxygen_concentration;
        molarFactor = 4.0;
    }
    else if(this->solutions.at(this->reactant).get_variablename() == hydrogen_molar_fraction)
    {
        
        //AssertThrow(false, ExcMessage("Hydrogen reaction not yet implemented in numerical ionomer agglomerate."));
        electrolyte->hydrogen_diffusivity(D_R_N);
        H_R_N = electrolyte->get_H_H2();
        tempReactantName = hydrogen_concentration;
        molarFactor = 2.0;
        
    }
    
    double D_eff = pow(epsilon_agg,1.5) * D_R_N;
    
    // Set the oxygen concentration and membrane potential at the boundary
    // Set the solid phase potential for the domain
    c_R = (this->solutions.at(this->reactant)[this->sol_index]*P)/H_R_N;
    phi_S = this->solutions.at(electronic_electrical_potential)[this->sol_index];
    phi_M = this->solutions.at(protonic_electrical_potential)[this->sol_index];
    
    
    double unit_size = 1;
    double agg_size = interface/unit_size;
    double volume = (4.0/3.0)*pi*pow(unit_size,3.0);
    double agg_volume = (4.0/3.0)*pi*pow(agg_size,3.0);
    
    //get the solution at the quadrature points
    std::vector<double> r_rate(1,0.0);
    
    
    std::vector<SolutionVariable> c_reactants;
    c_reactants.push_back(SolutionVariable (c_R,1,tempReactantName));
    SolutionVariable v_membrane(phi_M, 1, protonic_electrical_potential);
    SolutionVariable v_solid(phi_S, 1, electronic_electrical_potential);
    SolutionVariable temp_(this->solutions.at(temperature_of_REV)[this->sol_index], 1, temperature_of_REV);
    
    this->kinetics->set_reactant_concentrations(c_reactants);
    this->kinetics->set_electrolyte_potential(v_membrane);
    this->kinetics->set_solid_potential(v_solid);
    this->kinetics->set_temperature(temp_);
    this->kinetics->set_p_t(P);
    this->kinetics->current_density(r_rate);
    r_rate[0]/=c_R;
    
    std::map< VariableNames, std::vector<double> > derivatives;
    this->kinetics->set_derivative_flags(this->sol_names);
    this->kinetics->derivative_current(derivatives);
    
    double k_c = AV*r_rate[0] /(molarFactor*F);
    k_c *= volume/agg_volume;
    double dkc_dphi_s = AV*derivatives[electronic_electrical_potential][0] / (c_R*(molarFactor*F));
    double dkc_dphi_m = AV*derivatives[protonic_electrical_potential][0] / (c_R*(molarFactor*F));
    double E_agg = compute_Er(k_c, D_eff);
    double dE_agg = 0.0;
    
    //dRR_dxo2
    dI[0] = (molarFactor*F) * (P/H_R_N) * pow(1/(E_agg*k_c) + ((pow(r_agg,2)*delta_agg)/(3*(r_agg + delta_agg)*D_R_N)), -1);
    
    //dRR_dphi_m 
    dE_agg = compute_dEr(k_c, dkc_dphi_m, D_eff);
    dI[1] = (molarFactor*F) * c_R * pow(1/(E_agg*k_c) + ((pow(r_agg,2)*delta_agg)/(3*(r_agg + delta_agg)*D_R_N)), -2)
            * ((dE_agg*k_c + E_agg*dkc_dphi_m)/pow((E_agg*k_c),2.0));
    //dRR_dphi_s
    dE_agg = compute_dEr(k_c, dkc_dphi_s, D_eff);
    dI[2] = (molarFactor*F) * c_R * pow(1/(E_agg*k_c) + ((pow(r_agg,2)*delta_agg)/(3*(r_agg + delta_agg)*D_R_N)), -2)
            * ((dE_agg*k_c + E_agg*dkc_dphi_s)/pow((E_agg*k_c),2.0));
    
    
    for(unsigned int i=0;i<dI.size();++i)
        dI[i] *= agg_volume/volume;
    
    return dI;
}


//---------------------------------------------------------------------------
double
NAME::IonomerAgglomerateAnalytical::compute_Er (const double k_c, const double D)
{
    double phi_L = (r_agg / 3.0)*sqrt(k_c/D);
    double E_r = (1/phi_L)*(1/tanh(3*phi_L) - 1/(3*phi_L));  // Spherical agglomerate
    return E_r;
}


//---------------------------------------------------------------------------
double
NAME::IonomerAgglomerateAnalytical::compute_dEr (const double k_c, const double dk_c, const double D)
{
    //-- Compute phi_L
    double phi_L = (r_agg/3.0)*sqrt(k_c/D);
    double dphiL_dkc = (r_agg / 3.0)*pow(k_c/D,-0.5)*(0.5/D);
    
    //FcstUtilities::log << "phi_L = " << phi_L << " dphiL_dkc = " << dphiL_dkc << std::endl;
    //-- Dependance of E_r on phi_L
    double  dEr_dphiL = -1/(pow(phi_L,2.0)*tanh(3*phi_L)) - 3*(1-pow(tanh(3*phi_L),2.0))/(phi_L*pow(tanh(3*phi_L),2.0)) + 2.0/(3*pow(phi_L,3.0));
    
    double dE_r = dEr_dphiL*dphiL_dkc*dk_c;
    
    return dE_r;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
