//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: agglomerate_ionomer_sun.cc
//    - Description: Used to solve a system of equations representing a
//              spherical ionomer-filled agglomerate.
//    - Developers: M. Secanell and Phil Wardlaw
//    - $Id: agglomerate_ionomer_sun.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <microscale/agglomerate_ionomer_sun.h>

namespace NAME = FuelCellShop::MicroScale;

const std::string NAME::IonomerAgglomerateSun::concrete_name("IonomerAgglomerateSun");
NAME::IonomerAgglomerateSun const* NAME::IonomerAgglomerateSun::PROTOTYPE = new NAME::IonomerAgglomerateSun(concrete_name);


NAME::IonomerAgglomerateSun::IonomerAgglomerateSun(std::string concrete_name){
    //FcstUtilities::log <<" Register " + concrete_name +  " to FactoryMap"<<std::endl;
       this->get_mapFactory()->insert(
               std::pair<std::string, NAME::IonomerAgglomerateSun*>(
                       concrete_name, this));
}


//---------------------------------------------------------------------------
NAME::IonomerAgglomerateSun::IonomerAgglomerateSun()
{
	this->has_derivatives_ = true;
	checked_kinetics = false;
}

void
NAME::IonomerAgglomerateSun::check_kinetics(){

    bool error_state = false;
    std::string msg;

    //Check kinetics type
    if(typeid(*kinetics.get()).name() != typeid(FuelCellShop::Kinetics::TafelKinetics).name()
            and typeid(*kinetics.get()).name() != typeid(FuelCellShop::Kinetics::DualPathKinetics).name())
    {
        msg = "Incorrect kinetics model chosen for IonomerAgglomerateAnalytical";
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
NAME::IonomerAgglomerateSun::set_structure()
{
    CL_Properties props = this->layer->get_properties();

    electrolyte = this->layer->get_resource<FuelCellShop::Material::PolymerElectrolyteBase>();
    kinetics = this->layer->get_resource<FuelCellShop::Kinetics::BaseKinetics>();

    epsilon_V = props[CLPropNames::void_fraction];
    AV = props[CLPropNames::active_area_scaled]*(1.0 - epsilon_V); //Since Agglomerate_CL divides this by this factor at line 278
    r_agg *=1e-7;// convert to cm
    delta_agg *=1e-7; // convert to cm

    interface = 0;
}


//---------------------------------------------------------------------------
FuelCellShop::SolutionMap
NAME::IonomerAgglomerateSun::compute_current ()
{

    if(not checked_kinetics)
           check_kinetics();


    //True pressure is now available from the layer.
    //Was 0 at set structure time.
    P= this->layer->get_properties()[CLPropNames::pressure];


    double E_r = 0.0;
    double I = 0.0;
    if(this->solutions[this->reactant][this->sol_index]<= 0.0){
        SolutionMap sols;
        sols.push_back(SolutionVariable(0, 1, current_density));
        sols.push_back(SolutionVariable(0,1, CL_effectiveness));

        return sols; //To stabalize FEM solution
    }

        
    //Get the electrolyte effective properties for the current temperature
    electrolyte->set_T(this->solutions[temperature_of_REV][this->sol_index]);
    electrolyte->set_lambda(this->solutions[membrane_water_content][this->sol_index]);
    
    double molarFactor;
    if (this->solutions[this->reactant].get_variablename() == oxygen_molar_fraction)
    {
        electrolyte->oxygen_diffusivity(D_R_N);
        H_R_N = electrolyte->get_H_O2();
        tempReactantName = oxygen_concentration;
        molarFactor = 4.0;
    }
    else if(this->solutions[this->reactant].get_variablename() == hydrogen_molar_fraction)
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
    c_R = (this->solutions[this->reactant][this->sol_index]*P)/H_R_N;
    phi_S = this->solutions[electronic_electrical_potential][this->sol_index];
    phi_M = this->solutions[protonic_electrical_potential][this->sol_index];
    
    //get the solution at the quadrature points
    std::vector<double> r_rate(1,0.0);
    
    std::vector<SolutionVariable> c_reactants;
    c_reactants.push_back(SolutionVariable (c_R,1, tempReactantName));
    SolutionVariable v_membrane(phi_M, 1, protonic_electrical_potential);
    SolutionVariable v_solid(phi_S, 1, electronic_electrical_potential);
    SolutionVariable temp_(this->solutions[temperature_of_REV][this->sol_index], 1, temperature_of_REV);
    
    this->kinetics->set_reactant_concentrations(c_reactants);
    this->kinetics->set_electrolyte_potential(v_membrane);
    this->kinetics->set_solid_potential(v_solid);
    this->kinetics->set_temperature(temp_);
    this->kinetics->set_p_t(P);
    //this->kinetics->molar_reaction_rate(r_rate, "Oxygen molar fraction");
    this->kinetics->current_density(r_rate);
    r_rate[0]/=c_R;
    double k_c = AV*r_rate[0] / (molarFactor*F*(1-epsilon_V));
    double E_agg = compute_Er(k_c, D_eff);
       
    double a_agg = n_agg*4*pi*pow(r_agg + delta_agg, 2.0)*epsilon_V;
    
    I = (molarFactor*F) * c_R * pow(1/(E_agg*k_c*(1-epsilon_V)) + (((r_agg + delta_agg)*delta_agg)/(a_agg*r_agg*D_R_N)), -1);

    
    SolutionMap sols;
    sols.push_back(SolutionVariable(I/(1.0-epsilon_V), 1, current_density));//Since Agglomerate_CL multiplies by this factor at line 398


	//Find effectiveness by normalizing I by the current produced by a similar sphere with perfect reactant transport
	double Er = I/(r_rate[0]*c_R*(4/3)*pi*std::pow(r_agg, 3.0)); //r_rate[0] is multiplied by c_R to undo previous normalization
    sols.push_back(SolutionVariable(Er,1, CL_effectiveness));


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
NAME::IonomerAgglomerateSun::compute_derivative_current ()
{
    std::vector<double> dI(3);
    if(this->solutions[this->reactant][this->sol_index]< 0.0 < 0.0)
        return dI; // Quit on non-real values
        
    //Get the electrolyte effective properties for the current temperature
    electrolyte->set_T(this->solutions[temperature_of_REV][this->sol_index]);
    electrolyte->set_lambda(this->solutions[membrane_water_content][this->sol_index]);
    
    double molarFactor;
    if (this->solutions[this->reactant].get_variablename() == oxygen_molar_fraction)
    {
        electrolyte->oxygen_diffusivity(D_R_N);
        H_R_N = electrolyte->get_H_O2();
        tempReactantName = oxygen_concentration;
        molarFactor = 4.0;
    }
    else if(this->solutions[this->reactant].get_variablename() == hydrogen_molar_fraction)
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
    c_R = (this->solutions[this->reactant][this->sol_index]*P)/H_R_N;
    phi_S = this->solutions[electronic_electrical_potential][this->sol_index];
    phi_M = this->solutions[protonic_electrical_potential][this->sol_index];
    
    //get the solution at the quadrature points
    std::vector<double> r_rate(1,0.0);
    
    
    std::vector<SolutionVariable> c_reactants;
    c_reactants.push_back(SolutionVariable (c_R,1,tempReactantName));
    SolutionVariable v_membrane(phi_M, 1, protonic_electrical_potential);
    SolutionVariable v_solid(phi_S, 1, electronic_electrical_potential);
    SolutionVariable temp_(this->solutions[temperature_of_REV][this->sol_index], 1, temperature_of_REV);
    
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
    
    double k_c = AV*r_rate[0] /(molarFactor*F*(1-epsilon_V));
    double dkc_dphi_s = AV*derivatives[electronic_electrical_potential][0] / (c_R*(molarFactor*F*(1-epsilon_V)));
    double dkc_dphi_m = AV*derivatives[protonic_electrical_potential][0] / (c_R*(molarFactor*F*(1-epsilon_V)));
    double E_agg = compute_Er(k_c, D_eff);
    double dE_agg = 0.0;
    
    double a_agg = n_agg*4*pi*pow(r_agg+ delta_agg, 2.0)*epsilon_V;
    
    //Original
    //-- dRR_dxo2
    dI[0] = (molarFactor*F*P)/H_R_N*pow(1/(E_agg*k_c*(1-epsilon_V)) + ((r_agg + delta_agg)*delta_agg)/(a_agg*r_agg*D_R_N), -1);
    
    //-- dRR_dphi_m
    dE_agg = compute_dEr(k_c, dkc_dphi_m, D_eff);
    dI[1] = (molarFactor*F*c_R)*(-1)*pow(1/(E_agg*k_c*(1-epsilon_V)) + ((r_agg + delta_agg)*delta_agg)/(a_agg*r_agg*D_R_N), -2)
    *(-(dE_agg*k_c*(1-epsilon_V) + E_agg*dkc_dphi_m*(1-epsilon_V) )/pow(E_agg*k_c*(1-epsilon_V),2.0));
    
    //-- dRR_dphi_s
    dE_agg = compute_dEr(k_c, dkc_dphi_s, D_eff);
    dI[2] = (molarFactor*F*c_R)*(-1)*pow(1/(E_agg*k_c*(1-epsilon_V)) + ((r_agg + delta_agg)*delta_agg)/(a_agg*r_agg*D_R_N), -2)
    *(-(dE_agg*k_c*(1-epsilon_V) + E_agg*dkc_dphi_s*(1-epsilon_V) )/pow(E_agg*k_c*(1-epsilon_V),2.0));
    
    dI[0] /=(1.0-epsilon_V);dI[1] /=(1.0-epsilon_V);dI[2] /=(1.0-epsilon_V);

    return dI;
}


//---------------------------------------------------------------------------
double
NAME::IonomerAgglomerateSun::compute_Er (const double k_c, const double D)
{
    //Works the same as before
    double phi_L = (r_agg / 3.0)*sqrt(k_c/D);
    double E_r = (1/phi_L)*(1/tanh(3*phi_L) - 1/(3*phi_L));  // Spherical agglomerate
    return E_r;

}


//---------------------------------------------------------------------------
double
NAME::IonomerAgglomerateSun::compute_dEr (const double k_c, const double dk_c, const double D)
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
