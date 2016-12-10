/*
 * DT_test.h
 *
 *  Created on: Feb 11, 2014
 *      Author: wardlawp
 */


#include<DT_test.h>




void
DoubleTrapTest::setup(){



    FuelCellShop::Material::CatalystBase::declare_Catalyst_parameters(param);
    FuelCellShop::Material::PolymerElectrolyteBase::declare_PolymerElectrolyte_parameters(param);
    FuelCellShop::Kinetics::BaseKinetics::declare_Kinetics_parameters(param);

    param.enter_subsection("Materials");
    param.enter_subsection("Platinum");
    param.set("Method for kinetics parameters (ORR)","Double_trap");
    param.leave_subsection();
    param.leave_subsection();


    catalyst = FuelCellShop::Material::CatalystBase::create_Catalyst(param, "Platinum");
    electrolyte = FuelCellShop::Material::PolymerElectrolyteBase::create_PolymerElectrolyte(param, "Nafion");
    kinetics = FuelCellShop::Kinetics::BaseKinetics::create_Kinetics(param, "DoubleTrapKinetics");
    catalyst->set_reaction_kinetics(ORR);

    kinetics->set_catalyst(catalyst.get());
    kinetics->set_electrolyte(electrolyte.get());
    kinetics->set_p_t(101325.0);
    kinetics->set_temperature(SOL(353.0, 1, temperature_of_REV));
    kinetics->set_solid_potential(SOL(0.7, 1, electronic_electrical_potential));
    kinetics->set_electrolyte_potential(SOL(-0.05, 1, protonic_electrical_potential));


    SOLS reacts; reacts.push_back(SOL(1e-5, 1, oxygen_concentration));

    double c_h_ref = 0.001;
    reacts.push_back(SOL(c_h_ref,1,proton_concentration));

    kinetics->set_reactant_concentrations(reacts);



    kinetics->set_reaction_kinetics(ORR);

    std::vector<VariableNames> flag_names;
    flag_names.push_back(oxygen_molar_fraction);
    flag_names.push_back(protonic_electrical_potential);
    flag_names.push_back(electronic_electrical_potential);
    flag_names.push_back(proton_concentration);
    flag_names.push_back(temperature_of_REV);

    kinetics->set_derivative_flags(flag_names);

}

void
DoubleTrapTest::testCurrent(){

    std::vector<double> current;
    kinetics->current_density(current);

    double tolerance = 1e-4;
    double diff;
    double expected_curret = 0.0021476;

    diff = std::abs((current[0] -expected_curret)/current[0]);

    std::string msg = "Calculated current density value is incorrect. Value was " + std::to_string(current[0]);
    TEST_ASSERT_MSG(diff < tolerance, msg.c_str());
}



void
DoubleTrapTest::testDerivatives(){

    double perterb = 1e-8;

    double c_h_ref = 0.001;

    std::map< VariableNames, std::vector<double> > derivatives;
    kinetics->derivative_current(derivatives);

    //Calculate derivatives using central differencing

    std::map<VariableNames, std::vector<double>> forwardPterb;

    kinetics->set_solid_potential(SOL(0.7 + perterb, 1, electronic_electrical_potential));
    kinetics->current_density(forwardPterb[electronic_electrical_potential]);


    kinetics->set_solid_potential(SOL(0.7, 1, electronic_electrical_potential));
    kinetics->set_electrolyte_potential(SOL(-0.05 + perterb, 1, protonic_electrical_potential));
    kinetics->current_density(forwardPterb[protonic_electrical_potential]);


    kinetics->set_electrolyte_potential(SOL(-0.05, 1, protonic_electrical_potential));
    SOLS reacts; reacts.push_back(SOL(1e-5 + perterb, 1, oxygen_concentration));
    kinetics->set_reactant_concentrations(reacts);
    kinetics->current_density(forwardPterb[oxygen_concentration]);

    reacts.clear() ; reacts.push_back(SOL(1e-5, 1, oxygen_concentration));
    reacts.push_back(SOL(c_h_ref+perterb,1,proton_concentration));
    kinetics->set_reactant_concentrations(reacts);
    kinetics->current_density(forwardPterb[proton_concentration]);

    //Custom perturbation for temperature derivative
    kinetics->set_temperature(SOL(353.0 + 1.0, 1, temperature_of_REV));
    kinetics->current_density(forwardPterb[temperature_of_REV]);
    kinetics->set_temperature(SOL(353.0 , 1, temperature_of_REV));

    std::map<VariableNames, std::vector<double>> backPterb;

    reacts.clear() ; reacts.push_back(SOL(1e-5, 1, oxygen_concentration));
    reacts.push_back(SOL(c_h_ref,1,proton_concentration));
    kinetics->set_reactant_concentrations(reacts);
    kinetics->set_solid_potential(SOL(0.7 - perterb, 1, electronic_electrical_potential));
    kinetics->current_density(backPterb[electronic_electrical_potential]);


    kinetics->set_solid_potential(SOL(0.7, 1, electronic_electrical_potential));
    kinetics->set_electrolyte_potential(SOL(-0.05 - perterb, 1, protonic_electrical_potential));
    kinetics->current_density(backPterb[protonic_electrical_potential]);


    kinetics->set_electrolyte_potential(SOL(-0.05, 1, protonic_electrical_potential));
    reacts.clear(); reacts.push_back(SOL(1e-5 - perterb, 1, oxygen_concentration));
    kinetics->set_reactant_concentrations(reacts);
    kinetics->current_density(backPterb[oxygen_concentration]);

    reacts.clear() ; reacts.push_back(SOL(1e-5, 1, oxygen_concentration));
    reacts.push_back(SOL(c_h_ref-perterb,1,proton_concentration));
    kinetics->set_reactant_concentrations(reacts);
    kinetics->current_density(backPterb[proton_concentration]);

    //Custom perturbation for temperature derivative
    kinetics->set_temperature(SOL(353.0 - 1.0, 1, temperature_of_REV));
    kinetics->current_density(backPterb[temperature_of_REV]);







    std::map<VariableNames, double> numDerv;


    typedef std::map<VariableNames, std::vector<double>>::iterator it_type;
    for(it_type iterator = forwardPterb.begin(); iterator != forwardPterb.end(); iterator++) {

        numDerv[iterator->first] = (iterator->second[0] -backPterb[iterator->first][0])/(2*perterb);

    }



    //Custom perturbation for temperature derivative
    numDerv[temperature_of_REV] = (forwardPterb[temperature_of_REV][0] - backPterb[temperature_of_REV][0])/(2*1.0);

    double tolerance = 1e-4;
    double diff;

    diff = std::abs((derivatives[electronic_electrical_potential][0] -numDerv[electronic_electrical_potential])/numDerv[electronic_electrical_potential]);
    TEST_ASSERT_MSG(diff < tolerance, "electronic_electrical_potential Derivatives are incorrect");


    diff = std::abs((derivatives[protonic_electrical_potential][0] -numDerv[protonic_electrical_potential])/numDerv[protonic_electrical_potential]);
    TEST_ASSERT_MSG(diff < tolerance, "protonic_electrical_potential Derivatives are incorrect");

    diff = std::abs(((derivatives[oxygen_molar_fraction][0]* (3.1664e10/101325.0)) -numDerv[oxygen_concentration])/(numDerv[oxygen_concentration]));
    TEST_ASSERT_MSG(diff < tolerance, "oxygen Derivatives are incorrect");


    diff = std::abs((derivatives[proton_concentration][0] -numDerv[proton_concentration])/(numDerv[proton_concentration]));
    TEST_ASSERT_MSG(diff < tolerance, "proton_concentration Derivatives are incorrect");

    //Custom tolerance for temperature derivative
    diff = std::abs((derivatives[temperature_of_REV][0] -numDerv[temperature_of_REV])/numDerv[temperature_of_REV]);
    TEST_ASSERT_MSG(diff < 0.001, "temperature_of_REV Derivatives are incorrect");

    //Alternative method to check derivatives base on the fact that errors should be O(h^2)
    //Doesn't work
    //    TEST_ASSERT_DELTA_MSG(derivatives[electronic_electrical_potential][0],numDerv[protonic_electrical_potential],
    //            perterb, "electronic_electrical_potential Derivatives are incorrect" );
    //
    //
    //    TEST_ASSERT_DELTA_MSG(derivatives[protonic_electrical_potential][0],numDerv[protonic_electrical_potential],
    //                perterb, "protonic_electrical_potential Derivatives are incorrect" );
    //
    //    double val =derivatives[oxygen_molar_fraction][0]* (3.1664e10/101325.0);
    //    TEST_ASSERT_DELTA_MSG((derivatives[oxygen_molar_fraction][0]* (3.1664e10/101325.0)),numDerv[oxygen_concentration],
    //                    perterb, "oxygen Derivatives are incorrect" );



}

