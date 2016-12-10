//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: 
//    - Description: 
//    - Developers: 
//    - $Id: analytical_agglomerate_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------
#include <analytical_agglomerate_test.h>

using namespace FuelCellShop;

std::string parseDouble(double toParse){
	return boost::lexical_cast<std::string>(toParse);
}




void AnalyticalAgglomerateTest::testO2CurrentDensity()
{
    ParameterHandler param;
    //kinetics.declare_parameters(param);
    FuelCellShop::Layer::CatalystLayer<dim>::declare_CatalystLayer_parameters("CathodeCL", param);
    FuelCellShop::Material::PolymerElectrolyteBase::declare_PolymerElectrolyte_parameters(param);
    FuelCellShop::Material::CatalystSupportBase::declare_CatalystSupport_parameters(param);
    FuelCellShop::Material::CatalystBase::declare_Catalyst_parameters(param);
    OC.declare_parameters(param);

    param.enter_subsection("Fuel cell data");
        param.enter_subsection("CathodeCL");
        param.set("Material id", "4");
        param.set("Catalyst layer type", "MultiScaleCL");
        
            param.enter_subsection("ConventionalCL");
            param.set("Platinum loading on support (%wt)", "4:0.46");
            param.set("Platinum loading per unit volume (mg/cm3)", "4:400");
            param.set("Electrolyte loading (%wt)", "4:0.30");
            param.set("Active area [cm^2/cm^3]", "4:2.0e5");
            param.leave_subsection();
        
        param.enter_subsection("MultiScaleCL");
            param.enter_subsection("MicroScale");
            param.set("Microscale type", "IonomerAgglomerateAnalytical");
            param.leave_subsection();
        param.leave_subsection();
        
        param.leave_subsection();
    param.leave_subsection();


    CCL = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("CathodeCL", param);
    CCL->set_local_material_id(4);
    OC.initialize(param);

    CCL->set_reaction_kinetics(ORR);
    CCL->set_constant_solution(OC.get_pc_Pa(), total_pressure);
    CCL->set_constant_solution(OC.get_T(), temperature_of_REV);

    std::vector<FuelCellShop::SolutionVariable> sols;
    sols.push_back(FuelCellShop::SolutionVariable(0.1,1,oxygen_molar_fraction));
    sols.push_back(FuelCellShop::SolutionVariable(0.625,1,electronic_electrical_potential));
    sols.push_back(FuelCellShop::SolutionVariable(-0.1,1,protonic_electrical_potential));
    sols.push_back(FuelCellShop::SolutionVariable(8,1,membrane_water_content));
    CCL->set_solution(sols);

    std::vector<double> current(1, 0.0);

    auto t_start = std::chrono::high_resolution_clock::now();
    CCL->current_density(current);
    auto t_end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();


    //Ionomer analytical : 206.25 match within 5%
    double expected = 2113.66;
    double match = expected*0.05;
    
    std::string msg = "Current from model ionomer analytical does not match expected results! Value obtained is " + std::to_string(current[0]);
    TEST_ASSERT_DELTA_MSG(expected, current[0], match, msg.c_str());

    std::string time_err_msg = "The model took to long too solve. Time taken: " + FcstUtilities::number_to_string(time_taken) + " ms";
    TEST_ASSERT_MSG(time_taken < 15.0, time_err_msg.c_str() ); 

}

void AnalyticalAgglomerateTest::test02CurrentDerivative(){
    //std::vector<std::vector<double> >Dcurrent(3, std::vector<double>(1,0.0));

    std::vector<VariableNames> sol_names;
    sol_names.push_back(oxygen_molar_fraction);
    sol_names.push_back(protonic_electrical_potential);
    sol_names.push_back(electronic_electrical_potential);
    std::map< VariableNames, std::vector<double> > derivatives;

    CCL->set_derivative_flags(sol_names);
    CCL->derivative_current_density(derivatives);
    double dCurrentdO2 = derivatives[oxygen_molar_fraction][0];    //ionomer numerical : 2062.3 //ionomer analytical : 2062.4  //water numerical: 1885.3
    double dCurrentdPhi_m = derivatives[protonic_electrical_potential][0]; //ionomer numerical : 6082.1 //ionomer analytical : 6074.5  //water numerical: 6090.75
    double dCurrentdPhi_s = derivatives[electronic_electrical_potential][0]; //ionomer numerical :-6066   //ionomer analytical : -6074.5 //water numerical: -6027.7

    double expDCurrentDO2 = 21136;     double diffDCurrentDO2    = expDCurrentDO2*0.05;
    double expDCurrentPhi_m = 29393;  double diffDCurrentPhi_m  = expDCurrentPhi_m*0.05;
    double expDCurrentPhi_s = -29393;  double diffDCurrentPhi_s  = expDCurrentPhi_s*-0.05; //Difference has to be positive value

    std::string msg_1 = "Current derivative w.r.t. O2 is incorrect! Value obtained is " + std::to_string(dCurrentdO2);
    std::string msg_2 = "Current derivative w.r.t. Phi_m is incorrect! Value obtained is " + std::to_string(dCurrentdPhi_m);
    std::string msg_3 = "Current derivative w.r.t. Phi_s is incorrect! Value obtained is " + std::to_string(dCurrentdPhi_s);

    TEST_ASSERT_DELTA_MSG(dCurrentdO2, expDCurrentDO2, diffDCurrentDO2, msg_1.c_str());
    TEST_ASSERT_DELTA_MSG(dCurrentdPhi_m, expDCurrentPhi_m, diffDCurrentPhi_m, msg_2.c_str());
    TEST_ASSERT_DELTA_MSG(dCurrentdPhi_s, expDCurrentPhi_s, diffDCurrentPhi_s, msg_3.c_str());
}


void AnalyticalAgglomerateTest::testInvalidKineticsDT(){
    ParameterHandler param;
    FuelCellShop::Layer::CatalystLayer<dim>::declare_CatalystLayer_parameters("CathodeCL", param);
    FuelCellShop::Material::PolymerElectrolyteBase::declare_PolymerElectrolyte_parameters(param);
    FuelCellShop::Material::CatalystSupportBase::declare_CatalystSupport_parameters(param);
    FuelCellShop::Material::CatalystBase::declare_Catalyst_parameters(param);
    OC.declare_parameters(param);

    param.enter_subsection("Fuel cell data");
        param.enter_subsection("CathodeCL");
        
        param.set("Material id", "4");
        param.set("Kinetics type","DoubleTrapKinetics");
        param.set("Catalyst layer type", "MultiScaleCL");

            param.enter_subsection("ConventionalCL");
            param.set("Platinum loading on support (%wt)", "4:0.46");
            param.set("Platinum loading per unit volume (mg/cm3)", "4:400");
            param.set("Electrolyte loading (%wt)", "4:0.30");
            param.set("Active area [cm^2/cm^3]", "4:2.0e5");
            param.leave_subsection();
        
        param.enter_subsection("MultiScaleCL");
            param.enter_subsection("MicroScale");
            param.set("Microscale type", "IonomerAgglomerateAnalytical");
            param.leave_subsection();
        param.leave_subsection();
        
        param.leave_subsection();
    param.leave_subsection();

    CCL = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("CathodeCL", param);
    CCL->set_local_material_id(4);
    OC.initialize(param);

    CCL->set_reaction_kinetics(ORR);
    CCL->set_constant_solution(OC.get_pc_Pa(), total_pressure);
    CCL->set_constant_solution(OC.get_T(), temperature_of_REV);

    std::vector<FuelCellShop::SolutionVariable> sols;
    sols.push_back(FuelCellShop::SolutionVariable(0.1,1,oxygen_molar_fraction));
    sols.push_back(FuelCellShop::SolutionVariable(0.625,1,electronic_electrical_potential));
    sols.push_back(FuelCellShop::SolutionVariable(-0.1,1,protonic_electrical_potential));
    sols.push_back(FuelCellShop::SolutionVariable(8,1,membrane_water_content));
    CCL->set_solution(sols);

    std::vector<double> current(1, 0.0);

    TEST_THROWS_ANYTHING_MSG(CCL->current_density(current), "Running analytical agglomerate with DT kinetics should throw an exception");;
    
}
void AnalyticalAgglomerateTest::testInvalidKineticsORR(){
    ParameterHandler param;
        FuelCellShop::Layer::CatalystLayer<dim>::declare_CatalystLayer_parameters("CathodeCL", param);
        FuelCellShop::Material::PolymerElectrolyteBase::declare_PolymerElectrolyte_parameters(param);
        FuelCellShop::Material::CatalystSupportBase::declare_CatalystSupport_parameters(param);
        FuelCellShop::Material::CatalystBase::declare_Catalyst_parameters(param);
        OC.declare_parameters(param);

        param.enter_subsection("Fuel cell data");
			param.enter_subsection("CathodeCL");
               param.set("Material id", "4");
			param.enter_subsection("Materials");
				param.enter_subsection("Platinum");
				param.set("Oxygen reaction order (ORR)", "0.54");
				param.leave_subsection();
			param.leave_subsection();

			param.set("Catalyst layer type", "MultiScaleCL");

			param.enter_subsection("ConventionalCL");
			param.set("Platinum loading on support (%wt)", "4:0.46");
			param.set("Platinum loading per unit volume (mg/cm3)", "4:400");
			param.set("Electrolyte loading (%wt)", "4:0.30");
			param.set("Active area [cm^2/cm^3]", "4:2.0e5");
			param.leave_subsection();

			param.enter_subsection("MultiScaleCL");
				param.enter_subsection("MicroScale");
				param.set("Microscale type", "IonomerAgglomerateAnalytical");
				param.leave_subsection();
			param.leave_subsection();
			param.leave_subsection();
        param.leave_subsection();

        CCL = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("CathodeCL", param);
        //CCL->set_local_material_id(4);
        OC.initialize(param);

        CCL->set_reaction_kinetics(ORR);
        CCL->set_constant_solution(OC.get_pc_Pa(), total_pressure);
        CCL->set_constant_solution(OC.get_T(), temperature_of_REV);

        std::vector<FuelCellShop::SolutionVariable> sols;
        sols.push_back(FuelCellShop::SolutionVariable(0.1,1,oxygen_molar_fraction));
        sols.push_back(FuelCellShop::SolutionVariable(0.625,1,electronic_electrical_potential));
        sols.push_back(FuelCellShop::SolutionVariable(-0.1,1,protonic_electrical_potential));
        sols.push_back(FuelCellShop::SolutionVariable(8,1,membrane_water_content));
        CCL->set_solution(sols);

        std::vector<double> current(1, 0.0);

        TEST_THROWS_ANYTHING_MSG(CCL->current_density(current), "Running analytical agglomerate with ORR order not equal to 1.0 should throw an exception");;
}
