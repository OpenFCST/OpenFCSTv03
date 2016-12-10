/*
 * water_pore_agglomerate_test.cc
 *
 *  Created on: Feb 24, 2014
 *      Author: wardlawp
 */



#include<water_pore_agglomerate_test.h>


void WaterPoreAgglomerateTest::setup(){

    /*
    double AV = 0.0;
    double r = 50;
    double delta = 0.0;
    double epsilon_agg = 0.3;




    ParameterHandler param;

    agg.declare_parameters(param);


    param.enter_subsection("AgglomerateBase");
    param.set("Radius of the agglomerate [nm]", "50");
    param.set("Agglomerate porosity", "0.3");
    param.leave_subsection();

    agg.initialize(param);

     */

    //agg.set_P(101325);
    //agg.set_structure(support,catalyst,electrolyte,AV,r,delta,epsilon_agg);
    agg.set_structure();


}


void WaterPoreAgglomerateTest::testSolvePotentials(){
    std::map<VariableNames, FuelCellShop::SolutionVariable> sols;
    sols[temperature_of_REV] = FuelCellShop::SolutionVariable(353.0,1,temperature_of_REV);
    sols[oxygen_molar_fraction] = FuelCellShop::SolutionVariable(0.1,1,oxygen_molar_fraction);
    sols[electronic_electrical_potential] = FuelCellShop::SolutionVariable(0.5,1,electronic_electrical_potential);
    sols[protonic_electrical_potential] = FuelCellShop::SolutionVariable(0.0,1,protonic_electrical_potential);
    sols[proton_concentration] = FuelCellShop::SolutionVariable(1.2e-4,1,proton_concentration);
    agg.set_solution(sols,oxygen_molar_fraction,0);

    std::string error_msg;
    agg.initialize_problem(error_msg);


    //Call the solver
    TEST_ASSERT_MSG(agg.solveProtonPotentials(error_msg),error_msg.c_str());


    //Test solution
    TEST_ASSERT_MSG(not std::isnan(agg.PorePotential[agg.New][agg.M-1][0]), "Values of potential are not numbers!");

    double tolerance = 1e-4;
    TEST_ASSERT_DELTA_MSG(std::abs(agg.PorePotential[agg.New][0][0]), 0.322211, tolerance, "Potential value at point [0,0] is incorrect.");
    TEST_ASSERT_DELTA_MSG(std::abs(agg.PorePotential[agg.New][0][agg.N]), 0.483286, tolerance,"Potential value at point [0,N] is incorrect.");


    TEST_ASSERT_DELTA_MSG(std::abs(agg.PorePotential[agg.New][agg.M-1][0]), 0.496597, tolerance, "Potential value at point [M-1,0] is incorrect.");
    TEST_ASSERT_DELTA_MSG(std::abs(agg.PorePotential[agg.New][agg.M-1][agg.N]), 0.039768, tolerance, "Potential value at point [M-1,N] is incorrect.");





}



void WaterPoreAgglomerateTest::testSolveO2(){


    std::string error_msg;

    //Call the solver
    TEST_ASSERT_MSG(agg.solveO2(error_msg), error_msg.c_str());

    //Test solution
    TEST_ASSERT_MSG(not std::isnan(agg.cO2[agg.New][agg.M-1][0]), "Values of potential are not numbers!");

    double tolerance = 1e-7;
    TEST_ASSERT_DELTA_MSG(std::abs(agg.cO2[agg.New][0][0]), 4.58e-6, tolerance, "cO2 value at point [0,0] is incorrect.");
    TEST_ASSERT_DELTA_MSG(std::abs(agg.cO2[agg.New][0][agg.N]),  4.23e-6, tolerance,"cO2 value at point [0,N] is incorrect.");


    TEST_ASSERT_DELTA_MSG(std::abs(agg.cO2[agg.New][agg.M-1][0]), 0.8951, 1000*tolerance, "cO2 value at point [M-1,0] is incorrect.");
    TEST_ASSERT_DELTA_MSG(std::abs(agg.cO2[agg.New][agg.M-1][agg.N]), 0.7715, 1000*tolerance, "cO2 value at point [M-1,N] is incorrect.");





}


void WaterPoreAgglomerateTest::testCurrent(){

    double tolerance = 0.0001;
    double expected_eff = 0.1789;
    double eff;
    double expected_i;
    double i = agg.calculate_j(eff);

    std::string msg = "Incorrect effectivness value " + FcstUtilities::number_to_string(eff) + ", should be " + FcstUtilities::number_to_string(expected_eff) + ".";
    TEST_ASSERT_DELTA_MSG(eff, expected_eff, tolerance, msg.c_str());

    msg = FcstUtilities::number_to_string(i);
    TEST_ASSERT_MSG(false, msg.c_str());
    //TEST_ASSERT_DELTA_MSG(i, expected_i, tolerance, msg);


}
