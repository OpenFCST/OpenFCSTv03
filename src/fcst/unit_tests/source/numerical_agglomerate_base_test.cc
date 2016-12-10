#include <numerical_agglomerate_base_test.h>

void NumericalAgglomerateBaseTest::testAV(){
    agg.setAV(1);
    agg.interface = 1.0;
    agg.loadingWeigths.push_back(1); agg.loadingWeigths.push_back(2); agg.loadingWeigths.push_back(3);
    agg.setUpLoadings();

    std::string msg = FcstUtilities::number_to_string(agg.getAV(0.5));

    TEST_ASSERT_MSG(agg.getAV(0.5) == 0.75, msg.c_str());

}

void NumericalAgglomerateBaseTest::testAVSetup(){
    agg.setAV(1);
    agg.interface = 1.0;
    agg.loadingWeigths.push_back(1); agg.loadingWeigths.push_back(2); agg.loadingWeigths.push_back(3);
    agg.setUpLoadings();

    std::string msg = FcstUtilities::number_to_string(agg.actualLoadings[0]) + " "
            + FcstUtilities::number_to_string(agg.actualLoadings[1]) + " "
            + FcstUtilities::number_to_string(agg.actualLoadings[2]);

    TEST_ASSERT_MSG(((agg.actualLoadings[0] == 0.375) && (agg.actualLoadings[1] == 0.75) &&
            (agg.actualLoadings[2] == 1.125)), msg.c_str() );

}

void NumericalAgglomerateBaseTest::testSetSolution(){
    std::map<VariableNames,FuelCellShop::SolutionVariable> solutions;
    solutions[oxygen_molar_fraction] = FuelCellShop::SolutionVariable(0.5,2,oxygen_molar_fraction);

    //Check for setting solution without providing key variables (x_O2, phi_M, and phi_S)
    TEST_THROWS_MSG(agg.set_solution(solutions, oxygen_molar_fraction, 1), ExcMessage, "An exception for missing solution was expected.");


    solutions[protonic_electrical_potential] = FuelCellShop::SolutionVariable(0.6, 2, protonic_electrical_potential);
    solutions[electronic_electrical_potential] = FuelCellShop::SolutionVariable(0.6, 2, protonic_electrical_potential);

    std::vector<VariableNames> expected_sol_names;
    expected_sol_names.push_back(oxygen_molar_fraction);
    expected_sol_names.push_back(protonic_electrical_potential);
    expected_sol_names.push_back(electronic_electrical_potential);

    agg.set_solution(solutions,oxygen_molar_fraction, 1);

    //Check that the solution names are set up correctly
    TEST_ASSERT_MSG(agg.sol_names == expected_sol_names, "Solution names inside agglomerate are incorrect");

    //Check that the agglomerate wont accept a reactant it can't solve for.
    TEST_THROWS_MSG(agg.set_solution(solutions, oxygen_concentration, 1), ExcMessage, "Agglomerate should not except this reactant");

    //Check the solution data
    TEST_ASSERT_MSG(agg.solutions[protonic_electrical_potential][0] == 0.6, "Solution inside agglomerate is incorrect");
}
