/*
 * fcst_db_test.cc
 *
 *  Created on: Jun 21, 2013
 *      Author: wardlawp
 */

#include <solution_variable_test.h>


void
SolutionVariableTest::setup()
{
    sols.clear();
    //Constructor 1 (internally creates a vector)
    FuelCellShop::SolutionVariable phi_s(1.5, 8, electronic_electrical_potential);

    //Constructor 2 (copies a vector)
    std::vector<double> phiM_dummy(8, -0.2);
    FuelCellShop::SolutionVariable phi_m(phiM_dummy, protonic_electrical_potential);

    //Constructor 3 (copies a pointer)
    xO2_dummy = std::vector<double>(4, 0.5);
    FuelCellShop::SolutionVariable x_O2(&xO2_dummy, oxygen_molar_fraction);

    sols.push_back(phi_s);
    sols.push_back(phi_m);
    sols.push_back(x_O2);

}


void SolutionVariableTest::testAtOperator(){


    TEST_ASSERT_DELTA( sols.at(0)[0], 1.5, 1.e-6);
    TEST_ASSERT_DELTA( sols.at(1)[0], -0.2, 1.e-6);
    TEST_ASSERT_DELTA( sols.at(2)[2], 0.5, 1.e-6);

    //TEST_THROWS_MSG(sols.at(0)[8], ExcMessage, "An exception was expected for accessing a out of range element");     // If somebody wants to test this, Assert needs to be replaced by AssertThrow in SolutionVariable struct.


}


void SolutionVariableTest::testName(){

    TEST_ASSERT_MSG(sols.at(0).get_variablename() == electronic_electrical_potential, "Incorrect name returned");
    TEST_ASSERT_MSG(sols.at(1).get_variablename() == protonic_electrical_potential, "Incorrect name returned");
    TEST_ASSERT_MSG(sols.at(2).get_variablename() == oxygen_molar_fraction, "Incorrect name returned");

}


void SolutionVariableTest::testSize(){
    TEST_ASSERT_MSG(sols.at(0).size() == 8, "Incorrect size returned");
    TEST_ASSERT_MSG(sols.at(1).size() == 8, "Incorrect size returned");
    TEST_ASSERT_MSG(sols.at(2).size() == 4, "Incorrect size returned");
}


void SolutionVariableTest::testSolMap(){

    map.clear();

    //Add the sols created in SolutionVariableTest::setup()
    for (unsigned int i = 0; i < sols.size(); i++)
        map.push_back(sols[i]);

    //Check members have been correctly added to map
    TEST_ASSERT_DELTA( map.at(electronic_electrical_potential)[0], 1.5, 1.e-6);
    TEST_ASSERT_DELTA( map.at(protonic_electrical_potential)[0], -0.2, 1.e-6);
    TEST_ASSERT_DELTA( map.at(oxygen_molar_fraction)[2], 0.5, 1.e-6);


    //Test has() functionality
    TEST_ASSERT(map.has(electronic_electrical_potential));
    TEST_ASSERT(not map.has(temperature_of_REV));

}

