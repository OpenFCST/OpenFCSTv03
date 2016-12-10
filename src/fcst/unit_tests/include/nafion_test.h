//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: nafion_test.h
//    - Description: Unit testing class for Nafion
//    - Developers: Madhur Bhaiya
//    - Id: $Id: nafion_test.h 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

/**
 * A unit test class that tests the Nafion class. 
 * 
 */

#ifndef _FCST_Nafion_TESTSUITE
#define _FCST_Nafion_TESTSUITE

#include <cpptest.h>
#include <boost/lexical_cast.hpp>
#include <string.h>
#include <materials/nafion.h>


class NafionTest: public Test::Suite
{
public:
    NafionTest()
    {
        //Add a number of tests that will be called during Test::Suite.run()
        //Generic cases
        TEST_ADD(NafionTest::testSetT);
        TEST_ADD(NafionTest::testSetPT);
        TEST_ADD(NafionTest::testSetLambda);
        TEST_ADD(NafionTest::testSetTemperature);
        TEST_ADD(NafionTest::testSetMembraneWaterContent);
        TEST_ADD(NafionTest::testSetWaterMolarFraction);
        TEST_ADD(NafionTest::testProtonConductivitySpringer);
        TEST_ADD(NafionTest::testProtonConductivityNRE211);
        TEST_ADD(NafionTest::testProtonConductivityIden11);
        TEST_ADD(NafionTest::testProtonConductivityConstant);
        TEST_ADD(NafionTest::testWaterDiffusivityConstant);
        TEST_ADD(NafionTest::testWaterDiffusivitySpringer);
        TEST_ADD(NafionTest::testWaterDiffusivityMotupally);
        TEST_ADD(NafionTest::testElectroOsmoticDragConstant);
        TEST_ADD(NafionTest::testElectroOsmoticDragSpringer);
        TEST_ADD(NafionTest::testThermoOsmoticCoeffConstant);
        TEST_ADD(NafionTest::testThermoOsmoticCoeffKim);
        TEST_ADD(NafionTest::testSorptionEnthalpyConstant);
        TEST_ADD(NafionTest::testOxygenDiffusivity);
        TEST_ADD(NafionTest::testSorptionIsothermHinatsu);
        TEST_ADD(NafionTest::testSorptionIsothermLiu);

    }
protected:
    virtual void setup(){} // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down(){} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
private:
    FuelCellShop::Material::Nafion* electrolyte;    // Nafion material object which we are testing
    void testSetT();
    void testSetPT();
    void testSetLambda();
    
    void testSetTemperature();
    void testSetMembraneWaterContent();
    void testSetWaterMolarFraction();
    
    void testProtonConductivitySpringer();
    void testProtonConductivityNRE211();
    void testProtonConductivityIden11();
    void testProtonConductivityConstant();
    
    void testWaterDiffusivityConstant();
    void testWaterDiffusivitySpringer();
    void testWaterDiffusivityMotupally();
    
    void testElectroOsmoticDragConstant();
    void testElectroOsmoticDragSpringer();
    
    void testThermoOsmoticCoeffConstant();
    void testThermoOsmoticCoeffKim();

    void testSorptionEnthalpyConstant();
    
    void testOxygenDiffusivity();
    
    void testSorptionIsothermHinatsu();
    void testSorptionIsothermLiu();
};

#endif
