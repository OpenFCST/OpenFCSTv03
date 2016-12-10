//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: LiquidWater_test.h
//    - Description: Unit testing class for LiquidWater
//    - Developers: Madhur Bhaiya
//
//---------------------------------------------------------------------------

/**
 * A unit test class that tests the LiquidWater class. 
 * 
 */

#ifndef _FCST_LiquidWater_TESTSUITE
#define _FCST_LiquidWater_TESTSUITE

#include <cpptest.h>
#include <boost/lexical_cast.hpp>
#include <string.h>
#include <materials/PureLiquid.h>


class LiquidWaterTest: public Test::Suite
{
public:
    LiquidWaterTest()
    {
        //Add a number of tests that will be called during Test::Suite.run()
        //Generic cases
        TEST_ADD(LiquidWaterTest::testViscosity);
        TEST_ADD(LiquidWaterTest::testDerivViscosity);

    }
protected:
    virtual void setup(){} // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down(){} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
private:
    
    void testViscosity();
    void testDerivViscosity();
};

#endif
