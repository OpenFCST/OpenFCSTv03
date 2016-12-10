//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: puregas_oxygen_test.h
//    - Description: Unit testing class for Puregas Oxygen
//    - Developers: Jie Zhou
//    - $Id: puregas_oxygen_test.h 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

/**
 * A unit test class that tests the PureGas class. 
 * 
 */

#ifndef _FCST_PuregasOxygen_TESTSUITE
#define _FCST_PuregasOxygen_TESTSUITE

// STD
#include <cpptest.h>
#include <string.h>
#include <stdexcept>
#include <vector>
#include <iostream>

// FCST
#include <materials/PureGas.h>



class PuregasoxygenTest: public Test::Suite
{
public:
    PuregasoxygenTest()
    {
        
         TEST_ADD(PuregasoxygenTest::testmolarmass);
         TEST_ADD(PuregasoxygenTest::testcollision_diameter);
         TEST_ADD(PuregasoxygenTest::testeps_BY_k);
         TEST_ADD(PuregasoxygenTest::testPrandtl);
         TEST_ADD(PuregasoxygenTest::testID);
         TEST_ADD(PuregasoxygenTest::testchemical_formula);
         
         ///@name Service functions. EoS.
         
         TEST_ADD(PuregasoxygenTest::testpressure);
         TEST_ADD(PuregasoxygenTest::testpressuredensity);
         TEST_ADD(PuregasoxygenTest::testpressuredensitytemperature);
         TEST_ADD(PuregasoxygenTest::testDpressure_Ddensityconst);
         TEST_ADD(PuregasoxygenTest::testDpressure_Ddensityvector);
         TEST_ADD(PuregasoxygenTest::testDpressure_Dtemperatureconst);
         TEST_ADD(PuregasoxygenTest::testDpressure_Dtemperaturevector);
         
         ///@name Service functions. Sutherland dynamic viscosity.	
         
         TEST_ADD(PuregasoxygenTest::testSutherland_dynamic_viscosityconst);
         TEST_ADD(PuregasoxygenTest::testSutherland_dynamic_viscosityvector);
         TEST_ADD(PuregasoxygenTest::testDSutherland_dynamic_viscosity_Dtemperatureconst);
         TEST_ADD(PuregasoxygenTest::testDSutherland_dynamic_viscosity_Dtemperaturevector);
         
         ///@name Service functions. Chapman Enskog dynamic viscosity.
         
         TEST_ADD(PuregasoxygenTest::testChapmanEnskog_dynamic_viscosityconst);
         TEST_ADD(PuregasoxygenTest::testChapmanEnskog_dynamic_viscosityvector);
         TEST_ADD(PuregasoxygenTest::testDChapmanEnskog_dynamic_viscosity_Dtemperatureconst);
         
         ///@name Service functions. Chapman Enskog thermal conductivity.
         
         TEST_ADD(PuregasoxygenTest::testChapmanEnskog_thermal_conductivityconst);
         TEST_ADD(PuregasoxygenTest::testChapmanEnskog_thermal_conductivityvector);
         
         ///@name Service functions. Sutherland thermal conductivity.
         
         TEST_ADD(PuregasoxygenTest::testSutherland_thermal_conductivityconst);
         TEST_ADD(PuregasoxygenTest::testSutherland_thermal_conductivityvector);
         TEST_ADD(PuregasoxygenTest::testDSutherland_thermal_conductivity_Dtemperatureconst);
         
         ///@name Service functions. Molar enthalpy.
         
         TEST_ADD(PuregasoxygenTest::testmolar_enthalpyconst);
         TEST_ADD(PuregasoxygenTest::testmolar_enthalpyvector);
       
        
    }
protected:
    virtual void setup(); // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
private:
    FuelCellShop::Material::Oxygen oxy;
    void testmolarmass();
    void testcollision_diameter();
    void testeps_BY_k();
    void testPrandtl();
    void testID();
    void testchemical_formula();
    
    ///@name Service functions. EoS.
    
    void testpressure();
    void testpressuredensity();
    void testpressuredensitytemperature();
    void testDpressure_Ddensityconst();
    void testDpressure_Ddensityvector();
    void testDpressure_Dtemperatureconst();
    void testDpressure_Dtemperaturevector();
    
    ///@name Service functions. Sutherland dynamic viscosity.	
    
    void testSutherland_dynamic_viscosityconst();
    void testSutherland_dynamic_viscosityvector();
    void testDSutherland_dynamic_viscosity_Dtemperatureconst();
    void testDSutherland_dynamic_viscosity_Dtemperaturevector();
    
    ///@name Service functions. Chapman Enskog dynamic viscosity.
    
    void testChapmanEnskog_dynamic_viscosityconst();
    void testChapmanEnskog_dynamic_viscosityvector();
    void testDChapmanEnskog_dynamic_viscosity_Dtemperatureconst();
    
    ///@name Service functions. Chapman Enskog thermal conductivity.
    
    void testChapmanEnskog_thermal_conductivityconst();
    void testChapmanEnskog_thermal_conductivityvector();
    
    ///@name Service functions. Sutherland thermal conductivity.
    
    void testSutherland_thermal_conductivityconst();
    void testSutherland_thermal_conductivityvector();
    void testDSutherland_thermal_conductivity_Dtemperatureconst();
    
    ///@name Service functions. Molar enthalpy.
    
    void testmolar_enthalpyconst();
    void testmolar_enthalpyvector();
    
    
};

#endif










