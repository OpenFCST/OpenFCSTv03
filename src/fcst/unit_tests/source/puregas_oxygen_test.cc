//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: puregas_oxygen_test.cc
//    - Description: Unit testing class for Puregas Oxygen
//    - Developers: Jie Zhou
//    - $Id: puregas_oxygen_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include "puregas_oxygen_test.h"


void PuregasoxygenTest::setup()
{}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testmolarmass(){
    
    
    double answer = oxy.get_molar_mass();  
    double expectedAnswer = 31.999e-3;
    
    TEST_ASSERT_MSG(expectedAnswer == answer, "Molar mass failed!");
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testcollision_diameter(){
    
    
    double answer = oxy.get_collision_diameter();  
    double expectedAnswer = 3.4330000;
    
    TEST_ASSERT_MSG(expectedAnswer == answer, "collision_diameter failed!");
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testeps_BY_k(){
    
    
    double answer = oxy.get_eps_BY_k();  
    double expectedAnswer = 113.00000;
    
    TEST_ASSERT_MSG(expectedAnswer == answer, "eps_BY_k failed!");
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testPrandtl(){
    
    
    double answer = oxy.get_Prandtl();  
    double expectedAnswer = 0.7130000;
    
    TEST_ASSERT_MSG(expectedAnswer == answer, "Prandtl failed!");
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testID(){
    
    
    double answer = oxy.get_ID();  
    double expectedAnswer = 1;
    
    TEST_ASSERT_MSG(expectedAnswer == answer, "ID failed!");
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testchemical_formula(){
    
    
    std::string answer = oxy.get_chemical_formula();  
    std::string expectedAnswer = "O2";
    
    TEST_ASSERT_MSG(expectedAnswer == answer, "chemical_formula failed!");
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
///@name Service functions. EoS.
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void PuregasoxygenTest::testpressure(){
       
    double answer = oxy.get_pressure(1000.0,293.0);  
    double expectedAnswer = 76131673.4137942;

        TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-15, "pressure failed!");
		
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testpressuredensity(){
    
    std::vector<double> density(3);
    density[0] = 1.0;
    density[1] = 10.0;
    density[2] = 100.0;
    std::vector<double> pressure (density.size());
    double temperature = 293.0;
    oxy.get_pressure(density,temperature,pressure);
    std::vector<double> answer (pressure);
    std::vector<double> expectedAnswer (density.size());
    expectedAnswer[0] =  76131.6734137942;
    expectedAnswer[1] =  761316.734137942;
    expectedAnswer[2] =  7613167.34137942;
    
    for (int i=0; i<density.size(); i++)
    {
            TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-15, "testpressuredensity failed!");
    }
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testpressuredensitytemperature(){
    
    std::vector<double> density(3);
    density[0] = 1.0;
    density[1] = 10.0;
    density[2] = 100.0;
    
    std::vector<double> temperature(3);
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 313.0;
    
    std::vector<double> pressure (density.size());
    
    oxy.get_pressure(density,temperature,pressure);
    
    std::vector<double> answer (pressure);
    
    std::vector<double> expectedAnswer  (density.size());
    expectedAnswer[0] =  76131.6734137942;
    expectedAnswer[1] =  787300.240422513;
    expectedAnswer[2] =  8132837.46707085;
    
    for (int i=0; i<density.size(); i++)
    {
            TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-15, "testpressuredensitytemperature failed!");
    }
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testDpressure_Ddensityconst(){
    
    double answer = oxy.get_Dpressure_Ddensity(293.0);  
    double expectedAnswer = 76131.6734137942;
	
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-15, "testDpressure_Ddensityconst failed!");

}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testDpressure_Ddensityvector(){
    
    std::vector<double> temperature(3);
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 313.0;
    
    std::vector<double> dst (temperature.size());
    
    oxy.get_Dpressure_Ddensity(temperature,dst);
    
    std::vector<double> answer (dst);
    std::vector<double> expectedAnswer (temperature.size());
    expectedAnswer[0] =  76131.6734137942;
    expectedAnswer[1] =  78730.0240422513;
    expectedAnswer[2] =  81328.3746707085;
    
    for (int i=0; i<temperature.size(); i++)
    {
		TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-15, "testDpressure_Ddensityvector failed!");
    }
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testDpressure_Dtemperatureconst(){
    
    double answer = oxy.get_Dpressure_Dtemperature(1.0);  
    double expectedAnswer = 259.835062845714;
	
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-15, "testDpressure_Dtemperatureconst failed!");

}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testDpressure_Dtemperaturevector(){
    
    std::vector<double> density(3);
    density[0] = 1.0;
    density[1] = 10.0;
    density[2] = 100.0;
    
    std::vector<double> dst (density.size());
    
    oxy.get_Dpressure_Dtemperature(density,dst);
    
    std::vector<double> answer (dst);
    std::vector<double> expectedAnswer (density.size());
    expectedAnswer[0] =  259.835062845714;
    expectedAnswer[1] =  2598.35062845714;
    expectedAnswer[2] =  25983.5062845714;
    
    for (int i=0; i<density.size(); i++)
    {
		TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-15, "testDpressure_Dtemperaturevector failed!");
    }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
///@name Service functions. Sutherland dynamic viscosity.	
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void PuregasoxygenTest::testSutherland_dynamic_viscosityconst(){
    
    double answer = oxy.get_Sutherland_dynamic_viscosity(293.0);  
    double expectedAnswer = 0.0000202215570772966;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-15, "testSutherland_dynamic_viscosityconst failed!");

    
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testSutherland_dynamic_viscosityvector(){
    
    std::vector<double> temperature(3);
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 313.0;
    
    std::vector<double> dynamic_viscosity (temperature.size());
    
    oxy.get_Sutherland_dynamic_viscosity(temperature,dynamic_viscosity);
    
    std::vector<double> answer (dynamic_viscosity);
    std::vector<double> expectedAnswer (temperature.size());
    expectedAnswer[0] =  0.0000202215570772966;
    expectedAnswer[1] =  0.000020771025582657;
    expectedAnswer[2] =  0.0000213121018654091;
    
    for (int i=0; i<temperature.size(); i++)
    {
		TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-15, "testSutherland_dynamic_viscosityvector failed!");
    }
    
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testDSutherland_dynamic_viscosity_Dtemperatureconst(){
    
    double answer = oxy.get_DSutherland_dynamic_viscosity_Dtemperature(293.0);  
    double expectedAnswer = 0.0000000553767652774985;
	
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-15, "testDSutherland_dynamic_viscosity_Dtemperatureconst failed!");
    
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testDSutherland_dynamic_viscosity_Dtemperaturevector(){
    
    std::vector<double> temperature(3);
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 313.0;
    
    std::vector<double> dst (temperature.size());
    
    oxy.get_DSutherland_dynamic_viscosity_Dtemperature(temperature,dst);
    
    std::vector<double> answer (dst);
    std::vector<double> expectedAnswer (temperature.size());
    expectedAnswer[0] =  0.0000000553767652774985;
    expectedAnswer[1] =  0.0000000545221486627423;
    expectedAnswer[2] =  0.0000000536980783277445;
    
    for (int i=0; i<temperature.size(); i++)
    {
		TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-15, "testDSutherland_dynamic_viscosity_Dtemperaturevector failed!");
    }
    
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
///@name Service functions. Chapman Enskog dynamic viscosity.
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void PuregasoxygenTest::testChapmanEnskog_dynamic_viscosityconst(){
    
    double answer = oxy.get_ChapmanEnskog_dynamic_viscosity(293.0);  
    double expectedAnswer = 0.0000202581487149307;
	
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-15, "testChapmanEnskog_dynamic_viscosityconst failed!");   
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testChapmanEnskog_dynamic_viscosityvector(){
    
    std::vector<double> temperature(3);
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 313.0;
    
    std::vector<double> dynamic_viscosity (temperature.size());
    
    oxy.get_ChapmanEnskog_dynamic_viscosity(temperature,dynamic_viscosity);
    
    std::vector<double> answer (dynamic_viscosity);
    std::vector<double> expectedAnswer (temperature.size());
    expectedAnswer[0] =  0.0000202581487149307;
    expectedAnswer[1] =  0.0000208017199698285;
    expectedAnswer[2] =  0.0000213373588044888;
    
    for (int i=0; i<temperature.size(); i++)
    {
		TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-15, "testChapmanEnskog_dynamic_viscosityvector failed!");
    }
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testDChapmanEnskog_dynamic_viscosity_Dtemperatureconst(){
    
    double answer = oxy.get_DChapmanEnskog_dynamic_viscosity_Dtemperature(293.0);  
    double expectedAnswer = 0.0000000547662377;
	
    TEST_ASSERT_MSG(std::fabs((expectedAnswer-expectedAnswer)/expectedAnswer)<10e-15, "testDChapmanEnskog_dynamic_viscosity_Dtemperatureconst failed!");   
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
///@name Service functions. Chapman Enskog thermal conductivity.
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void PuregasoxygenTest::testChapmanEnskog_thermal_conductivityconst(){
    
    double answer = oxy.get_ChapmanEnskog_thermal_conductivity(293.0);  
    double expectedAnswer = 0.0019738518252423;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-15, "testChapmanEnskog_thermal_conductivityconst failed!");  
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testChapmanEnskog_thermal_conductivityvector(){
    
    std::vector<double> temperature(3);
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 313.0;
    
    std::vector<double> thermal_conductivity (temperature.size());
    
    oxy.get_ChapmanEnskog_thermal_conductivity(temperature,thermal_conductivity);
    
    std::vector<double> answer (thermal_conductivity);
    std::vector<double> expectedAnswer (temperature.size());
    expectedAnswer[0] =  0.0019738518252423;
    expectedAnswer[1] =  0.00202681466645387;
    expectedAnswer[2] =  0.00207900461264997;
    
    for (int i=0; i<temperature.size(); i++)
    {
		TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-15, "testChapmanEnskog_thermal_conductivityvector failed!");
    }    
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
///@name Service functions. Sutherland thermal conductivity.
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void PuregasoxygenTest::testSutherland_thermal_conductivityconst(){
    
    double answer = oxy.get_Sutherland_thermal_conductivity(293.0);  
    double expectedAnswer = 0.0260364152839446;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-15, "testSutherland_thermal_conductivityconst failed!");
    
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testSutherland_thermal_conductivityvector(){
    
    std::vector<double> temperature(3);
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 313.0;
    
    std::vector<double> thermal_conductivity (temperature.size());
    
    oxy.get_Sutherland_thermal_conductivity(temperature,thermal_conductivity);
    
    std::vector<double> answer (thermal_conductivity);
    std::vector<double> expectedAnswer (temperature.size());
    expectedAnswer[0] =  0.0260364152839446;
    expectedAnswer[1] =  0.0268120024388918;
    expectedAnswer[2] =  0.027581767789448;
    
    for (int i=0; i<temperature.size(); i++)
    {
		TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-15, "testSutherland_thermal_conductivityvector failed!");
    }
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testDSutherland_thermal_conductivity_Dtemperatureconst(){
    
    double answer = oxy.get_DSutherland_thermal_conductivity_Dtemperature(293.0);  
    double expectedAnswer = 0.000077862119578924;
	
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-15, "testDSutherland_thermal_conductivity_Dtemperatureconst failed!");

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
///@name Service functions. Molar enthalpy.
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void PuregasoxygenTest::testmolar_enthalpyconst(){
    
    double answer = oxy.get_molar_enthalpy(293.0);  
    double expectedAnswer = -146.879957088299;

    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-15, "testmolar_enthalpyconst failed!");
    
}

//---------------------------------------------------------------------------
void PuregasoxygenTest::testmolar_enthalpyvector(){
    
    std::vector<double> temperature(3);
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 313.0;
    
    std::vector<double> molar_enthalpy (temperature.size());
    
    oxy.get_molar_enthalpy(temperature,molar_enthalpy);
    
    std::vector<double> answer (molar_enthalpy);
    std::vector<double> expectedAnswer (temperature.size());
    expectedAnswer[0] =  -146.879957088299;
    expectedAnswer[1] =  147.25405066151;
    expectedAnswer[2] =  442.90747486407;
    
    for (int i=0; i<temperature.size(); i++)
    {
		TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-15, "testmolar_enthalpyvector failed!");
    }
    
}
