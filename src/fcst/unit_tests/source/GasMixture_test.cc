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
//    - Description: Unit testing class for GasMixture
//    - Developers: Jie Zhou and Marc Secanell
//    - $Id: GasMixture_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include "GasMixture_test.h"

namespace NAME = FuelCellShop::Material;

//---------------------------------------------------------------------------
void GasMixtureTest::setup()
{
    fluid.declare_parameters(param);
    OC.declare_parameters(param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Operating conditions");
        {
            param.set("Cathode pressure [Pa]",
                      "100000.0");
            param.set("Temperature cell [K]",
                      "293.0");  
        }
        param.leave_subsection();
        
        param.enter_subsection("Materials");
        {
            param.enter_subsection("fluid");
            {
                param.set("Isobaric fluid flow",
                          "true");
                param.set("Isothermal fluid flow",
                          "true");  
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    OC.initialize(param);
    fluid.initialize(param);
    fluid.set_temperature(OC.get_T());
    fluid.set_total_pressure(OC.get_pc_Pa());
}

//---------------------------------------------------------------------------
void GasMixtureTest::testChapmanEnskog_isobaric_diffusion_coefficient(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    double answer = fluid.get_ChapmanEnskog_isobaric_diffusion_coefficient(); 
    
    double expectedAnswer = 7.8840471941;
        
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-5, "testChapmanEnskog_isobaric_diffusion_coefficient failed!");
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_isobaric_diffusion_coefficient(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    std::vector<double> answer(3);
    
    fluid.get_ChapmanEnskog_isobaric_diffusion_coefficient(answer); 
    
    std::vector<double> expectedAnswer (answer.size());
    
    expectedAnswer[0] = 7.8840471941;
    expectedAnswer[1] = 7.8840471941;
    expectedAnswer[2] = 7.8840471941;
    
    for (int i=0; i<answer.size(); i++)
    {
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-5, "testget_ChapmanEnskog_isobaric_diffusion_coefficient failed!");
    }
}

//---------------------------------------------------------------------------
void GasMixtureTest::testtemp_ChapmanEnskog_isobaric_diffusion_coefficient(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    double answer = fluid.get_ChapmanEnskog_isobaric_diffusion_coefficient(303.0); 
    
    double expectedAnswer = 8.3512973681;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-5, "testtemp_ChapmanEnskog_isobaric_diffusion_coefficient failed!");
}

//---------------------------------------------------------------------------
void GasMixtureTest::testtemp_diff_ChapmanEnskog_isobaric_diffusion_coefficient(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    std::vector<double> temperature(3);
    
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 293.0;
    
    std::vector<double> answer(3);
    
    fluid.get_ChapmanEnskog_isobaric_diffusion_coefficient(temperature,answer); 
    
    std::vector<double> expectedAnswer (answer.size());
    
    expectedAnswer[0] = 7.8840471941;
    expectedAnswer[1] = 8.3512973681;
    expectedAnswer[2] = 7.8840471941;
    
    for (int i=0; i<answer.size(); i++)
    {
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-5, "testtemp_diff_ChapmanEnskog_isobaric_diffusion_coefficient failed!");
    }
}

///@name Service functions. Derivatives of Chapman Enskog isobaric diffusion coefficient. Binary gas mixture only.
//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    double answer = fluid.get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(293.0); 
    
    double expectedAnswer = 0.046204176679;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-5, "testget_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature failed!");
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperaturevector(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    std::vector<double> temperature(3);
    
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 293.0;
    
    std::vector<double> answer(3);
    
    fluid.get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(temperature,answer); 
    
    std::vector<double> expectedAnswer (answer.size());
    
    expectedAnswer[0] = 0.046204176679;
    expectedAnswer[1] = 0.047244210264;
    expectedAnswer[2] = 0.046204176679;
    
    for (int i=0; i<answer.size(); i++)
    {
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-5, "testget_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperaturevector failed!");
    }
}

///@name Service functions. Chapman Enskog diffusion coefficient. Binary gas mixture only.
//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficient(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    double answer = fluid.get_ChapmanEnskog_diffusion_coefficient(); 
    
    double expectedAnswer = 7.8840471941/100000;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-5, "testget_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature failed!");
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficientvector(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    std::vector<double> answer(3);
    
    fluid.get_ChapmanEnskog_diffusion_coefficient(answer); 
    
    std::vector<double> expectedAnswer (answer.size());
    
    expectedAnswer[0] = 7.8840471941/100000;
    expectedAnswer[1] = 7.8840471941/100000;
    expectedAnswer[2] = 7.8840471941/100000;
    
    for (int i=0; i<answer.size(); i++)
    {
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-5, "testget_ChapmanEnskog_diffusion_coefficientvector failed!");
    }
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficient_at_constant_pressure(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    double answer = fluid.get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure(293.0); 
    
    double expectedAnswer = 7.8840471941/100000;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-5, "testget_ChapmanEnskog_diffusion_coefficient_at_constant_pressure failed!");
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficient_at_constant_pressurevector(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    std::vector<double> temperature(3);
    
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 293.0;
    
    std::vector<double> answer(3);
    
    fluid.get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure(temperature,answer); 
    
    std::vector<double> expectedAnswer (answer.size());
    
    expectedAnswer[0] = 7.8840471941/100000;
    expectedAnswer[1] = 8.3512973681/100000;
    expectedAnswer[2] = 7.8840471941/100000;
    
    for (int i=0; i<answer.size(); i++)
    {
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-5, "testget_ChapmanEnskog_diffusion_coefficient_at_constant_pressurevector failed!");
    }
}

//---------------------------------------------------------------------------
///@name Service functions. Derivatives of Chapman Enskog diffusion coefficient. Binary gas mixture only.
//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_diffusion_coefficient_Dpressure(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    double answer = fluid.get_DChapmanEnskog_diffusion_coefficient_Dpressure(100000.0); 
    
    double expectedAnswer = -7.8840471941/100000/100000;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-5, "testget_DChapmanEnskog_diffusion_coefficient_Dpressure failed!");
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_diffusion_coefficient_Dpressurevector(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    std::vector<double> temperature(3);
    
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 293.0;
    
    std::vector<double> pressure(3);
    
    pressure[0] = 100000.0;
    pressure[1] = 200000.0;
    pressure[2] = 100000.0;
    
    std::vector<double> answer(3);
    
    fluid.get_DChapmanEnskog_diffusion_coefficient_Dpressure(pressure,temperature,answer); 
    
    std::vector<double> expectedAnswer (answer.size());
    
    expectedAnswer[0] = -7.8840471941/100000/100000;
    expectedAnswer[1] = -8.3512973681/200000/200000;
    expectedAnswer[2] = -7.8840471941/100000/100000;
    
    for (int i=0; i<answer.size(); i++)
    {
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-5, "testget_DChapmanEnskog_diffusion_coefficient_Dpressurevector failed!");
    }
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_diffusion_coefficient_Dtemperature(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    double answer = fluid.get_DChapmanEnskog_diffusion_coefficient_Dtemperature(293.0); 
    
    double expectedAnswer = 0.046204176679/100000;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-5, "testget_DChapmanEnskog_diffusion_coefficient_Dtemperature failed!");
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_diffusion_coefficient_Dtemperaturevector(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    std::vector<double> temperature(3);
    
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 293.0;
    
    std::vector<double> pressure(3);
    
    pressure[0] = 100000.0;
    pressure[1] = 200000.0;
    pressure[2] = 100000.0;
    
    std::vector<double> answer(3);
    
    fluid.get_DChapmanEnskog_diffusion_coefficient_Dtemperature(pressure,temperature,answer); 
    
    std::vector<double> expectedAnswer (answer.size());
    
    expectedAnswer[0] = 0.046204176679/100000;
    expectedAnswer[1] = 0.047244210264/200000;
    expectedAnswer[2] = 0.046204176679/100000;
    
    for (int i=0; i<answer.size(); i++)
    {
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-5, "testget_DChapmanEnskog_diffusion_coefficient_Dtemperaturevector failed!");
    }
}

//---------------------------------------------------------------------------
///@name Service functions. Chapman Enskog isobaric diffusion coefficients. Ternary and more complicated gas mixtures.
//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.clear();
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    gases.push_back(&nitrogen);
    
    fluid.set_gases(gases);
    
    Table< 2, double > answer(gases.size(), gases.size());
    
    answer = fluid.get_ChapmanEnskog_isobaric_diffusion_coefficients();
    
    Table< 2, double > expectedAnswer(gases.size(), gases.size());
    expectedAnswer(0,1) = 7.8840471941;
    expectedAnswer(1,0) = 7.8840471941;
    expectedAnswer(0,2) = 1.9944813200;
    expectedAnswer(2,0) = 1.9944813200;
    expectedAnswer(1,2) = 7.4629860961;
    expectedAnswer(2,1) = 7.4629860961;
    
    for(unsigned int i = 0; i < gases.size(); ++i)
    {
        for(unsigned int j = 0; j < gases.size(); ++j)
        {
            if (i != j)
            {
                TEST_ASSERT_MSG(std::fabs((answer(i,j)-expectedAnswer(i,j))/expectedAnswer(i,j))<10e-5, "testget_DChapmanEnskog_diffusion_coefficient_Dtemperature failed!");
            }
        }
        
    }
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table_vector(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.clear();
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    gases.push_back(&nitrogen);
    
    fluid.set_gases(gases);
    
    Table< 2, double > answer(gases.size(), gases.size());
    
    double temperature = 293.0;
    
    answer = fluid.get_ChapmanEnskog_isobaric_diffusion_coefficients(temperature);
    
    Table< 2, double > expectedAnswer(gases.size(), gases.size());
    expectedAnswer(0,1) = 7.8840471941;
    expectedAnswer(1,0) = 7.8840471941;
    expectedAnswer(0,2) = 1.9944813200;
    expectedAnswer(2,0) = 1.9944813200;
    expectedAnswer(1,2) = 7.4629860961;
    expectedAnswer(2,1) = 7.4629860961;
    
    for(unsigned int i = 0; i < gases.size(); ++i)
    {
        for(unsigned int j = 0; j < gases.size(); ++j)
        {
            if (i != j)
            {
                TEST_ASSERT_MSG(std::fabs((answer(i,j)-expectedAnswer(i,j))/expectedAnswer(i,j))<10e-5, "testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table_vector failed!");
            }
        }
        
    }
}

//---------------------------------------------------------------------------
///@name Service functions. Derivatives of Chapman Enskog isobaric diffusion coefficients. Ternary and more complicated gas mixtures.
//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature_table(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.clear();
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    gases.push_back(&nitrogen);
    
    fluid.set_gases(gases);
    
    std::vector< Table< 2, double > > answer (3);
    
    std::vector <double> atemperature(3);
    atemperature[0] = 293.0;
    atemperature[1] = 293.0;
    atemperature[2] = 293.0;
    
    fluid.get_DChapmanEnskog_diffusion_coefficients_Dtemperature(atemperature,answer);
    
    Table< 2, double > expectedAnswer(gases.size(), gases.size());
    expectedAnswer(0,1) = 0.046204176679/100000;
    expectedAnswer(1,0) = 0.046204176679/100000;
    expectedAnswer(0,2) = 0.012129471075/100000;
    expectedAnswer(2,0) = 0.012129471075/100000;
    expectedAnswer(1,2) = 0.043597154276/100000;
    expectedAnswer(2,1) = 0.043597154276/100000;
    
    unsigned int q = 2;
    
    for(unsigned int i = 0; i < gases.size(); ++i)
    {
        for(unsigned int j = 0; j < gases.size(); ++j)
        {
            if (i != j)
            {
                TEST_ASSERT_MSG(std::fabs((answer[q](i,j)-expectedAnswer(i,j))/expectedAnswer(i,j))<10e-5, "testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table_vector failed!");
            }
        }
        
    }
}

//---------------------------------------------------------------------------
///@name Service functions. Chapman Enskog diffusion coefficients. Ternary and more complicated gas mixtures.
//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficients_table(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.clear();
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    gases.push_back(&nitrogen);
    
    fluid.set_gases(gases);
    
    std::vector< Table< 2, double > > answer (3);
    
    fluid.get_ChapmanEnskog_diffusion_coefficients(answer);
    
    Table< 2, double > expectedAnswer(gases.size(), gases.size());
    expectedAnswer(0,1) = 7.8840471941/100000;
    expectedAnswer(1,0) = 7.8840471941/100000;
    expectedAnswer(0,2) = 1.9944813200/100000;
    expectedAnswer(2,0) = 1.9944813200/100000;
    expectedAnswer(1,2) = 7.4629860961/100000;
    expectedAnswer(2,1) = 7.4629860961/100000;
    
    
    unsigned int q = 2;
    
    for(unsigned int i = 0; i < gases.size(); ++i)
    {
        for(unsigned int j = 0; j < gases.size(); ++j)
        {
            if (i != j)
            {
                TEST_ASSERT_MSG(std::fabs((answer[q](i,j)-expectedAnswer(i,j))/expectedAnswer(i,j))<10e-5, "testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table_vector failed!");
            }
        }
        
    }
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficients_table_tp(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.clear();
    
    gases.push_back(&oxygen);
            
    gases.push_back(&hydrogen);
    
    gases.push_back(&nitrogen);
    
    fluid.set_gases(gases);
    
    std::vector< Table< 2, double > > answer (3);
    
    std::vector <double> tpressure(3);
    tpressure[0] = 500000.0;
    tpressure[1] = 500000.0;
    tpressure[2] = 500000.0;
    
    std::vector <double> atemperature(3);
    atemperature[0] = 500.0;
    atemperature[1] = 500.0;
    atemperature[2] = 500.0;
    
    fluid.get_ChapmanEnskog_diffusion_coefficients(tpressure,atemperature,answer);
    
    Table< 2, double > expectedAnswer(gases.size(), gases.size());
    expectedAnswer(0,1) = 19.52273379/500000;
    expectedAnswer(1,0) = 19.52273379/500000;
    expectedAnswer(0,2) = 5.0620398945/500000;
    expectedAnswer(2,0) = 5.0620398945/500000;
    expectedAnswer(1,2) = 18.4376329159709/500000;
    expectedAnswer(2,1) = 18.4376329159709/500000;
    
    unsigned int q = 2;
    
    for(unsigned int i = 0; i < gases.size(); ++i)
    {
        for(unsigned int j = 0; j < gases.size(); ++j)
        {
            if (i != j)
            {
                TEST_ASSERT_MSG(std::fabs((answer[q](i,j)-expectedAnswer(i,j))/expectedAnswer(i,j))<10e-5, "testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table_vector failed!");
            }
        }
    }
}

//---------------------------------------------------------------------------
///@name Service functions. Derivatives of Chapman Enskog diffusion coefficients. Ternary and more complicated gas mixtures.
//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_diffusion_coefficients_Dtemperature(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.clear();
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    gases.push_back(&nitrogen);
    
    fluid.set_gases(gases);
    
    std::vector< Table< 2, double > > answer (3);
    
    std::vector <double> tpressure(3);
    tpressure[0] = 500000.0;
    tpressure[1] = 500000.0;
    tpressure[2] = 500000.0;
    
    std::vector <double> atemperature(3);
    atemperature[0] = 500.0;
    atemperature[1] = 500.0;
    atemperature[2] = 500.0;
    
    fluid.get_DChapmanEnskog_diffusion_coefficients_Dtemperature(tpressure,atemperature,answer);
    
    Table< 2, double > expectedAnswer(gases.size(), gases.size());
    expectedAnswer(0,1) = 0.0655584367221128/500000;
    expectedAnswer(1,0) = 0.0655584367221128/500000;
    expectedAnswer(0,2) = 0.0173378420615506/500000;
    expectedAnswer(2,0) = 0.0173378420615506/500000;
    expectedAnswer(1,2) = 0.0617900564100835/500000;
    expectedAnswer(2,1) = 0.0617900564100835/500000;
    
    unsigned int q = 2;
    
    for(unsigned int i = 0; i < gases.size(); ++i)
    {
        for(unsigned int j = 0; j < gases.size(); ++j)
        {
            if (i != j)
            {
                TEST_ASSERT_MSG(std::fabs((answer[q](i,j)-expectedAnswer(i,j))/expectedAnswer(i,j))<10e-5, "testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table_vector failed!");
            }
        }
        
    }
}

//---------------------------------------------------------------------------
///@name Service functions. Calculate partial viscosity for different empirical equations
//---------------------------------------------------------------------------
void GasMixtureTest::testget_isothermal_nonisobaric_Wilke_partial_viscosity(){
    //Initialize vector of gas properties
    std::vector<double> pureComponentViscosity(3, 0.0);
    std::vector<double> gasDensity(3, 0.0);
    std::vector<double> gasMolarMass(3, 0.0);
    
    //Properties are at T=353.15 [K] and P=101325.0 [Pa]
    pureComponentViscosity[0] = 2.34147e-4; // [g/cm*s]; O2
    pureComponentViscosity[1] = 1.25395e-4; // [g/cm*s]; H2O
    pureComponentViscosity[2] = 1.98067e-4; // [g/cm*s]; N2
    
    gasDensity[0] = 1.10423e-3;  // [g/cm^3]; O2
    gasDensity[1] = 0.621667e-3; // [g/cm^3]; H2O
    gasDensity[2] = 0.966694e-3; // [g/cm^3]; N2
    
    gasMolarMass[0] = 0.031999;  // [g/mol]; O2
    gasMolarMass[1] = 0.018015;  // [g/mol]; H2O
    gasMolarMass[2] = 0.0280134; // [g/mol]; N2
    
    //Create vector of hand calculated solutions
    std::vector<double> expectedAnswer(3, 0.0);
    expectedAnswer[0] = 7.731091942e-5; // [g/cm*s]; O2
    expectedAnswer[1] = 4.256173901e-5; // [g/cm*s]; H2O
    expectedAnswer[2] = 6.642618664e-5; // [g/cm*s]; N2
    
    //Get solution from function to be compared with expectedAnswer
    std::vector<double> answer(3, 0.0);
    answer = fluid.get_isothermal_nonisobaric_Wilke_partial_viscosity(gasDensity, pureComponentViscosity, gasMolarMass, 1.0);
    
    //Check accuracy of function
    for(unsigned int i = 0; i < answer.size(); ++i)
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-10, "testget_isothermal_nonisobaric_Wilke_partial_viscosity failed!");
}

void GasMixtureTest::testget_Wilke_variation_partial_viscosity_wrt_density(){
    //Initialize vector of gas properties
    std::vector<double> porosity(1, 1.0); //One quadrature point
    std::vector<double> pureComponentViscosity(3, 0.0);
    std::vector<double> gasMolarMass(3, 0.0);
    
    std::vector< std::vector<double> > gasDensity(3, std::vector<double>(1, 0.0));
    
    std::vector< std::vector< std::vector<double> > > xi(1, std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0)));
    std::vector< std::vector< std::vector<double> > > gasDensityVariation(3, std::vector< std::vector<double> >(1, std::vector<double>(1, 0.0)));
    
    //Properties are at T=353.15 [K] and P=101325.0 [Pa]
    pureComponentViscosity[0] = 2.34147e-4; // [g/cm*s]; O2
    pureComponentViscosity[1] = 1.25395e-4; // [g/cm*s]; H2O
    pureComponentViscosity[2] = 1.98067e-4; // [g/cm*s]; N2
    
    gasDensity[0][0] = 1.10423e-3;  // [g/cm^3]; O2
    gasDensity[1][0] = 0.621667e-3; // [g/cm^3]; H2O
    gasDensity[2][0] = 0.966694e-3; // [g/cm^3]; N2
    
    gasDensityVariation[0][0][0] = 5.763e-3; // O2
    gasDensityVariation[1][0][0] = 3.947e-3; // H2O
    gasDensityVariation[2][0][0] = 7.251e-3; // N2
    
    gasMolarMass[0] = 0.031999;  // [g/mol]; O2
    gasMolarMass[1] = 0.018015;  // [g/mol]; H2O
    gasMolarMass[2] = 0.0280134; // [g/mol]; N2
    
    xi[0][0][0] = 1.0000; xi[0][0][1] = 1.0118; xi[0][0][2] = 1.0168;
    xi[0][1][0] = 0.9625; xi[0][1][1] = 1.0000; xi[0][1][2] = 0.9837;
    xi[0][2][0] = 0.9825; xi[0][2][1] = 0.9992; xi[0][2][2] = 1.0000;
    
    //Create vector of hand calculated solutions
    std::vector<double> expectedAnswer(3, 0.0);
    expectedAnswer[0] = -8.84139176e-5; // [g/cm*s]; O2
    expectedAnswer[1] = -6.54850635e-7; // [g/cm*s]; H2O
    expectedAnswer[2] = 7.558422036e-5; // [g/cm*s]; N2
    
    //Get solution from function to be compared with expectedAnswer
    std::vector< std::vector< std::vector<double> > > answer(3, std::vector< std::vector<double> >(1, std::vector<double>(1, 0.0)));
    fluid.get_Wilke_delta_partial_viscosity_wrt_density(xi, porosity, gasMolarMass, pureComponentViscosity, gasDensityVariation, gasDensity, answer);
    
    //Check accuracy of function
    for(unsigned int i = 0; i < answer.size(); ++i)
        TEST_ASSERT_MSG(std::fabs((answer[i][0][0]-expectedAnswer[i])/expectedAnswer[i])<10e-10, "testget_Wilke_variation_partial_viscosity_wrt_density failed!");
}

void GasMixtureTest::testget_isothermal_nonisobaric_OmegaKG_partial_viscosity(){
    //properties
    const double temperature = 373.0; // [K]
    
    //Set gases
    std::vector< NAME::PureGas* > gases;
    gases.push_back(&argon);
    gases.push_back(&neon);
    gases.push_back(&helium);
    fluid.set_gases(gases);
    const unsigned int numOfSpecies = gases.size();
        
    //Initialize vector of gas properties
    std::vector<double> gasDensity(numOfSpecies, 0.0);
    std::vector<double> collisionDiameter(numOfSpecies, 0.0);
    std::vector<double> gasMolarMass(numOfSpecies, 0.0);
    
    gasMolarMass[0] = 39.94800; // [g/mol]; Ar
    gasMolarMass[1] = 20.18300; // [g/mol]; Ne
    gasMolarMass[2] = 4.002602; // [g/mol]; He
    
    collisionDiameter[0] = 3.418e-10; // [meter]; Ar
    collisionDiameter[1] = 2.789e-10; // [meter]; Ne
    collisionDiameter[2] = 2.576e-10; // [meter]; He
    
    //We use density instead of molar fraction
    gasDensity[0] = 1.0e-3; // [g/cm^3]; Ar
    gasDensity[1] = 1e-3 * 0.5576 * gasMolarMass[1] / (0.2670 * gasMolarMass[0]); // [g/cm^3]; Ne
    gasDensity[2] = 1e-3 * 0.1754 * gasMolarMass[2] / (0.2670 * gasMolarMass[0]); // [g/cm^3]; He
    
    //Create vector of hand calculated solutions
    std::vector<double> expectedAnswer(numOfSpecies, 0.0);
    expectedAnswer[0] = 1.0692960568693e-4;   // [g/cm*s]; Ar
    expectedAnswer[1] = 1.89490636926156e-4;  // [g/cm*s]; Ne
    expectedAnswer[2] = 0.268908694753035e-4; // [g/cm*s]; He

    //Get solution from function to be compared with expectedAnswer
    std::vector<double> answer(numOfSpecies, 0.0);
    answer = fluid.get_isothermal_nonisobaric_OmegaKG_partial_viscosity(gasDensity, collisionDiameter, gasMolarMass, temperature, 1.0);
    
    //Check accuracy of function
    for(unsigned int i = 0; i < answer.size(); ++i)
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-10, "testget_isothermal_nonisobaric_OmegaKG_partial_viscosity failed!");
}

void GasMixtureTest::testget_OmegaKG_variation_partial_viscosity_wrt_density(){
    //Initialize vector of gas properties
    const double temperature = 373.0;
    
    //Set gases
    std::vector< NAME::PureGas* > gases;
    gases.push_back(&argon);
    gases.push_back(&neon);
    gases.push_back(&helium);
    fluid.set_gases(gases);
    const unsigned int numOfSpecies = 3;
    
    std::vector<double> porosity(1, 1.0); //One quadrature point
    std::vector<double> collisionDiameter(numOfSpecies, 0.0);
    std::vector<double> gasMolarMass(numOfSpecies, 0.0);
    
    std::vector< std::vector<double> > gasDensity(numOfSpecies, std::vector<double>(1, 0.0)); //For variation calc
    
    std::vector< std::vector< std::vector<double> > > omegaIntegralTable(1, std::vector< std::vector<double> >(numOfSpecies*2, std::vector<double>(numOfSpecies, 0.0)));
    std::vector< std::vector< std::vector<double> > > gasDensityVariation(numOfSpecies, std::vector< std::vector<double> >(1, std::vector<double>(1, 0.0)));
    
    std::vector< FullMatrix<double> > PInv;
    PInv.resize(1);
    PInv[0].reinit(numOfSpecies, numOfSpecies);
    
    gasMolarMass[0] = 39.94800; // [g/mol]; Ar
    gasMolarMass[1] = 20.18300; // [g/mol]; Ne
    gasMolarMass[2] = 4.002602; // [g/mol]; He
    
    collisionDiameter[0] = 3.418e-10; // [meter]; Ar
    collisionDiameter[1] = 2.789e-10; // [meter]; Ne
    collisionDiameter[2] = 2.576e-10; // [meter]; He
    
    gasDensity[0][0] = 1.0e-3; // [g/cm^3]; Ar
    gasDensity[1][0] = 1e-3 * 0.5576 * gasMolarMass[1] / (0.2670 * gasMolarMass[0]); // [g/cm^3]; Ne
    gasDensity[2][0] = 1e-3 * 0.1754 * gasMolarMass[2] / (0.2670 * gasMolarMass[0]); // [g/cm^3]; He
    
    //Omega_11 integrals                                                                   Omega_22 integrals
    omegaIntegralTable[0][0][0] = 0.0;                                                     omegaIntegralTable[0][3][0] = 1.1988325544050172232895870796352435978282060006442e-16;
    omegaIntegralTable[0][0][1] = 4.7715015527639230252974280435096692983945438856254e-17; omegaIntegralTable[0][3][1] = 1.0517198560579174580595809709537204335065278029129e-16;
    omegaIntegralTable[0][0][2] = 7.6639456808781214154357224378942582574101256907666e-17; omegaIntegralTable[0][3][2] = 1.6999700684017038379711014404616462875757041241318e-16;
    
    omegaIntegralTable[0][1][0] = 4.7715015527639230252974280435096692983945438856254e-17; omegaIntegralTable[0][4][0] = 1.0517198560579174580595809709537204335065278029129e-16;
    omegaIntegralTable[0][1][1] = 0.0;                                                     omegaIntegralTable[0][4][1] = 8.8444408472299732050485430193368488361097420294866e-17;
    omegaIntegralTable[0][1][2] = 5.8075461309260823998933506283164509454474575759608e-17; omegaIntegralTable[0][4][2] = 1.2973840884644438842046965679389411544654237484157e-16;
    
    omegaIntegralTable[0][2][0] = 7.6639456808781214154357224378942582574101256907666e-17; omegaIntegralTable[0][5][0] = 1.6999700684017038379711014404616462875757041241318e-16;
    omegaIntegralTable[0][2][1] = 5.8075461309260823998933506283164509454474575759608e-17; omegaIntegralTable[0][5][1] = 1.2973840884644438842046965679389411544654237484157e-16;
    omegaIntegralTable[0][2][2] = 0.0;                                                     omegaIntegralTable[0][5][2] = 1.4107299464170787764039664623566099723710237829905e-16;
    
    PInv[0](0, 0) = 9.3607780822693117158805192756787505459215026348829e-06;
    PInv[0](0, 1) = 1.2204026151289843300583965943206798954179248539731e-06;
    PInv[0](0, 2) = 1.1177987059654905623823053150908690689391278283438e-07;
    
    PInv[0](1, 0) = 1.2204026151289843300583965943206798954179248539731e-06;
    PInv[0](1, 1) = 1.7492177428501410327895543295184666021668817847967e-05;
    PInv[0](1, 2) = 2.3648364774809272789059667700672928702942954259925e-07;
    
    PInv[0](2, 0) = 1.1177987059654905623823053150908690689391278283438e-07;
    PInv[0](2, 1) = 2.3648364774809272789059667700672928702942954259925e-07;
    PInv[0](2, 2) = 2.3408234290101494305951380231478609061923634726554e-06;
    
    gasDensityVariation[0][0][0] = 5.763e-3; // Ar
    gasDensityVariation[1][0][0] = 3.947e-3; // Ne
    gasDensityVariation[2][0][0] = 7.251e-3; // He
    
    std::vector< std::vector< std::vector<double> > > answer(numOfSpecies, std::vector< std::vector<double> >(1, std::vector<double>(1, 0.0)));
    fluid.get_OmegaKG_delta_partial_viscosity_wrt_density(temperature, omegaIntegralTable, PInv, porosity, gasMolarMass, collisionDiameter, gasDensityVariation, gasDensity, answer);
    
    //Create vector of hand calculated solutions
    std::vector<double> expectedAnswer(numOfSpecies, 0.0);
    expectedAnswer[0] = 8.730384947107410e-4; // Ar
    expectedAnswer[1] = 2.1029079665357200e-3; // Ne
    expectedAnswer[2] = -2.5133608407666700e-3; // He
    
    //Check accuracy of function
    for(unsigned int i = 0; i < answer.size(); ++i)
        TEST_ASSERT_MSG(std::fabs((answer[i][0][0]-expectedAnswer[i])/expectedAnswer[i])<10e-10, "testget_OmegaKG_variation_partial_viscosity_wrt_density failed!");
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------