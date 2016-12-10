//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: LiquidWater_test.cc
//    - Description: Unit testing class for LiquidWater
//    - Developers: Madhur Bhaiya
//    - Id: $Id: LiquidWater_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include "LiquidWater_test.h"

// ------------------------------------------------------------
void
LiquidWaterTest::testViscosity()
{
    double mu = FuelCellShop::Material::LiquidWater::viscosity(353.);
    double expected_answer = 351.655607e-5;
    
    std::ostringstream streamOut;
    streamOut <<"The value of the mu is: "<<mu<<". The expected value is: "<<expected_answer<<std::endl;
    streamOut <<"testViscosity failed!"<<std::endl;
    
    TEST_ASSERT_DELTA_MSG(mu, 351.655607e-5, 1.e-6, streamOut.str().c_str());
}

// ------------------------------------------------------------
void
LiquidWaterTest::testDerivViscosity()
{
    double mu_1 = FuelCellShop::Material::LiquidWater::viscosity(353.);
    double delta = std::pow(10., -7);
    double mu_2 = FuelCellShop::Material::LiquidWater::viscosity(353.+ delta);
    double dmu_dT = FuelCellShop::Material::LiquidWater::deriv_viscosity(353.);
    
    TEST_ASSERT_DELTA_MSG(dmu_dT, (mu_2-mu_1)/delta, 1.e-6, "deriv_viscosity method failure!");
}