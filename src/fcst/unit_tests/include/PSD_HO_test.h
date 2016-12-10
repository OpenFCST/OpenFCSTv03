//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: PSD_test.h
//    - Description: Unit testing class for PSD
//    - Developers: Prafful Mangal
//    - Id: $Id: PSD_HO_test.h 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

/**
 * A unit test class that tests the PSD class. 
 * 
 */

#ifndef _FCST_PSD_HO_TESTSUITE
#define _FCST_PSD_HO_TESTSUITE

#include <cpptest.h>
#include <boost/lexical_cast.hpp>
#include <string.h>
#include <microscale/PSD_HO.h>
#include <layers/gas_diffusion_layer.h>
#include <layers/porous_layer.h>

class PSD_HO_Test: public Test::Suite
{
public:
    PSD_HO_Test()
    :
    psd_object("BasePSD")
    {
        //Add a number of tests that will be called during Test::Suite.run()
        //Generic cases
//         TEST_ADD(PSD_HO_Test::testcompute_rc_HO);
//         TEST_ADD(PSD_HO_Test::testcompute_k_sat_HO);
//         TEST_ADD(PSD_HO_Test::testcompute_sat_HO);
//         TEST_ADD(PSD_HO_Test::testcompute_k_L_HO);
//         TEST_ADD(PSD_HO_Test::testcompute_kr_L_HO);
//         TEST_ADD(PSD_HO_Test::testcompute_k_G_HO);
//         TEST_ADD(PSD_HO_Test::testcompute_kr_G_HO);
//         TEST_ADD(PSD_HO_Test::testcompute_interfacial_area_per_volume_HO);
//         TEST_ADD(PSD_HO_Test::testcompute_avg_rk_HO);
    }
protected:
    virtual void setup();
      // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down(){}; // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
private:
    FuelCellShop::MicroScale::HOPSD<dim> psd_object;    // PSD object which we are testing
    void testcompute_rc_HO();
//     void testcompute_k_sat_HO();
//     void testcompute_sat_HO();
//     void testcompute_k_L_HO();
//     void testcompute_kr_L_HO();
//     void testcompute_k_G_HO();
//     void testcompute_kr_G_HO();
//     void testcompute_interfacial_area_per_volume_HO();
//     void testcompute_avg_rk_HO();
};

#endif
