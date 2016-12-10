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
//    - Description: Unit testing class for PSD_HI
//    - Developers: 2009-13 by Marc Secanell, University of Alberta
//                  2013-14 by Jie Zhou, University of Alberta
//    - Id: $ 
//
//---------------------------------------------------------------------------

/**
 * A unit test class that tests the PSD class. 
 * 
 */

#ifndef _FCST_PSD_HI_TESTSUITE
#define _FCST_PSD_HI_TESTSUITE

#include <cpptest.h>
#include <boost/lexical_cast.hpp>
#include <string.h>
#include <microscale/PSD_HI.h>
#include <microscale/PSD_base.h>
#include <layers/gas_diffusion_layer.h>
#include <layers/porous_layer.h>


class PSD_HI_Test: public Test::Suite
{
public:
    PSD_HI_Test()
    :
    psd_object("BasePSD")
    {
        
        //Add a number of tests that will be called during Test::Suite.run()
        //Generic cases
       // TEST_ADD(PSD_HI_Test::testcompute_rc_HI);
        TEST_ADD(PSD_HI_Test::testcompute_k_sat_HI);
        TEST_ADD(PSD_HI_Test::testcompute_sat_HI);
//        TEST_ADD(PSD_HI_Test::testcompute_k_L_HI);
        TEST_ADD(PSD_HI_Test::testcompute_kr_L_HI);
//        TEST_ADD(PSD_HI_Test::testcompute_k_G_HI);
        TEST_ADD(PSD_HI_Test::testcompute_kr_G_HI);
        TEST_ADD(PSD_HI_Test::testcompute_interfacial_area_per_volume_HI);
//       TEST_ADD(PSD_HI_Test::testcompute_avg_rk_HI);
    }
protected:
    virtual void setup();
      // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down(){}; // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
private:
    FuelCellShop::MicroScale::DualPSD<dim>* psd_pointer;
    FuelCellShop::MicroScale::DualPSD<dim> psd_object;  // PSD object which we are testing
    void testcompute_rc_HI();
    void testcompute_k_sat_HI();
    void testcompute_sat_HI();
    void testcompute_k_L_HI();
    void testcompute_kr_L_HI();
    void testcompute_k_G_HI();
    void testcompute_kr_G_HI();
    void testcompute_interfacial_area_per_volume_HI();
   void testcompute_avg_rk_HI();
   
};

#endif
