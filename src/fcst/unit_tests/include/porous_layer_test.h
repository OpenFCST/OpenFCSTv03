//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: 
//    - Description: 
//    - Developers: 
//    - $Id: porous_layer_test.h 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------
#ifndef _Porous_Layer_Test
#define _Porous_Layer_Test

#include <cpptest.h>
#include <boost/lexical_cast.hpp>

#include <layers/porous_layer.h>
#include <layers/design_fibrous_GDL.h>
#include <catalyst_layer.h>
#include <application_core/fcst_variables.h>

/**
 * A unit test class that tests the FCST enumeration provided in system managment.
 *
 */
class PorousLayerTest: public Test::Suite
{
    public:
    PorousLayerTest()
    {
       //Add a number of tests that will be called during Test::Suite.run()
       //Generic cases
        TEST_ADD(PorousLayerTest::molecular_diffusion_test);
        TEST_ADD(PorousLayerTest::Knudsen_diffusion);
        TEST_ADD(PorousLayerTest::Knudsen_diffusion_derivatives);
        TEST_ADD(PorousLayerTest::diffusion_test);

    }
    protected:
        virtual void setup(); // setup resources... called before Test::Suite.run() ..not implemented for this test suite
        virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
        
        
    private:
        
        void molecular_diffusion_test();
        void Knudsen_diffusion();
        void Knudsen_diffusion_derivatives();
        void diffusion_test();
        
        boost::shared_ptr< FuelCellShop::Layer::GasDiffusionLayer< dim > > layer;
        FuelCellShop::Material::Oxygen oxygen;
        FuelCellShop::Material::Nitrogen nitrogen;
        
        
};

#endif
