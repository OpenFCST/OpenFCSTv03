// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: application_step3_test.h
// - Description: Test for application. This application solves Step-3 in the deal.II tutorial
// - Developers: Marc Secanell,    University of Alberta
// - $Id: application_step3_test.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _APPLICATION_STEP3_TESTSUITE
#define _APPLICATION_STEP3_TESTSUITE

#include <cpptest.h>
#include <applications/app_step3.h>

namespace FuelCell
{
    namespace UnitTest
    {
        
        class ApplicationStep3Test: public Test::Suite
        {
        public:
            ApplicationStep3Test()
            {
                //Add a number of tests that will be called during Test::Suite.run()
                //Generic cases
                TEST_ADD(ApplicationStep3Test::runApplication);
            }
        protected:
            virtual void setup()     {} // setup resources... called before Test::Suite.run() ..not implemented for this test suite
            virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
        private:
            //Generic cases
            void runApplication();
            
            
        };
    }
}

#endif
