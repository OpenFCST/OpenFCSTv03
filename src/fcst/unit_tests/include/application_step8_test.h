// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: application_step8_test.h
// - Description: Test for application. This application solves Step-8 in the deal.II tutorial
// - Developers: Marc Secanell,    University of Alberta
// - $Id: application_step8_test.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _APPLICATION_STEP8_TESTSUITE
#define _APPLICATION_STEP8_TESTSUITE

#include <cpptest.h>
#include <applications/app_step8.h>

namespace FuelCell
{
    namespace UnitTest
    {
        
        class ApplicationStep8Test: public Test::Suite
        {
        public:
            ApplicationStep8Test()
            {
                //Add a number of tests that will be called during Test::Suite.run()
                //Generic cases
                TEST_ADD(ApplicationStep8Test::runApplication);
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
