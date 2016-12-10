//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: utils_test.cc
//    - Description: Unit testing class for FCST Utility classes
//    - Developers: Philip Wardlaw
//    - Id: $Id: utils_test.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

/**
 * A unit test class for testing the functions in fcst_utilities.cc
 * 
 * @authors Phil Wardlaw and M. Secanell
 * 
 */

#ifndef _FCST_Utils_TESTSUITE
#define _FCST_Utils_TESTSUITE

#include <cpptest.h>
#include <boost/lexical_cast.hpp>
#include <string.h>
#include <utils/fcst_utilities.h>

class UtilsTest: public Test::Suite
{
public:
    UtilsTest()
    {
        //Add a number of tests that will be called during Test::Suite.run()
        //Generic cases
        TEST_ADD(UtilsTest::testIsNumber);
        TEST_ADD(UtilsTest::testModify_parameter_file_double);
        TEST_ADD(UtilsTest::testModify_parameter_file);
        TEST_ADD(UtilsTest::testModify_parameter_file_list);

    }
protected:
    virtual void setup(){} // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down(){} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
private:
    /**
     * 
     */
    void testIsNumber();
    /**
     * Check when there is only one value of type id:value
     */
    void testModify_parameter_file_double();
    /**
     * Check when there is only one value of type id:value
     */
    void testModify_parameter_file();
    /**
     * Check when there is a list of values
     */
    void testModify_parameter_file_list();
    
};

#endif
