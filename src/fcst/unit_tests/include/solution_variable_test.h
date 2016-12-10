/**
 * A unit test class that tests the FCSTdatabase class.
 *
 */



#ifndef _FCST_SolutionVariableTest_TESTSUITE
#define _FCST_SolutionVariableTest_TESTSUITE

#include <cpptest.h>
#include <utils/fcst_db.h>
#include <utils/fcst_utilities.h>
#include <stdexcept>
#include <layers/base_layer.h>

class SolutionVariableTest: public Test::Suite
{
    public:
    SolutionVariableTest()
    {
        //Add a number of tests that will be called during Test::Suite.run()
        //Generic cases
        TEST_ADD(SolutionVariableTest::testAtOperator);
        TEST_ADD(SolutionVariableTest::testName);
        TEST_ADD(SolutionVariableTest::testSize);
        TEST_ADD(SolutionVariableTest::testSolMap);


    }
    protected:
        virtual void setup(); // setup resources using SolutionVariables 3 different constructors
        virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
    private:
        std::vector< FuelCellShop::SolutionVariable> sols;
        FuelCellShop::SolutionMap map;
        std::vector<double> xO2_dummy;

        void testAtOperator();
        void testName();
        void testSize();

        void testSolMap();

};

#endif
