/**
 * A unit test class that tests the fcst_units class. 
 * 
 */



#ifndef _FCST_AnalyticalAgglomerate_TESTSUITE
#define _FCST_AnalyticalAgglomerate_TESTSUITE

#include <map>
#include <cpptest.h>
#include <string>
#include <layers/catalyst_layer.h>
#include <reactions/tafel_kinetics.h>
#include <utils/operating_conditions.h>
#include <boost/shared_ptr.hpp>
#include <chrono>
#include <utils/fcst_utilities.h>
#include <exception>

class AnalyticalAgglomerateTest: public Test::Suite
{
public:
	AnalyticalAgglomerateTest()
    {
	   //Add a number of tests that will be called during Test::Suite.run()
	   //Generic cases
        TEST_ADD(AnalyticalAgglomerateTest::testO2CurrentDensity);
        TEST_ADD(AnalyticalAgglomerateTest::test02CurrentDerivative);
        #ifndef _OPENMP // The following tests are removed for parallel execution
            TEST_ADD(AnalyticalAgglomerateTest::testInvalidKineticsDT);
            TEST_ADD(AnalyticalAgglomerateTest::testInvalidKineticsORR);
        #endif


    }
protected:
    virtual void setup() {} // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
private:


    void testO2CurrentDensity();
    void test02CurrentDerivative();
    void testInvalidKineticsDT();
    void testInvalidKineticsORR();

    FuelCell::OperatingConditions OC;
    boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > CCL;
};

#endif
