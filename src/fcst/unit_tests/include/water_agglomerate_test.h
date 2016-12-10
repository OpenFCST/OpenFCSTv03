/**
 * A unit test class that tests the fcst_units class. 
 * 
 */



#ifndef _FCST_WaterAgglomerate_TESTSUITE
#define _FCST_WaterAgglomerate_TESTSUITE

#include <map>
#include <layers/base_layer.h>
#include <cpptest.h>
#include <boost/lexical_cast.hpp>
#include <string>
#include <microscale/agglomerate_water_1D.h>
#include <layers/catalyst_layer.h>
#include <reactions/tafel_kinetics.h>
#include <utils/operating_conditions.h>
#include <boost/shared_ptr.hpp>

class WaterAgglomerateTest: public Test::Suite
{
public:
	WaterAgglomerateTest()
    {

        TEST_ADD(WaterAgglomerateTest::testO2CurrentDensity);
        TEST_ADD(WaterAgglomerateTest::test02CurrentDerivative);
    }
protected:
    virtual void setup() {} // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
private:


    void testO2CurrentDensity();
    void test02CurrentDerivative();

    //Layer
    boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > CCL;

    FuelCell::OperatingConditions OC;

};

#endif
