/**
 * A unit test class that tests the design_MPL class. 
 * 
 * \author Madhur Bhaiya
 * 
 */



#ifndef _FCST_DesignMPL_TESTSUITE
#define _FCST_DesignMPL_TESTSUITE

#include <cpptest.h>
#include <boost/lexical_cast.hpp>
#include <string.h>
#include <layers/design_MPL.h>
#include <layers/design_fibrous_GDL.h>



class DesignMPLTest: public Test::Suite
{
public:
    DesignMPLTest()
    {
        //Add a number of tests that will be called during Test::Suite.run()
        //Generic cases
        TEST_ADD(DesignMPLTest::testPCapillary);
        TEST_ADD(DesignMPLTest::testInterfacialArea);
        TEST_ADD(DesignMPLTest::testKumburLiquidPermeability);
    }
protected:
    virtual void setup() ; // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
    
private:
    
    boost::shared_ptr< FuelCellShop::Layer::GasDiffusionLayer<dim> > dummy; // somehow MPL alone is not working...
    boost::shared_ptr< FuelCellShop::Layer::MicroPorousLayer<dim> > layer;
    
    /** Parameter file */
    ParameterHandler param;
    
    void testPCapillary();
    void testInterfacialArea();
    void testKumburLiquidPermeability();


};

#endif
