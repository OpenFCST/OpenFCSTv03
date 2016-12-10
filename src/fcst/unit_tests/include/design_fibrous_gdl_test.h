/**
 * A unit test class that tests the design_fibrous_GDL class. 
 * 
 * \author Madhur Bhaiya
 * 
 */



#ifndef _FCST_DesignFibrousGDL_TESTSUITE
#define _FCST_DesignFibrousGDL_TESTSUITE

#include <cpptest.h>
#include <boost/lexical_cast.hpp>
#include <string.h>
#include <layers/design_fibrous_GDL.h>



class DesignFibrousGDLTest: public Test::Suite
{
public:
    DesignFibrousGDLTest()
    {
        //Add a number of tests that will be called during Test::Suite.run()
        //Generic cases
        TEST_ADD(DesignFibrousGDLTest::testPCapillary);
        TEST_ADD(DesignFibrousGDLTest::testInterfacialArea);
        TEST_ADD(DesignFibrousGDLTest::testKumburLiquidPermeability);
    }
protected:
    virtual void setup() ; // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
    
private:
    
    boost::shared_ptr< FuelCellShop::Layer::GasDiffusionLayer<dim> > layer;
    
    /** Parameter file */
    ParameterHandler param;
    
    void testPCapillary();
    void testInterfacialArea();
    void testKumburLiquidPermeability();


};

#endif
