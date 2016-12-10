/**
 * A unit test class that tests the conventional_CL class. 
 * 
 * \author Madhur Bhaiya
 * 
 */



#ifndef _FCST_ConventionalCLLayer_TESTSUITE
#define _FCST_ConventionalCLLayer_TESTSUITE

#include <cpptest.h>
#include <boost/lexical_cast.hpp>
#include <string.h>
#include <layers/conventional_CL.h>



class ConventionalCLTest: public Test::Suite
{
public:
    ConventionalCLTest()
    {
        //Add a number of tests that will be called during Test::Suite.run()
        //Generic cases
        TEST_ADD(ConventionalCLTest::testVolumeFraction);
        TEST_ADD(ConventionalCLTest::testICRatio);
        TEST_ADD(ConventionalCLTest::testNafionLoading);
        TEST_ADD(ConventionalCLTest::testPCapillary);
        TEST_ADD(ConventionalCLTest::testInterfacialArea);
        TEST_ADD(ConventionalCLTest::testKumburLiquidPermeability);
        TEST_ADD(ConventionalCLTest::testActiveArea);
        TEST_ADD(ConventionalCLTest::testPtLoadingSupport);
    }
protected:
    virtual void setup() ; // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
    
private:
    
    boost::shared_ptr< FuelCellShop::Layer::CatalystLayer<dim> > layer; //A child of Base Layer
    
    /** Parameter file */
    ParameterHandler param;
    
    void testVolumeFraction();
    void testICRatio();
    void testNafionLoading();
    void testPCapillary();
    void testInterfacialArea();
    void testKumburLiquidPermeability();
    void testActiveArea();
    void testPtLoadingSupport();

    std::vector<unsigned int> test_ids; //Material id's coresponding to graded sub layers

};

#endif
