/**
 * A unit test class that tests the fcst_units class. 
 * 
 */



#ifndef _FCST_Platinum_TESTSUITE
#define _FCST_Platinum_TESTSUITE

#include <cpptest.h>
#include <boost/lexical_cast.hpp>
#include <string.h>

#include <materials/catalyst_base.h>
#include <layers/dummy_CL.h>



class PlatinumTest: public Test::Suite
{
public:
    PlatinumTest()
    {
        //Add a number of tests that will be called during Test::Suite.run()
        //Generic cases
        TEST_ADD(PlatinumTest::testSetSolution);
        
        TEST_ADD(PlatinumTest::testAlphaAnodic);
        
        TEST_ADD(PlatinumTest::testReferenceConcentration);
        
        TEST_ADD(PlatinumTest::testNeyerlin);
        
        TEST_ADD(PlatinumTest::testParthasarathy);
        
        TEST_ADD(PlatinumTest::testParthasarathyHCD);
        
    }
protected:
    
    virtual void setup(); // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down(); // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
    
private:
    
    // Test setsolution
    void testSetSolution();
    
    // Test alpha anodic
    void testAlphaAnodic();
    
    // Test ReferenceConcentration
    void testReferenceConcentration();
    
    // Test Neyerlin parameters:
    void testNeyerlin();
    
    // Test Parthasarathy parameters:
    void testParthasarathy();
    
    // Test Parthasarathy parameters:
    void testParthasarathyHCD();
    
    /** Parameter file */ 
    ParameterHandler param;
    
    /** Pointer that stores the catalyst value */
    boost::shared_ptr<FuelCellShop::Material::CatalystBase > catalyst;
    
};

#endif
