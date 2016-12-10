/**
 * A unit test class that tests the fcst_units class. 
 * 
 */



#ifndef _FCST_MultiScaleCLLayer_TESTSUITE
#define _FCST_MultiScaleCLLayer_TESTSUITE

#include <cpptest.h>
#include <boost/lexical_cast.hpp>
#include <string.h>
#include <layers/catalyst_layer.h>
#include <reactions/tafel_kinetics.h>



class MultiScaleCLTest: public Test::Suite
{
public:
    MultiScaleCLTest()
    {
        //Add a number of tests that will be called during Test::Suite.run()
        //Generic cases
        TEST_ADD(MultiScaleCLTest::testSetSolution);
        
    }
protected:
    virtual void setup() ; // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
    
private:
    
    boost::shared_ptr< FuelCellShop::Layer::CatalystLayer<dim> > layer; //A child of Base Layer
    std::vector<double> xO2_dummy;
    /** Parameter file */
    ParameterHandler param;
    
    void testSetSolution();
    
    
};

#endif
