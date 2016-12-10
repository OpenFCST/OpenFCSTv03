/**
 * A unit test class that tests the fcst_units class. 
 * 
 */

#ifndef _FCST_NumericalAgglomerateBase_TESTSUITE
#define _FCST_NumericalAgglomerateBase_TESTSUITE


#include <cpptest.h>
#include <string>
#include <microscale/numerical_agglomerate_base.h>
#include <utils/fcst_utilities.h>
#include <application_core/fcst_variables.h>


typedef FuelCellShop::MicroScale::NumericalAgglomerateBase NumAggBase;
class DummyInstance: public NumAggBase
{
public:
    virtual FuelCellShop::SolutionMap compute_current(){
        return FuelCellShop::SolutionMap();
    }
    virtual std::string get_name(){
        return "test";
    }
    virtual double aux_volume_fraction(){
        return 0.0;
    }

protected:
    virtual void set_structure(){
        //Do nothing
    }
    virtual boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> create_replica(){
        //return boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> (new FuelCellShop::MicroScale::NumericalBaseInstance());
    }
    virtual double get_film_thickness(){
        return 0.0;
    }
    virtual double get_radius(){
        return 0.0;
    }


};




class NumericalAgglomerateBaseTest: public Test::Suite
{
public:
    NumericalAgglomerateBaseTest()
    {
	   //Add a number of tests that will be called during Test::Suite.run()
	   //Generic cases
        TEST_ADD(NumericalAgglomerateBaseTest::testAV);
        TEST_ADD(NumericalAgglomerateBaseTest::testAVSetup);
        TEST_ADD(NumericalAgglomerateBaseTest::testSetSolution);

    }
protected:
    virtual void setup() {} // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
private:
    DummyInstance  agg;

    void testAVSetup();
    void testAV();
    void testSetSolution(); //These functions really test the agglomerate base

};

#endif
