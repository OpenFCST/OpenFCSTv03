/**
 * A unit test class that tests the fcst_units class. 
 * 
 */



#ifndef _FCST_UNITS_TESTSUITE
#define _FCST_UNITS_TESTSUITE

#include <cpptest.h>
#include <utils/fcst_units.h>
#include <stdexcept>



class FcstsUnitsTestSuite: public Test::Suite
{
public:
    FcstsUnitsTestSuite()
    {
		//Add a number of tests that will be called during Test::Suite.run()
		//Generic cases
        TEST_ADD(FcstsUnitsTestSuite::perBigToSmallTest);
        TEST_ADD(FcstsUnitsTestSuite::bigToSmallTest);
        TEST_ADD(FcstsUnitsTestSuite::perSmallToBig);
        TEST_ADD(FcstsUnitsTestSuite::smallToBig);

        //specific Cases
		TEST_ADD(FcstsUnitsTestSuite::btuToKwh);
		TEST_ADD(FcstsUnitsTestSuite::kwhToBtu);

		//Expected error
		TEST_ADD(FcstsUnitsTestSuite::testExceptions);
    }
protected:
    virtual void setup()     {} // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
private:
	//Generic cases
    void perBigToSmallTest();
    void bigToSmallTest();
    void perSmallToBig();
    void smallToBig();

    //Specific cases
    void btuToKwh();
    void kwhToBtu();

    //Expected error
    void testExceptions();

};

#endif
