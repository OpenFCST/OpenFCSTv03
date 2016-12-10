/**
 * A unit test class that tests the FCST enumeration provided in system managment.
 *
 */



#ifndef _Enum_Test
#define _Enum_Test

#include <cpptest.h>
#include <application_core/system_management.h>

class EnumerationTest: public Test::Suite
{
    public:
    EnumerationTest()
    {
       //Add a number of tests that will be called during Test::Suite.run()
       //Generic cases
        TEST_ADD(EnumerationTest::testDefault);
        TEST_ADD(EnumerationTest::testComparisons);

    }
    protected:
        virtual void setup() {} // setup resources... called before Test::Suite.run() ..not implemented for this test suite
        virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
    private:

        void testDefault();
        void testComparisons();


};

#endif
