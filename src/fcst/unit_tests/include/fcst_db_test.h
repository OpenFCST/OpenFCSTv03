/**
 * A unit test class that tests the FCSTdatabase class.
 *
 */



#ifndef _FCST_FCSTdatabaseTest_TESTSUITE
#define _FCST_FCSTdatabaseTest_TESTSUITE

#include <cpptest.h>
#include <fcst_db.h>
#include <fcst_utilities.h>
#include <stdexcept>

class FCSTdatabaseTest: public Test::Suite
{
    public:
    FCSTdatabaseTest()
    {
       //Add a number of tests that will be called during Test::Suite.run()
       //Generic cases
        TEST_ADD(FCSTdatabaseTest::testCon);
        TEST_ADD(FCSTdatabaseTest::testRequest);
        TEST_ADD(FCSTdatabaseTest::testHasData);
        TEST_ADD(FCSTdatabaseTest::testGetData);
        TEST_ADD(FCSTdatabaseTest::testExceptions);
        TEST_ADD(FCSTdatabaseTest::testCommitData);
        TEST_ADD(FCSTdatabaseTest::testTolerance);
        TEST_ADD(FCSTdatabaseTest::testCreateDb);
        TEST_ADD(FCSTdatabaseTest::testValidColumnNames);
        TEST_ADD(FCSTdatabaseTest::testCompareOC);
        TEST_ADD(FCSTdatabaseTest::testCleanUp);


    }
    protected:
        virtual void setup(); // setup resources... called before Test::Suite.run() ..not implemented for this test suite
        virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
    private:
        FcstUtilities::FCSTdatabase testDb;
        std::string addr;
        void testCon();
        void testRequest();
        void testHasData();
        void testGetData();
        void testCommitData();
        void testExceptions();
        void testTolerance();
        void testCreateDb();
        void testValidColumnNames();
        void testCompareOC();
        void testCleanUp();

};

#endif
