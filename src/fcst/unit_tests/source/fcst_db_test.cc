/*
 * fcst_db_test.cc
 *
 *  Created on: Jun 21, 2013
 *      Author: wardlawp
 */

#include <fcst_db_test.h>


void FCSTdatabaseTest::setup(){
    testDb = FcstUtilities::FCSTdatabase();
    addr = FcstUtilities::find_fcst_root();
    addr += "databases/main_db";
}

void FCSTdatabaseTest::testCon(){

    std::string msg = "Error connecting to " + addr;
    TEST_ASSERT_MSG(testDb.connect() , "Error connecting to default db");
    TEST_ASSERT_MSG(testDb.disconnect() , "Error disconnecting from default db");
    TEST_ASSERT_MSG(testDb.connect(addr) , "Error disconnecting from db.");
    TEST_ASSERT_MSG(testDb.disconnect(), "Could not disconnect");
    TEST_ASSERT_MSG(testDb.connect("fail") == false , "FCSTdatabase connects to non-existent db!");




}

void FCSTdatabaseTest::testRequest(){


    testDb.connect(addr);

    std::vector<std::vector<std::string>> answer = testDb.request("SELECT OC1_NAME, OC1_DAT FROM HEAD WHERE MODEL_NAME = 'TEST';");
    std::vector<std::vector<std::string>> expected_answer;
    std::vector<std::string> first_line;
    first_line.push_back("OCV");
    first_line.push_back("1");
    expected_answer.push_back(first_line);

    TEST_ASSERT_MSG(answer == expected_answer, "Data from db interface does not match expected results.");

    testDb.disconnect();



}

void FCSTdatabaseTest::testHasData(){

    testDb.connect(addr);
    FcstUtilities::DatabaseOC OC;
    OC.add_param("OCV", 1);

    TEST_ASSERT_MSG(testDb.has_data("TEST", OC) == true, "According to db, data does not exist (but in reality it does).");

    OC.add_param("IAJWOAJSDOI", 2);

    TEST_ASSERT_MSG(testDb.has_data("TEST", OC) == false , "According to db, data exits (but in reality it doesn't).");

    testDb.disconnect();
}

void FCSTdatabaseTest::testExceptions(){
    //The following should throw an exception be cause we do not open the db before making a request
    TEST_THROWS(testDb.request("SELECT OC1_NAME, OC1_DAT FROM HEAD WHERE MODEL_NAME = 'TEST';"), std::exception);
}

void FCSTdatabaseTest::testGetData(){

    testDb.connect(addr);
    FcstUtilities::DatabaseOC OC;
    OC.add_param("OCV", 1);
    OC.add_param("lambda", 3);


    std::vector<std::vector<double>> expected_answer;
    std::vector<double> answer_first_line;
    answer_first_line.push_back(1);
    answer_first_line.push_back(3);
    expected_answer.push_back(answer_first_line);
    std::vector<double> answer_second_line;
    answer_second_line.push_back(2);
    answer_second_line.push_back(4);
    expected_answer.push_back(answer_second_line);

    std::vector<std::vector<double>> answer = testDb.get_data("TEST", OC);
    TEST_ASSERT_MSG(answer == expected_answer, "Data from table TEST does not match expected results.");

    testDb.disconnect();

}

void FCSTdatabaseTest::testCommitData(){
    testDb.connect(addr);

    std::string model_name = "CommitTest";
    FcstUtilities::DatabaseOC OC;
    OC.add_param("oc_name", 1);
    OC.add_param("fav_colour", "blue");

    std::vector<std::string> columnNames;
    columnNames.push_back("col1");
    columnNames.push_back("col2");

    std::vector<std::vector<double>> data;
    std::vector<double> answer_first_line;
    answer_first_line.push_back(1);
    answer_first_line.push_back(3);
    data.push_back(answer_first_line);
    std::vector<double> answer_second_line;
    answer_second_line.push_back(2);
    answer_second_line.push_back(4);
    data.push_back(answer_second_line);


    //Clear the data
    testDb.clear_data(model_name, OC);

    //Check that there is no data
    TEST_ASSERT_MSG(testDb.has_data(model_name, OC) == false , "According to db, data exits (but it should have been cleared).")

    //Add the data
    testDb.commit_data(model_name,OC, columnNames,data);




    //Check that the data is there
    TEST_ASSERT_MSG(testDb.has_data(model_name, OC) == true , "According to db, data does not exits (but it should).")


    //Pull the data and compare
    std::vector<std::vector<double>> answer = testDb.get_data(model_name, OC);
    TEST_ASSERT_MSG(answer == data, "Data from table CommitTest does not match expected results.");

    testDb.disconnect();
}

void FCSTdatabaseTest::testTolerance(){
    testDb.connect(addr);

    FcstUtilities::DatabaseOC OC;
    OC.add_param("OCV", 1.05);

    TEST_ASSERT_MSG(testDb.has_data("TEST", OC, 0.1) == true, "According to db, data does not exist (but in reality it does).");
    TEST_ASSERT_MSG(testDb.has_data("TEST", OC) == false , "According to db, data exits (but in reality it doesn't).");

    testDb.disconnect();


}

void FCSTdatabaseTest::testCreateDb(){
    std::string tempAddr = FcstUtilities::find_fcst_root();
    tempAddr += "databases/test_db";

    //Delete the file if it exists
    remove(tempAddr.c_str());

    //Create the new database by connecting with create_if_does_not_exist set to true
    testDb.connect(tempAddr, true);

    FcstUtilities::DatabaseOC OC;
    OC.add_param("OCV", 1);

    TEST_ASSERT_MSG(testDb.has_data("TEST", OC) == true, "According to db, data does not exist (but in reality it does).");

    OC.add_param("IAJWOAJSDOI", 2);

    TEST_ASSERT_MSG(testDb.has_data("TEST", OC) == false , "According to db, data exits (but in reality it doesn't).");

    testDb.disconnect();

    //Clean up
    remove(tempAddr.c_str());
}

void FCSTdatabaseTest::testValidColumnNames(){

    testDb.connect(addr);

    //Create the operating conditions for the submission
    std::string model_name = "FailTest";
    FcstUtilities::DatabaseOC OC;
    OC.add_param("abc", 1);

    //Create invalid column names
    std::vector<std::string> columnNames;
    columnNames.push_back("a+");

    //Create some data for submission
    std::vector<std::vector<double>> data;
    std::vector<double> answer_first_line;
    answer_first_line.push_back(1);
    data.push_back(answer_first_line);

    TEST_THROWS_MSG(testDb.commit_data(model_name,OC,columnNames,data), std::exception, "Invalid column name accepted");
    testDb.disconnect();


}

void FCSTdatabaseTest::testCompareOC(){
    FcstUtilities::DatabaseOC OC1;
    OC1.add_param("a", "b");
    OC1.add_param("c", 4);

    FcstUtilities::DatabaseOC OC2 = OC1; //copy of OC1

    FcstUtilities::DatabaseOC OC3; //Empty OC

    FcstUtilities::DatabaseOC OC4; //Similar to OC1 but different length
    OC4.add_param("a", "b");
    OC4.add_param("c", 4);
    OC4.add_param("d", 5);

    FcstUtilities::DatabaseOC OC5; //Identical to OC1
    OC5.add_param("a", "b");
    OC5.add_param("c", 4);

    FcstUtilities::DatabaseOC OC6; //Very similar to OC1
    OC6.add_param("a", "b");
    OC6.add_param("c", 4.1);

    FcstUtilities::DatabaseOC OC7; //Completly different to OC1
    OC7.add_param("a", "v");
    OC7.add_param("c", 4);

    TEST_ASSERT_MSG(OC1.compare(OC2, 0.00000001), "OCs should match");
    TEST_ASSERT_MSG(not OC1.compare(OC3, 0.00000001), "OCs should not match");
    TEST_ASSERT_MSG(not OC1.compare(OC4, 0.00000001), "OCs should not match");
    TEST_ASSERT_MSG(OC1.compare(OC5, 0.00000001), "OCs should match");
    TEST_ASSERT_MSG(not OC1.compare(OC6, 0.00000001), "OCs should not match");
    TEST_ASSERT_MSG(OC1.compare(OC6, 0.03), "OCs should match within this tolerance");
    TEST_ASSERT_MSG(not OC1.compare(OC7, 0.1), "OCs should not match");
}

void FCSTdatabaseTest::testCleanUp(){
    testDb.connect(addr);

    //Make a head entry, but don't commit corresponding data table
    std::string cmd = "INSERT INTO HEAD VALUES('cleanUpTest','joe','bloggs',";
    for(unsigned int i = 0; i < 9; i++)
        cmd +=  "'','',";
    cmd += "'incorrect_table_ref');";
    testDb.request(cmd);

    //Attempt cleanup
    testDb.cleanUp();

    //Check invalid data is no longer there
    std::vector<std::vector<std::string>> answer = testDb.request("SELECT * FROM HEAD WHERE TABLE_REF = 'incorrect_table_ref';");
    TEST_ASSERT_MSG(answer.size() ==0, "Cleanup failed");

}



