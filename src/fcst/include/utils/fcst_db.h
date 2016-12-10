//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: FCSTdatabase
//    - Description: A class for interfacing with SQL databases,
//          for the purpose of accessing and storing simulation results.
//    - Developers: Philip Wardlaw
//
//---------------------------------------------------------------------------

#ifndef FCSTDATABASE_H_
#define FCSTDATABASE_H_

#include <utils/fcst_utilities.h>
#include  <utils/logging.h>
//SQLITE3
#include <sqlite3.h>

//STL
#include <iostream>
#include <ctime>
#include <string>
#include <cstddef>        // std::size_t
#include <vector>
#include <stdexcept>
#include <fstream>
#include <chrono>
#include <thread>
#include <random>

//Deal.ii



//Boost
#include <boost/algorithm/string.hpp>

namespace FcstUtilities
{
   /**
    * This class is for storing a list of up to 10 parameters and
    * is used as a way of informing the FCSTdatabase class
    * about the model you wish to store/access.
    */
    class DatabaseOC;

    /**
     *
     * @brief This class is for interfacing with SQL databases, for the purpose of accessing and storing simulation results.
     *
     *
     * <h3> Usage details</h3>
     * This class provides a database interface for storing generic numerical data to disc.
     * It provides tolerance based selection of data based upon explicit operating conditions.
     *
     *
     * @code
     *
     * //Create the database interface and address
     * FcstUtilities::FCSTdatabase testDb;
     * std::string address = "/some/location"
     *
     * //Connect to that database
     * testDb.connect(address);
     *
     * //Alternatively if the following line will create the database if it does not exist
     * bool create_if_does_not_exit = true;
     * testDb.connect(address, create_if_does_not_exit);
     *
     * //We want to get some data from the database
     * //First we must describe the model and operating conditions
     * //for which we would like some data
     * std::string model_name = "test";
     * FcstUtilities::DatabaseOC OC;
     * OC.add_param("OCV", 1);
     * OC.add_param("lambda", 3);
     *
     * //Data is stored in tables, i.e. 2d array or a vector of vector of doubles
     * std::vector<std::vector<double>> data;
     *
     * //Check if the database had data matching our OC within a 10% tolerance
     * if (testDb.has_data(model_name, OC, 0.1)){
     *    //If it does retrieve this data
     *    data= testDb.get_data(model_name, OC, 0.1);
     * }
     *
     * //We wish to store data to the database
     * //First we create the relevant name and OC
     * FcstUtilities::DatabaseOC new_OC;
     * std::string new_model_name = "test";
     * new_OC.add_param("OCV", 4);
     * new_OC.add_param("lambda",12);
     *
     *
     * //Create some data we wish to save.
     * std::vector<std::vector<double>> data_to_save = create_some_data();
     *
     * //Column names correspond to the data we are saving,
     * //in this case our data table looks like this:
     * //
     * //   |_x_|_phi_M_|
     * //   |0.1|  0.1  |
     * //   |0.3|  0.2  |
     * //   |0.5|  0.4  |
     * //   |0.6|  0.4  |
     * //   |1.1|  0.5  |
     * //   |2.1|  0.7  |
     * //
     * //To submit the data we must also let the database know what the column names are.
     *
     * std::vector<std::string> columnNames;
     * columnNames.push_back("x");
     * columnNames.push_back("phi_M");
     *
     * //Now we are ready to store our data in the database.
     * testDb.commit_data(new_model_name, new_OC, columnNames, data_to_save);
     *
     * //It is good practice to always close the database connection
     * //once we are done submitting/retrieving data.
     * testDb.close()
     *
     * @endcode
     *
     *
     * @author Philip Wardlaw
     * @date 2014
     */
	class FCSTdatabase {
		public:
	        /**
	         * Constructor
	         */
            FCSTdatabase();
            /**
             * Destructor
             */
            virtual ~FCSTdatabase();

            /**
             * Function for connecting to SQL database. Takes database file path as argument.
             * Returns true if connection successful
             */
            bool connect(const std::string& db_path_, const bool& create_if_does_not_exist = false);

            /**
             * Function for connecting to default SQL database, located at fcst_root/main_db.
             * Returns true if connection successful
             */
            bool connect();

            /**
            * Function for disconnecting to SQL database.
            * Returns true if disconnection successful
            */
            bool disconnect();

            /**
            * Function for testing if db that you are connected to has data for a given model and operating conditions.
            * Takes model name and operating conditions as arguments. OC is an object of DatabaseOC type.
            *
            * Returns true if data exists.
            *
            * See Unit test "testHasData" for usage demo.
            */
            bool has_data(const std::string& model_name, const DatabaseOC& OC);

            /**
            * Function for testing if db that you are connected to has data for a given model and operating conditions.
            * Takes model name, operating conditions, and a tolerance as arguments. OC is an object of DatabaseOC type.
            * Differences in operating condition values will be compared with provided tolerance (0.0 - 1.0).
            * Returns true if data exists.
            *
            * See Unit test "testHasData" for usage demo.
            */
            bool has_data(const std::string& model_name, const DatabaseOC& OC, const double& tolerance);

            /**
            * Function for getting numerical information for a given model and operating conditions.
            * Takes model name and operating conditions as arguments.  OC is an object of DatabaseOC type.
            *
            * N.b: returns an potentially unordered (by rows) list of results. See third overloading of this function.
            *
            * Returns an empty vector(vector(double)) if no data exists.
            *
            * See Unit test "testGetData" for usage demo.
            */
            std::vector<std::vector<double>> get_data(const std::string& model_name, const DatabaseOC& OC);

            /**
            * Function for getting numerical information for a given model and operating conditions.
            * Takes model name, operating conditions, and a tolerance as arguments. OC is an object of DatabaseOC type.
            * Differences in operating condition values will be compared with provided tolerance (0.0 - 1.0).
            *
            * N.b: returns an potentially unordered (by rows) list of results. See third overloading of this function.
            *
            * Returns an empty vector(vector(double)) if no data exists.
            *
            * See Unit test "testGetData" for usage demo.
            */
            std::vector<std::vector<double>> get_data(const std::string& model_name, const DatabaseOC& OC, const double& tolerance);

            /**
            * Function for getting numerical information for a given model and operating conditions.
            * Takes model name, operating conditions, orderby string,and a tolerance as arguments. OC is an object of DatabaseOC type.
            * Differences in operating condition values will be compared with provided tolerance (0.0 - 1.0).
            * Orderby is the name of the column that the results will be ordered by before being returned.
            *
            * Returns an empty vector(vector(double)) if no data exists.
            *
            * See Unit test "testGetData" for usage demo.
            */
            std::vector<std::vector<double>> get_data(const std::string& model_name, const DatabaseOC& OC, const double& tolerance, const std::string& orderby);

            /**
            * Function for committing numerical information for a given model and operating conditions.
            * Takes model name, operating conditions, data column titles, and data as arguments.
            * OC is a 2d array (vector vector) in which for a given line the first element is the operating condition name and the second element is the value.
            * Column_names is a list of column that pertains to the numerical data stored in argument data.
            * Data is a 2d array (vector vector) of doubles containing the data we whish to store.
            *
            * Returns an empty vector(vector(double)) if no data exists.
            *
            * Throws exceptions if invalid column names are provided (only characters "123456789abcdefghijklmnopqrstuvwxyz " are permitted, i.e. alpha numeric and whitespace).
            *
            * Throws exception if data or column_name vectors are empty or do not match in dimension.
            *
            * See Unit test "testCommitData" for usage demo.
            */
            bool commit_data(const std::string& model_name, const DatabaseOC& OC, const std::vector<std::string>& column_names, const std::vector<std::vector<double>>& data);

            /**
            * Wrapper Function for performing direct SQL requests. Takes a SQL statement as argument.
            * Returns the results of the request in a 2d array (vector vector) of strings.
            *
            *
            * Note: If you don't know SQL then this is probably not the function for you, try another.
            *
            * See Unit test "testRequest" for usage demo.
            */
            std::vector<std::vector<std::string>> request(const std::string& sqlCmd);

            /**
            * Function for clearing data for a given model and operating conditions.
            * Takes model name and operating conditions as arguments. OC is a 2d array (vector vector) in which for
            * a given line the first element is the operating condition name and the second element is the value.
            *
            * Returns true if data is successfully deleted.
            *
            * See Unit test "testCommitData" for usage demo.
            */
            bool clear_data(const std::string& model_name, const DatabaseOC& OC);

            /**
             * Function which performs database cleanup, removing erroneous data occurring from SQL errors.
             *
             * Currently the function removes invalid HEAD entries (entries which are missing corresponding
             * data tables).
             */
            void cleanUp();

		private:
            /*
             * SQLite object for all our db interactions
             */
            sqlite3 *db;

            /*
             * Status variables
             */
            bool connected;
            std::string db_path;

            /*
             * Temporary variable for storing requested data
             */
            std::vector<std::vector<std::string>> temp;

            /*
             * Static call back function used for interfacing with sqlite3_exec
             * Required by public function request
             */
            static int callback(void *db_, int argc, char **argv, char **azColName){

                FCSTdatabase* db = reinterpret_cast<FCSTdatabase*>(db_);
                std::vector<std::string> line_temp;
                std::string str_temp;
                for(int i=0; i<argc; i++){
                    str_temp = argv[i] ? argv[i] : "NULL";
                    line_temp.push_back(std::string(str_temp));

                }
                db->temp.push_back(line_temp);
                return 0;
            }

            /*
             * Function for creating new entry in HEAD table
             * Required by public function commit_data
             */
            std::string make_new_head_entry(const std::string& model_name, const DatabaseOC& OC);

            /*
             * Function for creating new entry in HEAD table
             * Required by public function commit_data
             */
            bool create_table(const std::string& table_name, const std::vector<std::string>& column_names);

            /*
             * Function for filling an empty table with data
             * Required by public function commit_data
             */
            bool fill_empty_table(const std::string&  table_name, const std::vector<std::vector<double>>& data, const std::vector<std::string>& column_names);

            /*
             * Function for finding a the table name for a given model name and OC
             * Required by public function has_data and clear_data
             */
            std::string find_table(const std::string&  model_name, const DatabaseOC& OC);

            std::string find_table(const std::string&  model_name, const DatabaseOC& OC, const double& tolerance);


            bool find_matching_head_term(const std::vector<std::string>& OC, const std::vector<std::string>& head_line, const double& tolerance);


            /*
             * Private wrapper function for making requests that don't require an answer other than a yes or no
             * The main difference between this function and request is that here callback is simply 0, i.e. we do not read the returned data
             * Required by several public and private functions
             */
            bool request_no_callback(const std::string&  sqlCmd);

            /*
             * Variables associated with database locking.
             */
            const static int max_lock;
            unsigned int max_milli;
            double lock_time;

            /*
             * Private function that implements the all request functionality.
             * Request(sqlCmd) and request_no_callback functions rely on this function.
             * Standard variable lockCounter is used to recursively call this function in the event of a database lock (for reattempts)
             * Bool useCallBack determines if we read the returned data
             */

            bool request(const std::string&  sqlCmd, const bool& useCallBack,  int lockCounter = max_lock);


            /**
             * Function for creating a new default db.
             */
            bool create_new_db(const std::string& filePath);
    };


    /**
     * This class is for storing a list of up to 10 parameters and
     * is used as a way of informing the FCSTdatabase class
     * about the model you wish to store/access.
     */
    class DatabaseOC{

        //FCSTdatabase is allowed to access private members of this class
        friend class FCSTdatabase;

        public:
           /**
            * Constructor
            */
            DatabaseOC();

            /**
             * Function for adding a parameter (name and value pair) to the list
             * of up to 5 parameters. Returns false if unsuccessful.
             *
             * See Unit test "testHasData" for usage demo.
             */
            bool add_param(const std::string& name, const double& value);
            bool add_param(const std::string&, const std::string& );

            /**
             * Function for clearing parameter data so one can reuse the same object more than once.
             */
            inline void clear(){
                num_param = 0;
                param_data.clear();
            }


            /**
             * Function for performing a tolerance based comparison of OC object with another OC object.
             *
             * Returns false if the two OC objects do not match within a certain tolerance.
             *
             */
            bool compare(const DatabaseOC& other, const double& tolerance);

        private:


            int num_param;
            int max_param;
            std::vector<std::vector<std::string>> param_data;

    };

}



#endif /* FCSTDATABASE_H_ */
