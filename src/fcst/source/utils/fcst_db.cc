//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: fcst_db.cc
//    - Description: 
//    - Developers: P. Wardlaw
//
//---------------------------------------------------------------------------

#include <utils/fcst_db.h>

namespace NAME =  FcstUtilities;

//------------------------------------------------------------------------------------//
NAME::DatabaseOC::DatabaseOC(){
    num_param = 0;
    max_param = 10;
}

//-----------------------------------------------------------------------------------//
bool
NAME::DatabaseOC::add_param(const std::string&  name, const double& value){
    return  add_param(name, std::to_string(value));
}

//------------------------------------------------------------------------------------//
bool
NAME::DatabaseOC::add_param(const std::string&  name, const std::string&  value){
    bool answer = false;

    if(num_param < max_param){
        std::vector<std::string> line;
        line.push_back(name);
        line.push_back(value);
        param_data.push_back(line);
        answer = true;
        num_param++;
    }

    return answer;
}

//------------------------------------------------------------------------------------//
bool
NAME::DatabaseOC::compare(const DatabaseOC & other, const double & tolerance){

    bool answer = true;
    double diff;

    if(num_param == other.num_param){
        for(unsigned int i =0; i< num_param; i++){
            if(param_data[i][0] == other.param_data[i][0]){ //If the names match
                if(FcstUtilities::is_number(param_data[i][1])){
                    //Yes we are comparing numeric values

                    //Zero value case
                    if(std::atof(param_data[i][1].c_str()) == 0.0){
                        if (std::atof(other.param_data[i][1].c_str()) == 0.0)
                            diff = 0;
                    }
                    else
                        diff =((std::abs(std::atof(param_data[i][1].c_str()) - std::atof(other.param_data[i][1].c_str())))/std::abs(std::atof(param_data[i][1].c_str())));

                    //check if the value is equal within tolerance

                    if(diff > tolerance){
                        answer = false;
                        break;
                    }

                }
                else {
                    //We are comparing strings
                    if (param_data[i][1] != other.param_data[i][1]){
                        answer = false;
                        break;
                    }
                }
            }
            else{
                answer = false;
                break;
            }
        }
    }
    else
    {
        answer = false;
    }

    return answer;
}

//------------------------------------------------------------------------------------//
const int NAME::FCSTdatabase::max_lock= 40;
NAME::FCSTdatabase::FCSTdatabase() {
    connected = false;

    //In event of database lock, db interface may reattempt for up to max_lock*max_mill ~ 4 seconds
    max_milli = 100;
    lock_time = 0.0;

    //set default db path;
    db_path = FcstUtilities::find_fcst_root() + "databases/main_db";
}

//------------------------------------------------------------------------------------//
NAME::FCSTdatabase::~FCSTdatabase() {
    if (connected)
        disconnect();

    if(lock_time > 0.0){
        FcstUtilities::log << "FCSTdatabase: Db locks occurred, accumulating to " << lock_time << "[s] wasted time. \n";
        if(lock_time > 60.0)
            FcstUtilities::log << "FCSTdatabase: If locks cause a persistent problem, try using separate databases \n";
    }

}

//------------------------------------------------------------------------------------//
bool
NAME::FCSTdatabase::connect(){
    //db path is set to fcst_dir/main_db in constructor
    return connect(db_path, false);
}

//------------------------------------------------------------------------------------//
bool
NAME::FCSTdatabase::connect(const std::string&  db_path_, const bool& create_if_does_not_exist){
    db_path = db_path_;
    int status;

    if(connected){
        FcstUtilities::log << "FCSTdatabase: Database connection already established." << std::endl;
        return false;
    }

    std::ifstream ifile(db_path);

    if ((not ifile)){ //If the file does not exist
        if(create_if_does_not_exist){ //And we should make one
            if(not create_new_db(db_path_))
                return false; //We do not successfully create the db
        }
        else{
            return false; //It did not exist and we were not to create a new one
        }
    }




    status = sqlite3_open(db_path.c_str(), &db); //sqlite returns ok signal if it connects

    if(status == SQLITE_OK){
        connected =true;

        //enable WAL mode to increase performance when multiple applications are running
        std::string sqlCmd= "PRAGMA journal_mode=WAL;";
        request(sqlCmd);
        //To debug journal_mode print temp[0][0] now, if change of mode was successfull repsonse will be "wal". Make sure Sqlite3 version is 3.7 or greater.

    }
    else{
        FcstUtilities::log << "FCSTdatabase:Database not found; connection not established." << std::endl;
        connected = false;
    }


    //unable to connect or exception occured
    return connected;
}

//------------------------------------------------------------------------------------//
bool
NAME::FCSTdatabase::disconnect(){
    int rc;

    if (connected){
        rc = sqlite3_close(db); //sqlite returns true if it remains connected

        if(rc == SQLITE_OK)
           connected =false;
    }

    //unable to disconnect or exception occured
    return !connected; //logic: if we are not connected then we have disconnected
}

//------------------------------------------------------------------------------------//
bool
NAME::FCSTdatabase::has_data(const std::string&  model_name, const DatabaseOC& OC){
    bool answer = false;

    if (find_table(model_name, OC) != "NONE")
        answer = true; //we have found data for given args

    return answer;
}

//------------------------------------------------------------------------------------//
std::vector<std::vector<double>>
NAME::FCSTdatabase::get_data(const std::string&  model_name, const DatabaseOC& OC, const double& tolerance){

    return get_data(model_name, OC, tolerance, "");

}

//------------------------------------------------------------------------------------//
std::vector<std::vector<double>>
NAME::FCSTdatabase::get_data(const std::string&  model_name, const DatabaseOC& OC, const double& tolerance, const std::string&  orderby){

    std::vector<std::vector<double>> answer;
    std::vector<double> temp_line;
    //Check for data
    std::string table_name = find_table(model_name, OC, tolerance);
    if(table_name !="NONE"){
        //Get table name

        if(orderby != "")
            table_name += " ORDER BY " + orderby;
        //Get everything from the table
        std::vector<std::vector<std::string>> string_answer = request("SELECT * FROM " + table_name);

        //Change type from string to double
        for(int i = 0; i != string_answer.size(); i++){
            temp_line.clear();
            for(int j = 0; j != string_answer[i].size(); j++){

                try
                {
                    temp_line.push_back(std::atof(string_answer[i][j].c_str()));
                }
                catch (std::exception e)
                {
                    FcstUtilities::log << "FCSTdatabase:function get_data() cannot return value of double type./n";
                    FcstUtilities::log << e.what() << std::endl;
                    exit(1);
                }

            }
            answer.push_back(temp_line);
        }

    }
    return answer;

}

//------------------------------------------------------------------------------------//
std::vector<std::vector<double>>
NAME::FCSTdatabase::get_data(const std::string&  model_name, const DatabaseOC& OC){
    return get_data(model_name, OC, 0.00000001);
}

//------------------------------------------------------------------------------------//

bool
NAME::FCSTdatabase::request(const std::string&  sqlCmd, const bool& useCallBack, int lockCounter){


    bool answer = false;
    if(!connected){
        FcstUtilities::log << "FCSTdatabase:Database connection not established." << std::endl;;
        throw std::exception();
    }

    temp.clear(); //clear the temporary data

    int status;
    char *zErrMsg = 0;

    if(useCallBack)
        status = sqlite3_exec(db, sqlCmd.c_str(), callback, this, &zErrMsg);
    else
        status = sqlite3_exec(db, sqlCmd.c_str(), 0, this, &zErrMsg);

    if(status != SQLITE_OK ){


        if(status == SQLITE_ERROR){
            //For general errors attempt perform a cleanup
            cleanUp();
        }
        else if(status == SQLITE_BUSY){
            #ifdef DEBGUG
            FcstUtilities::log << "FCSTdatabase: SQL error: " << zErrMsg << ",error code:" << status << "\n" ;
            #endif

            if(lockCounter >0){
                //If a lock occurs call this function recursively after waiting a  random time

                #ifdef DEBGUG
                FcstUtilities::log << "FCSTdatabase: Reattempting request(" << lockCounter << ")\n";
                #endif

                std::srand(std::time(0));                                //Seed the random number gen
                unsigned int t = std::rand()%max_milli;                  //
                std::chrono::milliseconds dura(t);                       //up to max_milli milliseconds pause
                std::this_thread::sleep_for(dura);                       //wait for this time
                answer = request(sqlCmd, useCallBack, lockCounter -1);   //decrement lockCounter so that we only make a limited amount of attempts
                lock_time += double(t)/1000.0;
            }
            else{

                FcstUtilities::log << "FCSTdatabase: A unresolved database lock occurred. \n";
                return false;
            }

        }
        else{
            FcstUtilities::log << "FCSTdatabase: SQL error: " << zErrMsg << ",error code:" << status << "\n" ;
        }
    }
    else{
        #ifdef DEBUG
        if(lockCounter < max_lock) //10 is a hard coded default value for this function, see header definition of this function
            FcstUtilities::log << "FCSTdatabase: Success!\n";
        #endif

        answer = true;
    }

    sqlite3_free(zErrMsg);
    return answer;
}

//------------------------------------------------------------------------------------//
std::vector<std::vector<std::string>>
NAME::FCSTdatabase::request(const std::string&  sqlCmd){
    request(sqlCmd, true);
    return temp;
}

//------------------------------------------------------------------------------------//
bool
NAME::FCSTdatabase::request_no_callback(const std::string&  sqlCmd){
    //The main difference between this function and request is that here callback is simply 0, i.e. we do not read the returned data

    return request(sqlCmd, false);
}

//------------------------------------------------------------------------------------//

bool
NAME::FCSTdatabase::commit_data(const std::string&  model_name, const DatabaseOC& OC, const std::vector<std::string>& column_names, const std::vector<std::vector<double>>& data){
    bool answer = false;



    //Make some checks
    if((column_names.size() ==0) or (data.size() ==0)){
        std::string msg =  "FCSTdatabase:No data provided for " + model_name + " whilst trying to commit to db.";
        throw std::runtime_error(msg);
    }

    if(column_names.size() != data[0].size()){
        std::string msg = "FCSTdatabase: Column_names vector and data row vector must match in dimension." ;
        throw std::runtime_error(msg);
    }

    for(unsigned int a = 0; a < column_names.size(); a++){
        //Check that the column titles all have legal characters
        std::size_t found = column_names.at(a).find_first_not_of("0123456789abcdefghijklmnopqrstuvwxyz_ ");
        if (found!=std::string::npos){
            std::stringstream msg;
            msg << "FCSTdatabase: There is an illegal character in this column title '" << column_names.at(a)[found] << "'." ;
            throw std::runtime_error(msg.str());

        }
    }




    //Check model does not already exist for given operating conditions.
    if(has_data(model_name,OC)){
        FcstUtilities::log << "FCSTdatabase:Data for " + model_name + " already exists in db." << std::endl;
    }
    else{
        //insert new entry to head

        std::string table_name = make_new_head_entry(model_name,OC);

        if (table_name != "NONE"){
            if (create_table(table_name, column_names)){
                if (fill_empty_table(table_name, data, column_names))
                    answer = true;
            }

            if (answer == false){
                //Tidy up the table entry we made in the head since we failed to commit the data
                std::string  sqlCmd = "delete from head where table_ref ='" + table_name + "'";
                request_no_callback(sqlCmd);
            }
        }
    }

    return answer;




}

//------------------------------------------------------------------------------------//
std::string
NAME::FCSTdatabase::make_new_head_entry(const std::string&  model_name, const DatabaseOC& OC){
    std::string answer = "NONE"; //the answer is either the new table name or NONE

    //Find how many results there are for the model name

    //Counting head entry method
    std::vector<std::vector<std::string>> string_answer = request("SELECT COUNT(MODEL_NAME) FROM HEAD WHERE MODEL_NAME = '" + model_name + "'");
    int num = std::atoi(string_answer[0][0].c_str());
    num++;

    if(OC.param_data.size() == 0){
        FcstUtilities::log << "FCSTdatabase:Cannot commit entry with zero operating condition information." << std::endl;;
        return answer;
    }

    //Add your process id to the start of the string and the counting number to the end. This should ensure uniquness.
    answer =  model_name + std::to_string(num) + "_" + std::to_string(getpid());

    //Build  the command to create the new head entry
    std::string  sqlCmd = "INSERT INTO HEAD VALUES('" + model_name + "', '";
    int blank_OC_counter = 10; //TODO: refactor this hard coded variable to constructor

    for(int i =0; i < OC.param_data.size(); i++){
        sqlCmd += OC.param_data[i][0] + "', '" + OC.param_data[i][1]  + "','";
        blank_OC_counter--;
    }

    //we need to fill in all the blanks for OC
    for(int j =0; j < blank_OC_counter; j++){
        sqlCmd +=  "', '','";
    }
    //Finish the command string
    sqlCmd += "" + answer + "');";

    //Execute the command
    if(!request_no_callback(sqlCmd))
        answer = "NONE";
    return answer;
}

//------------------------------------------------------------------------------------//
bool
NAME::FCSTdatabase::create_table(const std::string&  table_name, const std::vector<std::string>& column_names){
    bool answer = false; //True if we successfully commit the table
    std::string  sqlCmd = "CREATE TABLE " + table_name + " (";

    if(column_names.size() >0){

        //Build the command to create the table
        sqlCmd += column_names[0] + " varchar(255) ";

        for (int i = 1; i < column_names.size(); i++)
            sqlCmd += ", " + column_names[i] + " varchar(255) ";

        sqlCmd += ");";

        if(request_no_callback(sqlCmd))
           answer = true;
        else
            FcstUtilities::log << "FCSTdatabase: Error creating table." << std::endl;
    }

    return answer;
}



//------------------------------------------------------------------------------------//
void
NAME::FCSTdatabase::cleanUp(){

    #ifdef DEBGUG
        FcstUtilities::log << "FCSTdatabase: Cleaning up database. \n";
    #endif

    //Remove invalid head entries (which have no corresponding tables)

    std::vector<std::vector<std::string>> HEAD_TABLES = request("SELECT TABLE_REF FROM HEAD");
    std::vector<std::vector<std::string>> ACTUAL_TABLES = request("SELECT name FROM sqlite_master WHERE type = 'table';");

    for(auto h: HEAD_TABLES){
        bool found = false;

        for(auto a: ACTUAL_TABLES){

            //flatten the strings
            boost::algorithm::to_lower(h[0]);
            boost::algorithm::to_lower(a[0]);

            if(a[0]  == h[0])
                found = true;

        }

        if(not found)
            request_no_callback("DELETE FROM HEAD WHERE TABLE_REF = '" + h[0] +"';");
    }


}
//------------------------------------------------------------------------------------//
bool
NAME::FCSTdatabase::fill_empty_table(const std::string&  table_name, const std::vector<std::vector<double>>& data, const std::vector<std::string>& column_names){
    bool answer = false;

    int columns = request("pragma table_info("+ table_name +");").size();
    int rows = std::atoi(request("SELECT COUNT(*) FROM "+ table_name +";")[0][0].c_str());
    std::string  sqlCmd;


    //Check Dimensions of table and data match & make sure table is empty
    if ((rows == 0) && (columns == data[0].size())){

        bool finished = false;
        int chunk_size = 480/data[0].size(); //ensure that number of selects is under 500 (SQLITE restriction)
        for(int c = 0; c < data.size() ; c+=chunk_size){
            //Create the first line of the sql command to add a row of data
            sqlCmd = "INSERT INTO "+  table_name + " SELECT ";

            for(int a = 0; a < data[c].size(); a++){
                           sqlCmd += "'" + FcstUtilities::number_to_string(data[c][a]) + "' as '" + column_names[a] + "',";
            }
            sqlCmd.resize(sqlCmd.length() -1); //remove the last ","

            //Add the data progressively
            for(int i = c+1; i < (c + chunk_size); i++){

                //Ensure chunk does not go out of bounds
                if (i == data.size()){
                    finished = true;
                    break;
                }


                sqlCmd += " UNION SELECT ";

                for(int j = 0; j < data[i].size(); j++){
                    sqlCmd += "'" + FcstUtilities::number_to_string(data[i][j]) + "',";
                }
                sqlCmd.resize(sqlCmd.length() -1); //remove the last ","


            }

            if (request_no_callback(sqlCmd))
                answer = true;

            if (finished){
                break;
            }
        }


    }



    return answer;
}

//------------------------------------------------------------------------------------//
bool
NAME::FCSTdatabase::clear_data(const std::string&  model_name, const DatabaseOC& OC){
    bool answer = false;

    std::string table_name = find_table(model_name, OC);
    if(table_name != "NONE"){
        //Clear head entry
        std::string sqlCmd = "DELETE FROM HEAD WHERE TABLE_REF ='"+table_name +"';";
        if (request_no_callback(sqlCmd))
            answer = true;


        //Drop the table
        sqlCmd = "DROP TABLE "+table_name +";";
        if (!request_no_callback(sqlCmd))
            answer = false; //if the above request is false then we must change answer from true to false
    }

    return answer;
}

//------------------------------------------------------------------------------------//
std::string
NAME::FCSTdatabase::find_table(const std::string&  model_name, const DatabaseOC& OC){
    return find_table(model_name, OC, 0.00000001);
}

//------------------------------------------------------------------------------------//
std::string
NAME::FCSTdatabase::find_table(const std::string&  model_name, const DatabaseOC& OC, const double& tolerance){

    //Get the relavent head data
    std::string sqlCmd = "SELECT * FROM HEAD WHERE MODEL_NAME = '" + model_name + "'";

    std::vector<std::vector<std::string>> head_info = request(sqlCmd);
    std::vector<std::string> head_line; // A line of the head info

    bool term_found;
    std::string table_ref;

    //for each line in head_info
    for(int i = 0; i != head_info.size(); i++) {

        head_line = head_info[i];
        //note the line's table_ref
        table_ref = head_line.back();

        //for each oc in OC
        for(int j = 0; j != OC.param_data.size(); j++) { //for all the OCs
           term_found = true;

            if (!find_matching_head_term(OC.param_data[j], head_line, tolerance))
            {
                term_found = false;
                break;
            }
        }

        if(term_found){ //We  found matching operating conditions
            return table_ref;
        }
    }

    return "NONE";
}

//------------------------------------------------------------------------------------//
bool
NAME::FCSTdatabase::find_matching_head_term(const std::vector<std::string>& OCparam, const std::vector<std::string>& head_line, const double& tolerance){

    double diff;
    //For all the elements in head_line
    for(int i =0; i < head_line.size() -1; i++){
        //If the OCparam name is equal to the element
        if (OCparam[0] == head_line[i]){

            //Is the data entry numeric?
            if(FcstUtilities::is_number(OCparam[1])){
                //Yes we are comparing numeric values

                //Zero value case
                if(std::atof(OCparam[1].c_str()) == 0.0){
                    if (std::atof(head_line[i+1].c_str()) == 0.0)
                        diff = 0;
                }
                else
                    diff =((std::abs(std::atof(OCparam[1].c_str()) - std::atof(head_line[i+1].c_str())))/std::abs(std::atof(OCparam[1].c_str())));

                //check if the value is equal within tolerance

                if(diff < tolerance) //if the values match within tolerance
                    return true;
            }
            else {
                //We are comparing strings
                if (OCparam[1] == head_line[i+1])
                    return true;
            }


        }
    }
    return false;
}

//------------------------------------------------------------------------------------//
bool
NAME::FCSTdatabase::has_data(const std::string&  model_name, const DatabaseOC& OC, const double& tolerance){
    bool answer = false;
    if (find_table(model_name, OC, tolerance) != "NONE")
        answer = true; //we have found data for given args

    return answer;
}

//------------------------------------------------------------------------------------//
bool
NAME::FCSTdatabase::create_new_db(const std::string& file_path){
    bool answer = false;



    FcstUtilities::log << "FCSTdatabase: Creating new database." << std::endl;

    //Open the database Connection (Creates a blank database)
    int status = sqlite3_open(file_path.c_str(), &db);

    if(status == SQLITE_OK){
        connected = true;
        //create the head table
        std::string sqlCMD = "CREATE TABLE 'head' ('Model_name' TEXT NOT NULL, 'OC1_name' TEXT NOT NULL, 'OC1_dat' TEXT NOT NULL, 'OC2_name' TEXT,"
                " 'OC2_dat' TEXT, 'OC3_name' TEXT, 'OC3_dat' TEXT, 'OC4_name' TEXT, 'OC4_dat' TEXT, 'OC5_name' TEXT, 'OC5_dat' TEXT, 'OC6_name' TEXT,"
                " 'OC6_dat' TEXT, 'OC7_name' TEXT, 'OC7_dat' TEXT, 'OC8_name' REAL, 'OC8_dat' TEXT, 'OC9_name' TEXT, 'OC9_dat' TEXT, 'OC10_name' TEXT,"
                " 'OC10_dat' TEXT, 'Table_REF' TEXT NOT NULL)";

        if(request_no_callback(sqlCMD)){ //True means head table was created

            //If the head table was created successfully create a test data entry
            DatabaseOC oc;
            oc.add_param("OCV", 1);
            oc.add_param("lambda", 3);

            std::vector<std::string> column_titles;
            column_titles.push_back("col1");
            column_titles.push_back("col2");

            std::vector<std::vector<double>> data;
            std::vector<double> row;
            row.push_back(1); row.push_back(3);
            data.push_back(row);
            row.clear(); row.push_back(2); row.push_back(4);
            data.push_back(row);


            if (commit_data("TEST",oc,column_titles, data))
                answer = true; //Finally we have created a template database.
        }

    }
    sqlite3_close(db);
    connected = false;


    return answer;
}

