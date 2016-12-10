//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: fcst_utilities.cc
//    - Description: Namespace carrying a lot of utility functions
//    - Developers: Madhur Bhaiya, Valentin N. Zingan
//
//---------------------------------------------------------------------------

#include "utils/fcst_utilities.h"

template <typename NumType>
NumType FcstUtilities::string_to_number(const std::string& str)
{
    try
    {
        return boost::lexical_cast<NumType>(str);
    }
    catch ( boost::bad_lexical_cast const& )
    {
        AssertThrow ( false, ExcWrongString(str) );
    }
}

//---------------------------------------------------------------------------
template <>
unsigned char FcstUtilities::string_to_number<unsigned char>(const std::string& str)
{  
    return (unsigned char)atoi(str.c_str());    
}
    
//---------------------------------------------------------------------------


bool FcstUtilities::is_number(const std::string& s)
{
    //Only allow for one encounter of a sign of decimal
    //Rules:
    //  The sign must be at the start
    //  There may only be one sign and one decimal place
    //  No other alpha characters are allowed

    bool sign_encountered = false;
    bool decimal_encountered = false;

    std::string::const_iterator it = s.begin();
    while (it != s.end()){

        if(!std::isdigit(*it)){
            if((*it == '-' or *it == '+') and not sign_encountered and it == s.begin()){
                sign_encountered = true;
            }
            else if(*it == '.' and not decimal_encountered){
                decimal_encountered = true;
            }
            else
                break;

        }

        ++it;
    }


    return !s.empty() && it == s.end();
}

//---------------------------------------------------------------------------

template <typename NumType>
std::string FcstUtilities::number_to_string(const NumType &num)
{
    try
    {
        return boost::lexical_cast<std::string>(num);
    }
    catch ( boost::bad_lexical_cast const& )
    {
        AssertThrow ( false, ExcInternalError());
    }
}

//---------------------------------------------------------------------------

template <typename NumType>
std::vector<NumType> FcstUtilities::string_to_number(const std::vector<std::string>& vec)
{
    std::vector<NumType> temp( vec.size() );
    for (unsigned int i=0; i<vec.size(); ++i)
        temp[i] = string_to_number<NumType>(vec[i]);

    return temp;
}

//---------------------------------------------------------------------------

template <typename KeyType, typename ValueType>
std::map< KeyType, ValueType > FcstUtilities::string_to_map(const std::vector<std::string>& vec)
{
    std::map< KeyType, ValueType > return_map;
    for (std::vector<std::string>::const_iterator p = vec.begin(); p != vec.end(); ++p)
    {
        //Checking if ":" exists in a string or not.
        std::string key_str = *p;
        AssertThrow ( key_str.find(":") != std::string::npos, dealii::StandardExceptions::ExcMessage("Incorrect map string (" + key_str + "); No '':'' found for FcstUtilities::string_to_map.") );
        key_str.erase( key_str.find(":"), std::string::npos );
        while ( (key_str.length() > 0) && (std::isspace(key_str[key_str.length()-1])) )
            key_str.erase( key_str.length()-1, 1 );
        KeyType key_num = string_to_number<KeyType>(key_str);

        //Checking if a particular key is already defined earlier or not.
        AssertThrow ( return_map.find(key_num) == return_map.end(), dealii::StandardExceptions::ExcMessage("Key value in <" + key_str + "> appears more than once in FcstUtilities::string_to_map.") );

        std::string value_str = *p;
        value_str.erase (0, value_str.find(":")+1);
        while ( (value_str.length() > 0) && (std::isspace(value_str[0])) )
            value_str.erase(0,1);
        ValueType value_num = string_to_number<ValueType>(value_str);

        return_map[key_num] = value_num;
    }

    return return_map;
}

//---------------------------------------------------------------------------

template< typename KeyType, typename ValueType >
std::map< KeyType, ValueType >
FcstUtilities::string_to_map(const std::string& name)
{
    const std::vector<std::string> vec = dealii::Utilities::split_string_list(name);

    std::map< KeyType, ValueType > return_map;
    for (std::vector<std::string>::const_iterator p = vec.begin(); p != vec.end(); ++p)
    {
        //Checking if ":" exists in a string or not.
        std::string key_str = *p;
        AssertThrow ( key_str.find(":") != std::string::npos, dealii::StandardExceptions::ExcMessage("Incorrect map string (" + key_str + "); No '':'' found for FcstUtilities::string_to_map.") );
        key_str.erase( key_str.find(":"), std::string::npos );
        while ( (key_str.length() > 0) && (std::isspace(key_str[key_str.length()-1])) )
            key_str.erase( key_str.length()-1, 1 );
        KeyType key_num = string_to_number<KeyType>(key_str);

        //Checking if a particular key is already defined earlier or not.
        AssertThrow ( return_map.find(key_num) == return_map.end(), dealii::StandardExceptions::ExcMessage("Key value in <" + key_str + "> appears more than once in FcstUtilities::string_to_map.") );

        std::string value_str = *p;
        value_str.erase (0, value_str.find(":")+1);
        while ( (value_str.length() > 0) && (std::isspace(value_str[0])) )
            value_str.erase(0,1);
        ValueType value_num = string_to_number<ValueType>(value_str);

        return_map[key_num] = value_num;
    }

    return return_map;
}

//---------------------------------------------------------------------------
int FcstUtilities::cellId_to_index(std::string& cell_id)
{
    boost::char_separator<char> sep{"_"};
    boost::tokenizer<boost::char_separator<char>> tokens(cell_id, sep);
    return boost::lexical_cast<int>(*(tokens.begin()));
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

template< typename KeyType, typename ValueType >
std::map< KeyType, std::vector<ValueType> >
FcstUtilities::split_mapvalue_list(const std::map< KeyType, std::string >& str_map, const char delimiter)
{
    std::map< KeyType, std::vector<ValueType> > return_map;

    for ( typename std::map< KeyType, std::string >::const_iterator iter = str_map.begin(); iter != str_map.end(); ++iter)
    {
        AssertThrow ((iter->second).find(delimiter) != std::string::npos, dealii::StandardExceptions::ExcMessage("Delimiter \"" + std::string(1, delimiter) + "\" should exist alteast once in the map value string."));
        return_map[iter->first] = string_to_number<ValueType>( dealii::Utilities::split_string_list(iter->second, delimiter) );
    }

    return return_map;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

std::string
FcstUtilities::find_fcst_root()
{

    #ifndef __gnu_linux__
        #error "Linux specific function in find_fcst_root(), function must be reimplemented for different OS. See http://stackoverflow.com/a/1024937"
    #endif

    std::string binaryDirname;
    const char *linkName = "/proc/self/exe";
    const size_t bufSize = 500;
    char dirNameBuffer[bufSize];
    const int ret = int(readlink(linkName, dirNameBuffer, bufSize - 1));

    std::string error_msg =  "Cannot find FCST Binary.";

    if (ret == -1) {
        AssertThrow( false , ExcMessage(error_msg) );
    }

    dirNameBuffer[ret] = 0;
    binaryDirname = dirNameBuffer;
    int string_position = binaryDirname.find("bin");
    if (string_position == std::string::npos){
        error_msg = "Cannot find FCST root. Executable has been moved from bin folder.";
        AssertThrow( false , ExcMessage(error_msg) );
    }

    return binaryDirname.substr(0, string_position);

}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void 
FcstUtilities::print_parameter_file_PRM(ParameterHandler& param, std::string path)
{
    std::filebuf fb;
    fb.open(path, std::ios::out);
    std::ostream out(&fb);
    param.print_parameters(out,
                            ParameterHandler::OutputStyle::Text);
    fb.close();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void
FcstUtilities::print_parameter_file_XML(ParameterHandler& param, std::string path)
{
    std::filebuf fb;
    fb.open(path, std::ios::out);
    std::ostream out(&fb);
    param.print_parameters(out,
                            ParameterHandler::OutputStyle::XML);
    fb.close();

    //Call python script to remove certain sections  from xml parameter file
    //run_python("ModifyXML.py", "default.xml");
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void
FcstUtilities::run_python(std::string script_name,std::string arg)
{
    //TODO: Implement error checking (check if script exists)
    std::string path= FcstUtilities::find_fcst_root();
    path += "python/" + script_name;
    std::string cmd = "python " + path +" " + arg;
    system(cmd.c_str());
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Mainly used in DakotaApplication
//---------------------------------------------------------------------------
template <typename NumType>
void
FcstUtilities::modify_parameter_file(const std::string name_design_var,
                                     const NumType value_design_var,
                                     ParameterHandler& param)
{
    // Parse the input parameter for subsections and parameter name
    std::vector< std::string > subsection_names;

    // Skip delimiters at beginning.
    std::string::size_type lastPos = name_design_var.find_first_not_of(">>", 0);
    // Find first "non-delimiter".
    std::string::size_type pos = name_design_var.find_first_of(">>", lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        subsection_names.push_back(name_design_var.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = name_design_var.find_first_not_of(">>", pos);
        // Find next "non-delimiter"
        pos = name_design_var.find_first_of(">>", lastPos);
    }

    // Enter appropriate subsection
    for (unsigned int i = 0; i < subsection_names.size() -1; i++)
    {
        //FcstUtilities::log<< "Entering " << subsection_names.at(i)<<std::endl;
        param.enter_subsection(subsection_names.at(i));
    }
    // Modify parameter value in the subsection.
    try
    {
        std::string parameter_name(subsection_names.at(subsection_names.size() -1));

        //FcstUtilities::log<< "Setting: " << parameter_name<<std::endl;


        // See if the parameter contains an Identifier (comma separated element)
        if (parameter_name.find(":") != std::string::npos)
        {
            // -- First get the parameter name
            std::string key_str(parameter_name);
            key_str.erase( key_str.find(":"), std::string::npos );
            std::string parameter_name_short = key_str;


            // -- Now get the key
            key_str.clear();
            key_str = parameter_name;
            key_str.erase(0, parameter_name.find(":") + 1);
            unsigned int key_num = string_to_number<unsigned int>(key_str);


            // -- Once I have the ID, I now need to get the comma separated list, modify
            //   the appropriate component and set the variable in the file.
            std::string current_value = param.get(parameter_name_short);


            // -- Check file has the appropriate format, i.e. id:value, id:value
            AssertThrow(current_value.find(":") != std::string::npos,
                        ExcMessage("Parameter should have format id:value,id:value -- If only a double do not use :"));

            std::map<unsigned int, NumType> current_value_map;
            current_value_map = FcstUtilities::string_to_map<unsigned int, NumType>(Utilities::split_string_list(current_value));

            // -- Modify value:
            current_value_map[key_num]=value_design_var;

            // -- Create new string:
            std::string strToReturn; //This is no longer on the heap

            for (auto iter = current_value_map.begin(); iter != current_value_map.end(); iter++) {
                strToReturn.append(std::to_string(iter->first)); //Not a method call
                strToReturn.append(":");
                strToReturn.append(std::to_string(iter->second));
                if (iter != (--current_value_map.end()))
                    strToReturn.append(",");
            }


            param.set(parameter_name_short, strToReturn);
        }
        else
        {
            param.set(parameter_name, value_design_var);
        }

    }
    catch(std::exception &exc)
    {
         FcstUtilities::log<< "Exception on processing: " << std::endl
                 << exc.what() << std::endl
                 << "Aborting!" << std::endl
                 << "----------------------------------------------------"<< std::endl;
        abort();
    }
    // Leave subsections
    for (unsigned int j = 0; j < subsection_names.size() -1; j++)
    {
        //FcstUtilities::log<< "Leaving subsection" << std::endl;
        param.leave_subsection();
    }
}

template <typename NumType>
void
FcstUtilities::modify_parameter_file(const std::vector<std::string> name_design_var,
                                     const std::vector<NumType> value_design_var,
                                     ParameterHandler& param)
{
    for(unsigned int i = 0; i < name_design_var.size(); ++i)
        modify_parameter_file(name_design_var[i], value_design_var[i], param);
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void
FcstUtilities::read_parameter_files(ParameterHandler&  param,
                                    const std::string& filename)
{
    if(!filename.empty())
    {
        std::size_t pos = filename.find_last_of(".");
        std::string ext = filename.substr(pos);

        if( ext.compare(".prm") == 0 )
        {
            FcstUtilities::log << "Parameters: " << filename << std::endl;
            param.read_input(filename,
                             true);
        }
        else if( ext.compare(".xml") == 0 )
        {
            std::filebuf fb;
            fb.open(filename,
                    std::ios::in);
            std::istream is(&fb);
            FcstUtilities::log << "Parameters: " << filename << std::endl;
            param.read_input_from_xml(is);
        }
        else
        {
            FcstUtilities::log << "Extension " << ext << " is not handled" << std::endl;
            AssertThrow( false , ExcNotImplemented() );
        }
    }
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
bool
FcstUtilities::file_exists(const std::string &file_name)
{
    return ( access( file_name.c_str(), F_OK ) != -1 );
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantations
//---------------------------------------------------------------------------
template unsigned int FcstUtilities::string_to_number(const std::string&);
template double FcstUtilities::string_to_number(const std::string&);
//---------------------------------------------------------------------------
template std::vector<unsigned int> FcstUtilities::string_to_number(const std::vector<std::string>&);
template std::vector<double> FcstUtilities::string_to_number(const std::vector<std::string>&);
//---------------------------------------------------------------------------
template std::map<unsigned int, double> FcstUtilities::string_to_map(const std::vector<std::string>&);
template std::map<unsigned int, double> FcstUtilities::string_to_map(const std::string&);

template std::map<std::string, double> FcstUtilities::string_to_map(const std::vector<std::string>&);
template std::map<std::string, double> FcstUtilities::string_to_map(const std::string&);

template std::map<unsigned int, std::string> FcstUtilities::string_to_map(const std::vector<std::string>&);
template std::map<unsigned int, std::string> FcstUtilities::string_to_map(const std::string&);

template std::map<unsigned char, double> FcstUtilities::string_to_map(const std::vector<std::string>&);
template std::map<unsigned char, double> FcstUtilities::string_to_map(const std::string&);
//---------------------------------------------------------------------------
template std::map<unsigned int, std::vector<double> > FcstUtilities::split_mapvalue_list(const std::map< unsigned int, std::string >&, const char);
//---------------------------------------------------------------------------
template std::string FcstUtilities::number_to_string(const double &num);
template std::string FcstUtilities::number_to_string(const int &num);
//---------------------------------------------------------------------------
template void FcstUtilities::modify_parameter_file(const std::string, const long int, ParameterHandler& );
template void FcstUtilities::modify_parameter_file(const std::string, const double, ParameterHandler& );
template void FcstUtilities::modify_parameter_file(const std::string, const bool, ParameterHandler& );
//template void FcstUtilities::modify_parameter_file(const std::string, const char , ParameterHandler& );
//---------------------------------------------------------------------------
template void FcstUtilities::modify_parameter_file(const std::vector<std::string>, const std::vector<long int>, ParameterHandler& );
template void FcstUtilities::modify_parameter_file(const std::vector<std::string>, const std::vector<double>, ParameterHandler& );
template void FcstUtilities::modify_parameter_file(const std::vector<std::string>, const std::vector<bool>, ParameterHandler& );
//template void FcstUtilities::modify_parameter_file(const std::vector<std::string>, const std::vector<char>, ParameterHandler& );
