//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: fcst_utilities.h
//    - Description: Namespace carrying a lot of utility functions
//    - Developers: Madhur Bhaiya, Valentin N. Zingan, Marc Secanell, Phil Wardlaw
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__FCST_UTILITIES_H
#define _FUELCELLSHOP__FCST_UTILITIES_H

//Include STL and boost libraries
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cctype>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <algorithm>
#include <boost/tokenizer.hpp>
#include <stdlib.h>


//Include deal.ii classes
#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>

#include <application_core/application_wrapper.h>

//Posix access function
#include <unistd.h>

using namespace dealii;

namespace FcstUtilities
{
    /**
     * Exception shown when a particular string can't be converted to Number type.
     */
    DeclException1 (ExcWrongString,
                    std::string,
                    << "Can't convert the string " << arg1
                    << " to the desired Number type");

    /**
     * Function to convert a std::string into Number template. Currently, it is templated only for
     * \p double and \p unsigned \p int types. But it can be easily extended for other numeric data types.
     *
     * Input string should be free from blank spaces at the beginning and end.
     *
     * \note There is a similar function in \p deal.ii \em Utilities namespace, but it has certain limitations
     * as it doesn't do effective error checking for wrong input arguments. For that purpose, a new function is written here using boost::lexical_cast.
     */
    template <typename NumType>
    NumType string_to_number(const std::string&);

    /**
     * Function to convert a std::vector< std::string > into std::vector < Number template >. Currently, it is templated only for
     * \p double and \p unsigned \p int types. But it can be easily extended for other numeric data types.
     *
     * Input strings should be free from blank spaces at the beginning and end.
     *
     * \note Before extending its capabilities to deal with other numeric data types, capability of #string_to_number should be extended.
     */
    template <typename NumType>
    std::vector<NumType> string_to_number(const std::vector<std::string>&);
    
    /**
     * Function to convert a std::string into an unsigned char such as types::material_id and types::boundary_id 
     */
    template <>
    unsigned char string_to_number<unsigned char>(const std::string& str);

    /**
     * Function to convert a std::vector< std::string > into std::map< KeyType, ValueType >.
     *
     * This function is used to convert a string like "2 : 3.86" into a \p std::map where Key is 2 and Value is 3.86.
     * This function is useful to convert a \p Patterns::Map (of \p ParameterHandler class of \p deal.ii) into a useful data map.
     *
     * \warning String should contain ":" to distinguish between key and corresponding value. Also a key should not be repeated in the
     * parameter strings. Also Input strings should be free from blank spaces at the beginning and end.
     *
     * \note Currently, it is templated only to have < \p KeyType, \p ValueType > of: < \p unsigned \p int, \p double >, < \p std::string, \p double >, < \p unsigned \p int, \p std::string >, 
     * < \p unsigned \p char, \p double >. But, the capabilities can be extended.
     */
    template <typename KeyType, typename ValueType>
    std::map< KeyType, ValueType > string_to_map(const std::vector<std::string>&);

    /**
     * This function is the same as the previous one.
     * The only difference is that this function takes
     * std::string instead of std::vector<std::string>.
     */
    template< typename KeyType, typename ValueType >
    std::map< KeyType, ValueType > string_to_map(const std::string& name);

    /**
     * This function determines whether a string is a "valid" representation of a number
     *
     * Examples:
     *  Valid: "2.2", "-33", "-1.332"
     *  Invalid: "1+23.3", "2.2.2", "Bob"
     *
     * See UtilsTest::testIsNumber in unit_test/source/utils_test.cc for more information.
     */
    
    int cellId_to_index(std::string& cell_id);

    /**
     * This function extracts to parent cell index from the dealii CellID.
     * This is necessary for correlating the maps of field data in the microstructure simulations.
     * Example:
     *  CellID: 511_0:8 would return 511 which is the index of the parent cell
     */

    bool is_number(const std::string& s);

    /**
     * This function takes a std::map<KeyType, std::string>, and std::string contains text separated by a \p delimiter.
     * It splits the std::string into std::vector<ValueType>, removing all the leading and trailing spaces.
     * It returns std::map< KeyType, std::vector<ValueType> >. The default value of the \p delimiter is \p semi-colon <b>";"</b>.
     * \note It is necessary to have the delimiter "atleast once" in the \p Value string.
     */
    template< typename KeyType, typename ValueType >
    std::map< KeyType, std::vector<ValueType> > split_mapvalue_list(const std::map< KeyType, std::string >&, const char delimiter = ';');
    /**
     * This function prepares
     * an PRM file with ALL
     * your data inside.
     */
    void print_parameter_file_PRM(ParameterHandler& param, std::string path = "default.prm");
	
    /**
     * This function prepares
     * an XML file with ALL
     * your data inside.
     */
    void print_parameter_file_XML(ParameterHandler& param, std::string path = "default.xml");

    /**
     * This function returns the address of the fcst root directory.
     * It presumes that you are executing from within the data folder - an exception will be thrown otherwise.
     */
    std::string find_fcst_root();


    /**
     * This function runs external python scripts. User must provide the script name and argument.
     * It is presumed that the script resides in the python folder.
     */
    void run_python(std::string script_name,std::string arg);

    /**
     * This routine is used to parse the input parameter \param name_design_var for
     * subsections and parameter name. This values are used to set the parameters
     * specified as design variables.
     *
     * The string should be of the form:
     * @code
     * set DV_0 = subsection_1>>subsection_2>>subsection_3>>parameter name:material_ID
     * @endcode
     *
     * This routine will go to the appropriate subsection in the ParameterHandler passed
     * as input and will add the value in \param value_design_var to the ParameterHandler file
     *
     * @author M. Secanell and P. Wardlaw, 2014
     *
     */
    template <typename ValType>
    void modify_parameter_file(const std::string name_design_var,
                               const ValType value_design_var,
                               ParameterHandler& param);
    
    template <typename ValType>
    void modify_parameter_file(const std::vector<std::string> name_design_var,
                               const std::vector<ValType> value_design_var,
                               ParameterHandler& param);
    
    /**
     * Function for casting from double to string. Uses boost lexical cast.
     */
    template <typename NumType>
    std::string number_to_string(const NumType &num);

    /**
     * This function reads parameter files written in both
     * prm and xml formats.
     */
    void read_parameter_files(ParameterHandler&  param,
                              const std::string& filename);

    /**
     * Check if a file exists.
     *
     * Adapted from http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
     *
     * @param file_name The file name or relative file path.
     */
     bool file_exists(const std::string &file_name);

} //FcstUtilities
#endif //_FUELCELLSHOP__FCST_UTILITIES_H
