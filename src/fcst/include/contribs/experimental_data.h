//---------------------------------------------------------------------------
//
// C++ Implementation: experimental_data
//
// Description: Used to read in an array of experimental data with column headers.
// 		For use with parameter estimation or generating polarization curves.
//
// Author: Marc Secanell <secanell@ualberta.ca>, Peter Dobson <pdobson@ualberta.ca> (C) 2009,2010
//
// Copyright: See COPYING file that comes with this distribution
//
// Created by: Marc Secanell, Peter Dobson (December 2009)
// Last updated: January, 2010
//
//---------------------------------------------------------------------------

#ifndef _FUELCELL__EXPERIMENTAL_DATA_H
#define _FUELCELL__EXPERIMENTAL_DATA_H

//appframe
#include <application_core/fe_vectors.h>
#include <application_core/application_data.h>
#include <application_core/event.h>
#include <application_core/application_wrapper.h>

//STL libraries:
#include<string>
#include<fstream>
#include<iostream>
#include<vector>

using namespace dealii;

namespace SIM
{
  /**
   * Description: Used to read in an array of experimental data with column headers.
   * For use with parameter estimation or generating polarization curves.Define the info you need
   */
  class ExperimentalData
  {
  public:

    /**
     * Constructor
     */
    ExperimentalData(std::string filename);
    /**
     * Destructor
     */
    ~ExperimentalData()
	{};

    /**
     * Function that parses the data names and compares them against a list of acceptable parameters
     */
    void parse_data_names() const;

    /**
     * Function that prints the data header names and data values
     */
    void print_data() const;
    /**
     * Function to extract data from the experimental_values array into a vector.
     * Removes the extracted data from the array, and the name (label) from the data_names vector
     */
    void extract_vector(std::string label, std::vector<double>& data_vector);

    /**
     * function that gets the string vector and value array
     */
    void get_experimental_values(std::vector<std::string>& name_OC, std::vector< std::vector<double> >& value_OC) const;

//     /**
//      * function that gets the string vector and value array
//      */
//     void store_OCV_values(std::vector<std::string>& name_OC, std::vector< std::vector<double> >& value_OC, std::vector<double>& current);

//     /**
//      * function that gets the string vector and value array
//      */
//     bool find_OCV(std::vector<double>& conditions, double& OCV);

   // Inline Functions
    /**
     * function that returns the number of data points (rows) in the data file
     */
    inline int get_num_rows() const
	{
	return experimental_values.size();
	}
    /**
     * function that returns the number of data points (rows) in the data file
     */
    inline int get_num_columns() const
	{
	return experimental_values[0].size();
	}


  private:
     /**
     * Function that reads the input file and populates the data header names vector
     * and data values array
     */
    void read_data_file();

    /**
     * array of experimental values used in NLS parameter study
     */
    std::vector< std::vector<double> > experimental_values;
    /**
     * names of operating conditions for experimental values used in NLS parameter study
     */
    std::vector<std::string> data_names;
    /**
     * String specifying the data file to be read
     */
    const std::string data_file;
    /**
     * vector of values holding the open circuit voltage for the data provided
     */
    std::vector<double> OCV_list;

    std::vector <std::vector<double> > OCV_conditions;

  };

}

#endif



