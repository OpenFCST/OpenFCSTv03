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
// Last updated: January 2010
//
//---------------------------------------------------------------------------

#include<contribs/experimental_data.h>

using namespace SIM;

//---------
ExperimentalData::ExperimentalData(std::string filename)
	: data_file(filename)
{
read_data_file();
}

void
ExperimentalData::read_data_file() 
{
	std::ifstream input;
	std::string line, word;
	double data;
	std::vector<double> data_line;
	
	input.open(data_file.c_str());
	if(!input.is_open())
	{
		FcstUtilities::log << "Cannot open experimental data file '"<< data_file.c_str() <<"\n";
		abort();
	}
	int first_line=0;
	while(getline(input,line))
	{
		std::istringstream stream(line);
		while (stream >> word)
		{
			if (first_line == 0)
			{
				//FcstUtilities::log << first_line << " " << word << " ";
				data_names.push_back(word);
			}
			else
			{
				{
					data = atof(word.c_str());
					//FcstUtilities::log << first_line << " " << data << " ";
					data_line.push_back(data);
				}
			}
		}
		if(first_line != 0 && !data_line.empty())
		{
			experimental_values.push_back(data_line);
			data_line.clear();
		}
		first_line++;
		parse_data_names();
	}
	input.clear();
	input.close();
}

//---------------------------------------------------------------------------
void
ExperimentalData::parse_data_names() const
{
	for(unsigned int i=0;i<data_names.size();i++)
	{
		if (	(data_names[i] != "T_cell") &&
			(data_names[i] != "RH_a") &&
			(data_names[i] != "RH_c") &&
			(data_names[i] != "P_a") &&
			(data_names[i] != "P_c") &&
			(data_names[i] != "V_cell") &&
			(data_names[i] != "Current"))
		{
		FcstUtilities::log << "Parameter #" << i+1 << ": " << data_names[i] << ", is not permited in the experimental data set"<< "\n";
		abort();
		}
	}

}

//---------------------------------------------------------------------------
void
ExperimentalData::print_data() const
{
	// testing output
	for(unsigned int i=0;i<data_names.size();i++)
	{
		FcstUtilities::log << data_names[i] << " ";
	}
	FcstUtilities::log << "\n";
	for(std::vector< std::vector<int> >::size_type i=0; i<experimental_values.size();i++)
	{
		for(std::vector<int>::size_type j=0;j<experimental_values[i].size();j++)
		{
			FcstUtilities::log << experimental_values[i][j] << " ";
		}
		FcstUtilities::log << "\n";
	}

}

//---------------------------------------------------------------------------
void
ExperimentalData::get_experimental_values(std::vector<std::string>& names, std::vector< std::vector<double> >& values) const
{
names = data_names;
values = experimental_values;
}

//---------------------------------------------------------------------------
void 
ExperimentalData::extract_vector(std::string label, std::vector<double>& data_vector)
{
	int rows = get_num_rows();
	int columns = get_num_columns();
	data_vector.resize(rows);
	bool data_exists = false;
	int data_index = 0;

	for(int i=0; i<columns;i++)
	{
		if(data_names[i] == label)
		{
			data_exists = true;
			data_index = i;
			for(int j=0;j<rows;j++)
				data_vector[j] = experimental_values[j][i];
		}
	}
	if(data_exists)
	{
		data_names.erase(data_names.begin()+data_index);
		for(int j=0;j<rows;j++)
		{
			experimental_values[j].erase(experimental_values[j].begin()+data_index);
		}
	}
	else
		FcstUtilities::log << "The vector name requested does not appear in the data" << std::endl;
}

// //---------------------------------------------------------------------------
// void
// ExperimentalData::store_OCV_values(std::vector<std::string>& OC_names, 
// 				std::vector< std::vector<double> >& OC_values,
// 				std::vector<double>& current)
// {
// double OCV_tol = 1e-4;
// bool found_OCV = false;
// bool found_V_cell = false;
// int V_cell_index;
// FcstUtilities::log << "Checking for the open circuit voltages in the data provided" << std::endl;
// 
// for(unsigned int j=0; j<OC_names.size(); ++j)
// {
// 	if (OC_names[j] == "V_cell")
// 	{
// 		found_V_cell = true;
// 		V_cell_index = j;
// 		break;
// 	}
// }
// if(!found_V_cell)
// {
// 	FcstUtilities::log << "There are no voltage values in the data set" << std::endl;
// 	return;
// }
// // Store the operating conditions for the OCV value
// //Store the value of the cell voltage for use in the application
// for(unsigned int i=0; i<current.size(); ++i)
// {
// 	if (current[i] < OCV_tol)
// 	{
// 		found_OCV = true;
// 		OCV_list.push_back(OC_values[i][V_cell_index]);
// 		OCV_conditions.push_back(OC_values[i]);
// 	}
// }
// if(!found_OCV)
// 	FcstUtilities::log << "The data provided does not contain a value corresponding to the open circuit voltage" << std::endl;
// else
// {
// 	FcstUtilities::log << "Found open circuit voltages: \n";
// 	for (unsigned int i=0; i<OCV_list.size();++i)
// 	{
// 		FcstUtilities::log << "V[" << i << "] = " << OCV_list[i] << " at: ";
// 		for (unsigned int j=0; j<OC_names.size();++j)
// 			FcstUtilities::log << OC_names[j] << " = " << OCV_conditions[i][j] << " ";
// 		FcstUtilities::log << std::endl;
// 	}	
// }
// }

// //---------------------------------------------------------------------------
// bool
// ExperimentalData::find_OCV(std::vector<double>& conditions, double& OCV)
// {
// for (unsigned int i=0; i<OCV_list.size(); ++i)
// {
// 	unsigned int j=0;
// 	for (;j<conditions.size();j++)
// 	{
// 		if(data_names[j] == "V_cell") continue; //Ignore V_cell values
// 		else if (conditions[j] != OCV_conditions[i][j]) break;
// 	}
// // If all the conditions matched, j reaches the end of conditions 
// if (j == conditions.size()) 
// {
// 	OCV = OCV_list[i];
// 	return true;
// }
// }
// 
// return false;
// }