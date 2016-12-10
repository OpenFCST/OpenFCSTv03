// ------------------------------------------------------------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2016 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: scaling.cc
// - Description: This class allows scaling of any of the vector valued equations by a constant to help with 
//                stability of finding a solution.
//
// - Developers: Chad Balen, University of Alberta
//
// ------------------------------------------------------------------------------------------------------------------------------

#include "utils/scaling.h"

//---------------------------------------------------------------------------

FuelCell::Scaling::Scaling(FuelCell::SystemManagement& sys_management)
:
system_management(&sys_management)
{}

//---------------------------------------------------------------------------
FuelCell::Scaling::~Scaling()
{}

//---------------------------------------------------------------------------
void
FuelCell::Scaling::declare_parameters (ParameterHandler &param) const
{
    param.enter_subsection("Equations");
    {
        param.declare_entry("Apply scaling",
                            "false",
                            Patterns::Bool(),
                            "Bool for determining if scaling will be applied by application.");
        
        param.declare_entry("Equation matrix scaling",
                            "",
                            Patterns::Map( Patterns::Anything(), Patterns::Double() ),
                            "Specify a comma separated list of the equation name followed by a colon and then the scaling factor you wish to use."
                            "Ex. Electron Transport Equation:1e-4; to scale Electron Transport Equation by a factor of 1e-4.");
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
void
FuelCell::Scaling::initialize (ParameterHandler& param)
{   
    param.enter_subsection("Equations");
    {
        applyScaling  = param.get_bool("Apply scaling");
        ScalingMap = FcstUtilities::string_to_map<std::string, double>( param.get("Equation matrix scaling") );
    }
    param.leave_subsection();
    
}

//---------------------------------------------------------------------------
void
FuelCell::Scaling::scale_cell_matrix(MatrixVector& cell_matrices)
{
    //Initialize
    std::vector< std::string > solution_names = this->system_management->get_solution_names ();
    std::map<std::string, double>::iterator it;
    
    //Loop through elements in map of equation names (LHSScalingMap) and scale equations matrix
    for(it = ScalingMap.begin(); it != ScalingMap.end(); it++) // it->first: equation name, it->second: equation scaling factor
        for(unsigned int i = 0; i < solution_names.size(); ++i)
            if (this->system_management->matrix_block_index_exists(it->first, solution_names[i]))
                cell_matrices[this->system_management->matrix_block_index(it->first, solution_names[i])].matrix *= it->second; // Scale solution block of equation
}

//---------------------------------------------------------------------------
void
FuelCell::Scaling::scale_cell_residual(FEVector& cell_res)
{
    //Initialize
    std::map<std::string, double>::iterator it;
    
    //Loop through elements in map of solution names (RHSScalingMap) and scale residual vector
    for(it = ScalingMap.begin(); it != ScalingMap.end(); it++) // it->first: equation name, it->second: equation scaling factor
        cell_res.block(this->system_management->equation_name_to_index(it->first) ) *= it->second; // Scale solution block of residual
}

//---------------------------------------------------------------------------
void
FuelCell::Scaling::scale_bdry_matrix(MatrixVector& local_bdry_matrices, MatrixVector& bdry_matrices)
{
    //Initialize
    std::vector< std::string > solution_names = this->system_management->get_solution_names ();
    std::vector< std::string > equation_names = this->system_management->get_equation_names ();
    std::map<std::string, double>::iterator it;
    
    //Loop through elements in map of equation names (LHSScalingMap) and scale equation boundaries matrix
    for(it = ScalingMap.begin(); it != ScalingMap.end(); ++it) // it->first: equation name, it->second: equation scaling factor
        for(unsigned int i = 0; i < solution_names.size(); ++i)
            if (this->system_management->matrix_block_index_exists(it->first, solution_names[i]))
                local_bdry_matrices[this->system_management->matrix_block_index(it->first, solution_names[i])].matrix *= it->second; // Scale solution block of equation
    
    for(unsigned int i = 0; i < equation_names.size(); ++i)
        for(unsigned int j = 0; j < solution_names.size(); ++j)
            if (this->system_management->matrix_block_index_exists(equation_names[i], solution_names[j]))
                bdry_matrices[this->system_management->matrix_block_index(equation_names[i], solution_names[j])].matrix.add(1.0, local_bdry_matrices[this->system_management->matrix_block_index(equation_names[i], solution_names[j])].matrix); //Add scaled component to solution matrix
}

//---------------------------------------------------------------------------
void
FuelCell::Scaling::scale_bdry_residual(FEVector& local_bdry_res, FEVector& bdry_res)
{
    //Initialize
    std::map<std::string, double>::iterator it;
    
    //Loop through elements in map of solution names (RHSScalingMap) and scale residual vector
    for(it = ScalingMap.begin(); it != ScalingMap.end(); it++) // it->first: solution name, it->second: solution scaling factor
        local_bdry_res.block(this->system_management->equation_name_to_index(it->first)) *= it->second;
    bdry_res += local_bdry_res; //Add scaled component to solution residual
}

//---------------------------------------------------------------------------
void
FuelCell::Scaling::create_local_bdry_matrix(MatrixVector& local_bdry_matrices, MatrixVector& bdry_matrices)
{
    local_bdry_matrices.resize(bdry_matrices.size());
    for (unsigned int i = 0; i<bdry_matrices.size(); i++)
        local_bdry_matrices[i].matrix = FullMatrix<double> (bdry_matrices[i].matrix.m(), bdry_matrices[i].matrix.n() );
}
