// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: equation_auxiliaries.cc
// - Description: This is a base class for all auxiliary openFCST equation functions.
// - Developers: Aslan Kosakian,        University of Alberta
//
// ----------------------------------------------------------------------------

#include "equations/equation_auxiliaries.h"

namespace DT = FuelCellShop::Equation::DebugTools;

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
DT::DebugOutput<dim>::DebugOutput(FuelCell::SystemManagement& sys_management, boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data)
:
system_management(&sys_management),
data(data)
{ }

// ---            ---
// --- Destructor ---
// ---            ---
template<int dim>
DT::DebugOutput<dim>::~DebugOutput()
{ }

// ---               ---
// --- output_matrix ---
// ---               ---

template<int dim>
void
DT::DebugOutput<dim>::output_matrix(FuelCell::ApplicationCore::MatrixVector                                  cell_matrices,
                                    const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                    FuelCellShop::Equation::VariableInfo                                     xi
                                   )
{
    FuelCell::ApplicationCore::FEVector matrix_by_rows;        
    
    matrix_by_rows.reinit(this->system_management->block_info->global);
    matrix_by_rows.block(xi.solution_index).reinit((cell_info.fe(xi.fetype_index)).dofs_per_cell*(cell_info.fe(xi.fetype_index)).dofs_per_cell);
    
    unsigned int index=0;
    
    //---------------LOOP over i -----------------------------------------------------------------
    for (unsigned int i=0; i < (cell_info.fe(xi.fetype_index)).dofs_per_cell; ++i)
    {
        //--------------LOOP(s) over j-------------------------------------------------------------
        
        //-----------Assembling Matrix for terms corresponding to "xi" BLOCK------------------------
        for (unsigned int j=0; j < (cell_info.fe(xi.fetype_index)).dofs_per_cell; ++j)
        {
            matrix_by_rows.block(xi.solution_index)(index)=cell_matrices[xi.block_index].matrix(i,j);    
            index++;
        }            
    }  
    
    this->cell_number_m++;    
    this->fname="cell_matrices_"+std::to_string(this->cell_number_m)+".dat";
    
    this->filestream.open(this->fname, std::ofstream::out | std::ofstream::trunc);
    
    for (unsigned int i=0; i < (matrix_by_rows).block(xi.solution_index).size(); i++)
    {
        this->filestream << matrix_by_rows.block(xi.solution_index)(i) << "\n";
    }     
    
    this->filestream.close();
}


// ---               ---
// --- output_vector ---
// ---               ---

template<int dim>
void
DT::DebugOutput<dim>::output_vector(FuelCell::ApplicationCore::FEVector& cell_res)
{
    this->cell_number++;
    
    this->fname="cell_res_"+std::to_string(this->cell_number)+".dat";
    
    this->filestream.open(this->fname, std::ofstream::out | std::ofstream::trunc);
    
    for (unsigned int i=0; i < cell_res.size(); i++)
    {
        this->filestream << cell_res(i) << "\n";
    }     
    
    this->filestream.close();        
}


/////////////////////////////
// EXPLICIT INSTANTIATIONS //
/////////////////////////////
template class DT::DebugOutput<deal_II_dimension>;