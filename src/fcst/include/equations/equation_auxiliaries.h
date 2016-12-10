// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: equation_auxiliaries.h
// - Description: This is a base class for all auxiliary openFCST equation functions.
// - Developers: Aslan Kosakian,        University of Alberta
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_EQUATION_EQUATION_AUXILIARIES_H_
#define _FCST_FUELCELLSHOP_EQUATION_EQUATION_AUXILIARIES_H_

#include <boost/shared_ptr.hpp>

#include <application_core/system_management.h>
#include <application_core/dof_application.h>
#include <application_core/initial_and_boundary_data.h>
#include <layers/base_layer.h>
#include <utils/fem_extras.h>
#include <utils/fcst_utilities.h>
#include <utils/fcst_constants.h>
#include <application_core/application_base.h>

using namespace dealii;
using namespace FuelCell::ApplicationCore;

namespace FuelCellShop
{
    namespace Equation
    {
        /**
         * This simple structure stores certain information
         * regarding a particular variable for the equation (all of them retrieved from #SystemManagement).
         *
         * The purpose of defining such a structure is to avoid
         * repeated use of #SystemManagement functions, (which use \p find methods),
         * hence improving code speed.
         *
         * This structure will normally be filled in #make_assemblers_generic_constant_data of
         * derived classes.
         *
         * \author Madhur Bhaiya, 2013
         */
        
        struct VariableInfo
        {
            /**
             * Index of the user-defined solution variable, retrieved from #SystemManagement.
             */
            unsigned int solution_index;
            
            /**
             * Block index of the matrix relating to the variable corresponding to an equation, retrieved from #SystemManagement.
             */
            unsigned int block_index;
            
            /**
             * Index corresponding to type of \p fevalue object used for this variable.
             * This information is also retrieved from #SystemManagement, from a group of actual \p FEVALUES objects (in \b AppFrame).
             */
            unsigned int fetype_index;
            
            /**
             * Boolean storing whether indices exist or not.
             * This boolean member is very useful for checks inside the derived equation classes, to
             * avoid errors.
             */
            bool indices_exist;
        };
        
        namespace DebugTools
        {
            
            template<int dim>
            class DebugOutput : public Subscriptor
            {
            public:
                
                /**
                 *      Constructor.
                 */
                DebugOutput(FuelCell::SystemManagement& sys_management, boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data);      
                               
                /**
                 *       Outputs a vector that contains cell matrix written row by row.
                 */
                void output_matrix(FuelCell::ApplicationCore::MatrixVector                                  cell_matrices,
                                   const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                   FuelCellShop::Equation::VariableInfo                                     xi
                                  );
                
                /**
                 *       Outputs a vector that contains cell residual (right hand side).
                 */
                void output_vector(FuelCell::ApplicationCore::FEVector& cell_res);
                
                /**
                 * Destructor.
                 */
                virtual ~DebugOutput();      
                
            private:           
                
                ///@name Debug output
                //@{                    
                /**
                 * Number of the cell. Used for RHS.
                 */
                int cell_number=0;    
                
                /**
                 * Number of the cell. Used for matrices.
                 */
                int cell_number_m=0;              
                
                /**
                 * File name.
                 */
                std::string fname="";
                
                /**
                 * Variable used for streaming to file.
                 */
                std::ofstream filestream;
                
                /**
                 * FEVectors object used for output.
                 */
                FEVectors out;
                
                /**
                 * Contains system matrix written row by row as column.
                 */
                FEVector matrix_by_rows;
                //@}            
              
                ///@name Generic data
                //@{                     
                /**
                 * Pointer to the external YourApplication<dim>::system_management object.
                 */
                FuelCell::SystemManagement* system_management;     
                
                /**
                 * Data object for the application data to be passed to the equation classes
                 */
                boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data;                
                //@}            
            };
        }
        
    }
}

#endif