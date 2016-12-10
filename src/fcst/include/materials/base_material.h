//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: base_material.h
//    - Description: Base material class
//    - Developers: M. Secanell, Madhur Bhaiya, Valentin N. Zingan
//    - $Id: base_material.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__BASE_MATERIAL__H
#define _FUELCELLSHOP__BASE_MATERIAL__H

// Include deal.II classes
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/subscriptor.h>

// Include STL
#include<cmath>
#include<iostream>

// Include OpenFCST routines:
#include <application_core/fcst_variables.h>
#include <application_core/system_management.h>

using namespace dealii;

namespace FuelCellShop
{
    namespace Material
    {
        /**
         * Virtual class used to provide the interface for all material classes.
         * 
         * No object of type BaseMaterial should ever be created. 
         * 
         * \todo create a proper tree for materials. I think it should go, material_base, then childs: electrolyte, catalyst support, catalyst (...)
         * 
         * @author Secanell
         * @author Bhaiya
         * @author Zingan
         * 
         */
        class BaseMaterial : public dealii::Subscriptor
        {
        public:
            
            ///@name Initalization
            //@{    
            /**
             * Set the names of FCST solution variables
             * with respect to which you would like
             * to compute the derivatives of material properties.
             */
            void set_derivative_flags(const std::vector<VariableNames>& flags)
            {
                derivative_flags = flags;
            }

            //@}
            
            ///@name Accessors and info
            //@{
            /**
             * Return the name of the layer.
             */
            inline const std::string& name_material() const
            {
                return name;
            }
            
            /**
             * This function prints out
             * the material properties.
             */
            virtual void print_material_properties() const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            }
            //@}
        protected:
            
            ///@name Constructors, destructor, and parameter initalization
            //@{    
            /**
             * Constructor.
             */
            BaseMaterial()
            :
            name("no_name")
            { }

            /**
             * Constructor.
             * 
             * \deprecated Material classes should not take any string. The string is not
             * used for anything.
             */
            BaseMaterial(const std::string& name)
            : name(name)
            { }

            /**
             * Destructor.
             */
            virtual ~BaseMaterial()
            { }
            
            /**
             * Declare parameters for a parameter file.
             * 
             * \warning This is a PureFunction and it does not declare anything, so please do not call this function in the children.
             */
            virtual void declare_parameters(ParameterHandler&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            }
            
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             * 
             * \warning This is a PureFunction and it does not initialize anything, so please do not call this function in the children.
             * 
             */
            virtual void initialize(ParameterHandler&)
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            }
            //@}
            
            ///@name Internal variables
            //@{            
            /** Name of the layer. This value is used as a header in a subsection for the parameter file
             * where all the information concerning this layer is held.*/
            const std::string name;
            
            /** Flags for derivatives: These flags are used to request derivatives of material properties.*/
            std::vector<VariableNames> derivative_flags;
            //@}
        };
        
    } // Material
    
}  // FuelCellShop
    
    #endif