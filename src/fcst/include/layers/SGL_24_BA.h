//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2012, 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: SGL_24_BA.h
//    - Description: Header file for a specific type of gas diffusion layer, i.e. SIGRACET 24 BA
//    - Developers: M. Secanell
//    - Id: $Id: SGL_24_BA.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__SGL_24_BA_H
#define _FUELCELLSHOP__SGL_24_BA_H

// FCST classes
#include <materials/base_material.h>
#include <layers/gas_diffusion_layer.h>

// Include deal.II classes
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_values.h>

//Include STL
#include<cmath>
#include<iostream>

using namespace dealii;

namespace FuelCellShop
{
    namespace Layer
    {
        /**
         * This class defines a SGL-24-BA GDL, for which effective transport properties are constant
         * 
         * @todo Effective electron and thermal properties are estimates. Better values need to be obtained from the literature
         * for accurate use of this class. The values for the electron and thermal conducitivities should be
         * added to the constructor (lines 26-34 in SGL_24_BA.cc)
         * 
         * @todo Premeability values also needed.
         */
        template <int dim>
        class SGL24BA : 
        public  GasDiffusionLayer<dim>
        {
        public:
            
            /**
             * Concrete name used for objects of this class. This name is used when
             * setting up the subsection where the data is stored in the input file.
             * 
             * The data will be store under
             * \code
             * subsection name_specified_in_constructor
             *    set Material id = 2
             *    set Gas diffusion layer type = DummyGDL # <-here I select the type of object of type GasDiffusionLayer
             *    subsection DummyGDL # <- this is the concrete_name for this class
             *       set all info relevant to this object
             *    end
             * end
             * \endcode
             */
            static const std::string concrete_name;
            
            ///@name Constructors, destructor, and initalization
            //@{    
            /** 
             * Replica Constructors
             * 
             * \warning For internal use only.
             * 
             * Constructor used only to create a prototype. Do not use
             * in general since this will not include the name of the section
             * in the parameter file you need. 
             * 
             */
            SGL24BA();
            
            /**
             * Constructor
             * 
             * \note Eventually, I would like to make this private.
             * 
             * \deprecated Use create_GasDiffusionLayer
             */
            SGL24BA(const std::string& name);
            
            /** Destructor */
            ~SGL24BA()
            {}          
            
            /**
             * 
             * \note Eventually, I would like to make this private.
             * 
             * \deprecated Use declare_all_GasDiffusionLayer_parameters
             */
            void declare_parameters (ParameterHandler &param) const
            {
                declare_parameters(this->name, param);                
            }
            
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            virtual void initialize (ParameterHandler &param);
            //@}
            
            ///@name Accessors and info
            //@{
            /**
             * Returns the effective diffusivity of the media.
             */
            void effective_gas_diffusivity(Table<2, Tensor<2,dim> >&) const;
            /**
             * Returns the effective electron conductivity of the media.
             */
            void effective_electron_conductivity(Tensor<2,dim>&) const;
            /**
             * Returns the effective thermal conductivity of the media.
             */
            void effective_thermal_conductivity(Tensor<2,dim>&) const;
            //@}
            
        private:
            ///@name Constructors, destructor, and initalization
            //@{               
            /**
             * Declare parameters for a parameter file.
             * 
             * Nothing to declare here other than the declaration of the base layer. 
             * However, the member function must be implemented so that the layer has the
             * same structure as other layers.
             * 
             */
            void declare_parameters (const std::string& name, 
                                     ParameterHandler &param) const;
            //@}
            ///@name Instance Delivery
            //@{
             /**
             * This member function is used to create an object of type gas diffusion layer
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > create_replica (const std::string &name)
            {
                return boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > (new FuelCellShop::Layer::SGL24BA<dim> (name));
            }
            /**
             * Create prototype for the layer
             */            
            static SGL24BA<dim> const* PROTOTYPE;
            //@}
            
            ///@name Internal variables
            //@{
            /** Tensor of effective porosity over tortuosity used for gas diffusivity */
            Tensor<2, dim>    porosity_over_tortuosity;
            /** Effective porosity over tortuosity used for gas diffusivity in X direction */
            double porosity_over_tortuosity_X;
            /** Effective porosity over tortuosity used for  gas diffusivity in Y direction */
            double porosity_over_tortuosity_Y;
            /** Effective porosity over tortuosity used for  gas diffusivity in Z direction */
            double porosity_over_tortuosity_Z;
            /** Tensor storing the effective electronic conductivity */
            Tensor<2, dim> electron_conductivity;
            /** Effective electronic conductivity in X direction */
            double electron_conductivity_X;
            /** Effective electronic conductivity in Y direction */
            double electron_conductivity_Y;
            /** Effective electronic conductivity in Z direction */
            double electron_conductivity_Z;
            /** Effective thermal conductivity (W/m-K) in X direction */
            Tensor<2, dim> thermal_conductivity;
            /** Effective thermal conductivity (W/m-K) in X direction */
            double thermal_conductivity_X;
            /** Effective thermal conductivity (W/m-K) in Y direction */
            double thermal_conductivity_Y;
            /** Effective thermal conductivity (W/m-K) in Z direction */
            double thermal_conductivity_Z;
            //@}
        };
    }
}

#endif
        