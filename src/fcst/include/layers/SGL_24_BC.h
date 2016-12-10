//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2012, 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: SGL_24_BC.h
//    - Description: Header file for a specific type of microporous layer, i.e. SIGRACET 24 BC
//    - Developers: Madhur Bhaiya
//    - Id: $Id: SGL_24_BC.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__SGL_24_BC_H
#define _FUELCELLSHOP__SGL_24_BC_H

// FCST classes
#include <materials/base_material.h>
#include <layers/micro_porous_layer.h>

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
         * This class defines a SGL-24-BC MPL, for which effective transport properties are constant
         * 
         * \warning This is a template class for different materials. Currently all effective properties
         * are made-up!!!
         * 
         * \author M. Secanell, 2013
         */
        template <int dim>
        class SGL24BC : 
        public  MicroPorousLayer<dim>
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
             * Replica Constructor 
             */
            SGL24BC();
            
            /**
             * Constructor
             * 
             * \note Eventually, I would like to make this private.
             * 
             * \deprecated Use create_MicroPorousLayer
             */
            SGL24BC(std::string name);
            
            /** Destructor */
            ~SGL24BC()
            {}
            
            /**
             * Declare parameters for a parameter file. 
             * 
             * \deprecated Use declare_all_MicroPorousLayer_parameters
             * 
             * \note Since there are two declare parameters, this one is hidden by the former, so it has to be
             * redefined otherwise it cannot be called :( 
             */
            void declare_parameters (ParameterHandler &param) const
            {
                declare_parameters(this->name, param);
            };
	    
	    void initialize (ParameterHandler &param)
	    {FuelCellShop::Layer::MicroPorousLayer<dim>::initialize(param);}
            //@}
            
            ///@name Accessors and info
            //@{
            /**
             * Returns the effective diffusivity of the media.
             */
            void effective_gas_diffusivity(Table<2, Tensor<2,dim> >&prop_eff) const;
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
            ///@name Declaration
            //@{    
            /**
             * Declare parameters for a parameter file.
             * 
             * Nothing to declare here other than the declaration of the base layer. 
             * However, the member function must be implemented so that the layer has the
             * same structure as other layers.
             * 
             */
            void declare_parameters (std::string name, ParameterHandler &param) const
            {   FuelCellShop::Layer::MicroPorousLayer<dim>::declare_parameters(name,param);}
            
            ///@name Instance Delivery (function)
            //@{
            /**
             * This member function is used to create an object of type micro porous layer
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Layer::MicroPorousLayer<dim> > create_replica (const std::string &name)
            {
                return boost::shared_ptr<FuelCellShop::Layer::MicroPorousLayer<dim> > (new FuelCellShop::Layer::SGL24BC<dim> (name));
            }   
            //@}
            
            ///@name Instance Delivery (Data member)
            /**
             * Create prototype for the layer
             */          
            static SGL24BC<dim> const* PROTOTYPE;
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
