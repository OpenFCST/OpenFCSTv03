//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: dummy_GDL.h
//    - Description: Implementation of a GDL class that setup us all properties from file
//    - Developers: M. Secanell
//    - $Id: dummy_GDL.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__DUMMY_GDL_H
#define _FUELCELLSHOP__DUMMY_GDL_H

// FCST classes
#include <utils/fcst_units.h>
#include <layers/gas_diffusion_layer.h>

using namespace dealii;

namespace FuelCellShop
{
    namespace Layer
    {
        /**
         * This class is used when we want to input the effective properties to the GDL directly, without taking
         * into account the structure of the GDL \\
         * UNDER DEVELOPMENT
         * @author M. Secanell, 2011-13
         */
        template <int dim>
        class DummyGDL : 
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

            
            ///@name Destructor, and initalization
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
            DummyGDL();
            
            /** Destructor */
            ~DummyGDL()
            {}

            /**
             * Declare parameters for a parameter file
             * 
             * \deprecated Use declare_all_GasDiffusionLayer_parameters
             * 
             */
            void declare_parameters (ParameterHandler &param) const
            {
                this->declare_parameters(this->name, param);
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
             * Compute effective diffusivity in the GDL. 
             * 
             * Reimplementation of parent class. See parent class for info.
             */
            virtual void effective_gas_diffusivity(std::vector< Tensor<2,dim> >& prop_eff_vec) const;

            /**
             * Compute the effective diffusivty in the GDL. This routine takes the 
             * gas diffusivity from FuelCellShop::BinaryDiffusion and transforms
             * it into an effective property taking into account the porosity and 
             * structure of the GDL
             * 
             * \deprecated Use #compute_gas_diffusion with appropriate gases and then #effective_gas_diffusivity
             */
            virtual void effective_gas_diffusivity(Table< 2, double> &) const;
                                           
            /**
             * Compute the effective diffusivty in the GDL. This routine takes the 
             * gas diffusivity from FuelCellShop::BinaryDiffusion and transforms
             * it into an effective property taking into account the porosity and 
             * structure of the GDL (Anisotropic case)
             * 
             * \deprecated Use #compute_gas_diffusion with appropriate gases and then #effective_gas_diffusivity
             */
            virtual void effective_gas_diffusivity(Table< 2, Tensor< 2, dim > > &) const;
            
            /**
             * Reimplementation of parent class. See parent class for info.
             */
            virtual void derivative_effective_gas_diffusivity(std::map< VariableNames, std::vector< Tensor<2,dim> > >&) const;
                                           
            /**
             * Compute the effective electron conductivity in the GDL
             */
            virtual void effective_electron_conductivity(double& ) const;
            /**
             * Compute the effective electron conductivity in the GDL
             */
            virtual void effective_electron_conductivity(Tensor<2,dim>& ) const;
            /**
             * Compute the effective thermal conductivity in the GDL
             */
            virtual void effective_thermal_conductivity(std::vector< Tensor<2,dim> >& prop_eff) const;
            /**
             * 
             */
            virtual void derivative_effective_thermal_conductivity(std::vector< Tensor<2,dim> >& dK_dT) const;
            //@}
        private:
            ///@name Constructors, destructor, and initalization
            //@{    
            /** 
             * Constructor
             * 
             * \deprecated Use create_GasDiffusionLayer
             */
            DummyGDL(std::string name);
            
            /**
             * Declare parameters for a parameter file.
             * 
             */
            virtual void declare_parameters (const std::string& name, 
                                             ParameterHandler &param) const;   
                                             
                                             
            //@}
            ///@name Instance Delivery (functions)
            //@{
            /**
             * This member function is used to create an object of type gas diffusion layer
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > create_replica (const std::string &name)
            {
                return boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > (new FuelCellShop::Layer::DummyGDL<dim> (name));
            } 
            //@}
            ///@name Instance Delivery (attributes)
            //@{
            /**
             * 
             */           
            static DummyGDL<dim> const* PROTOTYPE;
            //@}            
            ///@name Internal variables
            //@{
            /** Method used to compute the effective properties in the pores */
            std::string method_eff_property_pores;
            /** Diffusibility */
            std::vector<double> D_D0;
            /** Anisotropy ? */
            bool anisotropy;
            /** Porosity of the GDL */
            double porosity;
            /** Oxygen diffusion coefficient */
            std::vector<double> D_O2;
            /** Water vapour diffusion coefficient */
            std::vector<double> D_wv;
            /** Solid network conductivity */
            std::vector<double> sigma_e;
            /** Solid network thermal conductivity */
            std::vector<double> k_T;
            //@}
        };
        
    }
    
}
  
#endif