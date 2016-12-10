//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2012, 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: Dummy_MPL.h
//    - Description: Header file for a user defined MPL
//    - Developers: Marc Secanell
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__DUMMY_MPL_H
#define _FUELCELLSHOP__DUMMY_MPL_H

// FCST classes
#include <layers/micro_porous_layer.h>
#include <utils/fcst_units.h>

using namespace dealii;

namespace FuelCellShop
{
    namespace Layer
    {
        /**
         * This class defines an MPL with all effective properties given via input file.
         * 
         * 
         * \author M. Secanell, 2013
         */
        template <int dim>
        class DummyMPL : 
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
             *    set Microporous layer type = DummyMPL # <-here I select the type of object of type GasDiffusionLayer
             *    subsection DummyMPL # <- this is the concrete_name for this class
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
            DummyMPL();
                       
            /** Destructor */
            ~DummyMPL()
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
            
            void initialize (ParameterHandler &param);
            //@}
            
            ///@name Accessors and info
            //@{
            /**
             * Returns the effective diffusivity of the media.
             */
            void effective_gas_diffusivity(Table<2, Tensor<2,dim> >&prop_eff) const;
            
            void effective_gas_diffusivity(std::vector< Tensor<2,dim> >& prop_eff_vec) const;
            
            void derivative_effective_gas_diffusivity(std::map< VariableNames, std::vector< Tensor<2,dim> > >&dprop_eff)  const;
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
            /**
             * Constructor
             * 
             * \note Eventually, I would like to make this private.
             * 
             * \deprecated Use create_MicroPorousLayer
             */
            DummyMPL(std::string name);
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
            void declare_parameters (std::string name, ParameterHandler &param) const;
            
            ///@name Instance Delivery (function)
            //@{
            /**
             * This member function is used to create an object of type micro porous layer
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Layer::MicroPorousLayer<dim> > create_replica (const std::string &name)
            {
                return boost::shared_ptr<FuelCellShop::Layer::MicroPorousLayer<dim> > (new FuelCellShop::Layer::DummyMPL<dim> (name));
            }   
            //@}
            
            ///@name Instance Delivery (Data member)
            /**
             * Create prototype for the layer
             */          
            static DummyMPL<dim> const* PROTOTYPE;
            //@}
            
            ///@name Internal variables
            //@{

            //@}
            
        };
        
    }
}

#endif
