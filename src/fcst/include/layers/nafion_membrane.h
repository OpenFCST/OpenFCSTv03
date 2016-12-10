//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: nafion_membrane.h
//    - Description: Class representing Nafion membrane layer class - returning effective transport properties
//    - Developers: Madhur Bhaiya (2012-13)
//    - Id: $Id: nafion_membrane.h 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__NAFION_MEMBRANE_H
#define _FUELCELLSHOP__NAFION_MEMBRANE_H

// Include FCST classes
#include <layers/membrane_layer.h>

#include <boost/shared_ptr.hpp>

namespace FuelCellShop
{
    namespace Layer
    {
        /**
         * This class implements the necessary information for a Nafion membrane.
         * 
         * It contains information regarding the material_id used in the mesh to assign cells
         * with that belong to the membrane and provides an interface to compute
         * bulk membrane properties based on the FuelCellShop::Material::Nafion class
         *
         * <h3>Usage Details:</h3>
         * 
         * For usage information please see MembraneLayer documentation.
         * 
         * \author M. Secanell, 2009-13
         * \author M. Bhaiya, 2013 
         */
        template <int dim>
        class NafionMembrane : public MembraneLayer<dim>
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
             *    set Membrane type =  NafionMembrane # <-here I select the type of object of type GasDiffusionLayer
             *    subsection NafionMembrane # <- this is the concrete_name for this class
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
            NafionMembrane();
                        
            /**
            * Destructor
            */
            ~NafionMembrane();          
            
            //@}
            
            ///@name Accessors and info
            //@{
            /**
            * Compute the constant effective proton conductivity of the membrane.
            */
            virtual void effective_proton_conductivity(double&) const;
            /**
            * Compute the effective proton conductivity of the membrane at every quadrature point,
            * as a function of solution variables.
            */
            virtual void effective_proton_conductivity(std::vector<double>&) const;
            /**
            * Compute the derivative of effective proton conductivity of the membrane with respect to
            * the flags set by set_derivative_flags, at every quadrature point.
            */
            virtual void derivative_effective_proton_conductivity(std::map< VariableNames, std::vector<double> >&) const;
            
            /**
            * Compute the constant effective water diffusivity of the membrane.
            */
            virtual void effective_water_diffusivity(double&) const;
            /**
            * Compute the effective water diffusivity of the membrane at every quadrature point,
            * as a function of solution variables.
            */
            virtual void effective_water_diffusivity(std::vector<double>&) const;
            /**
            * Compute the derivative of effective water diffusivity of the membrane with respect to
            * the flags set by set_derivative_flags, at every quadrature point.
            */
            virtual void derivative_effective_water_diffusivity(std::map< VariableNames, std::vector<double> >&) const;

            /**
            * Compute the constant effective oxygen diffusivity of the membrane.
            */
            virtual void effective_oxygen_diffusivity(double&) const;
            /**
            * Compute the effective oxygen diffusivity of the membrane at every quadrature point,
            * as a function of Temperature.
            * \note This method should only be used when Temperature is one of the solution variables.
            */
            virtual void effective_oxygen_diffusivity(std::vector<double>&) const;
            /**
            * Compute the derivative of effective oxygen diffusivity of the membrane with respect to
            * the flags set by set_derivative_flags, at every quadrature point.
            * \note This method should only be used when Temperature is one of the solution variables.
            */
            virtual void derivative_effective_oxygen_diffusivity(std::map< VariableNames, std::vector<double> >&) const;
            
            /**
            * Compute the constant effective thermal conductivity of nafion membrane layer.
            */
            virtual void effective_thermal_conductivity(double& ) const;
            
            /**
            * Compute the effective thermal conductivity of membrane layer for isotropic case at all quadrature points.
            * Currently, we have only "Given" method for thermal conductivity, i.e. constant thermal conductivity.
            * \note Vector needs to be initialized with number of quadrature points, before passing to this function,
            * otherwise error is generated.
            */
            virtual void effective_thermal_conductivity(std::vector<double>& ) const;
            
            /**
             * Compute the effective thermal conductivity of membrane layer for anisotropic case at all quadrature points. 
             * Currently, we have only "Given" method for thermal conductivity, i.e. constant thermal conductivity.
             * \note Vector needs to be initialized with number of quadrature points, before passing to this function,
             * otherwise error is generated.
             */
            virtual void effective_thermal_conductivity(std::vector< Tensor<2,dim> >&) const;
            
            /**
            * Compute the derivative of the effective thermal conductivity in the membrane layer for isotropic case
            * at all quadrature points. Currently, we have only "Given" method for thermal conductivity, implying all
            * the derivatives are zero.
            * \note Vector needs to be initialized with number of quadrature points at inner level, and with number of
            * variables at outer level, before passing to this function, otherwise error is generated. 
            */
            virtual void derivative_effective_thermal_conductivity(std::vector< std::vector<double> >& ) const;
            
            /**
            * Compute the effective thermo-osmotic diffusivity of lambda (sorbed water),
            * at all quadrature points in the Nafion membrane layer.
            * \note Input vector needs to be initialized with number of quadrature points, before passing to this function.
            */
            virtual void effective_thermoosmotic_diffusivity(std::vector<double>& ) const;
            
            /**
            * Compute the derivative of the effective thermo-osmotic diffusivity of lambda (sorbed water) in the Nafion membrane layer
            * with respect to either the solution or design parameters. The parameters with respect to
            * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
            * \note Input vector needs to be initialized with number of quadrature points at inner level, and with
            * number of variable (derivative flags) at outer level, before passing to this function.
            */
            virtual void derivative_effective_thermoosmotic_diffusivity(std::map< VariableNames, std::vector<double> >& ) const;

            //@}
        
        protected:
            ///@name Construtor and initalization (Private)
            //@{
            /**
             * Constructor
             * 
             */
            NafionMembrane(std::string name);
            
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
                                 
            /**
            * Member function used to read in data and initialize the necessary data
            * to compute the coefficients.
            */
            void initialize (ParameterHandler &param);
            
            //@}
            ///@name Instance Delivery Functions
            //@{         
            /**
             * This member function is used to create an object of type gas diffusion layer
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Layer::MembraneLayer<dim> > create_replica (std::string &name)
            {
                return boost::shared_ptr<FuelCellShop::Layer::MembraneLayer<dim> > (new FuelCellShop::Layer::NafionMembrane<dim> (name));
            }            
            //@}                        
            ///@name Instance Delivery Data
            //@{         
            /**
             * Create prototype for the layer
             */            
            static NafionMembrane<dim> const* PROTOTYPE;
            //@}            
            
            ///@name Internal variables
            //@{           
            /**
            * String for storing method of computing thermal conductivity in the layer.
            */
            std::string method_thermal_conductivity;
            
            /**
            * Variable for storing thermal conductivity for isotropic case.
            */
            double thermal_conductivity;
            //@}            
        };
    }
}

#endif
