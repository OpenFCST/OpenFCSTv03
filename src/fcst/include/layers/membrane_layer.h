//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: membrane_layer.h
//    - Description: (base) class for membrane layers
//    - Developers: Marc Secanell (2013) and Madhur Bhaiya (2012-13)
//    - Id: $Id: membrane_layer.h 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__MEMBRANE_LAYER_H
#define _FUELCELLSHOP__MEMBRANE_LAYER_H

// Include FCST classes
#include <layers/base_layer.h>
#include <materials/polymer_electrolyte_material_base.h>
#include <materials/nafion.h>


namespace FuelCellShop
{
    namespace Layer
    {
        /**
         * Virtual class used to provide the interface for all MembraneLayer children. 
         * 
         * No object of type MembraneLayer should ever be created, instead this layer
         * is used to initialize pointers of type MembraneLayer. The class has a database of
         * children such that it will declare all necessary parameters for all children in the
         * input file, read the input file, create the appripriate children and return a pointer
         * to MembraneLayer with the children selected.
         * 
         * All public functions are virtual but the static functions used to declare parameters and
         * to initialize a pointer of MembraneLayer, i.e. declare_all_MembraneLayer_parameters,
         * set_all_MembraneLayer_parameters and create_MembraneLayer.
         * 
         * <h3>Usage Details:</h3>
         * 
         * In order to create a membrane layer within an application, the following steps need to be taken.
         * 
         * First, in the application .h file, create a pointer to a MembraneLayer object, i.e.
         * @code
         * boost::shared_ptr<FuelCellShop::Layer::MembraneLayer<dim> > PEM;
         * @endcode
         * 
         * This pointer object will be available anywhere inside the application. Because we do not want to
         * worry about deleting the pointer afterwards, we use a Boost pointer which has its own memory management
         * algorithms. See the <a href="http://www.boost.org/doc/libs/1_54_0/libs/smart_ptr/smart_ptr.htm">Boost website</a> for more information
         * 
         * Once the pointer is available, we need to do three things in the application
         * - Call declare_MembraneLayer_parameters in the application declare_parameters() member function. This 
         * member function is used to define all parameters that can be read from the input file
         * 
         * - Call create_MembraneLayer and initialize. The former member function will fill the pointer created above with the appropriate
         * membrane layer you have selected in the input file. Then, initialize reads from the input file all the data from the file in order
         * to setup the object.
         * 
         * The object is ready for use now.
         * 
         * Here is a code example from ??.cc:
         * @code
         * //--------- IN DECLARE_PARAMETERS ------------------------------------------------------
         * template <int dim>
         * void 
         * NAME::??<dim>::declare_parameters(ParameterHandler& param)
         * {
         *   (...)
         *   // Declare section on the input file where all info will be stored. In this case Fuel Cell Data > Membrane layer
         *   FuelCellShop::Layer::MembraneLayer<dim>::declare_MembraneLayer_parameters("Membrane layer", param);
         *   (...)
         * }
         * 
         * //--------- IN INITIALIZE ------------------------------------------------------
         * template <int dim>
         * void
         * NAME::App???<dim>::_initialize(ParameterHandler& param)
         * {   
         *  (...) 
         *   // Initialize layer classes:
         *    electrolyte.initialize(param);
         *   PEM = FuelCellShop::Layer::MembraneLayer<dim>::create_MembraneLayer("Membrane layer", param, &electrolyte);
         *   PEM->initialize(param);
         *   (...)
         * }
         * @endcode
         * 
         * 
         * \author M. Secanell, 2013
         * 
         */
        template <int dim>
        class MembraneLayer :
        public BaseLayer<dim>
        {
        public:
            ///@name Instance Delivery
            //@{         
            /**
             * Function used to declare all the data necessary in the parameter files for
             * all MembraneLayer children.
             * 
             * This member function should be used instead of declare_parameters() when we want
             * to use a MembraneLayer pointer that selects the type of PEM to run at runtime.
             * 
             * \param pem_section_name Name of the section that will encapuslate all the information about the PEM
             * \param param ParameterHandler object used to store all information about the simulation. Used
             * to read the parameter file.
             * 
             * The parameter file would look as follows:
             * 
             * @code
             * subsection Fuel cell data
             * (...)
             * subsection Membrane layer           #This is the name provided in gld_section_name
             * 
             *   set Membrane layer type = NafionMembrane #[ NafionMembrane ]
             *   set Material id = 2
             *
             *     subsection NafionMembrane         # This is the subsection for the first children
             *       set 
             *     end
             * end
             * (...)
             * end
             * @endcode 
             */
            static void declare_MembraneLayer_parameters (std::string pem_section_name, ParameterHandler &param)
            {
                for (typename FuelCellShop::Layer::MembraneLayer<dim>::_mapFactory::iterator iterator = FuelCellShop::Layer::MembraneLayer<dim>::get_mapFactory()->begin(); 
                     iterator != FuelCellShop::Layer::MembraneLayer<dim>::get_mapFactory()->end(); 
                iterator++)
                     {
                         iterator->second->declare_parameters(pem_section_name, param);
                     }   
            }
             
            /**
              * 
              * Function used to select the appropriate MembraneLayer type as specified in the ParameterHandler under
              * line 
              * @code 
              * set Membrane type = NafionMembrane 
              * @endcode
              * current options are [  ]
              * 
              * The class will read the appropriate section in the parameter file, i.e. the one with name \param pem_section_name ,
              * create an object of the desired type and return it.
              * 
              */
            static boost::shared_ptr<FuelCellShop::Layer::MembraneLayer<dim> > create_MembraneLayer (std::string pem_section_name, 
                                                                                                     ParameterHandler &param)
            {
                boost::shared_ptr<FuelCellShop::Layer::MembraneLayer<dim> > pointer;
                
                std::string concrete_name;
                param.enter_subsection("Fuel cell data");
                {
                    param.enter_subsection(pem_section_name);
                    {
                        concrete_name = param.get("Membrane layer type");
                        FcstUtilities::log << "name: "<<concrete_name.c_str()<<std::endl;
                    }
                    param.leave_subsection();
                }
                param.leave_subsection();
                
                typename FuelCellShop::Layer::MembraneLayer<dim>::_mapFactory::iterator iterator = FuelCellShop::Layer::MembraneLayer<dim>::get_mapFactory()->find(concrete_name);
                
                if (iterator != FuelCellShop::Layer::MembraneLayer<dim>::get_mapFactory()->end())
                {
                    if (iterator->second)
                    {
                        pointer = iterator->second->create_replica(pem_section_name);
                    }
                    else 
                    {
                        FcstUtilities::log<<"Pointer not initialized"<<std::endl;
                        abort();
                    }
                }
                                else
                 {
                     FcstUtilities::log<<"Concrete name in FuelCellShop::Layer::MembraneLayer<dim>::create_MembraneLayer does not exist"<<std::endl;
                     abort();
                 }
                 
                pointer->initialize(param);
                
                return pointer;
            }
            //@}

            ///@name Accessors and Initialization
            //@{
            /**
            * This member function returns a type_info object with the name of the 
            * base layer type the inherited class belongs to, i.e.
            * - GasDiffusionLayer
            * - MicroPorousLayer
            * - CatalystLayer
            * - MembraneLayer
            * - SolidLayer
            * - Channel
            * 
            * Note that this is necessary if we want to find out not the name of the actual class which can be obtain using
            * @code const std::type_info& name = typeid(*this) @endcode
            * but the name of the parent class.
            * 
            * @note Do not re-implement this class in children classes
            */
            const std::type_info& get_base_type() const
            {
                return typeid(MembraneLayer<dim>);
            }
                                    
            /**
             * Method used to set the variables for which you would like to compute the derivatives in the membrane layer. It takes vector of
             * #VariableNames as an input argument. It also sets the derivative flags in the electrolyte object of the catalyst layer.
             */
            virtual void set_derivative_flags(const std::vector<VariableNames>& flags)
            {
                this->derivative_flags = flags;
                this->electrolyte->set_derivative_flags(flags);
            }
            
            /**
             * Set those solution variables which are constant in the particular application. If the effective properties in the layer depend on other variables that are usually
             * part of the solution vector but are assumed to be constant in this simulation, the const solution value should be passed to the class using this member function. This 
             * method should be called in the initialization section of the application. This function takes value to be set as the first argument and the #VariableNames as 
             * second argument. For instance, it's required to store constant temperature value for an isothermal application, in that case this method can be used. \em e.g., in 
             * order to set temperature as \p 353.0 [\p Kelvin] in the layer, you can use the following code:
             * @code
             * // In the initialization section of the application.
             * layer.set_constant_solution(353.0, VariableNames::temperature_of_REV);
             * @endcode
             * 
             * If #temperature_of_REV is passed using this method, it also sets the temperature [\p Kelvin] values in the #electrolyte object.
             */
            virtual void set_constant_solution(const double& value, const VariableNames& name)
            {
                FuelCellShop::Layer::BaseLayer<dim>::set_constant_solution(value, name);
                
                if (name == temperature_of_REV)
                {
                    this->electrolyte->set_T(value);
                }
            }
            
            /** Method to provide access to the pointer to electrolyte material in the membrane layer. */
            FuelCellShop::Material::PolymerElectrolyteBase* get_electrolyte() const
            {
                return this->electrolyte.get();
            }
            
            //@}
            ///@name Effective property calculators
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
            * as a function of solution variables.
            */
            virtual void effective_oxygen_diffusivity(std::vector<double>&) const;
            /**
             * Compute the derivative of effective oxygen diffusivity of the membrane with respect to
             * the flags set by set_derivative_flags, at every quadrature point.
             */
            virtual void derivative_effective_oxygen_diffusivity(std::map< VariableNames, std::vector<double> >&) const;
            
            /**
             * Compute the constant effective thermal conductivity of membrane layer.
             */
            virtual void effective_thermal_conductivity(double&) const;
            /**
             * Compute the effective thermal conductivity of membrane layer for isotropic case at all quadrature points
             */
            virtual void effective_thermal_conductivity(std::vector<double>&) const;
            /**
             * Compute the effective thermal conductivity of membrane layer for anisotropic case at all quadrature points
             */
            virtual void effective_thermal_conductivity(std::vector< Tensor<2,dim> >&) const;
            /**
             * Compute the derivative of the effective thermal conductivity in the membrane layer for isotropic case
             * at all quadrature points
             */
            virtual void derivative_effective_thermal_conductivity(std::vector< std::vector<double> >&) const;
            /**
             * Compute the derivative of the effective thermal conductivity in the membrane layer for anisotropic case
             * at all quadrature points
             */
            virtual void derivative_effective_thermal_conductivity(std::vector< std::vector< Tensor<2,dim> > >&) const;
            
            /**
             * Compute the effective thermo-osmotic diffusivity of lambda (sorbed water),
             * at all quadrature points in the membrane layer.
             */
            virtual void effective_thermoosmotic_diffusivity(std::vector<double>& ) const;
            
            /**
             * Compute the derivative of the effective thermo-osmotic diffusivity of lambda (sorbed water) in the membrane layer
             * with respect to either the solution or design parameters. The parameters with respect to
             * which the derivatives are computed are setup in FuelCellShop::Layer::set_derivative_flags()
             */
            virtual void derivative_effective_thermoosmotic_diffusivity(std::map< VariableNames, std::vector<double> >& ) const;
            
            
            //@}
            
        protected:
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
            MembraneLayer();
            /**
            * Constructor
            * 
            * \note Eventually, I would like to make this private.
            * 
            * \deprecated Use create_MembraneLayer
            */
            MembraneLayer(std::string name);
            
            /**
            * Destructor
            */
            ~MembraneLayer();
            
            /**
            * Declare all necessary parameters in order to compute the coefficients
            * 
            * \deprecated Use declare_all_MembraneLayer_parameters
            */
            void declare_parameters(ParameterHandler &param) const
            {
                this->declare_parameters(this->name, param);
            }  
            
            /**
             * Declare parameters for a parameter file. 
             * 
             */            
            virtual void declare_parameters (const std::string& name, 
                                             ParameterHandler &param) const;           
                                   
            /**
            * Member function used to read in data and initialize the necessary data
            * to compute the coefficients.
            */
            void initialize (ParameterHandler &param);
            //@}
            ///@name Instance Delivery
            //@{         
             /** 
             * This object is used to store all objects of type GasDiffusionLayer. 
             * This information in then used in layer_generator.h in order to create the
             * correct object depending on the specified concrete type of layer selected
             * such as DummyGDL or SGL24BA.
             */
            typedef std::map< std::string, MembraneLayer<dim>* > _mapFactory;      
            
            /**
             * 
             */
            static _mapFactory * get_mapFactory()
            {
                static _mapFactory mapFactory;
                return &mapFactory;
            }  
            /**
             * This member function is used to create an object of type gas diffusion layer
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Layer::MembraneLayer<dim> > create_replica (std::string &name)
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name() << std::endl;
            }            
            //@}
                                                                          
            ///@name Internal variables
            //@{
            /**
             * Type of membrane
             */
             std::string electrolyte_type;
             
            /**
             * Pointer to the electrolyte object used in the constructor, and it is used 
             * to calculate the effective properties of the membrane layer.
             */
             boost::shared_ptr< FuelCellShop::Material::PolymerElectrolyteBase > electrolyte;    
            
            //@}
            
        };
        
    }
}

#endif
