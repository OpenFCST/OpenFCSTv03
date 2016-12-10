// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: base_layer.h
// - Description: This is a base class for all available FCST layers
// - Developers: Marc Secanell Gallart,    University of Alberta
//               Madhur Bhaiya, University of Alberta
//               Valentin N. Zingan, U of A
// - $Id: base_layer.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__BASE__LAYER_H
#define _FUELCELLSHOP__BASE__LAYER_H

// Include deal.II classes
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_values.h>

//Include STL
#include <cmath>
#include <iostream>

// Include OpenFCST routines:
#include <application_core/fcst_variables.h>
#include <application_core/system_management.h>
#include <utils/fcst_utilities.h>

using namespace dealii;

namespace FuelCellShop
{
    

    namespace Layer
    {
        /**
         * Virtual class used to characterize a generic layer interface. Note that his class does not
         * contain enough information to characterize any useful layer yet. Please see children.
         * 
         * The base layer is used to store an identification number for the layer and 
         * the geometry of the layer.
         * 
         * \author M. Secanell
         * \author M. Bhaiya
         * \author V. Zingan
         */
        template <int dim>
        class BaseLayer : public Subscriptor
        {
        public:
            ///@name Initialization
            //@{                      
            
            /**
             * Set the variables for which you would like to compute the derivatives. It takes a vector of #VariableNames as an input argument.
             */
            virtual void set_derivative_flags(const std::vector<VariableNames>& flags)
            {
                this->derivative_flags = flags;
            }
            
            /**
             * Member function used by some applications such as dummyGDL in order to know which value
             * to return. For other classes this class is not used.
             */
            void set_position(const std::vector<Point<dim> > &p)
            {
                point = p;
            }
            

            /**
             * Function for setting local material id, for unit testing purposes.
             */
            virtual void  set_local_material_id(const unsigned int& id){

            	AssertThrow(std::find(std::begin(material_ids), std::end(material_ids), id) != std::end(material_ids), ExcMessage("Material id not found in layer."));
            	local_material_id_ = id;
            }

            /**
             * Function for unsetting local material id, so that it isn't incorrectly used later
             * Once the key is "unset" to some invalid value, an error will be thrown if the key is requested again without being set.
             * The programmer should be concious of which material id they are using the layer for, therfore the must always set it before
             * using the layer. Manually unsetting the key is a good practise to ensure that you do not solve for the wrong layer id. 	   
             */
            inline void unset_local_material_id(){
				//Ideally a material id will be 0 to 255, so we assume the max integer value will be an invalid key
                local_material_id_ = std::numeric_limits<int>::max(); 
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
             */
            virtual void set_constant_solution(const double& value, const VariableNames& name)
            {
                constant_solutions[name] = value;
            }
            
            /**
             * If the effective properties in the layer depend on the solution, the solution for a given cell should be passed to
             * the class using this member function. It is used to set SolutionVariable structure inside the layer. This structure stores the solution variable values
             * at all quadrature points in the cell. For sample usage details, please see documentation of FuelCellShop::SolutionVariable structure.
             * 
             * Note, this function in the base layer sets the interface. It has to be reimplemented in respective child layer classes for respective uses.
             * 
             * \note Use only for solution variables.
             */
            virtual void set_solution(const std::vector< SolutionVariable >&)
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            }
            
            //@}
            
            ///@name Accessors and info
            //@{
                
            /**
             * Check if a given cell belongs to the catalyst layer and assign \param material_id to current_local_material_id_ of layers
             * that are graded.
             */
            bool belongs_to_material(const unsigned int material_id);
            
            
            /**
             * Return the name of the layer
             */
            inline const std::string& name_layer() const
            {
              return name;
            }

            /**
             * This member function return the name of the type of layer, i.e. 
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
             * @note Do not implement this class anywhere other than the following "base" classes:
             * - GasDiffusionLayer
             * - MicroPorousLayer
             * - CatalystLayer
             * - MembraneLayer
             * - SolidLayer
             * - Channel
             */
            
            virtual const std::type_info& get_base_type() const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            }    
            /**
             * This member function is a virtual class that can be used to output to screen information from the layer.
             *
             * This function should be re-implemented in each layer with data that is relevant in each case.
             */
            virtual void print_layer_properties() const;
            
            /**
             * This virtual class should be used for any derived class to be able to test the functionality of the class.
             */
            virtual bool test_layer()
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
                
                return false;
            }
            
            /**
             * Return the local material id of the layer
             */
            std::vector<unsigned int> get_material_ids()
            {
             return material_ids;
            }

            /**
             * Return the local material id of the layer, performs a check
             */
            inline unsigned int local_material_id() const
            {
                AssertThrow(std::find(std::begin(material_ids), std::end(material_ids), local_material_id_) != std::end(material_ids),
                        ExcMessage("The material id for layer " + name + " has not been set."));

                return local_material_id_;
            }
        
        protected:
            
            ///@name Constructors, destructor, and initalization
            //@{
            /**
             * Constructor
             * 
             * \warning For internal use only
             */
            BaseLayer()
            {};
    
            /**
             * Constructor
             */
            BaseLayer(const std::string& name);
            
            /**
             * Destructor
             */
            virtual ~BaseLayer();
            
            /**
             * Declare parameters for a parameter file
             */
            virtual void declare_parameters (const std::string &object_name, ParameterHandler &param) const
            {
                param.enter_subsection("Fuel cell data");
                {
                    param.enter_subsection(object_name);
                    {
                        param.declare_entry("Material id",
                                            "1",
                                            Patterns::List(Patterns::Integer(0, 255)),
                                            "Id number used to identify this layer");
                    }
                    param.leave_subsection();
                }
                param.leave_subsection();
            }
            
            /**
             * Declare parameters for a parameter file.
             * 
             * \deprecated
             */          
            virtual void declare_parameters (ParameterHandler &param) const
            {
                this->declare_parameters(this->name,param);
            }
                         
            /**
             * Member function used to change the values in the parameter file for a given list
             * of parameters.
             * - name_param should ideally contain the string as it would appear in the input file, i.e. "Thickness of the agglomerate radius"
             * - value_param contains the value that the variable should be set at
             * - param is the ParameterHandler that contains all the information that has been read from file
             * 
             * Note that this is a static function, therefore it requires as the first argument the string with the name
             * of the section you would like to create for the object of this class.
             */
            virtual void set_parameters(const std::string &object_name, 
                                       const std::vector<std::string>& name_dvar,
                                       const std::vector<double>& value_dvar,
                                       ParameterHandler &param)
            {
                param.enter_subsection("Fuel cell data");
                {
                    param.enter_subsection(object_name);
                    {
                        // ADD VARIABLES IF NEEDED
                        //param.set(name_var, value_param);
                    }
                    param.leave_subsection();
                }
                param.leave_subsection();
            };
            /**
             * Set parameters in parameter file
             * 
             * \deprecated Cannot be used if using layer_generator.h
             * 
             */
            virtual void set_parameters(const std::vector<std::string>& name_dvar,
                                        const std::vector<double>& value_dvar,
                                        ParameterHandler &param)
            {
                param.enter_subsection("Fuel cell data");
                {
                    param.enter_subsection(this->name);
                    {
                        // ADD VARIABLES IF NEEDED
                        //param.set(name_var, value_param);
                    }
                    param.leave_subsection();
                }
                param.leave_subsection();
            };
            
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            virtual void initialize (ParameterHandler &param);
            
            //@}
            
            ///@name Basic layer information
            //@{
            /**
             * Name of the layer. This value is used as a header in a subsection for the parameter file
             * where all the information concerning this layer is held
             */
            const std::string name;
            /** 
             * List of material IDs that belong to the layer. The parameter local_material_id_ is an integer in this list
             */
            std::vector<unsigned int> material_ids; 
            /** 
             * Coordinates of the point where we would like to compute the effective properties
             */
            std::vector<Point<dim> > point;
            /** 
             * Flags for derivatives: These flags are used to request derivatives
             */
            std::vector<VariableNames> derivative_flags;
            /**
             * Map storing values of solution variables constant in a particular application. 
             */
            std::map< VariableNames, double > constant_solutions;

            

        private:
            /**
             * Local material ID to select the appropriate properties. This value is used in graded layers which within
             * the given layer have several properties
             *
             * This variable is private, it should be accessed by acessor local_material_id_() which performs debug checks
             */
            unsigned int local_material_id_;

            //@}
        };
        
    } // Layer
    
}  // FuelCellShop
    
    #endif // _FUELCELLSHOP__GENERIC__LAYER_H
