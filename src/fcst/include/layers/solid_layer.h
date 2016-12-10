//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: solid_layer.h
//    - Description: Solid layer class
//    - Developers: M. Secanell and Jie Zhou
//    - $Id: solid_layer.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__SOLID_LAYER_H
#define _FUELCELLSHOP__SOLID_LAYER_H

// FCST classes
#include <materials/base_material.h>
#include <layers/base_layer.h>
#include <materials/carbon.h>
#include <materials/fiber_base.h>
#include <materials/carbon_fiber.h>
#include <materials/PureSolid.h>

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
         * This class is used to represent a solid layer. It takes only two arguments from the
         * input file, the material ID and the material name representing the mesh element number for the layer
         * and the type of material you would like the layer to return the properties
         * for respectively. Then, it creates a child object of PureSolid of the desired type and uses it to
         * compute all properties such as thermal conductivity.
         *  
         * <h3> Theory </h3>
         * This class is to compute the material properties at each quadrature points
         * base on the PureSolid class
         *
         * <h3> Input parameters </h3>
         * The input parameter section for this class looks like:
         * @code 
         * subsection Fuel cell data
         *   subsection solidlayer_section_name
         *     set Material id = 2
         *     set Solid material type = Graphite
         *   end
         * end
         * @endcode* 
         *
         * <h3> Usage details</h3>
         * 
         * //Create an object of SolidLayer
         * FuelCellShop::TemplateClass example; 
         * // Set necessary variables
         * temp = 358;
         * example.set_variable(temp);
         * // You can now request info from your class.
         * double temp = example.get_variable();
         * //Print to screen all properties
         * example.print_data();
         * 
         * \author J. Zhou and M. Secanell, 2014
         * 
         */
        template <int dim>
        class SolidLayer : 
        public  BaseLayer<dim>
        {
        public:
            
            ///@name Effective property calculators
            //@{
            /**
             * Compute the effective electron conductivity in the SolidLayer.
             * This class is used for both isotropic and anisotropic materials.
             */
             virtual void effective_electron_conductivity(Tensor<2,dim>& ) const;
            
            /**
             * Compute the effective electron conductivity in the SolidLayer
             * This class is used for both isotropic and anisotropic materials. It will return the conductivity
             * at each quadrature point based on the the solution provided in set_solution()
             */
             virtual void effective_electron_conductivity(std::vector<Tensor<2,dim>>& ) const;
            
            /**
             * Compute the derivative of the effective electron conductivity in the SolidLayer
             * with respect to temperature solutionvariable.
             * 
             * The map will contain the VariableNames and the derivative.
             */
            virtual void derivative_effective_electron_conductivity(std::map<VariableNames, Tensor<2,dim>>& ) const;
             
            /**
             * Compute the derivative of the effective electron conductivity in the SolidLayer
             * with respect to temperature solutionvariable.
             * 
             * The map will contain the VariableNames and the derivative.
             */
             virtual void derivative_effective_electron_conductivity(std::vector< std::map<VariableNames, Tensor<2,dim> > >& ) const;
            
            /**
             * Compute the effective thermal conductivity in the SolidLayer
             */
             virtual void effective_thermal_conductivity(Tensor<2,dim>& ) const;
            
            /**
             * Compute the effective thermal conductivity in the SolidLayer, dependent on various solution variables, eg: Temperature
             * 
             * Note it is a vector because it returns the solution at each quadrature point using set_solution to obtain the solution
             * at each location
             * 
             */
              virtual void effective_thermal_conductivity(std::vector< Tensor<2,dim> >&) const;
            
            /**
             * Compute the derivative of the effective thermal conductivity in the SolidLayer
             * with respect to temperature solutionvariable.
             */
            virtual void derivative_effective_thermal_conductivity(std::map<VariableNames, Tensor<2,dim> >& ) const;
            
            /**
             * Compute the derivative of the effective thermal conductivity in the SolidLayer
             * with respect to temperature solutionvariable.
             */
             virtual void derivative_effective_thermal_conductivity(std::vector< std::map<VariableNames, Tensor<2,dim> > >& ) const;
            
            /**
             * Member function used to set the solution [\p Kelvin] at every quadrature point
             * inside the cell. This function should particulary be used in the case of whether the solution solutionvariable is right or wrong.
             */
            inline void set_solution (std::map<VariableNames, SolutionVariable>& sol)
            {
                
                for (std::map<VariableNames, SolutionVariable>::iterator iterator = sol.begin();
                     iterator != sol.end();
                     iterator++)
                     {
                         Assert( iterator->second.get_variablename() == iterator->first, ExcMessage("Wrong solution variable passed in PureSolid::set_solution.") );
                     }
                     
                     this->Solutions = sol;
            }
            //@}
            
            ///@name Accessors and info
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
                return typeid(SolidLayer<dim>);
            }
            //@}
            
//             /**
//              * Class test
//              */
//             virtual void test_class();

        protected:
                     
            ///@name Constructors, destructor, and initialization
            //@{
            /**
             * Constructor
             */
            SolidLayer(const std::string &section_layer_name, ParameterHandler& param);
            
            /**
             * Destructor
             */
            ~SolidLayer();
            
            /**
              * 
              * Function used to select the appropriate SolidLayer type as specified in the ParameterHandler under
              * line 
              * @code 
              * subsection Fuel cell data
              *   subsection solidlayer_section_name
              *     set Solid material type = Graphite
              *   end
              * end
              * @endcode
              * current options are see FuelCellShop::Material::PureSolid
              * 
              * The class will read the appropriate section in the parameter file, i.e. the one with name \param gld_section_name ,
              * create an object of the desired type and return it.
              * 
              */
            void declare_parameters (const std::string& solidlayer_section_name,
                                     ParameterHandler &param) const;
                                     
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            void initialize (ParameterHandler &param) ;
            
            //@}
            
            ///@name Internal variables
            //@{
            
            /** Boundary ids */
            std::vector<unsigned int> boundary_ids;
            
            /** Data member that stores the solid material the layer is made of*/			
            boost::shared_ptr<FuelCellShop::Material::PureSolid>  solid;
            
            /** Name in the input file used to select the material */
            std::string concrete_solid_name;

            /** Solution Varaible at each quadrature point*/
            SolutionVariable T_vector;
            
            /** Solutions at every quadrature point inside the cell. 
             * for this mapping key is the Enumerator of VariableName_of_REV
             * the mapped value is of type SolutionVariable
             */
            std::map <VariableNames, SolutionVariable> Solutions;
            
            /*
             * Define an iterator for the mapping from VariableNames to SolutionVariable
             */
            std::map <VariableNames, SolutionVariable>::iterator solution_iterator;
            
            /*
             * Define an VariableNames of Enumerator of VariableName_of_REV
             */            
            VariableNames temperature_enum;
            
            //@}
        };
    }
}  // FuelCellShop

#endif // _FUELCELLSHOP__SOLID_LAYER_H
