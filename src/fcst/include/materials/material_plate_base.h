//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: material_plate_base.h
//    - Description: 
//    - Developers: M. Secanell and Madhur Bhaiya
//    - $Id: material_plate_base.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP_MATERIAL_PLATE_BASE__H
#define _FUELCELLSHOP_MATERIAL_PLATE_BASE__H

// Include deal.II classes
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

// Include FCST classes
#include <materials/base_material.h>


namespace FuelCellShop
{
    namespace Material
    { 
        /**
         * Base class for developing bipolar plate materials. Do not implement a material in this
         * class, this classs is only for setting up the interface.
         * This class returns:
         *  - thermal conductivity tensor and derivatives
         *  - electron (electrical) conductivity tensort and derivatives
         *  - Young's modulus tensor and derivatives
         *  - Poisson's ratio
         *  - expansion coefficient
         *  - surface properties (contact angle)
         * 
         * @todo2 Implement a liquid gas so that it can be passed to the member function that computes contact angles
         */
        class MaterialPlateBase 
        :
        public BaseMaterial
        {
        public:
            /**
             * Constructor
             */
            MaterialPlateBase(std::string name)
            : FuelCellShop::Material::BaseMaterial(name)
            {
                //implement routine  
            };
            /**
             * Destructor
             */
            ~MaterialPlateBase(){};
            /**
             * Declare parameters
             */
            void declare_parameters (ParameterHandler &param) const
            {
                FuelCellShop::Material::BaseMaterial::declare_parameters(param);
            };
            
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            void initialize (ParameterHandler &param)
            {
                FuelCellShop::Material::BaseMaterial::initialize(param);
            }
            
            /**
             * Member function to compute the electron conductivity (Isotropic properties). If the electron conductivity depends on the solution, then the solution is set in member function set_solution(std::vector< std::vector< double > * >, std::vector< std::string >)
             */
            virtual double get_electron_conductivity() const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
                
                return -1;
            };
            
            /**
             * Member function to compute the derivatives of the electron conductivity (Isotropic properties) with respect to the solution. The derivative is computed with respect to the variables that you have indicted using the member function set_derivative_flags(std::vector< std::string > &flags).
             */
            virtual void get_electron_conductivity_derivative(double &, std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            /**
             * Member function to compute the thermal conductivity (Isotropic properties). If the thermal conductivity depends on the solution, then the solution is set in member function set_solution(std::vector< std::vector< double > * >, std::vector< std::string >)
             */
            virtual double get_thermal_conductivity() const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
                
                return -1;
            };
            
            /**
             * Member function to compute the derivatives of the thermal conductivity (Isotropic properties) with respect to the solution. The derivative is computed with respect to the variables that you have indicted using the member function set_derivative_flags(std::vector< std::string > &flags).
             */
            virtual void get_thermal_conductivity_derivative(double &, std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            /**
             * Member function to compute the Young's modulus (Isotropic properties). If the thermal conductivity depends on the solution, then the solution is set in member function set_solution(std::vector< std::vector< double > * >, std::vector< std::string >)
             */
            virtual double get_youngs_modulus() const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
                
                return -1;
            };
            
            /**
             * Member function to compute the derivatives of the Young's modulus (Isotropic properties) with respect to the solution. The derivative is computed with respect to the variables that you have indicted using the member function set_derivative_flags(std::vector< std::string > &flags).
             */
            virtual void get_youngs_modulus_derivative(double &, std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            /**
             * Member function to compute the Poisson's ratio (Isotropic properties). If the thermal conductivity depends on the solution, then the solution is set in member function set_solution(std::vector< std::vector< double > * >, std::vector< std::string >)
             */
            virtual double get_poissons_ratio() const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
                
                return -1;
            };
            
            /**
             * Member function to compute the derivatives of the Poisson's ratio (Isotropic properties) with respect to the solution. The derivative is computed with respect to the variables that you have indicted using the member function set_derivative_flags(std::vector< std::string > &flags).
             */
            virtual void get_poissons_modulus_derivative(double &, std::vector<double>&) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            /**
             * Member function to compute the expansion coefficient (Isotropic properties). If the thermal conductivity depends on the solution, then the solution is set in member function set_solution(std::vector< std::vector< double > * >, std::vector< std::string >)
             */
            virtual double get_expansion_coefficient() const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
                
                return -1;
            };
            
            /**
             * Member function to compute the derivatives of the expansion coefficient (Isotropic properties) with respect to the solution. The derivative is computed with respect to the variables that you have indicted using the member function set_derivative_flags(std::vector< std::string > &flags).
             */
            virtual void get_expansion_coefficient_derivative(double &E, std::vector<double>& dE) const
            {
                const std::type_info& info = typeid(*this);
                FcstUtilities::log << "Pure function " << __FUNCTION__
                << " called in Class "
                << info.name()  << std::endl;
            };
            /**
             * 
             */
            protected:
                /** Variable storing electron conductivity */
                double electron_conductivity;
                /** Variable storing electron conductivity derivatives */
                std::vector<double> electron_conductivity_derivative;
                /** Variable storing thermal conductivity */
                double thermal_conductivity;
                /** Variable storing thermal conductivity derivatives */
                std::vector<double> thermal_conductivity_derivative;
                /** Variable storing Youngs modulus */
                double youngs_modulus;
                /** Variable storing Youngs modulus */
                std::vector<double> youngs_modulus_derivative;
                /** Variable storing Poisson ratio */
                double poissons_ratio;
                /** Variable storing Poisson modulus */
                std::vector<double> poissons_ratio_derivative;
                /** Variable storing expansion coefficient ratio */
                double expansion_coefficient;
                /** Variable storing expansion coefficient modulus */
                std::vector<double> expansion_coefficient_derivative;
        };
        
    }
}


#endif
