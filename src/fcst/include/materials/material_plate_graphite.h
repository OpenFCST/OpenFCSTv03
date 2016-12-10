//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: 
//    - Description: 
//    - Developers: M. Secanell and Madhur Bhaiya
//    - $Id: material_plate_graphite.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP_MATERIAL_PLATE_GRAPHITE__H
#define _FUELCELLSHOP_MATERIAL_PLATE_GRAPHITE__H

// Include deal.II classes
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

// Include FCST classes
#include <materials/material_plate_base.h>


namespace FuelCellShop
{
    namespace Material
    { 
        /**
         * Class to compute the properties of graphite used in bipolar plates.
         * class, this classs is only for setting up the interface.
         * This class returns:
         *  - thermal conductivity tensor and derivatives
         *  - electron (electrical) conductivity tensort and derivatives
         *  - Young's modulus tensor and derivatives
         *  - Poisson's ratio
         *  - expansion coefficient
         *  - surface properties (contact angle)
         * 
         * @TODO Implement a liquid gas so that it can be passed to the member function that computes contact angles
         */
        class MaterialPlateGraphite
        :
        public MaterialPlateBase
        {
        public:
            /**
             * Constructor
             */
            MaterialPlateGraphite(std::string name);
            /**
             * Destructor
             */
            ~MaterialPlateGraphite();
            /**
             * Declare parameters
             */
            void declare_parameters (ParameterHandler &param) const;
            
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            void initialize (ParameterHandler &param);
            
            /**
             * Member function to compute the electron conductivity (Isotropic properties). If the electron conductivity depends on the solution, then the solution is set in member function set_solution(std::vector< std::vector< double > * >, std::vector< std::string >)
             */
            double get_electron_conductivity() const;
            
            /**
             * Member function to compute the derivatives of the electron conductivity (Isotropic properties) with respect to the solution. The derivative is computed with respect to the variables that you have indicted using the member function set_derivative_flags(std::vector< std::string > &flags).
             */
            void get_electron_conductivity_derivative(double &, std::vector<double>&) const;
            /**
             * Member function to compute the thermal conductivity (Isotropic properties). If the thermal conductivity depends on the solution, then the solution is set in member function set_solution(std::vector< std::vector< double > * >, std::vector< std::string >)
             */
            //double get_thermal_conductivity() const;
            
            /**
             * Member function to compute the derivatives of the thermal conductivity (Isotropic properties) with respect to the solution. The derivative is computed with respect to the variables that you have indicted using the member function set_derivative_flags(std::vector< std::string > &flags).
             */
            //void get_thermal_conductivity_derivative(double &, std::vector<double>&) const;
            /**
             * Member function to compute the Young's modulus (Isotropic properties). If the thermal conductivity depends on the solution, then the solution is set in member function set_solution(std::vector< std::vector< double > * >, std::vector< std::string >)
             */
            //double get_youngs_modulus() const;
            
            /**
             * Member function to compute the derivatives of the Young's modulus (Isotropic properties) with respect to the solution. The derivative is computed with respect to the variables that you have indicted using the member function set_derivative_flags(std::vector< std::string > &flags).
             */
            //void get_youngs_modulus_derivative(double &, std::vector<double>&) const;
            /**
             * Member function to compute the Poisson's ratio (Isotropic properties). If the thermal conductivity depends on the solution, then the solution is set in member function set_solution(std::vector< std::vector< double > * >, std::vector< std::string >)
             */
            //double get_poissons_ratio() const;
            
            /**
             * Member function to compute the derivatives of the Poisson's ratio (Isotropic properties) with respect to the solution. The derivative is computed with respect to the variables that you have indicted using the member function set_derivative_flags(std::vector< std::string > &flags).
             */
            //void get_poissons_modulus_derivative(double &, std::vector<double>&) const;
            /**
             * Member function to compute the expansion coefficient (Isotropic properties). If the thermal conductivity depends on the solution, then the solution is set in member function set_solution(std::vector< std::vector< double > * >, std::vector< std::string >)
             */
            //double get_expansion_coefficient() const;
            
            /**
             * Member function to compute the derivatives of the expansion coefficient (Isotropic properties) with respect to the solution. The derivative is computed with respect to the variables that you have indicted using the member function set_derivative_flags(std::vector< std::string > &flags).
             */
            //void get_expansion_coefficient_derivative(double &E, std::vector<double>& dE) const;
        };
    }
}


#endif
