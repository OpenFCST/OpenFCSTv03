//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: carbon.h
//    - Description: Class characterizing a carbon black support
//    - Developers: M. Secanell and Madhur Bhaiya
//    - $Id: carbon.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__CARBON_H
#define _FUELCELLSHOP__CARBON_H

//Include STL
#include<cmath>
#include<iostream>

// FCST
#include <materials/catalyst_support_base.h>


namespace FuelCellShop
{
    namespace Material
    {
        /**
         * Class characterizing a carbon black support. This material is commonly used in catalyst layers
         * and micro-porous layers. 
         * 
         * There are many types of carbon black used in fuel cells such as
         * VulcanXC, Ketjen black, graphatized Ketjen black and HISPEC. Each one of these
         * carbons have different densities and bulk propeties. This class provides an example
         * for further development of other catalyst support objects. Note that in this case
         * the propeties of the material can be modified in the input file. This is
         * done in order to allow for investigation of different properties.
         * 
         * Parameters that can be specified in the input file for this class are as follows:
         * 
         * @code
         * subsection Fuel cell data
         *   subsection Materials
         *     subsection Carbon Black
         *       set Density [g/cm^3] = 2.0
         *       set Electrical conductivity [S/cm] = 88.84  # Electrical conductivity of a group of particles packed to have 0% porosity
         *       set Thermal conductivity [W/(cm-K)] = 1.0
         *     end
         *   end
         * end
         * 
         * @endcode
         * 
         * <h3>Usage Details:</h3>         
         * 
         * As with most routines you need to first declare_parameters then initialize. After this the class is ready for use.
         * 
         * See for example FuelCell::OperatingConditions class for more details.
         *
         * \author M. Secanell, 2011-13
         * \author M. Bhaiya, 2013
         * 
         */
        class CarbonBlack :
        public CatalystSupportBase
        {
        public:
            
            /** 
             * Name of the class. This name is used to select the layer.
             */
            static const std::string concrete_name;
            
            ///@name PROTOTYPE Constructor and destructor
            //@{
            /** PROTOTYE Constructor 
             * 
             * \warning For internal use only.
             */
            CarbonBlack(const bool);
            
                        /**  
             * Destructor  
             */
            ~CarbonBlack();
            
            
            //@}
            
            ///@name Information and accessors
            //@{
            /** Obtain the electrical conductivity [\p S/cm] of carbon black support material. */
            virtual double get_electrical_conductivity() const
            {
                return this->electrical_conductivity;
            }
            
            /** Obtain the thermal conductivity [\p W/\p (cm-K \p )] of carbon black support material. */
            virtual double get_thermal_conductivity() const
            {
                return this->thermal_conductivity;
            };

            /** Obtain the density [\p gm/cm^3] of carbon black support material. */
            virtual double get_density() const
            {
                return this->density;
            };
            //@}
            
        private:
            ///@name Constructors, destructor, and parameter initalization
            //@{
            /** Constructor 
             * The constructor initialize parameters using the default values. This is
             * so that if I do not want to call declare_parameters and initialize, I can
             * still use the routine with the hard coded values.
             * 
             */
            CarbonBlack();
                        
            /** Declare parameters.
             * 
             */
            virtual void declare_parameters(ParameterHandler &param) const;
            
            /**
             * Member function used to read in data and initialize the necessary data from the parameter file to compute the coefficients.
             */
            virtual void initialize (ParameterHandler &param);
            //@}
            
            ///@name Instance Delivery
            //@{
             /**
             * This member function is used to create an object of type carbon black material.
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Material::CatalystSupportBase > create_replica ()
            {
                return boost::shared_ptr<FuelCellShop::Material::CatalystSupportBase > (new FuelCellShop::Material::CarbonBlack ());
            }
            /**
             * Create prototype for the layer
             */            
            static CarbonBlack const* PROTOTYPE;
            //@}
            
        };
    }
}

#endif
