//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: PureLiquid.h
//    - Description: Material class for pure liquids.
//    - Developers: Madhur Bhaiya and Phil Wardlaw
//    - Id: $Id: PureLiquid.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__PURELIQUID__H
#define _FUELCELLSHOP__PURELIQUID__H

//Include STL
#include <cstdlib>
#include <vector>
#include <fstream>

#include <deal.II/base/parameter_handler.h>
#include <materials/base_material.h>
using namespace dealii;


namespace FuelCellShop

{
  namespace Material
  {

    /**
    * Virtual class used to describe different liquids and pure materials
    * for which viscority, diffusivity, etc. can be computed
    */
    class PureLiquid : public BaseMaterial
    {
    public:
          
        /** Molar mass in g/mol or kg/kmol */
        double M;
        /** dynamic viscosity at STP*/
        double Mu_0;
        
        /** Constructor */
        PureLiquid(){};
        /** Destructor */
        virtual ~PureLiquid(){};
        
        virtual char* get_name()     =0; 
        virtual char* get_formula()  =0;
    };
    //---------------------------------------------------------------------------
    
    /**
     * LiquidWater as a publicly derived class of Pureliquid.
     */
    class LiquidWater : public PureLiquid
    {
    public:

        LiquidWater();
        virtual ~LiquidWater();
        
        virtual char* get_name()     ; 
        virtual char* get_formula()  ; 
        
        /**
        * Declare all necessary parameters in order to compute the coefficients
        */
        void declare_parameters ( ParameterHandler &param ) const;

        /**
        * Member function used to read in data and initialize the necessary data
        * to compute the coefficients.
        */
        void initialize ( ParameterHandler &param );
        
        inline double get_DO2(){return oxygen_diffusion_coeff;};
        inline double get_DH() {return proton_diffusion_coeff;};
        inline double get_Relative_Permittivity(){return relative_permittivity;};
        
        /**
         * Return the density of liquid water in \f$ \left[ \frac{g}{cm^3} \right] \f$.
         */
        static double get_density() { return 1.0;}

        /**
         * Return the width of a water molecule in  \f$ \left[ M \right] \f$.
         */
        inline double get_molecular_width() {return molecular_width;}


        /**
         * Return the Henery's constant for oxygen in water \f$ \left[ Pa-m^3/mol \right] \f$.
         */
        inline double get_HO2(){return HenryO2;}

        /** 
         * Static function to return latent heat of vaporization of water \f$ \left[\frac{J}{mol}\right]\f$ as a function of input temperature \f$[K]\f$.
         * \b Ref: Simple Forumulas for Thermophysical Properties of Liquid Water for Heat Transfer Calculations (from 0C to 150C), C. O. Popiel and K. Wojtkowiak,
         * Heat Transfer Engineering, Vol.19, No.3, 1998.
         */
        static double latentVap_heat(const double&);
        
        /** 
         * Static function to return surface tension of liquid water \f$ \left[\frac{N}{m}\right]\f$ as a function of input temperature \f$[K]\f$.
         * \b Ref: International Tables of the surface tension of water N.B. Vargaftik,
         */
        
        static double surface_tension(const double&);
        /** 
         * Static function to return derivative of latent heat of vaporization of water with respect to temperature \f$ \left[\frac{J}{mol \cdot K}\right]\f$, as 
         * a function of input temperature \f$[K]\f$.
         * \b Ref: Simple Forumulas for Thermophysical Properties of Liquid Water for Heat Transfer Calculations (from 0C to 150C), C. O. Popiel and K. Wojtkowiak,
         * Heat Transfer Engineering, Vol.19, No.3, 1998.
         */
        static double deriv_latentVap_heat(const double&);
        
        /**
         * Static function to return viscosity of water \f$ \left[ \frac{g}{cm \cdot s}\right]\f$ as a function of input temperature \f$[K]\f$.
         * \b Ref: Fox R. W and McDonald A. T., Introduction to Fluid Mechanics. John Wiley & Sons, Inc., New York, fifth edition, 1998.
         */
        static double viscosity(const double&);
        /**
         * Static function return derivative of viscosity of water with respect to temperature, \f$ \left[\frac{g}{cm \cdot s \cdot K}\right] \f$, 
         * as a function of input temperature \f$[K]\f$.
         * \b Ref: Fox R. W and McDonald A. T., Introduction to Fluid Mechanics. John Wiley & Sons, Inc., New York, fifth edition, 1998.
         */
        static double deriv_viscosity(const double&);

        

    private:

        double oxygen_diffusion_coeff;
        double proton_diffusion_coeff;
        double relative_permittivity;
        double molecular_width;
        double HenryO2;


    };

  }
}

#endif
