// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: agglomerate_ionomer_analytical.h
// - Description: Used to solve a system of equations representing a spherical ionomer-filled agglomerate.
// - Developers: Peter Dobson <pdobson@ualberta.ca>
//               Marc Secanell Gallart, University of Alberta
// - $Id: agglomerate_ionomer_analytical.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef FUEL_CELL__IONOMER_AGGLOMERATE_ANALYTICAL__H
#define FUEL_CELL__IONOMER_AGGLOMERATE_ANALYTICAL__H


//------------------------------
// FUEL CELL DECLARATIONS
//-----------------------------

#include <microscale/agglomerate_base.h>
#include <reactions/tafel_kinetics.h>
#include <application_core/fcst_variables.h>
#include <reactions/dual_path_kinetics.h>
#include <materials/catalyst_base.h>


namespace FuelCellShop
{
    namespace MicroScale
    {
        /**
         * \brief Class that gives the analytical solution to an ionomer-filled agglomerate problem in 1D
         * 
         * This class implements the analytical ionomer-filled agglomerate first proposed in:
         * - M. Moore, P. Wardlaw, P. Dobson, R. Spiteri, J. B and M. Secanell, "", Journal of the Electrochemical Society
         * 
         * <h3> Theory </h3>
         * The volumetric current density is
         * 
         * \f[
         *    \nabla \cdot \vec{i} = \frac{1}{1-\epsilon_V} 4F \bar{V}_{agg} \frac{P_{O_2}}{H_{O_2,N}}\left[\frac{1}{E_r k_c} + \frac{\delta_{agg}r_{agg}^2}{3\left(r_{agg}+\delta_{agg}\right)D_{O_2,N}}\right]^{-1}
         * \f]
         * 
         * The effectiveness factor is obtained from the analytical solution of the oxygen transport equation on the agglomerate domain,
         * \f[
         * E_r = \frac{1}{\phi_L}\left( \frac{1}{\tanh(3\phi_L)} - \frac{1}{3\phi_L} \right)
         * \f]
         * where \f$ \phi_L \f$ is Thieles modulus, which characterizes the reaction-transport process for a given geometry.  
         * For a sphere, the characteristic length is \f$ \frac{r_{agg}}{3} \f$, so Thiele's modulus becomes 
         * \f[
         * \phi_L = \frac{r_{agg}}{3} \sqrt{\frac{k_c}{D^{eff}_{O_2}}}
         * \f]
         * 
         * The term \f$ \left(1-\varepsilon_V\right)\bar{V}_{agg} \f$ is an active area scaling factor.  
         * Typically, the active area for an electrode is given as the area per volume of catalyst layer. 
         * Since in the agglomerate model the platinum is only found in the core of the agglomerate, the active area 
         * has to be adjusted accordingly. Dividing by \f$ \left(1-\varepsilon_V\right) \f$, where \f$ \epsilon_V \f$ is 
         * the porosity of the electrode, gives the active area of Pt (\f$ cm^2_{Pt} \f$) per volume of agglomerate (\f$ cm^3_{agg} \f$). 
         * Then, dividing by \f$ \bar{V}_{agg} \f$ gives the active area of Pt (\f$ cm^2_{Pt} \f$) per volume of agglomerate core 
         * (\f$ cm^3_{agg, core} \f$). The variable $\bar{V}_{agg}$ is defined as
         * \f[
         * \bar{V}_{agg} = \frac{V_{agg}}{V_{tot}} = \frac{\frac{4\pi r_{agg}^3}{3}}{\frac{4\pi \left(r_{agg}+\delta_{agg}\right)^3}{3}} = \dfrac{r_{agg}^3}{\left(r_{agg}+\delta_{agg}\right)^3}
         * \f]
         * is a scaling factor determined as the ratio of the volume of the agglomerate core to the volume of the entire agglomerate.
         * 
         * \author P. Dobson, P. Wardlaw and M. Secanell 2009-13
         * 
         */
         
        class IonomerAgglomerateAnalytical : public AgglomerateBase
        {
        public:


            static const std::string concrete_name;
            
            /*
             * Set the composition and structure of the agglomerate
             */

            
            /**
             * Main function of the class used to compute the current over the whole agglomerate
             * at the local operating conditions
             */
            virtual SolutionMap compute_current ( );
            
            /**
             * Function to compute the derivative of the current density at the local operating conditions;
             */
            virtual std::vector<double> compute_derivative_current ();
            
            /**
             * Return name of class instance, i.e. concrete name.
             */
            virtual std::string get_name(){
                    return concrete_name;
                }

            /**
             * Returns extra contribution to volume of layer
             */
            virtual double aux_volume_fraction(){
                return 0;
            }



        protected:
            /*
             * Virtual function for returning film thickness in nano meters.
             *
             */
            virtual double get_film_thickness(){
                return delta_agg*1e7;
            }

            /*
             * Virtual function for returning agglomerate radius in nano meters.
             *
             */
            virtual double get_radius(){
                return r_agg*1e7;
            }

            /*
             * Set the composition and structure of the agglomerate
             */
            virtual void set_structure ();

            //name Instance Delivery (Prototype)
            static IonomerAgglomerateAnalytical const* PROTOTYPE;
            
            /**
             * This member function is used to create an object of MicroScaleBase
             */
            virtual boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> create_replica ()
            {
                return boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> (new FuelCellShop::MicroScale::IonomerAgglomerateAnalytical ());
            }
            
            /** Constructors */
            IonomerAgglomerateAnalytical ();
            IonomerAgglomerateAnalytical(std::string concrete_name);

            /*
             * Protected virtual member function for declaring parameters, pure
             * in MicroScaleBase, implemented here in IonomerAgglomerateAnalytical.
             * Calls parent AgglomerateBase
             */
            virtual void declare_parameters (ParameterHandler &param) const
            {
                AgglomerateBase::declare_parameters(param);
                param.enter_subsection(concrete_name);{
                    //No parameters at the moment
                }
                param.leave_subsection();
            }

            /*
             * Protected virtual member function for initializing parameters, pure
             * in MicroScaleBase, implemented here in IonomerAgglomerateAnalytical.
             * Calls parent AgglomerateBase
             */
            virtual void initialize (ParameterHandler &param)
            {
                AgglomerateBase::initialize(param);
                param.enter_subsection(concrete_name);{
                    //No parameters at the moment
                }
                param.leave_subsection();
            }

        private:
            
            /** Function to compute the effectiveness of the agglomerate core */
            double compute_Er (const double k_c, const double D);
            
            /** Function to compute the derivative of the effectiveness of the agglomerate core */
            double compute_dEr (const double k_c,const double dk_c, const double D);
            
            /* Private member function to check kinetics are appropriate for analytical formulation.
             * Throws exception in event of inappropriate kinetic conditions.
             */
            void check_kinetics();

            /*bool to monitor if we have checked kinetic conditions */
            bool checked_kinetics;


        };// class IonomerAgglomerate2
        
    } // namespace Layer
    
} // namespace FuelCellShop

#endif
