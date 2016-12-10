//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: platinum.h
//    - Description: Class representing Platinum material class
//    - Developers: Peter Dobson(2011), Madhur Bhaiya(2012-13) and M. Secanell(2013)
//    - $Id: platinum.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP_CATALYST_PLATINUM__H
#define _FUELCELLSHOP_CATALYST_PLATINUM__H

#include <materials/catalyst_base.h>

namespace FuelCellShop
{
    namespace Material
    {
        class Platinum
        :
        public CatalystBase
        {
        public:
            /** 
             * Name of the class. This name is used to select the layer.
             */
            static const std::string concrete_name;
            ///@name Constructors, destructor, and parameter initalization
            //@{
            /**
             * Constructor
             */
            Platinum(std::string name = "Platinum");
            /** PROTOTYE Constructor 
             * 
             * \warning For internal use only.
             */
            Platinum(const bool);
            
            /**
             * Destructor
             */
            ~Platinum();
            
            /**
             * 
             * \todo1 THIS SHOULD BE PRIVATE
             * 
             * Declare parameters.
             * 
             * The parameters for the Platinum class will usually be declared inside the layer classes.
             * 
             * The parameters that are declared in the input file are as follows:
             * 
             * @code
             * (...) <- Other sections here
             *   subsection Materials
             *     subsection name  <- This value is specified during the constructor. By default Platinum
             *       set Density [g/cm^3] = 21.5
             *       //############ Oxygen Reduction Reaction Parameters ############
             *       set Anodic transfer coefficient (ORR) = 0.5
             *       set Cathodic transfer coefficient (ORR) = 1.0
             *       set Method for kinetics parameters (ORR) = Parthasarathy # Options are: Given|Parthasarathy|Parthasarathy_hcd|Double_trap|Neyerlin
             *       set Given Open Cell Voltage (ORR) [V] = 1.229
             *       set Given Open Cell Voltage (HOR) [V] = 0.0
             *       set Reference exchange current density (ORR) [uA/cm2] = 2.707e-2
             *       set Reference oxygen concentration (ORR) = 0.725e-5
             *       set Reference proton concentration (ORR) = 1.818e-3
             *       set Oxygen reaction order (ORR) = 1.0
             *       set Proton reaction order (ORR) = 1.0
             *       //############ Hydrogen Oxidation Reaction Parameters ############
             *       set Anodic transfer coefficient (HOR) = 0.5
             *       set Cathodic transfer coefficient (HOR) = 0.5
             *       set Reference exchange current density (HOR) [uA/cm2] = 1e6
             *       set Reference hydrogen concentration (HOR) = 5.64e-5
             *       set Hydrogen reaction order (HOR) = 0.25
             *     end
             *   end
             * (...)
             * @endcode
             */
            virtual void declare_parameters (ParameterHandler &param) const;
                        
            /**
             * Member function used to read in data and initialize the necessary data from the parameter file
             * to compute the coefficients.
             */
            virtual void initialize (ParameterHandler &param);
            //@}            
            ///@name Kinetic Parameter Accessors
            //@{
            /**
             * Return anodic transfer coefficient for the reaction specified using #set_reaction_kinetics method.
             */
            virtual void alpha_anodic(double&) const;
            
            /**
             * Return cathodic transfer coefficient for the reaction specified using #set_reaction_kinetics method.
             */
            virtual void alpha_cathodic(double&) const;
            
            /**
             * Compute the exchange current density [\p A/cm^2] for the reaction specified using #set_reaction_kinetics method. It takes temperature [\p Kelvin]
             * as an input argument by reference and returns the exchange current density [\p A/cm^2].
             */
            virtual double exchange_current_density(const double&) const;
            
            /**
             * Compute the derivative of exchange current density [\p A/cm^2] w.r.t temperature [\p Kelvin] for the reaction specified using #set_reaction_kinetics method. It takes
             * temperature [\p Kelvin] as an input argument by reference and returns the derivative.
             */
            virtual double derivative_exchange_current_density(const double&) const;
                                                 
            /**
             * Compute the reference concentration for the reaction specified using #set_reaction_kinetics method.
             * The reaction might depend on more than one species, therefore this class returns a std::map of reference concentrations referenced using #VariableNames 
             * as \p Key. This map is returned by reference (second argument). It takes a vector of #VariableNames as first argument, corresponding 
             * to those solution variables for whose reference concentrations are required.
             */
            virtual void reference_concentration(const std::vector<VariableNames>&,
                                                 std::map<VariableNames, double>& ) const;
            
            /**
             * Compute the theroretical cell voltage [\p Volts] as a function of temperature, for the reaction specified using #set_reaction_kinetics
             * method. It takes temperature [\p Kelvin] as an input argument and returns the voltage.
             */
            virtual double voltage_cell_th(const double&) const;
            
            
            /**
             * Compute the derivative of theoretical cell voltage [\p Volts] w.r.t temperature [\p Kelvin], for the reaction specified using 
             * #set_reaction_kinetics method. It takes temperature [\p Kelvin] as an input argument and returns the derivative.
             */
            virtual double dvoltage_cell_th_dT(const double&) const;
                                        
            /**
             * Compute the reaction order of the electrochemical reaction with respect to each species involved in the reaction specified
             * using #set_reaction_kinetics method. The reaction might depend on more than one species, therefore this class returns a std::map of reaction orders referenced using #VariableNames 
             * as \p Key. This map is returned by reference (second argument). It takes a vector of #VariableNames as first argument, corresponding 
             * to those solution variables for whose reaction orders are required.
             */
            virtual void reaction_order(const std::vector<VariableNames>&,
                                        std::map<VariableNames, double>&) const;
            
            /**
             * Member function used to specify the reaction for which the kinetic parameters are needed. Method will fail if any other std::string is passed, 
             * other than the ones mentioned below:
             * - \b "ORR"
             * - \b "HOR"
             */
            virtual void set_reaction_kinetics(const ReactionNames name)
            {
                if (check_reaction_implementation(name))
                    name_reaction_kinetics = name;
                else
                    Assert(false, ExcNotImplemented());
            }
            //@}
                                        
        protected:
            
            /**
             * Check that whether a particular reaction is implemented in the class or not. It takes a std::string as an input argument, against
             * which the check is done. 
             * <b>For Developers:</b> As more reactions are implemented in the class, this method should be extended.
             */
            virtual bool check_reaction_implementation(const ReactionNames name) const
            {
                return ((name == ORR) || (name == HOR));
            }
            
        private:
            ///@name Instance Delivery (Member function to create replica)
            //@{
             /**
             * This member function is used to create an object of type gas diffusion layer
             * 
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Material::CatalystBase > create_replica ()
            {
                return boost::shared_ptr<FuelCellShop::Material::CatalystBase > (new FuelCellShop::Material::Platinum ());
            }
            //@}
            ///@name Instance Delivery (PROTOTYPE)
            //@{
            /**
             * Create prototype for the layer
             */            
            static Platinum const* PROTOTYPE;
            //@}
            
            ///@name Internal Data
            //@{
            /** Anodic transfer coefficient (ORR), given in the parameter file. */
            double alpha_a_ORR;
            
            /** Cathodic transfer coefficient (ORR), given in the parameter file. */
            double alpha_c_ORR;

            /** Given value of open cell voltage (OCV) for ORR, given  in the parameter file. */
            double given_OCV_ORR;
            
            /** Given value of open cell voltage (OCV) for HOR, given in the parameter file. */
            double given_OCV_HOR;
            
            /** Reference exchange current density [\f$ \mu \f$ \p A/cm^2] for ORR, given in the parameter file. */
            double i_0_ref_ORR;
            
            /** Reference oxygen concentration [\p mol/cm^3] for ORR, given in the parameter file. */
            double c_O2_ref_ORR;
            
            /** Reference proton concentration [\p mol/cm^3] for ORR, given in the parameter file. */
            double c_H_ref_ORR;
            
            /** Reaction order corresponding to oxygen gas for ORR, given in the parameter file. */
            double gamma_O2_ORR;
            
            /** Reaction order corresponding to protons for ORR, given in the parameter file. */
            double gamma_H_ORR;
            
            /** Anodic transfer coefficient (HOR), given in the parameter file. */
            double alpha_a_HOR;
            
            /** Cathodic transfer coefficient (HOR), given in the parameter file. */
            double alpha_c_HOR;
            
            /** Reference exchange current density [\f$ \mu \f$ \p A/cm^2] for HOR, given in the parameter file. */
            double i_0_ref_HOR;
            
            /** Reference hydrogen concentration [\p mol/cm^3] for HOR, given in the parameter file. */
            double c_H2_ref_HOR;
            
            /** Reaction order corresponding to hydrogen gas for HOR, given in the parameter file. */
            double gamma_H2_HOR;
            //@}
        };
        
    }//namespace Material
    
}//namespace FuelCellShop

#endif
