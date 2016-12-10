//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: reaction_heat.h
//    - Description: Header file for class computing heat source terms due to electrochemical reaction
//    - Developers: Madhur Bhaiya
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__REACTION_HEAT_H
#define _FUELCELLSHOP__REACTION_HEAT_H

// FCST classes
#include <reactions/base_kinetics.h>
#include <materials/PureLiquid.h>
#include <utils/fcst_constants.h>

// Include deal.II classes
#include <deal.II/base/subscriptor.h>

using namespace dealii;

namespace FuelCellShop
{
    namespace Equation
    {
        /**
        * This class is used to compute non-linear reaction heat source terms in a catalyst layer.
        * There are basically two sources of heat generation for an electro-chemical reaction. These are:
        * <b>Reversible Heat</b> and <b> Irreversible Heat</b>
        * 
        * Besides this, there is an additional heat sink, accounting for
        * complete evaporation of liquid water formed inside the <b>Cathode catalyst layer</b> \p(ORR \p). This is done
        * to maintain consistency with single phase formulation.
        * 
        * The \f$ \mathbf{S_{thermal}} \f$ can be written in following three components for \p ORR:
        * 
        * \f$ \qquad \mathbf{S_{thermal}} =
        *     \begin{cases}
        *     -j \left( \phi_s - \phi_m - E_{cathode,eq}(T) \right) \quad - \quad \text{Irreversible heat} \\ ~ \\
        *     \frac{-j}{2F} T \Delta S(T) \times f_{ORR} \quad - \quad \text{Reversible Heat} \\ ~ \\
        *     \frac{-j}{2F} h_{LV}(T) \quad - \quad \text{Evaporation of liquid water formed}
        *     \end{cases}
        * \f$
        * 
        * Similarly, the \f$ \mathbf{S_{thermal}} \f$ can be written for \p HOR:
        * 
        * \f$ \qquad \mathbf{S_{thermal}} = 
        *       \begin{cases} 
        *       j \left( \phi_s - \phi_m - E_{anode,eq}(T) \right) \quad - \quad \text{Irreversible heat} \\ ~ \\
        *       \frac{-j}{2F} T \Delta S(T) \times f_{HOR} \quad - \quad \text{Reversible Heat} \\ ~ \\
        *       \end{cases}
        * \f$
        * 
        * where,
        * - \f$ j \f$ is current density [\p A/cm^3], retrieved from kinetics object.
        * - \f$ \Delta S \f$ is Entropy change for overall reaction forming liquid water product and is a function of Temperature( \f$ T \f$ ).
        * - \f$ h_{LV} \f$ is Latent heat of vaporization of water and is a function of Temperature( \f$ T \f$ ).
        * - \f$ f_{ORR} \f$ and \f$ f_{HOR} \f$ are the fractions of overall entropic heat released, corresponding to half-cell reactions of \p ORR and \p HOR respectively. They are 
        * defined such that: \f$ f_{ORR} + f_{HOR} = 1.0 \f$.
        * 
        * @remarks Following points to be noted:
        * - Overpotential (\f$ \eta \f$) is negative for \em ORR (Cathode side). Hence, a negative sign is included in the formulation.
        * - Similarly, entropy change (\f$ \Delta S \f$) for an exothermic reaction (\em ORR) is negative, hence a negative sign.
        * - There is no water formation in \em HOR (Anode side), hence no water vaporization inside the anode catalyst layer.
        * - This class computes heat source in \p W/cm^3 and its derivatives.
        * 
        * \note 
        * - This class requires atleast these three solution variables, \a viz., \f$ \phi_s \f$ (<b>Electronic electric potential</b>),
        * \f$ \phi_m \f$ (<b>Protonic electric potential</b>), and, \f$ T \f$ (<b>Temperature</b>) to be set using #set_solid_potential, #set_electrolyte_potential and 
        * #set_temperature methods respectively, in order to compute heat source terms (#heat_source) and their derivatives (#derivative_heat_source).
        * - It is necessary to call #initialize_factors method to set the corresponding heat source flags before computing heat source terms or their derivatives.
        * - Also #kinetics object is required to be set using #set_kinetics method in the initialization section of application/equation.
        * 
        * \author Madhur Bhaiya, 2012-2013
        */
        class ReactionHeat
        {
        public:
            ///@name Constructor, Destructor and Initialization
            //@{
            /**
            * Constructor
            */
            ReactionHeat();
            
            /**
            * Destructor
            */
            ~ReactionHeat();
            
            /** Function to set Kinetics for enabling calculation source terms. */
            void set_kinetics(FuelCellShop::Kinetics::BaseKinetics* kin)
            {
                kinetics = kin;
            }
            //@}
            
            ///@name Initialization methods for computation
            //@{
            /**
            * Set the electrolyte phase potential. The potential should be in Volts.
            * 
            */
            void set_electrolyte_potential(const SolutionVariable& phi)
            {
                Assert( phi.get_variablename() == protonic_electrical_potential, ExcMessage("Wrong solution variable passed in BaseKinetics::set_electrolyte_potential.") );
                phi_m = phi;
            }
            /**
            * Set the solid phase potential. The potential should be in Volts.
            */
            void set_solid_potential(const SolutionVariable& phi)
            {
                Assert( phi.get_variablename() == electronic_electrical_potential, ExcMessage("Wrong solution variable passed in BaseKinetics::set_solid_potential.") );
                phi_s = phi;
            }
            /**
            * Set temperature. The temperature should be in Kelvin.
            */
            void set_temperature(const SolutionVariable& temperature)
            {
                Assert( temperature.get_variablename() == temperature_of_REV, ExcMessage("Wrong solution variable passed in BaseKinetics::set_temperature") );
                T = temperature;
            }
            /**
            * Set the variables for which we would like to compute the derivatives.
            */
            void set_derivative_flags(const std::vector<VariableNames>& flags)
            {
                derivative_flags = flags;
            }
            
            /**
             * Method to initialize the flags and factors corresponding to various source term components in the application. This method should normally be initialization 
             * of the application/equation class. It takes following arguments:
             * - \b flag_irrev_ORR is for \em Irreversible heat generation in \p ORR, defaulted to \b TRUE.
             * - \b flag_irrev_HOR is for \em Irreversible heat generation in \p HOR, defaulted to \b TRUE.
             * - \b flag_rev_heat is for \em Reversible heat generation in the overall reaction forming liquid water product, defaulted to \b TRUE.
             * - \b rev_heat_ORR_coef is for fraction of reversible heat released in the half-cell reaction of \p ORR, defaulted to \b 1.0 (\em i.e., complete reversible (entropic) heat is released in the \p ORR).
             * - \b flag_single_phase_ORR is for complete vaporization of liquid water formed during \p ORR, defaulted to \b TRUE.
             */
            void initialize_factors(const bool& flag_irrev_ORR = true,
                                    const bool& flag_irrev_HOR = true,
                                    const bool& flag_rev_heat = true,
                                    const double& rev_heat_ORR_coef = 1.0,
                                    const bool& flag_single_phase_ORR = true)
            {
                factor_irrev_ORR    = flag_irrev_ORR             ?   -1.0                                                  : 0.0;
                factor_rev_ORR      = flag_rev_heat              ?   ((-1.0*rev_heat_ORR_coef)/(2.0*Constants::F()))       : 0.0;
                factor_vap_ORR      = flag_single_phase_ORR      ?   ((-1.0)/(2.0*Constants::F()))                         : 0.0;
                
                factor_irrev_HOR    = flag_irrev_HOR             ?   1.0                                                   : 0.0;
                factor_rev_HOR      = flag_rev_heat              ?   ((rev_heat_ORR_coef-1.0)/(2.0*Constants::F()))        : 0.0;
                
                factors_initialized = true;
            }
            //@}
            
            ///@name Source term computation methods
            //@{
            /**
            * Member function that computes the non-linear thermal source terms, as a first argument by reference. This function takes the current [\p A/cm^3] computed from 
            * the catalyst layer object, as a \p const second argument.
            */
            void heat_source (std::vector<double>& heat,
                              const std::vector<double>& current) const;
            
            /**
            * Member function that computes the derivative of the source terms with respect to the solution variables/parameters set using #set_derivative_flags method. This 
            * function takes the current [\p A/cm^3] and its derivatives computed from the catalyst layer object, as \p const third and second arguments respectively. It is also 
            * necessary to ensure that the #derivative_flags set inside this object should atleast be there inside the derivative current map. The derivative of heat source is returned 
            * by reference in the first argument as a map.
            */
            void derivative_heat_source (std::map< VariableNames, std::vector<double> >& heat_derived,
                                         const std::map< VariableNames, std::vector<double> >& current_derived,
                                         const std::vector<double>& current) const;
            //@}

        protected:
            ///@name Internal data and methods
            //@{
            /**
            * Flags for derivatives: These flags are used to request derivatives which are 
            * computed using the derivative_heat_source function.
            */
            std::vector<VariableNames> derivative_flags;
            
            /** Store pointer to the kinetics object */
            FuelCellShop::Kinetics::BaseKinetics* kinetics;
            
            /** Struct storing a pointer to the solution vector for the electrolyte potential */
            SolutionVariable phi_m;
            /** Struct stroing a pointer to the solution vector for the electronic/solid potential */
            SolutionVariable phi_s;
            /** Struct stroing a pointer to the solution vector for the temperature */
            SolutionVariable T;
            
            /** Function to compute entropy change of reaction [\p J/\p (mol-K\p )] as a function of temperature [\p K].
            * Ref: G. Lewis and M. Randall. International Critical Tables, volum 7. McGraw Hill, New York, 1930.
            */
            double entropy_rxn(const double& Temp) const
            {
                double deltaS = ((8.0*(1.0+log(Temp))) - 92.84) * 4.184; // J/mol-K
                return deltaS;  
            }
            
            /** Function to compute derivative of entropy change of reaction with respect to temperature [\p K]. */
            double deriv_entropy_rxn(const double& Temp) const
            {
                double dDs_dT = (4.184*8.0)/Temp;
                return dDs_dT;
            }
            
            /** Factor for irreversible heating in ORR. */
            double factor_irrev_ORR;
            
            /** Factor for irreversible heating in HOR. */
            double factor_irrev_HOR;
            
            /** Factor for reversible heating in ORR. */
            double factor_rev_ORR;
            
            /** Factor for reversible heating in HOR. */
            double factor_rev_HOR;
            
            /** Factor for water vaporization heat sink in ORR. */
            double factor_vap_ORR;
            
            /** Flag to check whether factors are initialized or not. */
            bool factors_initialized;
            //@}
        };
    } //Thermal 
} //FuelCellShop
#endif //_FUELCELLSHOP__REACTION_HEAT_H