//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: tafel_kinetics.h
//    - Description: Header file for Tafel Kinetics model class.
//    - Developers: M. Moore, M. Bhaiya, M. Secanell and V. Zingan
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__TAFEL_KINETICS_H
#define _FUELCELLSHOP__TAFEL_KINETICS_H

// Include openFCST routines:
#include <reactions/base_kinetics.h>
#include <utils/fcst_constants.h>

// Include deal.II classes
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_values.h>

//Include STL
#include <cmath>
#include <iostream>

using namespace dealii;

namespace FuelCellShop
{
    namespace Kinetics
    {
        /**
         * This class defines a simple Tafel kinetic model. The model is based on the following
         * equation
         * \f[
         *   i = i_0 \left( \frac{c_R}{c^{ref}_R} \right)^\gamma exp \left[\frac{-\alpha_c RT}{nF} \left(\phi_s - \phi_m - E_0 \right)\right]
         * \f]
         * where
         * \f$ c_R \f$ and \f$ c^{ref}_R \f$ are the reactant concentration and reference concentration respectively, \f$ \gamma \f$ is the
         * reaction order,  \f$ \alpha_c \f$ is the cathodic transfer coefficient (usually between .5 and 1 for elementary (single electron transfer) reactions),
         * \f$ \phi_s \f$ and \f$ \phi_m \f$ are the solid and electrolyte potentials and \f$ E_0 \f$ is the open cell voltage for the
         * electrochemical reaction.
         *
         * In the parameter file there are no parameters to specify since the kinetics parameters are given by the
         * catalyst used.
         *
         * <h3>Usage Details:</h3>
         *
         *
         *
         *
         * \author M. Moore, 2011-12
         * \author M. Bhaiya, 2012-13
         * \author M. Secanell, 2006-13
         * \author V. Zingan, 2012-2014
         *
         */
        class TafelKinetics
        :
        public BaseKinetics
        {
        public:

            ///@name Constructor, destructor and initialization
            //@{
            /**
             * Constructor
             */
            TafelKinetics();

            /**
             * Constructor for PROTOTYPE
             * \warning For internal use only
             */
            TafelKinetics(const bool);

             /**
             * Destructor
             */
            ~TafelKinetics();

            /**
             * Declare parameters for a parameter file.
             * \remarks This method current doesn't implement anything.
             */
            virtual void declare_parameters(ParameterHandler& param) const;

            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             * \remarks This method currently doesn't implement anything.
             */
            virtual void initialize(ParameterHandler& param);

            //@}

            ///@name Instance Delivery (Public variables)
            //@{
            /**
             * Concrete name used for objects of this class. This name is used when
             * setting up the subsection where the data is stored in the input file.
             *
             * The data will be stored under
             * \code
             * subsection concrete_name
             *
             * end
             * \endcode
             */
            static const std::string concrete_name;
            //@}

            ///@name Computational Methods
            //@{
            /**
             * Member function that computes the current, at every quadrature point in the cell.
             */
            virtual void current_density (std::vector<double> &);

            /**
             * Function to return the derivative of the current density w.r.t solution variables. It returns a map of vectors containing
             * derivative w.r.t solution variables / design variables set using #set_derivative_flags method. Each vector can be accessed by using \p Key
             * of the map, which correpsonds to the #VariableNames (solution/design variable). This method takes input map by reference, hence the map is needed
             * to be created at application/equation level with default arguments and passed inside this method.
             */
            virtual void derivative_current (std::map< VariableNames, std::vector<double> > &);

            //@}

        protected:
            ///@name Instance Delivery
            //@{
             /**
             * This member function is used to create an object of type gas diffusion layer
             *
             * \warning This class MUST be redeclared in every child.
             */
            virtual boost::shared_ptr<FuelCellShop::Kinetics::BaseKinetics > create_replica ()
            {
                return boost::shared_ptr<FuelCellShop::Kinetics::BaseKinetics > (new FuelCellShop::Kinetics::TafelKinetics ( ));
            }
            /**
             * Create prototype for the layer
             */
            static TafelKinetics const* PROTOTYPE;
            //@}

            ///@name Kinetics parameters
            //@{
            /**
             * Method used to initialize reference concentrations, reactions orders and cathodic transfer coefficient for the reaction,
             * and number of quadrature points in the cell.
             */
            virtual void init_kin_param()
            {
                Assert( !kin_param_initialized, ExcInternalError() );
                Assert( catalyst != NULL, ExcMessage("Catalyst object not initialized in the TafelKinetics object.") );
                Assert( phi_m.is_initialized() && phi_s.is_initialized() && T.is_initialized(), ExcMessage("Either phi_m/phi_s/T is not set in the TafelKinetics object.") );
                Assert( reactants_map.size()>0, ExcMessage("At least one reactant should be set in the TafelKinetics object, using set_reactant_concentrations method.") );

                std::vector<VariableNames> names;
                for ( std::map< VariableNames, SolutionVariable >::const_iterator iter=reactants_map.begin(); iter!=reactants_map.end(); ++iter )
                    names.push_back(iter->first);

                catalyst->reference_concentration(names, ref_conc);
                catalyst->reaction_order(names, gamma);
                catalyst->alpha_cathodic(alpha_c);

                kin_param_initialized = true;
            }
            /**
             * If any of the oxygen concentrations are negative, make the ORR reaction negative
             * 
             * Note: This is used to stabilize the code at low O2 concentrations
             */
            inline double negative_concentration_correction(double concentration, double gamma, bool derivative = true) const
            {
                double correction(1.0);

                if (concentration < 0.0)
                {
                    if ( !derivative )
                        correction = -1.0;
                    else if (derivative && (gamma == 1.0 || gamma == 3.0 || gamma == 5.0) )
                        correction = 1.0; // <-- Note: 1.0 works well with gamma = 1.0 since (gamma-1) is zero, i.e. sign of derivative independent of concentration.
                    else if (derivative && ( (gamma != 1.0) || (gamma != 3.0) || (gamma != 5.0) ) )
                        correction = -1.0;
                }

                return correction;
            }
            /**
             * Cathodic transfer coefficient
             */
            double alpha_c;

            /**
             * Map of reference concentrations
             */
            std::map< VariableNames, double > ref_conc;

            /**
             * Map of reaction orders
             */
            std::map< VariableNames, double > gamma;
            //@}

        };

    } //Kinetics

} //FuelCellShop

#endif //_FUELCELLSHOP__TAFEL_KINETICS_H