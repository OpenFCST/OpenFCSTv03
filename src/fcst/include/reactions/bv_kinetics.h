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
//    - Description: Header file for Butler-Volmer Kinetics model class.
//    - Developers: M. Secanell, M. Moore and M. Bhaiya
//    - $Id: bv_kinetics.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__BV_KINETICS_H
#define _FUELCELLSHOP__BV_KINETICS_H

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
         * This class implements a Butler-Volmer kinetic model. The model is based on the following
         * equation 
         * \f[
         *   i = i_0 \left( \frac{c_R}{c^{ref}_R} \right)^\gamma \left[ exp \left(\frac{\alpha_a RT}{nF} (\phi_s - \phi_m - E_0)\right) - exp \left(\frac{-\alpha_c RT}{nF} (\phi_s - \phi_m - E_0)\right) \right]
         * \f]
         * where
         * \f$ c_R \f$ and \f$ c^{ref}_R \f$ are the reactant concentration and reference concentration respectively, \f$ \gamma \f$ is the
         * reaction order,  \f$ \alpha_a \f$ is the anodic transfer coefficient, \f$ \alpha_c \f$ is the cathodic transfer coefficient, 
         * \f$ \phi_s \f$ and \f$ \phi_m \f$ are the solid and electrolyte potentials and \f$ E_0 \f$ is the open cell voltage for the
         * electrochemical reaction.
         * 
         * In the parameter file there are no parameters to specify since the kinetics parameters are given by the
         * catalyst used.
         * 
         *
         * <h3>Usage Details:</h3>
         * 
         * 
         * 
         * 
         * \author M. Moore, 2011-12
         * \author M. Bhaiya, 2013
         * \author M. Secanell, 2006-13
         */
        class ButlerVolmerKinetics 
        : 
        public BaseKinetics
        {
        public:
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
            
            /**
             * Constructor
             */
            ButlerVolmerKinetics();
            
            /**
             * Constructor for PROTOTYPE
             * \warning For internal use only
             */
            ButlerVolmerKinetics(const bool);
            
            /**
             * Destructor
             */
            ~ButlerVolmerKinetics();
            
            /**
             * Declare parameters for a parameter file.
             * \remarks This method current doesn't implement anything.
             */
            virtual void declare_parameters(ParameterHandler&) const{};
            
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             * \remarks This method currently doesn't implement anything.
             */
            virtual void initialize(ParameterHandler&){};

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
                return boost::shared_ptr<FuelCellShop::Kinetics::BaseKinetics > (new FuelCellShop::Kinetics::ButlerVolmerKinetics ( ));
            }
            /**
             * Create prototype for the layer
             */            
            static ButlerVolmerKinetics const* PROTOTYPE;
            //@}
            
            ///@name Kinetics parameters
            //@{
            /**
             * Method used to initialize reference concentrations, reactions orders and cathodic/anodic transfer coefficient for the reaction, 
             * and number of quadrature points in the cell.
             */
            virtual void init_kin_param()
            {
                Assert( !kin_param_initialized, ExcInternalError() );
                Assert( catalyst != NULL, ExcMessage("Catalyst object not initialized in the ButlerVolmerKinetics object.") );
                Assert( phi_m.is_initialized() && phi_s.is_initialized() && T.is_initialized(), ExcMessage("Either phi_m/phi_s/T is not set in the ButlerVolmerKinetics object.") );
                Assert( reactants_map.size()>0, ExcMessage("Atleast one reactant should be set in the ButlerVolmerKinetics object, using set_reactant_concentrations method.") );
                
                std::vector<VariableNames> names;
                for ( std::map< VariableNames, SolutionVariable >::const_iterator iter=reactants_map.begin(); iter!=reactants_map.end(); ++iter )
                    names.push_back(iter->first);
                
                catalyst->reference_concentration(names, ref_conc);
                catalyst->reaction_order(names, gamma);
                catalyst->alpha_cathodic(alpha_c);
                catalyst->alpha_anodic(alpha_a);
                
                kin_param_initialized = true;
            }
            
            /**
             * Cathodic transfer coefficient
             */
            double alpha_c;
            
            /**
             * Anodic transfer coefficient
             */
            double alpha_a;
            
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

#endif //_FUELCELLSHOP__BV_KINETICS_H