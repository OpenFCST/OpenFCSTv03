//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: dual_path_kinetics.h
//    - Description:  Dual Path Kinetics model for hydrogen oxidation reaction
//    - Developers: M. Secanell, M. Moore and Madhur Bhaiya
//    - $Id: dual_path_kinetics.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__DUAL_PATH_KINETICS_H
#define _FUELCELLSHOP__DUAL_PATH_KINETICS_H

// Include openFCST routines:
#include <reactions/base_kinetics.h>

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
    namespace Kinetics
    {
        /**
         * This class will contain the implementation of the dual path kinetic kinetic model as developed by
         * Wang et al and described in the following paper:
         * 
         * J.X. Wang, T.E. Springer, R.R. Adzic, Dual-pathway kinetic equation for the
         * hydrogen oxidation reaction on pt electrodes, Journal of the Electrochemical
         * Society 153~(9) (2006) A1732--A1740.
         * 
         * In order to use this class, the user needs to provide the following in set_solution
         * 
         * 
         * \author M. Moore, 2011-12
         * \author M. Bhaiya, 2013
         * \author M. Secanell, 2006-13
         * 
         */
        class DualPathKinetics : 
        public  BaseKinetics
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
            
            ///@name Constructor, destructor and initialization
            //@{ 
            /**
             * Constructor. It also initializes the default values for kinetics parameters.
             */
            DualPathKinetics();
            
            /**
             * Constructor
             * 
             * \warning For internal use only
             */
            DualPathKinetics(const bool);
            
            /**
             * Destructor
             */
            ~DualPathKinetics();
            
            /**
             * Declare parameters for a parameter file
             */
            virtual void declare_parameters(ParameterHandler &param) const;
            
            /**
             * Member function used to read in data and initialize the necessary data
             * to compute the coefficients.
             */
            virtual void initialize(ParameterHandler &param);
            //@}
            
            ///@name Initialization
            //@{
            /**
             * Member function used to set the reaction name in the Dual path kinetics object. It will return
             * error if any string other than "HOR" is passed as an input argument.
             */
            virtual void set_reaction_kinetics(const ReactionNames & name)
            {
                if (name == HOR)
                    name_reaction_kinetics = name;

                else
                {
                    const std::type_info& info = typeid(*this);
                    FcstUtilities::log << "Only HOR reaction is to be implemented in " << __FUNCTION__ << " called in Class " << info.name()  << std::endl;
                    exit(1);
                }
            }
            //@} 
            
            ///@name Computational methods
            //@{
            /**
             * Member function that computes the current, at every quadrature point in the cell.
             * 
             * \warning This routine can only be called after potentials and concentrations have been initialized using set_* functions. See parent class.
             */
            virtual void current_density (std::vector<double>&);
            
            /**
             * Function to return the derivative of the current density w.r.t solution variables. It returns a map of vectors containing
             * derivative w.r.t solution variables / design variables set using #set_derivative_flags method. Each vector can be accessed by using \p Key
             * of the map, which correpsonds to the #VariableNames (solution/design variable). This method takes input map by reference, hence the map is needed
             * to be created at application/equation level with default arguments and passed inside this method.
             * 
             * \warning Requesting derivatives with respect to molar fraction is deprecated. Please use concentrations at electrolyte|catalyst interface.
             * 
             */
            virtual void derivative_current (std::map< VariableNames, std::vector<double> >&);
            
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
                return boost::shared_ptr<FuelCellShop::Kinetics::BaseKinetics > (new FuelCellShop::Kinetics::DualPathKinetics ( ));
            }
            /**
             * Create prototype for the layer
             */            
            static DualPathKinetics const* PROTOTYPE;
            //@}

            ///@name Kinetics parameters
            //@{
            
            /**
             * Method used to initialize reference concentration of hydrogen for the reaction, and number of quadrature points in the cell.
             */
            virtual void init_kin_param()
            {
                Assert( !kin_param_initialized, ExcInternalError() );
                Assert( catalyst != NULL, ExcMessage("Catalyst object not initialized in the DualPathKinetics object.") );
                Assert( catalyst->get_reaction_name() == HOR, ExcMessage("Catalyst object in the DualPathKinetics not set to HOR reaction name.") );
                Assert( phi_m.is_initialized() && phi_s.is_initialized() && T.is_initialized(), ExcMessage("Either phi_m/phi_s/T is not set in the DualPathKinetics object.") );
                Assert( reactants_map.find(hydrogen_concentration) != reactants_map.end(), ExcMessage("Hydrogen concentration is not set in the DualPathKinetics object.") );

                std::vector<VariableNames> names(1, hydrogen_concentration);
                std::map< VariableNames, double > cref_map;
                catalyst->reference_concentration(names, cref_map);
                ref_conc_H2 = cref_map[hydrogen_concentration];
                
                kin_param_initialized = true;
            }
            
            /** Reference concentration for hydrogen, \f$ H_2 \f$. */
            double ref_conc_H2;
            
            /** TV exchange current density [\p A/cm^2] */
            double j_0T;
            
            /** HV exchange current density [\p A/cm^2] */
            double j_0H;
            
            /** Potential range constant */
            double potential_constant;
            
            /** Reference potential */
            double ref_potential;
            
            //@}
            
            
        };
        
    } //Kinetics 
    
} //FuelCellShop

#endif //_FUELCELLSHOP__DUAL_PATH_KINETICS_H
