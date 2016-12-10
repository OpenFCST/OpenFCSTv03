//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2010-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: double_trap_kinetics.h
//    - Description: DoubleTrap Kinetics implements the double trap model proposed by
//      Wang et al. and re-developed by Moore et al.
//    - Developers: M. Moore, M. Bhaiya, M. Secanell, P. Wardlaw
//    - $Id: double_trap_kinetics.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__DOUBLETRAP_KINETICS_H
#define _FUELCELLSHOP__DOUBLETRAP_KINETICS_H

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
         * This class contains the implementation of the double trap kinetic model as described in the 
         * following paper:
         * 
         * M. Moore, A. Putz and M. Secanell. Investigation of the ORR Using the Double-Trap Intrinsic Kinetic Model, 
         * Journal of the Electrochemical Society 160(6): F670-F681, 2013. doi: 10.1149/2.123306jes
         * 
         * A more detailed deriviation of the paper is given in:
         * 
         * J.X. Wang, F.A. Uribe, T. E. Springer, J. Zhangc, and R. R. Adzica. Intrinsic kinetic equation
         * for oxygen reduction reaction in acidic media: the double tafel slope and fuel cell 
         * applications. Faraday Discuss., 140:347–362, 2008.
         *
         * This model uses the activation energies of four basic, one electron transfer, reactions that 
         * comprise the ORR and the adsorption energies of two intermediates to find the current produced 
         * in the catalyst layer. 
         * 
         * NOTE: To match Parthasarathy's data this value should be multiplied by 5.2 since this is the roughness factor of his electrode (M. Secanell, 2015)
         * 
         * 
         * \author M. Moore, 2010-13
         * \author M. Bhaiya, 2013
         * \author M. Secanell, 2013
         * \author P. Wardlaw, 2014
         *
         */
        class DoubleTrapKinetics : 
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
            * Replica Constructor
            * 
            * \warning For internal use only
            */
            DoubleTrapKinetics(const bool);
            /**
            * Constructor. It also initializes the default values for kinetics parameters from Wang's 2008 paper.
            */
            DoubleTrapKinetics();
            /**
            * Destructor
            */
            ~DoubleTrapKinetics();
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
        
            ///@name Computational methods
            //@{
            /**
             * Member function that computes the current [\p A/cm^2], at every quadrature point in the cell.
             */
            virtual void current_density (std::vector<double>& );

            /**
             * Function to return the derivative of the current density w.r.t solution variables. It returns a map of vectors containing
             * derivative w.r.t solution variables / design variables set using #set_derivative_flags method. Each vector can be accessed by using \p Key
             * of the map, which correpsonds to the #VariableNames (solution/design variable). This method takes input map by reference, hence the map is needed
             * to be created at application/equation level with default arguments and passed inside this method.
             */
            virtual void derivative_current (std::map< VariableNames, std::vector<double> >& );
            //@}
        
            ///@name Free energies, species coverage computation and accessor methods
            //@{           
            /** Returns the value of the coverage of the OH intermediate. #compute_energies should be called before using this method. */
            inline void OH_coverage(std::vector<double>& temp) const
            { 
                temp = theta_OH;
            }
            
            /** Returns the value of the coverage of the O intermediate. #compute_energies should be called before using this method. */
            inline void O_coverage(std::vector<double>& temp) const
            { 
                temp = theta_O;
            }
            
            /** Returns the rate  of the forward DA step. #compute_energies should be called before using this method.*/
            inline void G_DA_rate(std::vector<double>& temp) const
            {
                double C;
                if (reactants_map.at(oxygen_concentration)[0] < 0.0)
                    C = 1.0e-12;
                else
                    C = pow(((reactants_map.at(oxygen_concentration)[0])/c_ref_oxygen),0.5);

                std::vector<double> t;

                for(unsigned int i =0; i < intrinsic_rates[0].size(); i++)
                    t.push_back(C*intrinsic_rates[0][i]);

                temp = t;//free_energies[0];
            }
            
            /** Returns the rate  of the backward DA step. #compute_energies should be called before using this method. */
            inline void g_da_rate(std::vector<double>& temp) const
            {
                std::vector<double> t;

                for(unsigned int i =0; i < intrinsic_rates[1].size(); i++)
                    t.push_back(intrinsic_rates[1][i]*theta_O[i]);

                temp = t;//free_energies[1];
            }
            
            /** Returns the rate  of the forward RA step. #compute_energies should be called before using this method. */
            inline void G_RA_rate(std::vector<double>& temp) const
            {
                double C;
                if (reactants_map.at(oxygen_concentration)[0] < 0.0)
                    C = 1.0e-12;
                else
                    C = pow(((reactants_map.at(oxygen_concentration)[0])/c_ref_oxygen),0.5);

                double Ch;
                if(reactants_map.count(proton_concentration)){
                    Ch=(reactants_map.at(proton_concentration)[0])/c_ref_protons;
                }
                else
                    Ch = 1.0;


                std::vector<double> t;
                for(unsigned int i =0; i < intrinsic_rates[2].size(); i++)
                    t.push_back(C*Ch*intrinsic_rates[2][i]);

                temp = t;//free_energies[2];
            }
            
            /** Returns the rate  of the backward RA step. #compute_energies should be called before using this method. */
            inline void g_ra_rate(std::vector<double>& temp) const
            {
                std::vector<double> t;
                for(unsigned int i =0; i < intrinsic_rates[3].size(); i++)
                    t.push_back(theta_OH[i]*intrinsic_rates[3][i]);//free_energies[3];

                temp = t;
            }
            
            /** Returns the rate  of the forward RT step. #compute_energies should be called before using this method. */
            inline void G_RT_rate(std::vector<double>& temp) const
            {

                double Ch;
                if(reactants_map.count(proton_concentration)){
                    Ch=(reactants_map.at(proton_concentration)[0])/c_ref_protons;
                }
                else
                    Ch = 1.0;

                std::vector<double> t;
                for(unsigned int i =0; i < intrinsic_rates[4].size(); i++)
                    t.push_back(theta_O[i]*Ch*intrinsic_rates[4][i]);//free_energies[4];

                temp = t;
            }
            
            /** Returns the rate  of the backward RT step. #compute_energies should be called before using this method. */
            inline void g_rt_rate(std::vector<double>& temp) const
            {
                std::vector<double> t;
                for(unsigned int i =0; i < intrinsic_rates[5].size(); i++)
                    t.push_back(theta_OH[i]*intrinsic_rates[5][i]);

                temp = t;//free_energies[5];
            }
            
            /** Returns the rate  of the forward RD step. #compute_energies should be called before using this method. */
            inline void G_RD_rate(std::vector<double>& temp) const
            {
                double Ch;
                if(reactants_map.count(proton_concentration)){
                    Ch=(reactants_map.at(proton_concentration)[0])/c_ref_protons;
                }
                else
                    Ch = 1.0;

                std::vector<double> t;
                for(unsigned int i =0; i < intrinsic_rates[6].size(); i++)
                    t.push_back(Ch*theta_OH[i]*intrinsic_rates[6][i]);//free_energies[6];

                temp = t;
            }
            
            /** Returns the rate  of the backward RD step. #compute_energies should be called before using this method. */
            inline void g_rd_rate(std::vector<double>& temp) const
            {
                std::vector<double> t;
                for(unsigned int i =0; i < intrinsic_rates[5].size(); i++)
                    t.push_back(intrinsic_rates[7][i]*(1.0 -theta_OH[i] -theta_O[i]));

                temp = t;//free_energies[7];
            }
            
            /**
            * This function will be used to compute the activation energies of each of the reactions
            * (\f$ \Delta G^*_i \f$), accounting for the overpotential and, in the case of the dissociative and reductive 
            * adsorptions, the oxygen concentration. The function will also compute the instrinsic rates \em i.e. \f$ g_i \f$. Before calling this method, 
            * various solution setting methods should be called accordingly.
            */
            void compute_energies();
            //@}
        
            ///@name Initialization
            //@{
            /**
            * Member function used to set the reaction name in the DoubleTrap kinetics object. It will return
            * error if any string other than "ORR" is passed as an input argument.
            */
            virtual void set_reaction_kinetics(const ReactionNames & name)
            {
                if (name == ORR)
                {
                    name_reaction_kinetics = name;
                }
                else
                {
                    const std::type_info& info = typeid(*this);
                    FcstUtilities::log << "Only ORR reaction is to be implemented in " << __FUNCTION__ << " called in Class " << info.name()  << std::endl;
                    exit(1);
                }
            }
            //@}


            virtual bool has_coverage(const VariableNames& type){

                if((type == VariableNames::OH_coverage)or (type == VariableNames::O_coverage))
                    return true;

                return false;
            }

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
                return boost::shared_ptr<FuelCellShop::Kinetics::BaseKinetics > (new FuelCellShop::Kinetics::DoubleTrapKinetics () );
            }
            /**
            * Create prototype for the layer
            */            
            static DoubleTrapKinetics const* PROTOTYPE;
            //@}
        
        private:
            ///@name Kinetics parameters and internal computational methods
            //@{
            /**
             * If any of the oxygen concentrations are negative, make the ORR reaction negative
             * 
             * Note: This is used to stabilize the code at low O2 concentrations
             */
            inline double negative_concentration_correction(unsigned int index) const
            {
                double correction(1.0);
                
                if (reactants_map.at(oxygen_concentration)[index] <= 0.0)
                    correction = 0.0;
                
                return correction;
            }

            /** 
             * Compute the oxygen concentration ratio, i.e.,
             * \f[
             * \frac{c_{O_2}}{c^{ref}_{O_2}}^{n}
             * \f]
             * where \f$ n \f$ is the variable exponential.
             */
            inline double oxygen_concentration_ratio(unsigned int index) const
            {
                double ratio;
                double exponential = 0.5;

                if (reactants_map.at(oxygen_concentration)[index] < 0.0)
                    ratio = 0.0;
                else
                    ratio = pow((reactants_map.at(oxygen_concentration)[index]/c_ref_oxygen), exponential);
                
                return ratio;
            }
            
            /** 
             * If needed, compute the proton concentration ratio, i.e.,
             * \f[
             * \frac{c_{H^{+}}}{c^{ref}_{H^+}}^{n}
             * \f]
             * where \f$ n \f$ is the variable exponential.
             */
            inline double proton_concentration_ratio(unsigned int index) const
            {
                double ratio = 1.0;
                
                if(reactants_map.count(proton_concentration))
                    ratio =(reactants_map.at(proton_concentration)[index])/c_ref_protons;
                
                return ratio;
            }
            /**
             * Method used to initialize reference concentration of oxygen for the reaction, and number of quadrature points in the cell.
             */
            virtual void init_kin_param()
            {
                Assert( !kin_param_initialized, ExcInternalError() );
                Assert( catalyst != NULL, ExcMessage("Catalyst object not initialized in the DoubleTrapKinetics object.") );
                Assert( catalyst->get_reaction_name() == ORR, ExcMessage("Catalyst object in the DoubleTrapKinetics not set to ORR reaction name.") );
                Assert( phi_m.is_initialized() && phi_s.is_initialized() && T.is_initialized(), ExcMessage("Either phi_m/phi_s/T is not set in the DoubleTrapKinetics object.") );
                Assert( reactants_map.find(oxygen_concentration) != reactants_map.end(), ExcMessage("Oxygen concentration is not set in the DoubleTrapKinetics object.") );

                std::vector<VariableNames> names(1, oxygen_concentration);
                names.push_back(proton_concentration);
                std::map< VariableNames, double > cref_map;
                catalyst->reference_concentration(names, cref_map);
                
                c_ref_oxygen = cref_map[oxygen_concentration];
                c_ref_protons = cref_map[proton_concentration];
                kin_param_initialized = true;
            }
        
            /** Function to compute the derivative of the PtOH coverage w.r.t to the solution **/
            void d_OH_coverage_du (std::map< VariableNames, std::vector<double> >& ) const;
            
            /** Function to compute the derivative of the PtO coverage w.r.t to the solution **/
            void d_O_coverage_du (std::map< VariableNames, std::vector<double> >& ) const;
            
            /** Function that will compute the derivative of the intrinsic rates wrt to
            * the requested solution variable.
            */
            void drate_du(std::vector<std::vector<double> >&, const VariableNames& ) const;
                       
            /** 
            * Values of the free activation energies for each reaction (\f$ \Delta G^*_i \f$).
            * Because the free activation energies change with oxygen and overpotential
            * we will need one for each quadrature point.
            * The numbering system is as follows:
            * 0 = G_DA, 1 = G_da
            * 2 = G_RA, 3 = G_ra
            * 4 = G_RT, 5 = G_rt
            * 6 = G_RD, 7 = G_rd
            * where the lower case subscript indicates a backward reaction.
            */
            std::vector<std::vector<double> > free_energies;
            
            /** 
            * Values of the equilibrium free activation energies for each reaction (\f$ \Delta G^{*0}_i
            * \f$). 
            * Note that the backward energies are given a lower case "subscript". These are the values 
            * at zero overpotential and are defined in the constructor, with the values taken from the 
            * following article:
            * Jia X. Wang, Junliang Zhang, and Radoslav R. Adzic. Double-trap kinetic equation for the 
            * oxygen reduction reaction on pt(111) in acidic media. 
            * The Journal of Physical Chemistry A, 111(49):12702–12710, 2007. PMID:18052309.
            */
            double G_DA_0, G_RA_0, G_RT_0, G_RD_0;
            
            /** 
            * Values of the free adsorption energies for each intermediate (\f$ \Delta G^{*0}_{OH} \f$
            * and \f$ \Delta G^{*0}_{O} \f$). 
            * These are the values at zero overpotential and are defined in the constructor, with the 
            * values taken from the following article:
            * Jia X. Wang, Junliang Zhang, and Radoslav R. Adzic. Double-trap kinetic equation for the 
            * oxygen reduction reaction on pt(111) in acidic media. 
            * The Journal of Physical Chemistry A, 111(49):12702–12710, 2007. PMID:18052309..
            */
            double G_O, G_OH;
            
            /** 
            * Values of the intrinsic rates for each reaction (\f$ g_i \f$). Note that
            * the backward energies are given a lower case "subscript".
            */
            std::vector<std::vector<double> > intrinsic_rates;
            
            /** Stores the value of the coverage of the OH and O intermediates.*/ 
            std::vector<double> theta_OH, theta_O;
            
            /** 
            * prefactor is the reference prefactor that sets the scale for the current produced. 
            * Given as 1000A/cm^2 in DoubleTrap paper.
            */
            double prefactor;
            
            /** Stores the reference concentration of oxygen */
            double c_ref_oxygen;
            
            /** Stores the reference concentration of protons */
            double c_ref_protons;

            /** Stores the transfer coefficient */
            double alpha;
            //@}
        };        
    } //Kinetics 
} //FuelCellShop

#endif //_FUELCELLSHOP__DOUBLETRAP_KINETICS_H
