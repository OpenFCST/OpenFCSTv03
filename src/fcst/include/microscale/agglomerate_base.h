// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: AgglomerateBase
// - Description: Base class for agglomerate objects, implements MicroScaleBase
// - Developers: Philip Wardlaw, Peter Dobson, Michael Moore, Marc Secanell, University of Alberta
// - $Id: agglomerate_base.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__LAYER__AGGLOMERATE_BASE_H
#define _FUELCELLSHOP__LAYER__AGGLOMERATE_BASE_H

//FCST materials and kinetics classes
#include <materials/catalyst_base.h>
#include <materials/polymer_electrolyte_material_base.h>
#include <reactions/base_kinetics.h>
#include <reactions/double_trap_kinetics.h>
#include <reactions/tafel_kinetics.h>
#include <reactions/dual_path_kinetics.h>

//Micro scale base class
#include <microscale/micro_scale_base.h>

//STD
#include <random>
#include <limits>       // std::numeric_limits


namespace FuelCellShop
{
    namespace MicroScale
    {
      
      
       /**
         *
         * @brief This class implements methods for calculating geometric parameters of a spherical agglomerate surrounded by a thin ionomer film.
         *
         * @note This class is used by Agglomerate base (which is used by most agglomerates) and ICCP. Some virtual methods are reimplemented in hybrid agglomerate.
	 * 
         * <h3> Usage details</h3>
         * Child class should set inherited members r_agg, delta_agg and epsilon_agg.
	 * Then base class can use member functions to calculate agglomerate geometry.
         *
         * @author Philip Wardlaw, Peter Dobson, Marc Secanell 
         * @date 2015
         */
      class SphericalAgglomerateGeometry
        {

        protected:

	    double pi;
	    
            /*
             * Type defs to make accessing layer info a bit easier
             */

            typedef FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>::Properties CLPropNames;
            typedef std::map<CLPropNames, double> CL_Properties;


            /*
             * Default Constructor
             */
            SphericalAgglomerateGeometry()  
	    {
	      pi = Constants::Pi();
	      
	    };

            /**
             * Member function to compute the thickness of the agglomerate thin film based on the radius and structure
             *
             * For a given ionomer volume fraction, \f$ \epsilon_N \f$; size of the agglomerate, \f$ r_{agg} \f$; and,
             * porosity inside the agglomerate, \f$ \epsilon_{agg} \f$, the necessary ionomer thin film thickness
             * is computed such that the correct volume fractions are obtained.
             *
             * A Newton loop is used in order to obtain the correct value for the thickness of the agglomerate film, \f$ \delta_{agg} \f$
             *
             * @note #compute_thickness_agg and #compute_epsilon_agg are mutually exclusive and the selection
             * between using one method or the other depends on the inpu parameter:
             * @code
             *      set Constant agglomerate parameter [Thickness | Porosity] = Thickness
             *      set Agglomerate porosity = 0.25               #specify if Porosity is selected in variable above
             *      set Thickness of the agglomerate film [nm] = 80               #specify if Thickness is selected in property above
             * @endcode
             *
             */
            virtual double compute_thickness_agg (FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>* layer);

            /**
             * Compute the volume fraction of nafion inside the agglomerate given the
             * catalyst layer ionomer volume fraction , \f$ \epsilon_N \f$; size of the agglomerate, \f$ r_{agg} \f$; and,
             * porosity inside the agglomerate, \f$ \epsilon_{agg} \f$, the necessary ionomer thin film thickness
             * is computed such that the correct volume fractions are obtained.
             *
             * A Newton loop is used in order to obtain the correct value for the agglomerate porosity, \f$ \epsilon_{agg} \f$.
             *
             * @note #compute_thickness_agg and #compute_epsilon_agg are mutually exclusive and the selection
             * between using one method or the other depends on the inpu parameter:
             * @code
             *      set Constant agglomerate parameter [Thickness | Porosity] = Thickness
             *      set Agglomerate porosity = 0.25               #specify if Porosity is selected in variable above
             *      set Thickness of the agglomerate film [nm] = 80               #specify if Thickness is selected in property above
             * @endcode
             */
            virtual double compute_epsilon_agg (FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>* layer);

            /**
             * Function used to compute the amount of electrolyte in the catalyst layer. This is
             * given by the agglomerate sturcture and the amount of solid phase. In particular,
             * \f[
             * \epsilon_N^{cat} = n \frac{4.0}{3.0}3.1416({r_agg*1e-7})^{3.0}*(\epsilon_agg-1) + (r_{agg}1e-4 + \delta_{agg}1e-4)^{3.0};
             * \f]
             * where n is computed by another function. For details on this equation and other see, M.Secanell et al., "",...
             *
             */
            virtual double compute_epsilon_N ( const double delta_agg, const double thickness_agg ) const;

            /**
             * Function to compute
             *      \f[
             *      \frac{\partial \epsilon_N^{cat}}{\partial thickness_{agg}}
             *      \f]
             * This function should only be used to calculate the thin film thickness using
             * Newton's method.  Physically, the thickness of the thin film should not have an effect
             * on the volume fractions defined in ConventionalCL
             */
            virtual  double compute_depsilonN_dthickness ( const double thickness_agg ) const;

            /**
             * Function to compute
             *      \f[
             *      \frac{\partial \epsilon_N^{cat}}{\partial \epsilon_{agg}}
             *      \f]
             * This function should only be used to calculate the thin film thickness using
             * Newton's method.  Physically, the thickness of the thin film should not have an effect
             * on the volume fractions defined in ConventionalCL
             */
            virtual double compute_depsilonN_depsilon_agg (FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>* layer) const;

            /**
             * Member function to compute the number of agglomerates. The number of agglomerates is computed differently depending
             * on the type of agglomerate and depends on the agglomerate radius, amount of ionomer per agglomerate and the
             * solid phase volume fraction in the CL.
             *
             * For a spherical agglomerate
             * \f[
             * n_{agg} = \frac{\epsilon_S}{\frac{4}{3} \pi r_{agg}^3 (1 - \epsilon_{agg})}
             * \f]
             */
            virtual double compute_n(FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>* layer) const;

            /**
             * Compute the derivative of n with respect to the porosity of the agglomerate.
             *
             * For a spherical agglomerate:
             * \f[
             * \frac{\partial n_{agg}}{\partial \epsilon_agg} = \frac{\epsilon_S}{\frac{4}{3} \pi r_{agg}^3 (1 - \epsilon_{agg})^2}
             * \f]
             */
            virtual double compute_dn_depsilon_agg(FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>* layer) const;


            /*
             * Radius of the agglomerate, typically [nm], but units may depend
             * on non abstract child
             */
            double r_agg;

            /*
             *  Thickness of the ionomer thin film, typically [nm], but units may depend
             * on non abstract child
             */
            double delta_agg;

            /*
             *  Porosity % of the agglomerate core
             */
            double epsilon_agg;


            /*
             * Number of agglomerates in CL domain
             */
            double n_agg;


        };
      
      
      
      
      
      
      
      
      
      
        /**
         *
         * @brief The base class for agglomerate objects in OpenFCST.
         *
         *
         * This base class implements some of the necessary interfaces specified by
         * MicroScaleBase, as well as providing spherical agglomerate specific
         * data and functionality.
         *
         * <h3> Input parameters </h3>
         * LIST OF INPUT PARAMETERS FOR THE CLASS.
         * @code
         * ...
         * subsection MicroScale
         *   subsection AgglomerateBase
         *     set Radius of the agglomerate [nm] = 200
         *     set Constant agglomerate parameter [Thickness | Porosity] = Porosity #Variable used to select if thickness or porosity should be constant
         *     set Thickness of the agglomerate film [nm] = 15
         *     set Agglomerate porosity = 0.25
         *   end
         * end
         * ...
         * @endcode
         *
         * <h3> Usage details</h3>
         * This class may be inherited by agglomerate object in order to gain
         * useful member data and methods. Use MicroScaleBase classes instance delivery
         * system for creating child.
         *
         * @author Philip Wardlaw
         * @date 2014
         */
        class AgglomerateBase: public MicroScaleBase, public SphericalAgglomerateGeometry
        {


        public:

            /*
             * Destructor
             */
            virtual ~AgglomerateBase(){};

            /**
             * Print out key agglomerate information (name, radius, film thickness, porosity).
             */
            virtual void print_properties();

            ///@name Initalization
            //@{
            /**
             * Function for setting the solution map(reactant concentration, phi_s, phi_m, etc.).
             * First argument provide the map of solution variables. The second argument provide
             * the name of the primary reactant. The final argument is the index of the solution
             * map that the micro scale object should solve for, i.e. if the solution map
	     * contains arrays of quadrature solutions then the index corresponds to the quadrature
	     * point we are solving for
             *
             * This function should be called immediatly before <b>compute_current</b> or
             * <b>compute_derivative_current</b>.
             *
             */
            virtual void set_solution(const std::map<VariableNames,SolutionVariable>&,const VariableNames&, const int&);
            //@}

            /**
             * Function to compute the derivative of the current density at the local operating conditions.
             * Returns current density derivatives in (A/cm^3)/(unit of the changing quantity)
             *

             * Solves for solution variables set by the last call to <b>set_solution</b>.
             *
             * <h3> Usage details</h3>
             * Call <b>has_derivatives</b> to check if it is OK to call this function.
             *
             */
            virtual std::vector<double> compute_derivative_current ()
            {
                Assert(false, ExcPureFunctionCalled());
                return std::vector<double>(sol_names.size(), 0.0);
            }


            /**
             * Returns true if the class instance can calculate
             * current density derivatives.
             *
             * <h3> Usage details</h3>
             * Call this function if unsure whether or not to call <b>compute_derivative_current</b>.
             *
             */
            virtual bool has_derivatives()
            {
                return has_derivatives_;
            }




            /*
             *  Agglomerate will deep copy it's pointers (kinetics, materials), use this for parallel execution.
             *
             * \warning This following function determines reaction name based on type of kinetics model
             * \warning which may not necessarily be correct and only will work for currently
             * \warning implemented kinetics models.
             * \warning A more robust mothod of determining the reaction name at this scope is required.
             * \warning Use of new kineitics models, or existing kinetics models for different reactions
             * \warning will result in bugs.
             */
            virtual void make_thread_safe(ParameterHandler &param, unsigned int thread_index);


        protected:

            /*
             * Method for setting up film thickness or porosity depending on parameters
             * called in initialize() or adjust_poly adjust_polydisperse_structure().
             */
            void _initialize_film_porosity();



            /*
             * Type defs to make accessing layer info a bit easier
             */

            typedef FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>::Properties CLPropNames;
            typedef std::map<CLPropNames, double> CL_Properties;


            /*
             * Constructor, sets up some constants
             */
            AgglomerateBase()
            {
                //Constants
                R = Constants::R();
                F = Constants::F();
                permittivity_0 = Constants::E0() *1e-2; // cm
            }


            /*
             * Protected virtual member function for declaring parameters, pure
             * in MicroScaleBase, implemented here in AgglomerateBase.
             * Parameter structure is hierarchical, therefore children should
             * call their parent's declare parameter function from their own declare
             * parameter function.
             */
            virtual void declare_parameters (ParameterHandler &param) const;

            /*
             * Protected virtual member function for initializing parameters, pure
             * in MicroScaleBase, implemented here in AgglomerateBase.
             * Parameter structure is hierarchical, therefore children should
             * call their parent's initialize function from their own declare
             * parameter function.
             */
            virtual void initialize (ParameterHandler &param);

            /*
             * Pure virtual function for returning film thickness in nano meters.
             *
             * Must be implemented in child, where informed unit conversion will occur.
             */
            virtual double get_film_thickness() = 0;

            /*
             * Pure virtual function for returning agglomerate radius in nano meters.
             *
             * Must be implemented in child, where informed unit conversion will occur.
             */
            virtual double get_radius() = 0;


            /**
             * Boost shared pointer to catalyst object.
             *
             * Used to determine kinetic properties of catalyst
             */
            boost::shared_ptr<FuelCellShop::Material::CatalystBase > catalyst;

            /**
             * Boost shared pointer to electrolyte object.
             *
             * Used to determine local effective properties.
             */
            boost::shared_ptr<FuelCellShop::Material::PolymerElectrolyteBase> electrolyte;

            /**
             * Boost shared pointer to kinetics object.
             *
             * Used to determine reaction rates.
             */
            boost::shared_ptr<FuelCellShop::Kinetics::BaseKinetics> kinetics;

            /**
             * Member data for storing solutions.
             */
            std::map<VariableNames, SolutionVariable> solutions;
            std::vector<VariableNames> sol_names;
            VariableNames reactant,tempReactantName;
            int sol_index;


            /**
             * Constants
             */
            double permittivity_0; //permittivity of free space (scaled in constructor)
            double F; //Faraday's Constant
            double R; //Gas constant (molar)

            /*
             * Pressure of domain in pascals
             */
            double P;

            /*
             * Active area, [cm^2/cm^3], Scaled to solid fraction by MSCL,
             * scaled to agglomerate core by non abstract agglomerate classes
             */
            double AV;

            /*
             * Concentration of primary reactant at the boundary, [mol/cm^3]
             */
            double c_R;

            /*
             * Boolean determining whether or not the agglomerate object returns derivatives.
             */
            bool has_derivatives_;

            /*
             * Concentration of protons at the boundary, [mol/cm^3]
             */
            double c_H;

            /*
             *  Electrolyte (Membrane) phase potential at the boundary, [V]
             */
            double phi_M;

            /*
             * Solid phase potential through the agglomerate, [V]
             */
            double phi_S;

            /*
             * Value of the boundary between the thin film and agglomerate domain
             *
             * Dimensionless (0.0 to 1.0)
             */
            double interface;

            /*
             * Oxygen Diffusion, [cm^2/s]
             */
            double D_R_N;

            /**
             *  Henry's Constant for primary reactant in nafion, [Pa cm^3/mol]
             */
            double H_R_N;

            /*
             * Variable used to select if thickness or porosity should be constant
             */
            std::string fixed_agg_variable;


        };
    }
}

#endif
