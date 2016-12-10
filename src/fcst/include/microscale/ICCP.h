//---------------------------------------------------------------------------
// C++ Interface: agglomerate_ionomer_1D.h
//
// Description: Class that solves solid carbon particle, with Pt surface loading,
//              surrounded by ionomer thin film
//
// Authors: Philip Wardlaw,
//          University of Alberta.
//
// Copyright: See COPYING file that comes with this distribution
//
//---------------------------------------------------------------------------

#include<microscale/micro_scale_base.h>

//FCST materials and kinetics classes
#include <materials/catalyst_base.h>
#include <materials/polymer_electrolyte_material_base.h>
#include <reactions/base_kinetics.h>
#include <reactions/double_trap_kinetics.h>
#include <reactions/tafel_kinetics.h>
#include <reactions/dual_path_kinetics.h>
#include <microscale/agglomerate_base.h> //For SphericalAgglomerateGeometry


#ifndef ICCO_H_
#define ICCO_H_


namespace FuelCellShop
{
    namespace MicroScale
    {

    /**
     * \brief Class that solves solid carbon particle, with Pt surface loading, surrounded by ionomer thin film
     *
     * This class is a very basic representation of the CL micro structure, intended to analyize the CL/MEA's
     * sensitivity to micro scale models for various operating conditions and material params.
     *
     * <h3> Input parameters </h3>
     * LIST OF INPUT PARAMETERS FOR THE CLASS.
     * @code
     * ...
     * subsection MicroScale
     *   subsection ICCP
     *      set Radius [nm] = 100.0 #Inner core radius
     *      set Non Equilibrium BC Rate constant = 0.13 #Non equilibrium Reaction rate coefficient
     *      set Use non equilibrium BC = false #Use non-equilibrium BC as described by Suzukhi et al.
     *   end
     * end
     * @endcode
     * 
     * \warning Treatment of pressure using get_properties from layer is incorrect for non-isobaric cases!!!!!!
     * 
     * @note Film thickness is calculated automatically from layer properties.
     * 
     */


    class ICCP: public MicroScaleBase, public SphericalAgglomerateGeometry {

        public:

            static const std::string concrete_name;


            /**
             * Function for setting the solution map(reactant concentration, phi_s, phi_m, etc.).
             * First argument provide the map of SolutionVariable. The second argument provide
             * the name of the primary reactant. The final argument is the index of the solution
             * map that the micro scale object should solve for.
             *
             * This function should be called immediatly before <b>compute_current</b> or
             * <b>compute_derivative_current</b>.
             *
             */
            virtual void set_solution(const std::map<VariableNames,SolutionVariable>& sols,const VariableNames& react, const int& index);

            /**
             * Function used to compute the current density produced by the micro structure.
             * Effectiveness is returned by reference, and defined by child class.
             *
             * Solves for solution variables set by the last call to <b>set_solution</b>.
             */
            virtual SolutionMap compute_current ( );

            /**
             * Returns true if the class instance can calculate
             * current density derivatives. In this case it will return false.
             *
             */
            virtual bool has_derivatives(){
                //numerical for the moment
                return false;
            }

            /**
             * Return name of class instance, i.e. concrete name.
             */
            virtual std::string get_name(){
                return concrete_name;
            }

            /**
             * MicroScale object may have extra contribution to volume of layer, e.g. water.
             * In this case it is zero.
             */
            virtual double aux_volume_fraction(){
                return 0;
            }


            /**
             * Print out key micro-structural dimensions, defined by child.
             */
            virtual void print_properties();

            virtual void make_thread_safe(ParameterHandler &param, unsigned int thread_index);


        protected:
	  
	    /*
             * Method for setting up film thickness or porosity depending on parameters
             * called in initialize() or adjust_poly adjust_polydisperse_structure().
             */
            void _initialize_film_porosity();


            /** Convenient typdef for getting properties */
            typedef FuelCellShop::Layer::MultiScaleCL<deal_II_dimension> CLPropNames;

            /**Default Constructor*/
            ICCP():
            F(Constants::F()),
            pi(Constants::Pi()){
                non_eq_BC = false;
            }

            /**Factory map registration constructor */
            ICCP(std::string name);


            virtual void set_structure ();
            /*
             * Protected pure virtual member function for declaring parameters, must be implemented
             * by all children. Parameter structure is hierarchical, therefore children
             * should call their parent's declare parameter function from their own declare
             * parameter function.
             */
            virtual void declare_parameters (ParameterHandler &param) const;

            /*
             * Protected pure virtual member function for initializing parameters, must be implemented
             * by all children. Parameter structure is hierarchical, therefore children
             * should call their parent's initialize function from their own initialize function.
             */
            virtual void initialize (ParameterHandler &param);

            /**
             * This member function is used to create an object of type MicroScaleBase
             */
            virtual boost::shared_ptr<MicroScaleBase> create_replica (){
                return boost::shared_ptr<MicroScaleBase>( new ICCP());
            }
            static ICCP const* PROTOTYPE;

        private:

            /**
             * Compute the residual functions describing oxygen concentration at the particle surface.
             *
             * Equilibrium BC:  0 = c_{O_2, \; f|c} - c_{O_2, \; g|f} - \frac{\delta_{agg} }{r_{agg}(r_{agg} +\delta_{agg}) } \frac{ j(c_{O_2, \; f|c}) \cdot A_s}{4 n F \pi D_{O_2,N}}  
             * Nonequilibrium BC (substitute the following term) :  c_{O_2, \; g|f} =   c_{O_2, \; g|f}^{eq} - \frac{j(c_{O_2, \; f|c}) \cdot A_s}{4 n F \cdot \pi (r_{agg} +\delta_{agg})^2 k_o}
             *
             * See equations 2.111 and 2.120 of Philip Wardlaw's Thesis
             */
            double residual(const double & c_inner, const double & c_outer);


            //Stored solutions
            std::vector<SolutionVariable> reactants;
            SolutionVariable proton_pot, electron_pot;

            //Transport factors
            double ActiveArea;      //Scaled to agglomerate surface
            bool non_eq_BC;         //If true then use. See residual()
            double non_eq_BC_coeff;
            double P;               //Pressure
            double H_R_N;            //Henry's number for reactant in Nafion
            double D_R_N;            // diffusion of reactant in bulk nafion

	    //Reaction factors
	    double molarNumerator;
	    VariableNames reactant;
	    
	    
            //Pointers to assests obtained from the CL
            boost::shared_ptr<FuelCellShop::Kinetics::BaseKinetics> kinetics;
            boost::shared_ptr<FuelCellShop::Material::PolymerElectrolyteBase> electrolyte;
            boost::shared_ptr<FuelCellShop::Material::CatalystBase> catalyst;


  
            //Some constants
            const double F, pi;
        };
    }

}



#endif /* ICCO_H_ */
