//---------------------------------------------------------------------------
// C++ Interface: agglomerate_water_1D.h
//
// Description: Used to solve a system of equations representing a
// 			spherical water-filled agglomerate.
//
// Author: Peter Dobson <pdobson@ualberta.ca>, (C) 2011
// 	       Philip Wardlaw,
//         University of Alberta
//
// Copyright: See COPYING file that comes with this distribution
//
//---------------------------------------------------------------------------

#ifndef FUEL_CELL__WATER_AGGLOMERATE_1D__H
#define FUEL_CELL__WATER_AGGLOMERATE_1D__H

//------------------------------
// STD LIBRARY DECLARATIONS
//------------------------------


//------------------------------
// DEAL.II DECLARATIONS
//------------------------------


//------------------------------
// FUEL CELL DECLARATIONS
//------------------------------

#include <microscale/numerical_agglomerate_base.h>
#include <contribs/DAE_wrapper.h>

//Liquid water
#include <materials/PureLiquid.h>


namespace FuelCellShop
{
namespace MicroScale
{

/**
 * \brief Class that solves a water-filled agglomerate problem in 1D
 *
 * This class uses the DAE solver interface to solve the problem using COLDAE
 * (included in the contrib folder) to solve a boundary value problem.
 * <h3> Input parameters </h3>
 * LIST OF INPUT PARAMETERS FOR THE CLASS.
 * @code
 * ...
 * subsection MicroScale
 *   subsection WaterAgglomerateNumerical
 *      set Relative core charge factor = 0.0 #Factor of ionic charge concentration, fraction of film sulphonic concentration.
 *   end
 *   subsection Materials
 *     subsection Water
 *        set Oxygen diffusion coefficient [cm^2/s]= 9.19e-5
 *        set Proton diffusion coefficient [cm^2/s]= 9.2e-5
 *        set Relative permittivity=60
 *     end
 *   end
 * end
 * @endcode
 */
class WaterAgglomerate : public NumericalAgglomerateBase, public FuelCell::ApplicationCore::DAEWrapper
{
public:



    static const std::string concrete_name;



    /**
    * Main function of the class used to compute the current over the whole agglomerate
    * at the local operating conditions
    */
    virtual SolutionMap compute_current ( );


    /**
     * Return name of class instance, i.e. concrete name.
     */
    virtual std::string get_name(){
            return concrete_name;
    }

    /**
     * MicroScale object may have extra contribution to volume of layer, e.g. water.
     *
     * The water filled model contributes to the porosity of the CL due to the presence of liquid
     * water in the agglomerate primary pores.
     */
    virtual double aux_volume_fraction(){


        typedef FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>::Properties Props;
        std::map<Props, double> p = this->layer->get_properties();

        return p[Props::solid_fraction]*epsilon_agg/(1-epsilon_agg);

    }




protected:

    /**
     * Function which implements a continuation based on the tolerance
     */
    int cont_tolerance (double start_tol, double end_tol);

    /**
     * Function which implements a continuation based on lambda -
     * modifying the source terms to scale back the equations
     */
    int cont_cd ();

    /** Setup the variables in the problem required by the DAE Solver */
    void setup_DAE_solver ();


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

    /** Set the composition and structure of the agglomerate */
    void set_structure ( );

    /** Constructors */
    WaterAgglomerate ( int verbose = 1 );
    WaterAgglomerate (std::string concrete_name);


    /*
     * Protected virtual member function for declaring parameters
     *
     *
     *
     */
    virtual void declare_parameters (ParameterHandler &param) const
    {
        NumericalAgglomerateBase::declare_parameters(param);
        water.declare_parameters(param);
        param.enter_subsection(concrete_name);{

            param.declare_entry("Relative core charge factor", "0.0", Patterns::Double(0.0,1.0),
                    "Factor dictating the amount of negative charged ions existing in the "
                    "water filled core. Fraction of film sulphonic concentration.");
        }
        param.leave_subsection();

    }

    /*
     * Protected virtual member function for initializing parameters
     */
    virtual void initialize (ParameterHandler &param)
    {
        NumericalAgglomerateBase::initialize(param);
        water.initialize(param);
        param.enter_subsection(concrete_name);{

            core_charge_factor = param.get_double("Relative core charge factor");
        }
        param.leave_subsection();


    }


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
    double compute_thickness_agg ();

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
    virtual double compute_epsilon_agg ();


    ///@name Instance Delivery
    //@{
    static WaterAgglomerate const* PROTOTYPE;

    /**
     * This member function is used to create an object of type MicroScaleBase
     */
    virtual boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> create_replica ()
    {
        return boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> (new FuelCellShop::MicroScale::WaterAgglomerate());
    }
    //@}


    /**
    * Define the DAE function.  In this case,it is simply a system of ODES.
    * COLDAE allows for the system to be defined as a system of 2 mixed-order ODEs.
    * See COLDAE.f for additional information about how to define fsub.
    * In particular, see COLDAE.f for the meaning of z and y.
    * Note that because a BVP is solved, y is not used.
    */
    void fsub ( double &, double [], double [], double [] );

    static void fsub_wrapper ( double &, double [], double [], double [] );

    /**
    * The Jacobian of fsub.  Until we decide on how to get AD support in FCST,
    * we must enter the Jacobian in manually.
    * Note that a one-dimensional array must be passed back to COLDAE.
    * However, it is easier to define a two-dimensional array as the matrix (see COLDAE.f).
    * After, use c_to_for_matrix to convert it to the correct one-dimensional array.
    */
    void dfsub ( double &, double [], double [], double [] );

    static void dfsub_wrapper ( double &, double [], double [], double [] );

    /**
    * Define the boundary conditions.
    *
    * There are 4 boundary conditions.  Note that i refers to the ith boundary condition.
    * See COLDAE.f
    */
    void gsub ( int &, double [], double & );

    static void gsub_wrapper ( int &, double [], double & );

    /**
    * The derivatives of the boundary conditions.
    *
    * See COLDAE.f
    */
    void dgsub ( int &, double [], double [] );

    static void dgsub_wrapper ( int &, double [], double [] );

    /**
    * The initial guess.
    *
    * This is optional, but a good idea for this problem.
    * If we do not provide this, CODAE will use a constant of 0.0 for an initial guess.
    */
    void guess ( double &, double [], double [], double [] );

    static void guess_wrapper ( double &, double [], double [], double [] );



    /** Object to store information on water, necessary for water filled agglomerate **/
    FuelCellShop::Material::LiquidWater water;

    // Transport Parameters
    /** Proton Diffusivity in Nafion and water */
    double D_H_NAF;
    double D_H_Water;

    /** Oxygen Diffusivity in Nafion and water */
    double D_O2_N;
    double DO2_Water;

    /** Henry's Constant for oxygen in nafion*/
    double HO2N;

    /*Permitivity and charge characteristics of the core*/
    double rel_permittivity_water;
    double rel_permittivity_Naf;
    double core_charge_factor;

    std::vector<std::string> derivative_flags;

    /** absolute tolerance in proton concentration */
    double H_tol;

    //Continuation variables
    double lambda;
    double lambda2;
    

};

}//namespace Layer
}//namespace FuelCellShop

#endif
