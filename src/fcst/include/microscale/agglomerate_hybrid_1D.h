//---------------------------------------------------------------------------
// C++ Interface: agglomerate_hybrid_1D.h
//
// Description: Used to solve a system of equations representing a
// 			spherical hybrid agglomerate.
//
// Author: Philip Wardlaw 2014,
//         University of Alberta
//
// Copyright: See COPYING file that comes with this distribution
//
//---------------------------------------------------------------------------

#ifndef FUEL_CELL__HYBRID_AGGLOMERATE_1D__H
#define FUEL_CELL__HYBRID_AGGLOMERATE_1D__H


//------------------------------
// Inherits everything necessary from water filled agglomerate
//------------------------------

#include <microscale/agglomerate_water_1D.h>


namespace FuelCellShop
{
namespace MicroScale
{

/**
 * \brief Class that solves a Hybrid (water/ionomer filled) agglomerate problem in 1D
 *
 * This class uses the DAE solver interface to solve the problem using COLDAE
 * (included in the contrib folder) to solve a boundary value problem.
 *
 */
class HybridAgglomerate : public WaterAgglomerate
{
public:



    static const std::string concrete_name;

    /**
     * Return name of class instance, i.e. concrete name.
     */
    virtual std::string get_name(){
        return concrete_name;
    }

    /**
    * Functions directly inherited from parent
    *
    * -double compute_current (  double &E_r)
    * -virtual std::string get_name()
    *
    */




    /**
     * MicroScale object may have extra contribution to volume of layer, e.g. water.
     *
     * The assumption of a hybrid core means that this contribution is different
     * than in the pure water filled case, thus it is reimplemented here.
     */
    virtual double aux_volume_fraction(){

        return WaterAgglomerate::aux_volume_fraction()*(1-hybrid_core_volume_fraction);

    }

    /**
     * Print out key agglomerate information (name, radius, film thickness, porosity).
     */
    virtual void print_properties();


protected:

    /**
     * Functions directly inherited from parent
     *
     * -int cont_tolerance (double start_tol, double end_tol)
     * -int cont_cd ()
     * -void setup_DAE_solver ()
     * -virtual double get_film_thickness()
     * -virtual double get_radius()
     * -void gsub ( int &, double [], double & )
     * -static void gsub_wrapper ( int &, double [], double & )
     * -void dgsub ( int &, double [], double [] )
     * -static void dgsub_wrapper ( int &, double [], double [] )
     * -void guess ( double &, double [], double [], double [] )
     * -static void guess_wrapper ( double &, double [], double [], double [] )
     */




    /** Set the composition and structure of the agglomerate
     *
     *  Structure is different to the pure water filled case, hence it is
     *  reimplemented here.
     *
     */
    void set_structure ( );

    /** Constructors */
    HybridAgglomerate ( int verbose = 1 );
    HybridAgglomerate (std::string concrete_name);


    /*
     * Protected virtual member function for declaring parameters
     */
    virtual void declare_parameters (ParameterHandler &param) const
    {
        NumericalAgglomerateBase::declare_parameters(param);
        water.declare_parameters(param);

        param.enter_subsection(concrete_name);{
            param.declare_entry("Hybrid core fraction", "0.01",
                    Patterns::Double(0.0,0.75),
                    "Fraction of the agglomerate core filled with ionomer.");

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
            hybrid_core_volume_fraction = param.get_double("Hybrid core fraction");
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
    static HybridAgglomerate const* PROTOTYPE;

    /**
     * This member function is used to create an object of type MicroScaleBase
     */
    virtual boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> create_replica ()
    {
        return boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> (new FuelCellShop::MicroScale::HybridAgglomerate());
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




    /** Object inheritted from water filled agglomerate
     *
     * -FuelCellShop::Material::LiquidWater water
     *
     * Transport Parameters
     * Proton Diffusivity in Nafion and water
     * -double D_H_NAF
     * -double D_H_Water
     *
     * Oxygen Diffusivity in Nafion and water
     * - double D_O2_N
     * -double DO2_Water
     *
     * Henry's Constant for oxygen in nafion
     * -double HO2N
     * -double rel_permittivity_water
     * -double rel_permittivity_Naf
     * -std::vector<std::string> derivative_flags
     *
     * absolute tolerance in proton concentration
     * double H_tol
     **/

private:

    double hybrid_interface;
    double hybrid_core_volume_fraction;
};

}//namespace Layer
}//namespace FuelCellShop

#endif
