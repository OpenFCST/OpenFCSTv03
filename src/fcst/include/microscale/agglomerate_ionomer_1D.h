//---------------------------------------------------------------------------
// C++ Interface: agglomerate_ionomer_1D.h
//
// Description: Used to solve a system of equations representing a
// 			spherical ionomer-filled agglomerate.
//
// Authors: Jason Boisvert <jjb701@mail.usask.ca>, (C) 2011
// 			University of Saskatchewan,
//          Peter Dobson <pdobson@ualberta.ca>,
//          Philip Wardlaw,
//          Michael Moore,
//		    University of Alberta.

// Copyright: See COPYING file that comes with this distribution
//
//---------------------------------------------------------------------------

#ifndef FUEL_CELL__IONOMER_AGGLOMERATE_NUMERICAL_1D__H
#define FUEL_CELL__IONOMER_AGGLOMERATE_NUMERICAL_1D__H


//------------------------------
// FUEL CELL DECLARATIONS
//-----------------------------
#include <microscale/numerical_agglomerate_base.h>
#include <contribs/DAE_wrapper.h>




namespace FuelCellShop
{
namespace MicroScale
{

/**
 * \brief Class that solves an ionomer-filled agglomerate problem in 1D
 *
 * This class uses the DAE solver interface to solve the problem using COLDAE
 * (included in the contrib folder) to solve a boundary value problem.
 * <h3> Input parameters </h3>
 * LIST OF INPUT PARAMETERS FOR THE CLASS.
 * @code
 * ...
 * subsection MicroScale
 *   subsection IonomerAgglomerateNumerical
 *      set Proton Conductivity factor = 1.0 #Factor modifying ionomer proton conductivity within the agglomerate and film
 *      set Non Equilibrium BC Rate constant = 0.13 #Non equilibrium Reaction rate coefficient
 *      set Use non equilibrium BC = false #Use non-equilibrium BC as described by Suzukhi et al.
 *   end
 * end
 * @endcode
 */

class IonomerAgglomerate : public NumericalAgglomerateBase, public FuelCell::ApplicationCore::DAEWrapper
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
     * The ionomer filled agglomerate model does not modify macro scale porosity
     */
    virtual double aux_volume_fraction(){
        return 0;
    }

protected:

    /** Setup the variables in the problem required by the DAE Solver */
    void setup_DAE_solver ();

    /**
     * Function which implements a continuation based on the tolerance
     */
    int cont_tolerance (double start_tol, double end_tol);

    /*
     * Virtual function for returning film thickness in nano meters.
     *
     */
    virtual double get_film_thickness(){
        return delta_agg*1e9;
    }

    /*
     * Virtual function for returning agglomerate radius in nano meters.
     *
     */
    virtual double get_radius(){
        return r_agg*1e9;
    }

    /** Set the composition and structure of the agglomerate */
    virtual void set_structure ();

    ///@name Instance Delivery
    //@{
    static IonomerAgglomerate const* PROTOTYPE;

    /**
     * This member function is used to create an object of type MicroScaleBase
     */
    virtual boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> create_replica ()
    {
        return boost::shared_ptr<FuelCellShop::MicroScale::MicroScaleBase> (new FuelCellShop::MicroScale::IonomerAgglomerate ());
    }
    //@}

    /** Constructors */
    IonomerAgglomerate ( int verbose = 1 );
    IonomerAgglomerate(std::string name);

    /*
     * Protected virtual member function for declaring parameters
     */
    virtual void declare_parameters (ParameterHandler &param) const
    {
        NumericalAgglomerateBase::declare_parameters(param);
        param.enter_subsection(concrete_name);{
            param.declare_entry("Proton Conductivity factor", "1.0",
                    Patterns::Double());
            param.declare_entry("Non Equilibrium BC Rate constant", "1.0", Patterns::Double());
            param.declare_entry("Use non equilibrium BC", "false",
                    Patterns::Bool());
        }
        param.leave_subsection();




    }

    /*
     * Protected virtual member function for initializing parameters
     */
    virtual void initialize (ParameterHandler &param)
    {
        NumericalAgglomerateBase::initialize(param);
        param.enter_subsection(concrete_name);{
            cond_factor = param.get_double("Proton Conductivity factor");
            k_O2 = param.get_double("Non Equilibrium BC Rate constant");
            non_equil_bc = param.get_bool("Use non equilibrium BC");
        }
        param.leave_subsection();


    }

private:



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


    /** Used to modify the proton conductivity in agglomerate */
    double cond_factor;

    /** Rate constant for non-equilibrium boundary condition*/
    double k_O2;

    /** Bool to determine if the non-equilibrium boundary condition is used */
    bool non_equil_bc;

    // Transport Parameters
    /** Proton conductivity */
    double sigma_p;
    double lambda; //used for continuity
    bool uFD; //use finte differences 
    double endTol; //set tolerance for final solution
	double prev_effectiveness;

};

}
}

#endif
