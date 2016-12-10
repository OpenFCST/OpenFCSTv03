//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: app_read_mesh.cc
//    - Description:
//    - Developers: ?????????
//    - $Id: app_read_mesh.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELL__APP_READ_MESH__H
#define _FUELCELL__APP_READ_MESH__H

//-- OpenFCST
#include <application_core/optimization_block_matrix_application.h>

// Use namespace of deal.II
using namespace dealii;

namespace FuelCell
{
/**
 * This namespace is used for all those auxiliary classes that are used by AppCathode.
 * These are mainly classes that inherit Function<dim> and that are necessary in order to call
 * some subroutine from deal.II. For example InitialSolution is created in order to use the
 * deal.II class VectorInterpolate which in turn is used to set up the initial solution to the problem.
 */
namespace InitialSolution
{
/**
 * This class is used when solving the problem using Newton's method to provide an initial solution.
 * This function is called in VectorTools::interpolate(..,..,InitialSolution<dim> marc,...)
 * It provides a solution that satisfies Dirichlet boundaries and has a gradient.
 */
template <int dim>
//class InitialSolution
class AppReadMeshIC
            :
            public Function<dim>
{
public:
    /**
     * Constructor
     */
    AppReadMeshIC (std::vector<std::string> );
    /**
     * Destructor
     */
    ~AppReadMeshIC ();
    /**
     * This is the member function that computes the value of the initial
     * solution for a given point.
     */
    void vector_value ( const Point<dim> &p,
                        Vector<double> &v ) const;

    double value (const Point<dim> &/*p*/, const unsigned int) const;

    /**
     * Function to set the component names of the initial solution.
     */
    void set_solution_names ( std::vector<std::string> names )
    {
        component_names = names;
    }

private:

    /**
     * List of solution variables
     */
    std::vector<std::string> component_names;
};
} //end namespace InitialSolution


namespace Application
{
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
 *
 * WHAT DOES THIS APPLICATION DO?
 *
 * @author ???
 */
template <int dim>
class AppReadMesh
            :
            public FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>
{
public:

    /**
     * Constructor.
     * @note the pointer data is initialized to boost::shared_ptr<> (), this means that
     * the pointer is empty and when we do data.get() it will return 0. This is good because at ApplicationBase
     * constructor an ApplicationData will be constructed.
     */
    AppReadMesh ( boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data =
                     boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> () );

    /**
     * Destructor
     */
    ~AppReadMesh();
    /**
     * Declare all parameters that are needed for:
     *   - the computation of the equation coefficients
     *   - the control of the linear system solution
     *   - ...
     */
    virtual void declare_parameters ( ParameterHandler& param );

    /**
     * Set up how many equations are needed and
     * read in parameters for the parameter handler in order to initialize data
     */
    void _initialize ( ParameterHandler& param );

    /**
     * Call the other initialize routines from the inherited classes
     */
    virtual void initialize ( ParameterHandler& param );
    
    /**
     * Initialize nonlinear solution
     */
    virtual void initialize_solution (FuelCell::ApplicationCore::FEVector& initial_guess,
                                      std::shared_ptr<Function<dim> > initial_function = std::shared_ptr<Function<dim> >());

    /**
     * Integration of the rhs of the equations. Here we loop over the quadrature
     * points and over degrees of freedom in order to compute the right
     * hand side for each cell
     * This routine depends on the problem at hand and is called by residual() in DoF_Handler
     * class
     * @note This function is called residual because in the case of nonlinear systems
     * the rhs is equivalent to the residual
     */
    void cell_residual ( FuelCell::ApplicationCore::FEVector&,
                                 const typename DoFApplication<dim>::CellInfo& )
	 {}

    /**
     * Compute the value of all objective function and constraints
     */
    void cell_responses ( std::vector<double>& ,
                                  const typename DoFApplication<dim>::CellInfo&,
                                  const FuelCell::ApplicationCore::FEVector& )
	 {}

    /**
     * This class is used to evaluate all responses that do not require looping over cells.
     * An example of one of this types of constraints is the solid volume fraction.
     */
    void global_responses ( std::vector<double>&,
                                    const FuelCell::ApplicationCore::FEVector& )
	 {}


    /**
     * Reimplementation of the routine in the base class BaseApplication in namespace AppFrame so
     * that the right labels are outputed and so that I can compute and output the source term.
     */
    virtual void data_out ( const std::string &basename,
                            const FuelCell::ApplicationCore::FEVectors &src );

    /**
     * Reimplementation of the pure virtual solve function. Function does nothing for this application.
     *
     * @note: Does nothing.
     */
    virtual void solve(FEVector& dst, const FEVectors& src){
        //Do nothing
    }

protected:
    /**
     * Structure where we store the problem we want to solve.
     * Each vector component contains a string with the
     * name of the equation we want to solve
     * Then, the number of components is equation_names.size()
     */
    std::vector<std::string> equation_names;
    /**
     * Structure where we store the name of each component
     * in our problem. The component names are stored in the same
     * way as they are stored in the solution.
     */
    std::vector<std::string> component_names;

	/** Stores the design variable names so that the name can be appended to the .vtk file name. */
	std::vector<std::string> design_var;

	/** Stores the values of the design variables so that the number can be appended to the .vtk file name. */
	std::vector<double> design_var_value;

};

}
}

#endif
