//---------------------------------------------------------------------------
// C++ Interface: DAE_wrapper.cc
//
// Description: Class used as a wrapper for all applications that will be
// 				solved using the DAE_solver.
//
// Author: Peter Dobson <pdobson@ualberta.ca>, (C) 2011
//		University of Alberta
//
// Copyright: See COPYING file that comes with this distribution
//
//---------------------------------------------------------------------------

#ifndef FUEL_CELL__DAE_WRAPPER__H
#define FUEL_CELL__DAE_WRAPPER__H

//------------------------------
// STD DECLARATIONS
//------------------------------
#include <iostream>


#include <boost/shared_ptr.hpp>
//------------------------------
// DEAL.II DECLARATIONS
//------------------------------
#include <deal.II/base/parameter_handler.h>

//------------------------------
// ALGLIB DECLARATIONS
//------------------------------
#include <integration.h>

//------------------------------
// FUEL CELL DECLARATIONS
//------------------------------
#include <contribs/DAE_solver.h>
#include <utils/logging.h>

/** \file */

// Must be very careful with this... it is now 'declared' as a global variable within FCST
// The 'definition' of the pointer in included in the .cc file.
// The implementation (assignment) of an object to the pointer should be within the derived class
#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_thread_num() 0
#endif
extern void *ptr_DAE_object[20];

using namespace dealii;
using namespace alglib;

namespace FuelCell
{
namespace ApplicationCore
{
	/**
	* \brief This class is to be used as a wrapper for the functions needed in the DAESolver
	*
	* All the functions defined within this class are pure and therefore must be reimplemented
	* 	in a derived class.  The functions here are define a Differential Algebraic Equation,
	* 	more specifically in this case, a system of ODE's.
	*
	* The derived classes must have the specific implementation of the functions and
	* 	static wrapper functions that point to the implementation.  This is so that
	* 	the DAESolver can be used in a generic form and interface with COLDAE.
	*
	* The static wrapper functions should take the form:
	*
	* 	\code static void fsub_wrapper (double &, double [], double [], double [])
	* 			DerivedClass *ptr_to_class = (DerivedClass*) ptr_DAE_object;
	*			ptr_to_class->fsub(x, z, y, f);
	*	\endcode
	*
	* 	where ptr_DAE_object is a void globally defined variable that can take any form.
	* 	This must be used with caution.
	*
	* \author Peter Dobson
	*/
	class DAEWrapper
	{
	public:

	/** Constructor */
	DAEWrapper();

	/** Destructor */
	~DAEWrapper(){};

	/**
	 * Member function that integrates a solution between a lower and upper bound.
	 *
	 * The function takes the values of the bounds, and the function evaluated at the quadrature points
	 *
	 * The integral is computed using gauss quadrature.
	 */
	double integrate(double lb, double ub, std::vector<double>& W, std::vector<double>& F);

	/**
	 * Function that obtains the gaussian quadrature points and weights.
	 *
	 * Number of points and weights based on based on the order of the polynomial, given by n_colloc + mm[i].
	 *
	 * The points returned in @param X are transformed to be points between @param lb and @param ub.
	 */
	void get_quadrature_points (double lb, double ub, std::vector<double>& X, std::vector<double>& W, FuelCell::ApplicationCore::DAESolver* prob);

	/** Setup the variables in the problem required by the DAE Solver */
	virtual void setup_DAE_solver () = 0;

	/**
	* Define the DAE function.  In this case,it is simply a system of ODES.
	* COLDAE allows for the system to be defined as a system of 2 mixed-order ODEs.
	* See COLDAE.f for additional information about how to define fsub.
	* In particular, see COLDAE.f for the meaning of z and y.
	* Note that because a BVP is solved, y is not used.
	*/
	virtual void fsub (double &, double [], double [], double []) = 0;


	/**
	* The Jacobian of fsub.  Until we decide on how to get AD support in FCST,
	* we must enter the Jacobian in manually.
	* Note that a one-dimensional array must be passed back to COLDAE.
	* However, it is easier to define a two-dimensional array as the matrix (see COLDAE.f).
	* After, use c_to_for_matrix to convert it to the correct one-dimensional array.
	*/
	virtual void dfsub (double &, double [], double [], double []) = 0;

	/**
	* Define the boundary conditions.
	*
	* There are 4 boundary conditions.  Note that i refers to the ith boundary condition.
	* See COLDAE.f
	*/
	virtual void gsub (int &, double [], double &) = 0;

	/**
	* The derivatives of the boundary conditions.
	*
	* See COLDAE.f
	*/
	virtual void dgsub (int &, double [], double []) = 0;

	/**
	* The initial guess.
	*
	* This is optional, but a good idea for this problem.
	* If we do not provide this, CODAE will use a constant of 0.0 for an initial guess.
	*/
	virtual void guess (double &, double [], double [], double []) = 0;

	/** Set the verbosity variable (controls output to screen) */
	inline void verbosity(int i)
	{n_output = i;}

	/**
	 * Indicates error in the solve function
	 *
	 * Prints error message and aborts the program.
	 */
	void DAE_Error(int flag);

	void clear_memory()
	{
		delete [] mm;
		delete [] zeta;
		delete [] fixpnt;
		delete [] tol;
		delete [] ltol;
		delete prob;
		if (ptr_DAE_object[omp_get_thread_num()] != NULL)
		{
			ptr_DAE_object[omp_get_thread_num()] = NULL;
			//delete ptr_DAE_object;
		}
	}

	protected:

	/** Number of mesh points */
	int n_mesh;

	/** Array of mesh points */
	double *mesh;

	/** Number of collocation points */
	int n_colloc;

	/** Output integer variable */
	int n_output;

	/** number of PDEs */
	int n_comp;

	/** number of Algebraic constraints */
	int n_y;

	/** array of integers storing the order of each PDE */
	int *mm;

	/** Integer representing the total number of variables
	 * 	given by \f$ \sum mm[i] \f$
	 */
	int m_star;

	/** Left boundary point */
	double boundary_0;

	/** Right boundary point */
	double boundary_1;

	/** Array of boundary points */
	double *zeta;

	/** DAE problem solver object */
	FuelCell::ApplicationCore::DAESolver *prob;

	/** Array of fixed points on the mesh */
	double *fixpnt;

	int *ltol;

	double *tol;

	/** Convert from centimetres to metres. */
	double cm_to_m;

	/** Convert from centimetres squared to metres squared. */
	double cm2_to_m2;

	/** Convert from centimetres cubed to metres cubed. */
	double cm3_to_m3;
	};

}
}
#endif
