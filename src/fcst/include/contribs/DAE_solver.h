//---------------------------------------------------------------------------
//    $Id: DAE_solver.h 2605 2014-08-15 03:36:44Z secanell $
//
//    Copyright (C) 2011 by Jason Boisvert
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__appframe__DAE_solver_h
#define __deal2__appframe__DAE_solver_h

#include <stdio.h>
#include <stdexcept>

#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_thread_num() 0
#endif

//Load COLDAE subroutines
extern "C" {
        void coldae_(int &, int &, int [], double &, double &,
        	double [], int [], int [], double [], double [], int [],
        	double [], int &,
        	void (*)(double &, double [], double [], double []),
        	void (*)(double &, double [], double [], double []),
        	void (*)(int &, double [], double &),
        	void (*)(int &, double [], double []),
        	void (*)(double &, double [], double [], double []));
	void appsln_ (double &, double [], double [], double [], int []);
}



namespace FuelCell
{
namespace ApplicationCore
{

	//Types
	typedef void (*fsub_ptr)(double &, double [], double [], double []);
	typedef void (*dfsub_ptr)(double &, double [], double [], double []);
	typedef void (*gsub_ptr)(int &, double [], double &);
	typedef void (*dgsub_ptr)(int &, double [], double []);
	typedef void (*guess_ptr)(double &, double [], double [], double []);

	//functions

	/**
	* A dummy guess function to be provided to COLDAE when a user wishes to
	* provide none.
	* @param x is a value between a <= x <= b.
	* @param z is the initial guess for z(u(x)).
	* @param y is the initial guess for y(x).
	* @param df is array that contains the
	* \f[ m_{i} \f]
	* derivative of u(x)
	*/
	void DAE_dummy_guess (double &x, double z[], double y[], double df[]);
	/**
	* Converts a FORTRAN 2D array to a C/C++ 2D array.
	* @param rows is the number of rows in the 2D array.
	* @param cols is the number of columns in the 2D array.
	* @param fmat is the 2D FORTRAN matrix.
	* @param cmat is the returned 2D C/C++ matrix.
	*/
	void for_to_c_matrix(int rows, int cols, double *fmat, double **cmat);
	/**
	* Converts a C/C++ 2D array to a Fortran 2D array.
	* @param rows is the number of rows in the 2D array.
	* @param cols is the number of columns in the 2D array.
	* @param cmat is the 2D C/C++ matrix.
	* @param fmat is the returned 2D FORTRAN matrix.
	*/
	void c_to_for_matrix(int rows, int cols, double ** cmat, double *fmat);



	/**
	* This class provides an interface to the Fortran 77 code COLDAE.
	* COLDAE solves multi-point boundary-value DAEs for a system of
	* mixed-order ODEs
	*\f[ u_{i}^{(m_{i})} =  f_{i} ( x; z(u(x)), y(x) ), \quad    i = 1, ... ,\tt{m\_comp}, \f]
        * and algebraic constraints
	*\f[ 0   =  f_{i}( x; z(u(x)), y(x) ), \quad   i = \tt{m\_comp}+1,...,\tt{m\_comp}+\tt{ny}, \f]
	* for
	* \f[ a < x < b. \f]
	* The DAE is subjected to the boundary conditions
	* \f[ g_{j}  ( \zeta_{j}; z(u(\zeta_{j})) ) = 0, \quad   j = 1, ... ,m^{*}, \f]
	* where
	* \f[ a \leq \zeta_{1} \leq \zeta_{2} \leq \dots \leq \zeta_{m^{*}} \leq b, \f]
	* \f[ m^{*} = \sum^{\tt{m\_comp}}_{i=1} m_{i}. \f]
	* The exact solution is represented by
	* \f[ z(u(x)) = ( u_{1}^{(1)}(x), u_{1}^{(m_{1}-1)}(x), \dots, u_{\tt{m\_comp}}^{m_{(\tt{m\_comp}-1)}}(x), \f]
	* where
	* \f[ u_{i}^{(m_{i})} \f]
	*is the ith derivative of
	* \f[ u_{i}. \f]
	*/
	class DAESolver
	{
		public:

		/** Constructor
		* @param m_comp is the number of ODEs defined in fsub.
		* @param ny is the number of algebraic constraints defined in fsub.
		* @param m is an array of size m_comp that contains the orders for each ODE defined in fsub.
		* @param a is the left-most boundary condition.
		* @param b is the right-most boundary condition.
		* @param zeta is an array that holds the location, in the problem domain, of each of the boundary conditions defined in gsub.
		*	The ith element in zeta is associated with the ith boundary condition defined in gsub.
		* @param fsub is a user-supplied function that defines the ODEs and algebraic constraints.
		* @param dfsub is a user-supplied function that defines the Jacobian matrix of fsub.
		* @param gsub is a user-supplied function that defines the boundary conditions.
		* @param dgsub is a user-supplied function that defines the partial derivatives of gsub.
		* @param guess is a optional user-supplied function that defines an initial guess.
		* @note ny should equal zero if solving a boundary-value ODE.
		*/
		DAESolver(int m_comp, int ny, int m[], double a, double b, double zeta[],
		void (*fsub)(double &, double [], double [], double []),
		void (*dfsub)(double &, double [], double [], double []),
		void (*gsub)(int &, double [], double &),
		void (*dgsub)(int &, double [], double []),
		void (*guess)(double &, double [], double [], double [])=NULL);

		/**
		 * Destructor - clear all data
		 */
		~DAESolver() {clear_mem();}

		/** Set tolerances for solution components.
		* @param ltol_size is the size of both arrays ltol and tol.
		* @param ltol is an array such that the ltol[i]=l specifies that tol[i] is assigned to the lth component of z(u).
		* @param tol is an array of tolerances.
		* @note Must be used before a call to DAE_solve().
		*/
		void set_tolerance(int ltol_size = 0, int *ltol  = NULL, double *tol = NULL);
		/** Indicate if the problem is linear.
		* @note Must be used before a call to DAE_solve().  If not, COLDAE treats the problem as a
		* 	nonlinear problem.
		*/
		void set_linear(void);
		/** Sets number of collocation points to be used in each subinterval.
		* @param pnts is the number of collocation points to be used in each subinterval.
		* @note Must be used before a call to DAE_solve().
		*/
		void set_collocation_points(int pnts);
		/** Sets the  initial mesh size.
		* @param pnts is the initial mesh size.
		* @note Must be used before a call to DAE_solve().
		*/
		void set_initial_mesh_size(int pnts);
		/** Sets the size of the integer array used by COLDAE.
		* @param ispace_size is the size of integer array.
		* @note Must be used before a call to DAE_solve().
		*/
		void set_integer_space(int ispace_size = 0);
		/** Sets the size and location of double array used by COLDAE.
		* @param fspace_size is the size of double array.
		* @param fspace is an optional parameter to pass a double array from a previous run.
		* @note Must be used before a call to DAE_solve().
		*/
		void set_float_space(int fspace_size= 0, double *fspace = NULL);
		/** Set output level for COLDAE.
		* @param level is the desired output level where
			level = -1 is for full output, level = 0 is for selected output, and level=1 is for no output.
		* @note Must be used before a call to DAE_solve().
		*/
		void set_output(int level);
		/** Set the fixed points in the mesh.
		* @param fixpnt_size is the number of fixed points.
		* @fixpnt is an array that contains the location of the fixed points.
		* @note Must be used before a call to DAE_solve().
		*/
		void set_fixpnts(int fixpnt_size= 1, double *fixpnt = NULL);


		/**
		 * Set the solver control parameter.
		 * @param control is the integer defined in ipar(10)
		 * -1 if the first relax factor is RSTART
		 * 0 if the problem is regular
		 * 1 if the newton iterations are not to be damped
		 * 2 if we are to return immediately upon  (a) two successive nonconvergences,
		 * or  (b) after obtaining error estimate for the first time.
		 */
		void set_solver_control(int control);


		/** Set index of DAE.
		* @param index is the index value of the DAE such that index =0,1,2.
		* @note Must be used before a call to DAE_solve().
		*/
		void set_DAE_index(int index);



		/** Allows for simple continuation
		* @note Should only be ran after atleast one use of DAE_solve
		*/
		void use_simple_cont(void);

		/** Get the order of ODE's. */
		inline int get_ODE_order(void) {return *ODEs_Orders;}


		/** Solves the DAE.
		* @return A flag indicating success is returned where 1 indicates a normal return,
		* 0 indicates the collocation matrix is singular,
		* -1 indicates the expected number of subintervals exceeds storage specifications,
		* -2 indicates the nonlinear iteration has not converged, and
		* -3 indicates there is an input data error.
		 */
		int DAE_solve(void);

		/** Return a solution for the DAE at point x.
		* @param x is a value between a <= x <= b.
		* @param z is an array that contains the solution z(u(z)).
		* @param y is an array that contains the solution y(x).
		*/
		void DAE_solution(double x, double z[], double y[]);

		/** Gets the integer array used by COLDAE.
		* @return A pointer to the integer array used by COLDAE is
		*	is returned.
		*/
		int *return_integer_space(void);

		/** Gets the double array used by COLDAE.
		* @return A pointer to the double array used by COLDAE is
		*	returned.
		*/
		double *return_float_space(void);

		/** Get the number of mesh points in the final mesh.
		* @return The size of the final mesh is returned.
		* @note This function should only be used after a call
		*	to DAE_solve().
		*/
		inline int get_size_final_mesh(void) { return (this->ispace[0] +1); }

		/** Returns a copy of the final mesh.
		* @param mesh is a pointer to an array used to store
		*	a copy of the mesh.
		* @note This function should only be used after a call
		*	to DAE_solve().
		*/
		void get_copy_final_mesh (double *mesh);

		/** Overloads delete operator */
// 		void operator delete (void *p);

		inline void use_sol_as_guess() { this->ispace[8] = 3; }

		protected:
		/** Set the ODE function.
		* @param fsub is a user-supplied function that
		*	defines the ODEs and  algebraic constraints.
		*/
		void set_fsub(void (*fsub)(double &, double [], double [], double []));
		/** Sets the Jacobian of the ODE function.
		* @param dfsub is a user-supplied function that
			defines the Jacobian matrix of fsub.
		*/
		void set_dfsub(void (*dfsub)(double &, double [], double [], double []));
		/** Sets the  boundary condition function.
		* @param gsub is a user-supplied function that
		*	defines the boundary conditions.
		*/
		void set_gsub(void (*)(int &, double [], double &));
		/** Sets the partial derivative of the boundary condition function
		* @param dgsub is a user-supplied function that defines
			the partial derivatives of gsub.
		*/
		void set_dgsub(void (*)(int &, double [], double []));
		/** Sets the  initial-guess function
		* @param guess is a optional user-supplied function
		*	that defines an initial guess.
		*/
		void set_guess(void (*)(double &, double [], double [], double []));


		/** * Pointer to DAE function */
		void (*DAE_fsub) (double &, double [], double [], double []);
		/** Pointer to wrapper for jacobian of fsub */
		void (*DAE_dfsub) (double &, double [], double [], double []);
		/**  Pointer to boundary condition function */
		void (*DAE_gsub) (int &, double [], double &);
		/** * Pointer to jacobian of boundary condition function */
		void (*DAE_dgsub) (int &, double [], double []);
		/** Pointer to geuss function */
		void (*DAE_guess) (double &, double [], double [], double []);

		/** Sets the problem size.
		* @param num_ODEs is the number of ODE.
		* @param num_Alg_Const is the number of algebraic constraints.
		* @param ODEs_Orders is an array that contains the orders of each
		*	of the ODEs defined in fsub.
		*/
		void set_prob_size(int num_ODEs, int num_Alg_Const, int *ODEs_Orders);
		/** Sets the left and right-most boundary points.
		* @param a is the left-most boundary point.
		* @param b is the right-most boundary point.
		*/
		void set_boundary_points(int a, int b);
		/** Sets additional boundary points
		* @param zeta_size is the size of the array zeta.
		* @param zeta is an array that holds boundary points other
			than a and b.
		*/
		void set_side_conditions(int zeta_size, double *zeta );
		/** Creates an array that contains various information about how to
			solve the DAE.*/
		void set_ipar(void);

		/** Deletes dynamic memory */
		void clear_mem(void);


		// user-supllied subroutines flags
		/** Boolian for DAE function */
		bool have_fsub;
		/** Boolian for jacobian of fsub wrapper function */
		bool have_dfsub;
		/**  Boolian for boundary condition function */
		bool have_gsub;
		/** Boolian for jacobian of boundary condition function */
		bool have_dgsub;
		/** Boolian for guess function */
		bool have_guess;


		/**  Array of boundary locations */
		double *zeta;
		/**  Size of zeta array */
		int zeta_size;
		/**  Set zeta flag */
		bool set_zeta;


		/** Array used to hold location of tolerances */
		int *ltol;
		/** Set ltol flag */
		bool set_ltol;
		/**  Size of ltol array */
		int ltol_size;

		/**  Array of tolerances */
		double *tol;
		/**  tol flag*/
		bool set_tol;

		/**  Array of info passed to COLDAE */
		int *ipar;

		/**  Linear flag */
		bool linear;

		/**  Number of collocation points */
		int collpnts;
		/**  Set collocation flag */
		bool set_collpnts;

		/**  Size of initial mesh*/
		int intialmeshsize;
		/**  Set initial mesh flag */
		bool set_intialmeshsize;

		/**  Integer array used by COLDAE */
		int *ispace;
		/** size of ispace */
		int ispace_size;
		/** ispace flag */
		bool set_ispace;

		/**  Double array used by COLDAE */
		double *fspace;
		/**  Size of fspace */
		int fspace_size;
		/**  Set fspace flag */
		bool set_fspace;
		/** USer old fspace for continuation */
		bool use_old_fspace;

		/**  Output level */
		int output_level;

		/**   Array of fixed points */
		double *fixpnt;
		/**  fixpnt flag*/
		bool set_fixpnt;
		/**  Size of fixpnt array*/
		int fixpnt_size;

		/** solver control flag */
		bool set_solvercontrol;
		/** solver control parameter (see ipar(10))*/
		int solvercontrol;

		/**  DAE index */
		int DAE_index;

		/** Number of ODEs */
		int num_ODEs;
		/**  Number of algebraic constants*/
		int num_Alg_Const;
		/**  A pointer to an array that contains the orders
			of each of the ODEs. */
		int *ODEs_Orders;
		/**  Leftmost boundary point */
		double a;
		/**  Rightmost boundary point */
		double b;
		/** use simple continuation */
		bool use_cont;
	};


}

}

#endif
