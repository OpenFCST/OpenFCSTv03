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

#include "DAE_wrapper.h"
#include <cmath>



namespace NAME = FuelCell::ApplicationCore;




// Must be very careful with this... it is now declared as a global variable within FCST
// These are the definitions of the pointers required for the DAE_wrapper functions
void *ptr_DAE_object[20]; //If you modify this please also correct DAESolver:DAE_Solve

//---------------------------------------------------------------------------
NAME::DAEWrapper::DAEWrapper ()
{
	cm_to_m = 0.01;
	cm2_to_m2 =  cm_to_m*cm_to_m;
	cm3_to_m3 = cm2_to_m2*cm_to_m;
}

//---------------------------------------------------------------------------
double
NAME::DAEWrapper::integrate (double lb, double ub, std::vector<double>& W, std::vector<double>& F)
{
// Make sure the two vectors of weights and function evaluations are of the same size.
// Assert (W.size() == F.size(), ExcDimensionMismatch(W.size(), F.size()));

double I=0.;
for (unsigned int i=0; i<W.size(); ++i)
	I += W[i]*F[i];

return 0.5*(ub-lb) * I;
}

//---------------------------------------------------------------------------
void
NAME::DAEWrapper::get_quadrature_points (double lb, double ub, std::vector<double>& X, std::vector<double>& W, FuelCell::ApplicationCore::DAESolver* prob)
{
// Find the order of the polynomial defined given by the solution
int orders = prob->get_ODE_order();
int poly_order = n_colloc + orders;

// Determine the optimal number of gauss points required for exactness
double np = ceil((poly_order - 1) / 2.0);
int N = np;
X.resize(N);
W.resize(N);

ae_int_t info;
real_1d_array x;
x.setlength(N);
real_1d_array w;
w.setlength(N);

//get quadrature points and weights
alglib::gqgenerategausslegendre(N, info, x, w);

if (info == 1)
{
	for (int i=0; i<N; ++i)
	{
		X[i] = 0.5 *(ub + lb + (ub-lb)*x[i]); // Transforms x from [-1,1] to [lb,ub]
		W[i] = w[i];
	}
}
}

//---------------------------------------------------------------------------
void
NAME::DAEWrapper::DAE_Error(int flag)
	{
		FcstUtilities::log << "FAILED TO OBTAIN A SOLUTION\n"
				<< "DAE_solve() returned code: " << flag << std::endl;
		FcstUtilities::log << "= 1 for normal return\n"
				<< "= 0 if the collocation matrix is singular.\n"
				<< "=-1 if the expected no. of subintervals exceeds storage specifications.\n"
				<< "=-2 if the nonlinear iteration has not converged.\n"
				<< "=-3 if there is an input data error.\n" << std::endl;
		//abort();
	}
