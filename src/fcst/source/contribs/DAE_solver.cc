//---------------------------------------------------------------------------
//    $Id: DAE_solver.cc 2605 2014-08-15 03:36:44Z secanell $
//
//    Copyright (C) 2011 by Jason Boisvert, University of Saskatchewan
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include "DAE_solver.h"
#include <iostream>


using namespace FuelCell::ApplicationCore;

void FuelCell::ApplicationCore::DAE_dummy_guess (double &/*x*/, double /*z*/[], double /*y*/[], double /*df*/[])
{
}
//------------------------------------------------------------------------------

void FuelCell::ApplicationCore::for_to_c_matrix(int rows, int cols, double *fmat, double **cmat)
{

	int k=0;
	for (int i=0; i < cols; i++)
	{
		for (int j=0; j < rows; j++)
		{
			cmat[j][i] = fmat[k++];
		}
	}	

}

//------------------------------------------------------------------------------

void FuelCell::ApplicationCore::c_to_for_matrix(int rows, int cols, double ** cmat, double *fmat)
{

	int k=0;
	for (int i=0; i < cols; i++)
	{
		for (int j=0; j < rows; j++)
		{
			fmat[k++] = cmat[j][i];
		}
	}

}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
DAESolver::DAESolver (int m_comp, int ny, int m[], double a, double b, double zeta[],
	void (*fsub)(double &, double [], double [], double []),
	void (*dfsub)(double &, double [], double [], double []),
	void (*gsub)(int &, double [], double &),
	void (*dgsub)(int &, double [], double []),
	void (*guess)(double &, double [], double [], double [])):
	set_zeta(false),
	set_ltol(false),
	set_tol(false),
	linear(false),
	set_collpnts(false),
	set_intialmeshsize(false),
	set_ispace(false),
	set_fspace(false),
	use_old_fspace(false),
	output_level(1),
	set_fixpnt(false),
	set_solvercontrol(false),
	DAE_index(0),
	use_cont(false)

	//Assign Parameters
{
	
	//Set info for DAE
	this->set_prob_size(m_comp,ny,m);
	
	//Set the boundary points
	this->set_boundary_points(a,b);

	//set side conditions
	int mcomp=0; // sum of orders
	for(int i = 0; i < this->num_ODEs; i ++)
	{
			mcomp += this->ODEs_Orders[i];
	}
	this->set_side_conditions(mcomp,zeta);

	//Set the user supplied functions
	this->set_fsub(fsub);
	this->set_dfsub(dfsub);
	this->set_gsub(gsub);
	this->set_dgsub(dgsub);
	this->set_guess(guess);
}
//------------------------------------------------------------------------------

void DAESolver::set_prob_size(int num_ODEs, int num_Alg_Const, int * ODEs_Orders)
{
	this->num_ODEs = num_ODEs;
	this->num_Alg_Const = num_Alg_Const;
	this->ODEs_Orders = ODEs_Orders;
}
//------------------------------------------------------------------------------

void DAESolver::set_boundary_points(int a, int b)
{
	this->a=a;
	this->b=b;
}
//------------------------------------------------------------------------------

void DAESolver::set_side_conditions(int zeta_size, double *zeta)
{
	
		this->zeta=zeta;
		this->set_zeta=true;
		this->zeta_size = zeta_size;
	
	
}
//------------------------------------------------------------------------------

void DAESolver::set_fsub(void (*fsub)(double &, double [], double [], double []))
{
	this->DAE_fsub = fsub;
}
//------------------------------------------------------------------------------

void DAESolver::set_dfsub(void (*dfsub)(double &, double [], double [], double []))
{
	this->DAE_dfsub=dfsub;
}
//------------------------------------------------------------------------------

void DAESolver::set_gsub(void (*gsub)(int &, double [], double &))
{
	this->DAE_gsub=gsub;
}
//------------------------------------------------------------------------------

void DAESolver::set_dgsub(void (*dgsub)(int &, double [], double []))
{
	this->DAE_dgsub=dgsub;
}
//------------------------------------------------------------------------------

void DAESolver::set_guess(void (*guess)(double &, double [], double [], double []))
{
	if (guess == NULL) 
	{
		this->have_guess = false;
		this->DAE_guess = &DAE_dummy_guess;
	}
	else
	{
		this->have_guess = true;
		this->DAE_guess=guess;
	}
	
}
//------------------------------------------------------------------------------

void DAESolver::set_tolerance(int ltol_size, int *ltol, double *tol)
{

	// attempt to get user supplied tolerances
	if( ltol != NULL && tol != NULL )
	{
		this->ltol = ltol;
		this->set_ltol = true;
		this->tol = tol;
		this->set_tol = true;
		this->ltol_size = ltol_size;
	}
	else
	{
		int mcomp=0; // sum of orders
		for(int i = 0; i < this->num_ODEs; i ++)
		{
			mcomp += this->ODEs_Orders[i];
		}
		int *tmp_ltol = new int[mcomp];
		double *tmp_tol = new double[mcomp];
		//populate ltol
		int k=0;
		for(int i = 0; i < this->num_ODEs; i ++)
		{
			for (int j=0; j < this->ODEs_Orders[i]; j++)			
			{
				tmp_ltol[k] = k+1;
				tmp_tol[k] = 1E-6;
				k++;
			}
		} 
		this->ltol = tmp_ltol;
		this->set_ltol = true;
		this->tol = tmp_tol;
		this->set_tol = true;
		this->ltol_size=mcomp;
		
		delete [] tmp_ltol;
		delete [] tmp_tol;
	}
	
}
//------------------------------------------------------------------------------

void DAESolver::set_linear(void)
{
	this->linear=true;
}

//------------------------------------------------------------------------------

void DAESolver::set_collocation_points(int pnts)
{
	this->collpnts = pnts;
	this->set_collpnts = true;
}

//------------------------------------------------------------------------------

void DAESolver::set_initial_mesh_size(int pnts)
{
	this->intialmeshsize = pnts;
	this->set_intialmeshsize = true;
}

//------------------------------------------------------------------------------

void DAESolver::set_integer_space(int ispace_size)
{
	//Ensure we do not have more than one copy of fspace in memory
	if(this->set_ispace ==  true)
	{
		//FcstUtilities::log << "Cleaning ispace" << std::endl;
		delete [] this->ispace;
		this->set_ispace = false;
		this->ispace_size = 0;
	} 

	if(ispace_size > 0)
	{
		this->ispace_size = ispace_size;
	}
	else
	{
		int mcomp = 0;
		for(int i = 0; i < this->num_ODEs; i ++) mcomp += this->ODEs_Orders[i];
		this->ispace_size = mcomp * 1000000;
	}
	int *tmp_ispace = new int[this->ispace_size];
	this->ispace = tmp_ispace;
	this->set_ispace = true;
}

//------------------------------------------------------------------------------

void DAESolver::set_float_space(int fspace_size, double *fspace)
{

	double *tmp_fspace;

	//Ensure we do not have more than one copy of fspace in memory
	if(this->set_fspace ==  true)
	{
		//FcstUtilities::log << "Cleaning fspace" << std::endl;
		delete [] this->fspace;
		this->set_fspace = false;
		this->fspace_size = 0;
	} 

	if (fspace != NULL)
	{
		this->use_old_fspace=true;
		tmp_fspace = fspace;
	}
	else if(fspace_size > 0)
	{
		this->fspace_size = fspace_size;
		tmp_fspace = new double[this->fspace_size];
	}
	else
	{
		int mcomp = 0;
		for(int i = 0; i < this->num_ODEs; i ++) mcomp += this->ODEs_Orders[i];
		this->fspace_size = mcomp * 10000000;
		tmp_fspace = new double[this->fspace_size];
		
	}
	this->fspace = tmp_fspace;
	this->set_fspace = true;
	
}

//------------------------------------------------------------------------------

void DAESolver::set_output(int level)
{
	if (level >= -1 && level <= 1) this->output_level = level;
	else level = 0;
}

//------------------------------------------------------------------------------

void DAESolver::set_fixpnts(int fixpnt_size, double *fixpnt)
{
	if(fixpnt != NULL)
	{
		this->fixpnt = fixpnt;
		this->fixpnt_size = fixpnt_size;
	}
	else
	{
		double *tmp_fixpnt = new double [1];
		this->fixpnt = tmp_fixpnt;
		this->fixpnt_size = 0;
		//delete [] tmp_fixpnt;
	}
	this->set_fixpnt=true; 
}


//------------------------------------------------------------------------------

void DAESolver::set_solver_control(int control)
{
	this->solvercontrol = control;
	this->set_solvercontrol = true;
}

//------------------------------------------------------------------------------

void DAESolver::set_DAE_index(int index)
{
	if (index >= 0 && index <= 2) this->DAE_index = index;
	else this->DAE_index = 0;
}

//------------------------------------------------------------------------------

void DAESolver::use_simple_cont(void)
{
	this->ipar[8]=2;
	this->ipar[2]= this->get_size_final_mesh() -1;
	this->use_cont=true;
}


//------------------------------------------------------------------------------

void DAESolver::set_ipar(void)
{
	int *tmp_ipar = new int[12]; 
	//initialize
	for (int i=0; i < 12; i++) tmp_ipar[i] = 0;
	//linear
	if (this->linear == false) tmp_ipar[0] = 1;
	//collocation points
	if (this->set_collpnts == true) tmp_ipar[1] = this->collpnts;
	// initial mesh size
	if (this->set_intialmeshsize == true) tmp_ipar[2] = this->intialmeshsize;
	// ltol size
	tmp_ipar[3] = this->ltol_size;
	//fspace size
	tmp_ipar[4] = this->fspace_size;
	// ipsace size
	tmp_ipar[5] = this->ispace_size;
	// output control
	tmp_ipar[6] = this->output_level;
	// continuation
	if (this->use_old_fspace)
	{
		tmp_ipar[7] = 2;
		tmp_ipar[8] = 2;
	}
	else
	{
		tmp_ipar[7] = 0;
		if (this->have_guess == false) tmp_ipar[8] = 0;
		else tmp_ipar[8] = 1;
	}
	// solver control
	if(this->set_solvercontrol == true) tmp_ipar[9] = this->solvercontrol;
	//fixpoints	
	tmp_ipar[10] = this->fixpnt_size;
	//dae index
	tmp_ipar[11] = this->DAE_index;
	this->ipar = tmp_ipar;
	
	//delete [] tmp_ipar;
}

//------------------------------------------------------------------------------

int DAESolver::DAE_solve(void)
{
	if (this->set_ltol == false && this->set_tol == false) this->set_tolerance();
	if (this->set_ispace == false) this->set_integer_space();
	if (this->set_fspace == false) this->set_float_space();
	if (this->set_fixpnt == false) this->set_fixpnts();
	if (this->use_cont == false) this->set_ipar();

	//get parameters and arrays
	int ncomp = this->num_ODEs;
	int ny = this->num_Alg_Const;
	int *m = this->ODEs_Orders;
	double aleft = this->a;
	double bright = this->b;
	double *zeta = this->zeta;
	int *ipar = this->ipar;
	int *ltol = this->ltol;
	double *tol = this->tol;
	double *fixpnt = this->fixpnt;
	int *ispace = this->ispace;
	double *fspace = this->fspace;
	int iflag=0;

	fsub_ptr fsub = this->DAE_fsub;
	gsub_ptr gsub = this->DAE_gsub;
	dgsub_ptr dgsub = this->DAE_dgsub;
	guess_ptr guess = this->DAE_guess;
	dfsub_ptr dfsub = this->DAE_dfsub;

	if (omp_get_thread_num() > 19)
		std::out_of_range("ptr_DAE_object is out of range, modify to run more threads");

	coldae_(ncomp, ny, m, aleft, bright, zeta, ipar, ltol, tol,                  
                fixpnt, ispace, fspace, iflag, fsub,                         
                 dfsub, gsub, dgsub, guess); 

             
    return iflag;
	
	
}

//------------------------------------------------------------------------------

void DAESolver::DAE_solution(double x, double z[], double y[])
{
	int *ispace = this->ispace;
	double *fspace = this->fspace;
	appsln_(x,z,y,fspace,ispace);
}

//------------------------------------------------------------------------------

int *DAESolver::return_integer_space(void)
{
	return (this->ispace);
}

//------------------------------------------------------------------------------

double *DAESolver::return_float_space(void)
{
	return (this->fspace);
}

//------------------------------------------------------------------------------
void DAESolver::get_copy_final_mesh (double *mesh)
{
	int npts = this->get_size_final_mesh();
	for (int i=0; i < npts; i++)
	{
		mesh[i] = this->fspace[i];
	}
} 

//------------------------------------------------------------------------------

void  DAESolver::clear_mem(void)
{
	delete [] this->ispace;
	delete [] this->fspace;
	delete [] this->ipar;
}

// -----------------------------------------------------------------------------

// //------------------------------------------------------------------------------
// 
// void DAESolver::operator delete (void *p)
// {
// 	DAESolver *tmp = (DAESolver *)p;
// 	tmp->clear_mem();
// 	
// }

//------------------------------------------------------------------------------


