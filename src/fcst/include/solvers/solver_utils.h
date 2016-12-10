//---------------------------------------------------------------------------
//    $Id: solver_utils.h 2605 2014-08-15 03:36:44Z secanell $
//
//    Copyright (C) 2008 by Marc Secanell
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#ifndef _SOLVER_UTILS__H
#define _SOLVER_UTILS__H

//-- deal.II 
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/petsc_parallel_block_sparse_matrix.h>

//-- OpenFCST
#include <application_core/application_wrapper.h>

// Include STD classes
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

// Use namespace of deal.II
using namespace dealii;

/**
 * This namespace is used to include routines that are used
 * in the solve() routine of several applications such as
 * AppCathode and AppPemfc.
 *
 * @author Marc Secanell, 2008
 */
namespace SolverUtils
{
    /**
     * Check that the diagonal of the matrix does not have any zeros. If it does,
     * the code issues a warning on optimized mode or throws an exception on debug mode
     */
    void check_diagonal(const BlockSparseMatrix<double>& A);

    /**
     * Output diagonal elements of the stiffness matrix to the screen
     */
    void output_diagonal(const BlockSparseMatrix<double>& A);

    /**
     * Print the diagonal elements of a matrix to a file named diag_matrix.dat
     *
     */
    void print_diagonal(const BlockSparseMatrix<double>& A,
            const std::string& file = std::string("diag_matrix.dat"));
    /**
     * This member function is used to make sure that the BlockSpareMatrix
     * has no zeros in the diagonal. If it has zeros, then, the diagonal is filled
     * with the average value between the largest and the smallest number in the
     * diagonal.
     *
     * In most problems you will need to use this class to remove zeros in the diagonal
     * for the equations that are not physical in certain domains. For example, when solving
     * the catalyst layer, the membrane potential inside the GDL is NOT a physical quantity.
     * Instead of adding a small number during assembly which would affect the solution,
     * we add zeros. Then, we modify the diagonal matrix so that the value is zero everywhere
     * on that domain and the fluxes at the boundaries between domains are not accounted for.
     *
     * PETScWrappers::MPI::SparseMatrix matrix implementation alsos included. It is very
     * similar to the BlockSparseMatrix<double> implementation.
     *
     * @note PETScWrappers::MPI::SparseMatrix implementation
     * requires that the matrix is writable, i.e. compressed
     *
     * @note PETScWrappers::MPI::SparseMatrix implementation compresses matrix before returning
     */
    void repair_diagonal (BlockSparseMatrix<double>& A);
#ifdef OPENFCST_WITH_PETSC
    void repair_diagonal (PETScWrappers::MPI::SparseMatrix & A);
#endif
    /**
     * This member function is used to make sure that the BlockSpareMatrix
     * has no zeros in the diagonal. If it has zeros, then we rewrite the equation
     * so that we ensure the solution at that node has a value of zero. This is done
     * by passing the solution vector to the function and setting the RHS to this value.
     * The zero in the system matrix is set to 1, so that we have \f$ 1(-\delta u) = \delta u \f$
     */
    void repair_diagonal (BlockSparseMatrix<double>& A, FuelCell::ApplicationCore::FEVector& ,  const FuelCell::ApplicationCore::FEVector& );


};


#endif //_SOLVER_UTILS__H
