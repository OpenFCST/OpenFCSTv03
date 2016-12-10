// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: linear_solvers.h
// - Description: This namespace contains various linear solvers and preconditioners
// - Developers: Valentin N. Zingan, University of Alberta
//               Marc Secanell Gallart, University of Alberta
//
// ----------------------------------------------------------------------------

#ifndef _FCST_LINEAR_SOLVERS_H_
#define _FCST_LINEAR_SOLVERS_H_

//-- dealII
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_matrix_array.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_ilu.h>

using namespace dealii;

/**
 * This namespace contains
 * various linear solvers
 * and preconditioners.
 *
 * Linear solvers:
 *
 * - SparseDirectUMFPACK solver,
 * - GMRES solver,
 * - Schur complement based solver.
 *
 * Preconditioners:
 *
 * - ILU preconditioner.
 *
 * \author Valentin N. Zingan, 2013
 * \author Marc Secanell Gallart, 2013
 */

namespace LinearSolvers
{
    
    /**
     * Memory management objects.
     */
    static GrowingVectorMemory< Vector<double> > vector_pool;
    static GrowingVectorMemory< BlockVector<double> > block_vector_pool;
    
    /**
     * This class implements
     * an interface to the sparse direct solver UMFPACK,
     * see the link below
     * http://www.cise.ufl.edu/research/sparse/umfpack/.
     *
     * \author Valentin N. Zingan, 2013
     */
    
    class SparseDirectUMFPACKSolver
    {
    public:
        
        ///@name Constructors, destructor, and initialization
        //@{
        
        /**
         * Constructor.
         */
        SparseDirectUMFPACKSolver() { FcstUtilities::log << "Solving the linear system using UMFPACK" << std::endl; }
        
        /**
         * Destructor.
         */
        ~SparseDirectUMFPACKSolver() { }
        
        //@}
        
        ///@name Solve function
        //@{
        
        /**
         * This function solves
         * the linear system Ax = b
         * where
         *
         * - \p A = \p matrix,
         * - \p x = \p solution,
         * - \p b = \p right_hand_side.
         */
        template<typename MATRIX, typename VECTOR>
        void solve(const MATRIX& matrix,
                   VECTOR&       solution,
                   const VECTOR& right_hand_side) const
                   {
                       Vector<double> tmp;
                       tmp = right_hand_side;
                       
                       SparseDirectUMFPACK solver;
                       solver.initialize(matrix);
                       solver.solve(tmp);
                       
                       solution = tmp;
                   }
                   
                   //@}
        
    };
    
    /**
     * This class implements
     * GMRES solver.
     *
     * \author Valentin N. Zingan, 2013
     * \author Marc Secanell Gallart, 2013
     */
    
    class GMRESSolver
    {
    public:
        
        ///@name Constructors, destructor, and initialization
        //@{
        
        /**
         * Constructor.
         */
        GMRESSolver() { FcstUtilities::log << "Solving the linear system using GMRES" << std::endl; }
        
        /**
         * Destructor.
         */
        ~GMRESSolver() { }
        
        //@}
        
        ///@name Solve function
        //@{
        
        /**
         * This function solves the linear system Ax = b where
         * - \p A = \p matrix,
         * - \p x = \p solution,
         * - \p b = \p right_hand_side.
         *
         * Other data includes
         * - \p n_iter is the max number of GMRES iterations,
         * - \p tolerance is the tolerance of GMRES solver,
         * - \p preconditioner is the preconditioner used.
         */
        template<typename MATRIX, typename VECTOR, typename PRECONDITIONER>
        void solve(const MATRIX&         matrix,
                   VECTOR&               solution,
                   const VECTOR&         right_hand_side,
                   const unsigned int&   n_iter,
                   const double&         tolerance,
                   const PRECONDITIONER& preconditioner) const
        {
            SolverControl solver_control(n_iter, tolerance);
            
            typename SolverGMRES< VECTOR >::AdditionalData additional_data(100, false, true);
            
            SolverGMRES< VECTOR > solver(solver_control, additional_data);

            solver.solve(matrix,
                         solution,
                         right_hand_side,
                         preconditioner);
            
        }
        
        /**
         * This routine is used when solver_control is an object already initialized.
         * 
         * This function solves the linear system Ax = b
         * where
         *
         * - \p A = \p matrix,
         * - \p x = \p solution,
         * - \p b = \p right_hand_side.
         * - \p solver_control
         */
        template<typename MATRIX, typename VECTOR, typename PRECONDITIONER>
        void solve(SolverControl& solver_control,
                   const MATRIX&         matrix,
                   VECTOR&               solution,
                   const VECTOR&         right_hand_side,
                   const PRECONDITIONER& preconditioner) const
        {
            SolverGMRES< VECTOR > solver(solver_control);
            
            solver.solve(matrix,
                         solution,
                         right_hand_side,
                         preconditioner);
            
        }
        //@}
                              
    };
    
    /**
     * This class implements
     * ILU preconditioner.
     *
     * \author Valentin N. Zingan, 2013
     * \author Marc Secanell Gallart, 2013
     */
    
    class ILUPreconditioner
    {
    public:
        
        ///@name Constructors, destructor, and initialization
        //@{
        
        /**
         * Constructor.
         */
        ILUPreconditioner(const BlockSparseMatrix<double>& matrix)
        :
        preconditioner_matrix(matrix.n_block_rows()),
        preconditioner_pointer(matrix.n_block_rows()),
        preconditioner(matrix.n_block_rows())
        {
            FcstUtilities::log << "ILUPreconditioner is used" << std::endl;
                        
            for(unsigned int i = 0; i < matrix.n_block_rows(); ++i)
            {
                preconditioner_matrix[i].initialize(matrix.block(i,i),
                                                    SparseILU<double>::AdditionalData());
                preconditioner_pointer[i].set_memory(&vector_pool);
                preconditioner_pointer[i] = &preconditioner_matrix[i];
            }
            
            for(unsigned int i = 0; i < matrix.n_block_rows(); ++i)
                for(unsigned int j = 0; j < matrix.n_block_cols(); ++j)
                {
                    if( i == j )
                        preconditioner.enter(preconditioner_pointer[i], i, i);
                    else
                        preconditioner.enter(matrix.block(i,j), i, j);
                }
        }
        
        /**
         * Destructor.
         */
        ~ILUPreconditioner() { }
        
        //@}
        
        ///@name Data
        //@{
        
        /**
         * Preconditioner.
         */
        BlockTrianglePrecondition<double> preconditioner;
        
        
        /**
         * Preconditioner matrix.
         */
        std::vector< SparseILU<double> > preconditioner_matrix;
        
        /**
         * Preconditioner pointer.
         */
        std::vector< PointerMatrixAux< SparseILU<double>, Vector<double> > > preconditioner_pointer;
        
        
        //@}
        
    };
    
} // LinearSolvers

#endif
