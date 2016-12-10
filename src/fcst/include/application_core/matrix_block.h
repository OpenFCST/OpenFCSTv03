// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2009 by Guido Kanschat
// Copyright (C) 2006-2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: matrix_block.h
// - Description: This class implements MATRIX object wrapper
// - Developers: Guido Kanschat,     Texas A&M University
//               Valentin N. Zingan, University of Alberta
//               Marc Secanell,      University of Alberta
// - Id: $Id: matrix_block.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _FUEL_CELL_APPLICATION_CORE_MATRIX_BLOCK_H_
#define _FUEL_CELL_APPLICATION_CORE_MATRIX_BLOCK_H_

#include <deal.II/lac/full_matrix.h>

using namespace dealii;

namespace FuelCell
{
namespace ApplicationCore
{

/**
 * A wrapper around a MATRIX object, storing the coordinates in the
 * global block matrix as well.
 *
 * @author Guido Kanschat
 */

template<typename MATRIX>
struct MatrixBlock
{
       /**
        * The matrix itself.
        */
       MATRIX matrix;

       /**
        * Row coordinate
        * in the global block matrix.
        */
       unsigned int row;

       /**
        * Column coordinate
        * in the global block matrix.
        */
       unsigned int column;



       // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



       /**
        * Constructor.
        */
       MatrixBlock()
       :
       row(static_cast<unsigned int>(-1)),
       column(static_cast<unsigned int>(-1))
       { }



       /**
        * Constructor.
        */
       MatrixBlock(unsigned int row,
                   unsigned int column)
       :
       row(row),
       column(column)
       { }



       /**
        * Copy constructor.
        */
       MatrixBlock(const MatrixBlock<MATRIX>& M)
       :
       matrix(M.matrix),
       row(M.row),
       column(M.column)
       { }
};

/**
 * The matrix vector used in the mesh loops.
 */
typedef std::vector< MatrixBlock< FullMatrix<double> > > MatrixVector;

/**
 * The std::vector of dealii::Vectors used in the mesh loops.
 */
typedef std::vector< Vector<double> > VectorVector;

} // ApplicationCore

} // FuelCell

#endif