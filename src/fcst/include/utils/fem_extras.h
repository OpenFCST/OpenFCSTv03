// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: fem_extras.h
// - Description: This namespace contains FEM related methods
// - Developers: Valentin N. Zingan, University of Alberta
//
// ----------------------------------------------------------------------------

#ifndef _FSCT_FEM_EXTRAS_H_
#define _FSCT_FEM_EXTRAS_H_

#ifndef dimension
#if   deal_II_dimension == 1
#define _1D_
#elif deal_II_dimension == 2
#define _2D_
#elif deal_II_dimension == 3
#define _3D_
#endif
#endif

#include <vector>

#include <deal.II/base/point.h>
#include <deal.II/base/exceptions.h>

using namespace dealii;

/**
 * This namespace contains
 * FEM related methods which are not
 * implemented in the standard deal.II
 * library.
 *
 * \author Valentin N. Zingan, 2013
 */

namespace FemExtras
{

  /**
   * This function computes
   * the tangential vectors \p dst
   * in the quadrature points of a \p dim - \p 1 dimensional face
   * based on the normal vectors \p src previously computed
   * in the same quadrature points of the same face.
   */
  template<int dim>
  void get_tangential_vectors(std::vector< std::vector< Point<dim> > >& dst,
                              const std::vector< Point<dim> >&          src);

} // FemExtras

#endif