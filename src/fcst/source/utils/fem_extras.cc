// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: fem_extras.cc
// - Description: This namespace contains FEM related methods
// - Developers: Valentin N. Zingan, University of Alberta
//
// ----------------------------------------------------------------------------

#include "utils/fem_extras.h"

// ---                        ---
// --- get_tangential_vectors ---
// ---                        ---

template<int dim>
void
FemExtras::get_tangential_vectors(std::vector< std::vector< Point<dim> > >& dst,
                                  const std::vector< Point<dim> >&          src)
{
  AssertThrow( dim != 1 , ExcInternalError() );
  AssertThrow( dst.size() == dim - 1 , ExcDimensionMismatch(dst.size(), dim - 1) );

  for(unsigned int alpha = 0; alpha < dim - 1; ++alpha)
  {
         AssertThrow( dst[alpha].size() == src.size() , ExcDimensionMismatch(dst[alpha].size(), src.size()) );
  }

  #ifdef _2D_

  for(unsigned int q = 0; q < src.size(); ++q)
  {
         dst[0][q](0) = - src[q](1);
         dst[0][q](1) =   src[q](0);
  }

  #endif

  #ifdef _3D_

  for(unsigned int q = 0; q < src.size(); ++q)
  {
         double tau11;
         double tau12;
         double tau13;

         double tau21;
         double tau22;
         double tau23;

         const double n1 = src[q](0);
         const double n2 = src[q](1);
         const double n3 = src[q](2);

         if( std::fabs(n1) >= 0.5 || std::fabs(n2) >= 0.5 )
         {
                const double norm12 = std::sqrt(n1*n1 + n2*n2);

                tau11 =   n2 / norm12;
                tau12 = - n1 / norm12;
                tau13 =   0.0;

                tau21 = - tau12 * n3;
                tau22 =   tau11 * n3;
                tau23 =   tau12 * n1 - tau11 * n2;
         }
         else
         {
                const double norm23 = std::sqrt(n2*n2 + n3*n3);

                tau11 =   0.0;
                tau12 = - n3 / norm23;
                tau13 =   n2 / norm23;

                tau21 =   tau13 * n2 - tau12 * n3;
                tau22 = - tau13 * n1;
                tau23 =   tau12 * n1;
         }

         dst[0][q](0) = tau11;
         dst[0][q](1) = tau12;
         dst[0][q](2) = tau13;

         dst[1][q](0) = tau21;
         dst[1][q](1) = tau22;
         dst[1][q](2) = tau23;
  }

  #endif
}

       /////////////////////////////
       /////////////////////////////
       // EXPLICIT INSTANTIATIONS //
       /////////////////////////////
       /////////////////////////////

// ---                        ---
// --- get_tangential_vectors ---
// ---                        ---

template void FemExtras::get_tangential_vectors<deal_II_dimension>(std::vector< std::vector< Point<deal_II_dimension> > >&, const std::vector< Point<deal_II_dimension> >&);