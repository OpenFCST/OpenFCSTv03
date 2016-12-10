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
// - Class: fe_vectors.cc
// - Description: This class implements data type
//                used in function calls of applications
// - Developers: Guido Kanschat,     Texas A&M University
//               Valentin N. Zingan, University of Alberta
//               Marc Secanell,      University of Alberta
// - Id: $Id: fe_vectors.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <application_core/fe_vectors.h>

namespace NAME = FuelCell::ApplicationCore;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NAME::FEVectors::FEVectors()
:
is_constant(false)
{ }

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::FEVectors::add_scalar(double&            s,
                            const std::string& name)
{
  Assert(!is_constant, ExcNotImplemented());
  scalars.push_back(&s);
  scalar_names.push_back(name);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::FEVectors::add_vector(NAME::FEVector&    v,
                            const std::string& name)
{
  Assert(!is_constant, ExcNotImplemented());
  vectors.push_back(&v);
  vector_names.push_back(name);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::FEVectors::add_vector(const NAME::FEVector& v,
                            const std::string&    name)
{
  Assert(!is_constant, ExcNotImplemented());
  vectors.push_back(const_cast<NAME::FEVector*> (&v));
  vector_names.push_back(name);
//   is_constant = true;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned int
NAME::FEVectors::count_vector(const std::string& name) const{
    return std::count(vector_names.begin(),vector_names.end(), name);
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::FEVectors::merge(FEVectors& other)
{
  bool constness = is_constant;

  for(unsigned int i = 0; i < other.n_scalars(); ++i)
    add_scalar(*other.scalars[i], other.scalar_names[i]);

  for(unsigned int i = 0; i < other.n_vectors(); ++i)
    add_vector(*other.vectors[i], other.vector_names[i]);

  is_constant = constness | other.is_constant;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::FEVectors::merge(const FEVectors& other)
{
  for(unsigned int i = 0; i < other.n_scalars(); ++i)
    add_scalar(*other.scalars[i], other.scalar_names[i]);

  for(unsigned int i = 0; i < other.n_vectors(); ++i)
    add_vector(*other.vectors[i], other.vector_names[i]);

//   is_constant = true;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned int
NAME::FEVectors::n_scalars() const
{
  return scalars.size();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned int
NAME::FEVectors::n_vectors() const
{
  return vectors.size();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double&
NAME::FEVectors::scalar(unsigned int i)
{
  Assert(!is_constant, ExcNotImplemented());
  Assert(i < n_scalars(), ExcIndexRange(i,0,n_scalars()));
  return *scalars[i];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const double&
NAME::FEVectors::scalar(unsigned int i) const
{
  Assert(i < n_scalars(), ExcIndexRange(i,0,n_scalars()));
  return *scalars[i];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NAME::FEVector&
NAME::FEVectors::vector(unsigned int i)
{
  Assert(!is_constant, ExcNotImplemented());
  Assert(i < n_vectors(), ExcIndexRange(i,0,n_vectors()));
  return *vectors[i];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const NAME::FEVector&
NAME::FEVectors::vector(unsigned int i) const
{
  Assert(i < n_vectors(), ExcIndexRange(i,0,n_vectors()));
  return *vectors[i];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const std::string&
NAME::FEVectors::scalar_name(unsigned int i) const
{
  Assert(i < n_scalars(), ExcIndexRange(i,0,n_scalars()));
  return scalar_names[i];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const std::string&
NAME::FEVectors::vector_name(unsigned int i) const
{
  Assert(i < n_vectors(), ExcIndexRange(i,0,n_vectors()));
  return vector_names[i];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned int
NAME::FEVectors::find_scalar(const std::string& name) const
{
  const std::vector<std::string>::const_iterator
    iter = std::find(scalar_names.begin(), scalar_names.end(), name);

  if(iter == scalar_names.end())
    return static_cast<unsigned int>(-1);

  return iter - scalar_names.begin();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned int
NAME::FEVectors::find_vector(const std::string& name) const
{
  const std::vector<std::string>::const_iterator
    iter = std::find(vector_names.begin(), vector_names.end(), name);

  if(iter == vector_names.end())
    return static_cast<unsigned int>(-1);

  return iter - vector_names.begin();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename OUT>
void
NAME::FEVectors::list(OUT&         out,
                      unsigned int verbosity) const
{
  if(verbosity == 0)
    return;
  if(verbosity == 1)
  {
        out << std::endl << "Scalars:";
        for(unsigned int i = 0; i < n_scalars(); ++i)
          out << ' ' << '\"' << scalar_names[i] << '\"';

        out << "Vectors:";
        for(unsigned int i = 0; i < n_vectors(); ++i)
          out << ' ' << '\"' << vector_names[i] << '\"';
        out << std::endl;
  }
  else
  {
        for(unsigned int i = 0; i < n_scalars(); ++i)
        {
          out << "Scalar-" << i << ' ' << scalar_names[i] << '=' << *scalars[i];
          out << std::endl;
        }
        for(unsigned int i = 0; i < n_vectors(); ++i)
        {
          out << "Vector-" << i << ' ' << vector_names[i];
          for(unsigned int b = 0; b < vectors[i]->n_blocks(); ++b)
            out << ' ' << vectors[i]->block(b).size();
          out << " Norm: " << vectors[i]->l2_norm();
          out << std::endl;
        }
  }
}
