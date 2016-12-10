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
// - Class: fe_vectors.h
// - Description: This class implements data type
//                used in function calls of applications
// - Developers: Guido Kanschat,     Texas A&M University
//               Valentin N. Zingan, University of Alberta
//               Marc Secanell,      University of Alberta
// - Id: $Id: fe_vectors.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _FUEL_CELL_APPLICATION_CORE_FE_VECTORS_H_
#define _FUEL_CELL_APPLICATION_CORE_FE_VECTORS_H_

#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/block_vector.h>

#include <algorithm>
#include <fstream>

using namespace dealii;

namespace FuelCell
{
namespace ApplicationCore
{

/**
 * The vector class used by applications.
 */
typedef BlockVector<double> FEVector;

/**
 * The data type used in function calls of Application.
 *
 * This class is a collection of pointers to BlockVector objects
 * representing finite element functions and additional parameters
 * that ought to be handed down to applications.
 *
 * Any class inheriting Application using its own data should
 * create a new FEVectors object and add its own data first, then
 * merge() with the data handed down from the enclosing Application.
 *
 * This policy ensures that the first vectors available at a certain
 * point are always the ones of the next enclosing Application.
 *
 * @author Guido Kanschat
 */

class FEVectors : public Subscriptor
{
public:

  /**
   * Constructor.
   */
  FEVectors();

  ///@name Adding members.
  //@{

  /**
   * Add a new scalar to the end
   * of the collection.
   */
  void add_scalar(double&            s,
                  const std::string& name);

  /**
   * Add a new vector to the end
   * of the collection.
   */
  void add_vector(FEVector&          v,
                  const std::string& name);

  /**
   * Add a new constant vector.
   */
  void add_vector(const FEVector&    v,
                  const std::string& name);

  /**
   * Merge the data of another
   * FEVectors to the end of this
   * object.
   */
  void merge(FEVectors& other);

  /**
   * Merge the data of another
   * FEVectors to the end of this
   * object.
   *
   * After this operation, all
   * data in this object will be
   * treated as const.
   */
  void merge(const FEVectors& other);

  //@}

  ///@name Accessors
  //@{

  /**
   * Number of stored scalars.
   */
  unsigned int n_scalars() const;

  /**
   * Number of stored vectors.
   */
  unsigned int n_vectors() const;

  /**
   * Access to a scalar stored.
   */
  double& scalar(unsigned int i);

  /**
   * Read-only access to a scalar stored.
   */
  const double& scalar(unsigned int i) const;

  /**
   * Access to a vector stored.
   */
  FEVector& vector(unsigned int i);

  /**
   * Read-only access to a vector stored.
   */
  const FEVector& vector(unsigned int i) const;

  /**
   * Name of a scalar.
   */
  const std::string& scalar_name(unsigned int i) const;

  /**
   * Name of a vector.
   */
  const std::string& vector_name(unsigned int i) const;

  /**
   * Find index of a named scalar.
   */
  unsigned int find_scalar(const std::string& name) const;

  /**
   * Find index of a named vector.
   */
  unsigned int find_vector(const std::string& name) const;

  /**
   * Counts the number of instance a vector by name.
   */
  unsigned int count_vector(const std::string& name) const;

  /**
   * List names of both stored scalars and vectors.
   * verbosity = 0 - nothing printed
   * verbosity = 1 - names pronted
   * verbosity > 1 - numbers, names, values, and L2-norms printed
   */
  template<typename OUT>
  void list(OUT&         out,
            unsigned int verbosity = 1) const;

  //@}

  /**
   * Exception indicating that a
   * function expected a vector
   * to have a certain name in this position, but
   * FEVectors had a different
   * name in that position.
   */
  DeclException2(ExcNameMismatch,
                 int,
                 std::string,
                 << "Name at position " << arg1 << " is not equal to " << arg2);

private:

  /**
   * True if the object is to be treated constant.
   */
  bool is_constant;

  /**
   * The scalars stored.
   */
  std::vector<double*> scalars;

  /**
   * The names of the scalars.
   */
  std::vector<std::string> scalar_names;

  /**
   * The vectors stored.
   */
  std::vector<FEVector*> vectors;

  /**
   * The names of the vectors.
   */
  std::vector<std::string> vector_names;
};

} // ApplicationCore

} // FuelCell

#endif