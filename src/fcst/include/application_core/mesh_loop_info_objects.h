// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2009 by Guido Kanschat
// Copyright (C) 2009-2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: mesh_loop_info_objects.h
// - Description: Objects used for looping over mesh
// - Developers: Guido Kanschat,     Texas A&M University
//               Valentin N. Zingan, University of Alberta
//
// ----------------------------------------------------------------------------

#ifndef _FUEL_CELL_APPLICATION_CORE_MESH_LOOP_INFO_OBJECTS_H_
#define _FUEL_CELL_APPLICATION_CORE_MESH_LOOP_INFO_OBJECTS_H_

//-- deal.II
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/vector_slice.h>
#include <deal.II/lac/block_indices.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_tools.h>

//-- OpenFCST
#include <application_core/fe_vectors.h>

//-- boost libraries:
#include <boost/shared_ptr.hpp>

//----
namespace dealii
{
       template<int,int> class DoFHandler;
}

using namespace dealii;

/**
 * A collection of functions and classes for the mesh loops that are
 * an ubiquitous part of each finite element program.
 *
 * The collection implements:
 *
 * - @p BlockInfo,
 * - @p DoFInfo,
 * - @p IntegrationInfo.
 *
 * They can be considered as an extension of mesh iterators
 * in the sense that they are updated for each mesh entity
 * and provide additional data on them.
 *
 * @author Guido Kanschat
 * @author Valentin N. Zingan
 */

namespace FuelCell
{
namespace ApplicationCore
{

/**
 * Helper functions computing the desired data in each quadrature point of a mesh entity by calling
 * @p FEValuesBase::get_function_values(), @p FEValuesBase::get_function_grads(),
 * and @p FEValuesBase::get_function_hessians().
 *
 * @param fe_values: The @p FEValues object.
 *
 * @param fe_vector: The global finite element function in the form of a global
 * nodal @p FEVector.
 *
 * @param local_dof_indices: The local DoF indices associated with the current mesh entity.
 *
 * @param first_index: The first index in @p local_dof_indices to be used.
 *
 * @param n_indices: The number of indices in @p local_dof_indices to be used.
 *
 * @param result: The result.
 * 
 * The access to @p ith component in @p qth quadrature point is @p result[i][q].
 *
 * @author Guido Kanschat
 */

// +++ IMPLEMENTATION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, typename TYPE>
inline
void
fill_data(const FEValuesBase<dim>&          fe_values,
          const FEVector&                   fe_vector,
          const std::vector<unsigned int>&  local_dof_indices,
          unsigned int                      first_index,
          unsigned int                      n_indices,
          std::vector< std::vector<TYPE> >& result)
{
       Assert(false, ExcNotImplemented());
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
inline
void
fill_data(const FEValuesBase<dim>&            fe_values,
          const FEVector&                     fe_vector,
          const std::vector<unsigned int>&    local_dof_indices,
          unsigned int                        first_index,
          unsigned int                        n_indices,
          std::vector< std::vector<double> >& result)
{
    
    VectorSlice<std::vector< std::vector<double> > > aux(result);
    fe_values.get_function_values(fe_vector,
                                  make_slice(local_dof_indices, first_index, n_indices),
                                  aux,
                                  true);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
inline
void
fill_data(const FEValuesBase<dim>&                     fe_values,
          const FEVector&                              fe_vector,
          const std::vector<unsigned int>&             local_dof_indices,
          unsigned int                                 first_index,
          unsigned int                                 n_indices,
          std::vector< std::vector< Tensor<1,dim> > >& result)
{
    VectorSlice< std::vector< std::vector< Tensor<1,dim> > > > aux(result);
    fe_values.get_function_gradients(fe_vector,
                                     make_slice(local_dof_indices, first_index, n_indices),
                                     aux,
                                     true);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
inline
void
fill_data(const FEValuesBase<dim>&                     fe_values,
          const FEVector&                              fe_vector,
          const std::vector<unsigned int>&             local_dof_indices,
          unsigned int                                 first_index,
          unsigned int                                 n_indices,
          std::vector< std::vector< Tensor<2,dim> > >& result)
{
       fe_values.get_function_hessians(fe_vector,
                                       make_slice(local_dof_indices, first_index, n_indices),
                                       result,
                                       true);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/**
 * A small structure collecting the different @p BlockIndices
 * of @p FEVector vectors (for instance, solution)
 * involved in the computations.
 *
 * @author Guido Kanschat
 */

struct BlockInfo : public Subscriptor
{
    /**
     * The block structure of the global @p FEVector solution.
     */
    BlockIndices global;
    
    /**
     * The block structure of a local @p FEVector solution.
     */
    BlockIndices local;
    
    /**
     * The multilevel block structure of the global @p FEVector solution.
     */
    std::vector<BlockIndices> levels;
    
    /**
     * A vector of base elements.
     *
     * @note index -> value = # block -> # fe in @p FESystem object.
     */
    std::vector<unsigned int> base_element;
    
    /**
     * A vector containing the internal renumbering of degrees of freedom from the
     * standard ordering on a cell to a block wise ordering on the same cell.
     *
     * @note index -> value = standard # dof -> block wise # dof.
     * @note if @p local_renumbering.size() = 0: the
     * standard internal renumbering of degrees of freedom
     * is activated.
     */
    std::vector<unsigned int> local_renumbering;
    
    /**
     * Fill the structure globally with values describing the @p DoFHandler.
     */
    template<int dim, int spacedim>
    void initialize(const DoFHandler<dim, spacedim>& dof_handler)
    {
        const FiniteElement<dim, spacedim>& fe = dof_handler.get_fe();
        std::vector<unsigned int> sizes(fe.n_blocks());
        DoFTools::count_dofs_per_block(dof_handler, sizes);
        global.reinit(sizes);
    }
    
    /**
     * Fill the structure locally with values describing the @p DoFHandler.
     */
    template<int dim, int spacedim>
    void initialize_local(const DoFHandler<dim, spacedim>& dof_handler)
    {
        const FiniteElement<dim, spacedim>& fe = dof_handler.get_fe();
        std::vector<unsigned int> sizes(fe.n_blocks());
        
        base_element.resize(fe.n_blocks());
        
        for(unsigned int i = 0; i < base_element.size(); ++i)
            base_element[i] = fe.block_to_base_index(i).first;
        
        local_renumbering.resize(fe.n_dofs_per_cell());
        FETools::compute_block_renumbering(fe,
                                           local_renumbering,
                                           sizes,
                                           false);
        local.reinit(sizes);
    }   
};


/**
 * Very basic info class only containing information on geometry and
 * degrees of freedom on a mesh entity.
 *
 * The information in this class is usually used by mesh loops.
 *
 * This class operates in two different modes:
 *
 * - @p MODE1: block wise local renumbering mode is activated if @p BlockInfo::local_renumbering.size() > 0
 *   by calling both @p BlockInfo::initialize() and @p BlockInfo::initialize_local() functions.
 *
 * - @p MODE2: standard deal.ii local renumbering mode is activated if @p BlockInfo::local_renumbering.size() = 0
 *   by only calling @p BlockInfo::initialize() function.
 *
 * @note @p MODE1 is used in @p FCST software by default.
 *
 * The @p BlockInfo object is stored here as a pointer. Therefore, if the
 * local block structure changes, for instance because of the mesh refinement,
 * the @p DoFInfo class will automatically use this new structure.
 *
 * @author Guido Kanschat
 * @author Valentin N. Zingan
 */
template<int dim, int spacedim = dim>
class DoFInfo
{
public:

///@name Constructors, destructor, and initialization
//@{

  /**
   * Constructor.
   */
  DoFInfo(const BlockInfo& block_info);

  /**
   * @deprecated Constructor.
   * The first argument is ignored.
   */
  DoFInfo(const FEVectors&,
          const BlockInfo& block_info);

//@}

///@name Reinitialization
//@{

  /**
   * Set the current cell and
   * fill @p indices.
   */
  template<typename DHCellIterator>
  void reinit(const DHCellIterator& c);

  /**
   * Set the current cell and face and
   * fill @p indices.
   */
  template<typename DHCellIterator, typename DHFaceIterator>
  void reinit(const DHCellIterator& c,
              const DHFaceIterator& f,
              const unsigned int    fn);

  /**
   * Set the current cell, face, and subface and
   * fill @p indices.
   */
  template<typename DHCellIterator, typename DHFaceIterator>
  void reinit(const DHCellIterator& c,
              const DHFaceIterator& f,
              const unsigned int    fn,
              const unsigned int    sn);

//@}

///@name Mesh iterators
//@{

  /**
   * The current @p DoFHandler<dim> cell.
   */
  typename DoFHandler<dim>::cell_iterator dof_cell;

  /**
   * The current @p DoFHandler<dim> active cell.
   */
  typename DoFHandler<dim>::active_cell_iterator dof_active_cell;

  /**
   * The current @p DoFHandler<dim> face.
   */
  typename DoFHandler<dim>::face_iterator dof_face;

  /**
   * The current @p Triangulation<dim> cell.
   */
  typename Triangulation<dim>::cell_iterator cell;

  /**
   * The current @p Triangulation<dim> face.
   */
  typename Triangulation<dim>::face_iterator face;

  /**
   * The number of the current
   * face on the current cell.
   *
   * This number is
   * @p static_cast<unsigned int>(-1)
   * if the @p DoFInfo object is
   * initialized with a cell only.
   */
  unsigned int face_number;

  /**
   * The number of the current
   * subface on the current face
   * of the current cell.
   *
   * This number is
   * @p static_cast<unsigned int>(-1)
   * if the @p DoFInfo object is
   * initialized with a cell and face only.
   */
  unsigned int sub_number;

//@}

///@name Other data
//@{

  /**
   * The local dof indices associated
   * with the current cell.
   */
  std::vector<unsigned int> indices;

  /**
   * The block structure of the problem at hand.
   */
  SmartPointer<const BlockInfo> block_info;

//@}

private:

  /**
   * Fill @p indices.
   */
  void get_indices(const typename DoFHandler<dim, spacedim>::cell_iterator c);

  /**
   * Auxiliary vector.
   */
  std::vector<unsigned int> indices_org;
};

// +++ IMPLEMENTATION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, int spacedim>
inline
DoFInfo<dim, spacedim>::DoFInfo(const BlockInfo& block_info)
:
block_info(&block_info,
           typeid(*this).name())
{ }

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, int spacedim>
inline
DoFInfo<dim, spacedim>::DoFInfo(const FEVectors&,
                                const BlockInfo& block_info)
:
block_info(&block_info,
           typeid(*this).name())
{ }

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, int spacedim>
template<typename DHCellIterator>
inline
void
DoFInfo<dim, spacedim>::reinit(const DHCellIterator& c)
{
       get_indices(c);

       dof_cell        = c;
       dof_active_cell = c;
       cell            = static_cast<typename Triangulation<dim>::cell_iterator> (c);

       face_number = static_cast<unsigned int>(-1);
       sub_number  = static_cast<unsigned int>(-1);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, int spacedim>
template<typename DHCellIterator, typename DHFaceIterator>
inline
void
DoFInfo<dim, spacedim>::reinit(const DHCellIterator& c,
                               const DHFaceIterator& f,
                               const unsigned int    fn)
{
       if(   cell.state() != IteratorState::valid || cell != static_cast<typename Triangulation<dim>::cell_iterator> (c)   )
       {
         get_indices(c);
       }

       dof_cell        = c;
       dof_active_cell = c;
       cell            = static_cast<typename Triangulation<dim>::cell_iterator> (c);

       dof_face = f;
       face     = static_cast<typename Triangulation<dim>::face_iterator> (f);

       face_number = fn;
       sub_number  = static_cast<unsigned int>(-1);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, int spacedim>
template<typename DHCellIterator, typename DHFaceIterator>
inline
void
DoFInfo<dim, spacedim>::reinit(const DHCellIterator& c,
                               const DHFaceIterator& f,
                               const unsigned int    fn,
                               const unsigned int    sn)
{
       if(   cell.state() != IteratorState::valid || cell != static_cast<typename Triangulation<dim>::cell_iterator> (c)   )
       {
         get_indices(c);
       }

       dof_cell        = c;
       dof_active_cell = c;
       cell            = static_cast<typename Triangulation<dim>::cell_iterator> (c);

       dof_face = f;
       face     = static_cast<typename Triangulation<dim>::face_iterator> (f);

       face_number = fn;
       sub_number  = sn;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, int spacedim>
inline
void
DoFInfo<dim, spacedim>::get_indices(const typename DoFHandler<dim, spacedim>::cell_iterator c)
{
       indices.resize(c->get_fe().dofs_per_cell);

       if(block_info->local_renumbering.size() == 0)
       {
         c->get_dof_indices(indices);
       }
       else
       {
         indices_org.resize(c->get_fe().dofs_per_cell);
         c->get_dof_indices(indices_org);

         for(unsigned int i = 0; i < indices.size(); ++i)
           indices[this->block_info->local_renumbering[i]] = indices_org[i];
       }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/**
 * This class is created for the objects handed to the mesh loops.
 *
 * The objects of this class contain an object of type @p FEVALUESBASE
 * which is
 * either @p FEValuesBase     (base for cell)
 * or     @p FEFaceValuesBase (base for face and subface).
 *
 * At the same time,
 * the actual type @p FEVALUES for @p FEValuesBase is @p FEValues (cell)
 * and
 * the actual type @p FEVALUES for @p FEFaceValuesBase is either @p FEFaceValues (face) or @p FESubfaceValues (subface).
 *
 * @note @p FEVALUESBASE is a template parameter of the class itself.
 * @note @p FEVALUES     is a template parameter of the @p initialize() function.
 *
 * Additionally, this class includes containers to store the @p values,
 * @p gradients and @p hessians of finite element functions
 * stored in @p global_data in the quadrature points of a mesh entity.
 *
 * @author Guido Kanschat
 * @author Valentin N. Zingan
 */

template<int dim, typename FEVALUESBASE>
class IntegrationInfo : public DoFInfo<dim>
{
public:

///@name Constructors, destructor, and initialization
//@{

  /**
   * Constructor.
   */
  IntegrationInfo(const BlockInfo& block_info);

  /**
   * Constructor.
   */
  IntegrationInfo(const FEVectors& data,
                  const BlockInfo& block_info);

  /**
   * Build internal structures @p fevalv and allocate memory for
   * @p values,
   * @p gradients,
   * @p hessians.
   *
   * @param FE_VALUES is an object of the actual type @p FEVALUES
   * which is either @p FEValues (cell) or @p FEFaceValues (face)
   * or @p FESubfaceValues (subface).
   * @note @p FE_VALUES is really not needed here and only serves to
   * define the actual type @p FEVALUES.
   *
   * @param fe is the finite element (for instance, @p FESystem object)
   * stored in the @p DoFHandler. It is used in the constructor of @p FE_VALUES.
   *
   * @param mapping is the @p Mapping object used to map the mesh cells.
   * It is used in the constructor of @p FE_VALUES.
   *
   * @param quadrature is a quadrature formula. It is used in the constructor
   * of @p FE_VALUES.
   *
   * @param flags are the @p UpdateFlags. They are used in the constructor
   * of @p FE_VALUES.
   */
  template<typename FEVALUES>
  void initialize(const FEVALUES*                                 FE_VALUES,
                  const FiniteElement<dim>&                       fe,
                  const Mapping<dim>&                             mapping,
                  const Quadrature<FEVALUES::integral_dimension>& quadrature,
                  const UpdateFlags                               flags);

  /**
   * Initialize @p global_data.
   */
  void initialize_data(const FEVectors& data);

//@}

///@name Reinitialization
//@{

  /**
   * Reinitialize internal data for use on a cell.
   */
  template<typename DHCellIterator>
  void reinit(const DHCellIterator& c);

  /**
   * Reinitialize internal data for use on a face of a cell.
   */
  template<typename DHCellIterator, typename DHFaceIterator>
  void reinit(const DHCellIterator& c,
              const DHFaceIterator& f,
              const unsigned int    fn);

  /**
   * Reinitialize internal data for use on a subface of a face of a cell.
   */
  template<typename DHCellIterator, typename DHFaceIterator>
  void reinit(const DHCellIterator& c,
              const DHFaceIterator& f,
              const unsigned int    fn,
              const unsigned int    sn);

  /**
   * Use the finite element functions in @p global_data and fill the vectors
   * @p values,
   * @p gradients,
   * @p hessians.
   */
  template<typename TYPE>
  void fill_local_data(std::vector< std::vector< std::vector<TYPE> > >& data,
                       bool split_fevalues) const;

//@}

///@name Accessors and info
//@{

  /**
   * Access to a single actual @p FEVALUES object.
   */
  const FEVALUESBASE& fe() const;

  /**
   * Access to a group of actual @p FEVALUES objects.
   */
  const FEVALUESBASE& fe(unsigned int i) const;

  /**
   * Access to @p fe_val_unsplit.
   */
  const FEVALUESBASE& get_fe_val_unsplit() const;

  /**
   * Resize @p fevalv to 0.
   */
  void clear();

//@}

///@name DATA
//@{

  /**
   * @p true if we assemble for multigrid.
   */
  bool multigrid;

  /**
   * The smart pointer to the @p FEVectors object called @p global_data.
   */
  SmartPointer<const FEVectors> global_data;

  /**
   * The vector containing the values of finite element
   * functions in the quadrature points.
   *
   * Example: @p nth solution, @p mth component, value at the @p qth quadrature point = @p values[n][m][q]
   */
  std::vector< std::vector< std::vector<double> > > values;

  /**
   * The vector containing the gradients of finite element
   * functions in the quadrature points.
   *
   * Example: @p nth solution, @p mth component, gradient at the @p qth quadrature point = @p gradients[n][m][q]
   */
  std::vector< std::vector< std::vector< Tensor<1,dim> > > > gradients;

  /**
   * @deprecated Reference to gradients for compatibility.
   */
  std::vector< std::vector< std::vector< Tensor<1,dim> > > >& derivatives;

  /**
   * The vector containing the hessians of finite element
   * functions in the quadrature points.
   *
   * Example: @p nth solution, @p mth component, hessian at the @p qth quadrature point = @p hessians[n][m][q]
   */
  std::vector< std::vector< std::vector< Tensor<2,dim> > > > hessians;

//@}

private:

///@name DATA
//@{

  /**
   * The vector of smart pointers to @p FEVALUESBASE objects.
   */
  std::vector< boost::shared_ptr<FEVALUESBASE> > fevalv;

  /**
   * In the @p MODE1 (or in the main mode), the class splits
   * an @p FEVALUESBASE object into several sub-objects and fills @p fevalv.
   * However, for the variety of problems, there is still a need
   * to keep the unsplit @p FEVALUESBASE object as well.
   */
  boost::shared_ptr<FEVALUESBASE> fe_val_unsplit;

//@}
};

// +++ IMPLEMENTATION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, typename FEVALUESBASE>
inline
IntegrationInfo<dim, FEVALUESBASE>::IntegrationInfo(const BlockInfo& block_info)
:
DoFInfo<dim>(block_info),
multigrid(false),
global_data(0, typeid(*this).name()),
derivatives(gradients),
fevalv(0)
{ }

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, typename FEVALUESBASE>
inline
IntegrationInfo<dim, FEVALUESBASE>::IntegrationInfo(const FEVectors& data,
                                                    const BlockInfo& block_info)
:
DoFInfo<dim>(block_info),
multigrid(false),
global_data(&data),
derivatives(gradients),
fevalv(0)
{ }

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, typename FEVALUESBASE>
template<typename FEVALUES>
inline
void
IntegrationInfo<dim, FEVALUESBASE>::initialize(const FEVALUES*,
                                               const FiniteElement<dim>&                       fe,
                                               const Mapping<dim>&                             mapping,
                                               const Quadrature<FEVALUES::integral_dimension>& quadrature,
                                               const UpdateFlags                               flags)
{
       if( this->block_info->local_renumbering.size() != 0 ) // MODE1
       {
              fevalv.resize(fe.n_base_elements());
              for(unsigned int i = 0; i < fevalv.size(); ++i)
              {
                fevalv[i] = boost::shared_ptr<FEVALUESBASE>( new FEVALUES(mapping,
                                                                          fe.base_element(i),
                                                                          quadrature,
                                                                          flags)
                                                           );
              }

              // --- fe_val_unsplit initialization ---

              fe_val_unsplit = boost::shared_ptr<FEVALUESBASE>( new FEVALUES(mapping,
                                                                             fe,
                                                                             quadrature,
                                                                             flags)
                                                              );
       }
       else // MODE2
       {
              fevalv.resize(1);
              fevalv[0] = boost::shared_ptr<FEVALUESBASE>( new FEVALUES(mapping,
                                                                        fe,
                                                                        quadrature,
                                                                        flags)
                                                         );
       }

       values.resize(global_data->n_vectors(),
                     std::vector< std::vector<double> >(fe.n_components(),
                                                        std::vector<double>(quadrature.size())
                                                       )
                    );
       gradients.resize(global_data->n_vectors(),
                        std::vector< std::vector< Tensor<1,dim> > >(fe.n_components(),
                                                                    std::vector< Tensor<1,dim> >(quadrature.size())
                                                                   )
                       );
       hessians.resize(global_data->n_vectors(),
                       std::vector< std::vector< Tensor<2,dim> > >(fe.n_components(),
                                                                   std::vector< Tensor<2,dim> >(quadrature.size())
                                                                  )
                      );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, typename FEVALUESBASE>
inline
void
IntegrationInfo<dim, FEVALUESBASE>::initialize_data(const FEVectors& data)
{
       global_data = &data;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, typename FEVALUESBASE>
template<typename DHCellIterator>
inline
void
IntegrationInfo<dim, FEVALUESBASE>::reinit(const DHCellIterator& c)
{
       DoFInfo<dim>::reinit(c);

       for(unsigned int i = 0; i < fevalv.size(); ++i)
       {
         FEVALUESBASE&  fe_values_base = *fevalv[i];
         FEValues<dim>& fe_values      = dynamic_cast<FEValues<dim>&>(fe_values_base);
         fe_values.reinit(this->cell);
       }

       FEVALUESBASE&  fe_values_base_unsplit = *fe_val_unsplit;
       FEValues<dim>& fe_values_unsplit      = dynamic_cast<FEValues<dim>&>(fe_values_base_unsplit);
       fe_values_unsplit.reinit(this->dof_active_cell);

       const bool split_fevalues = this->block_info->local_renumbering.size() != 0;
       fill_local_data(values,
                       split_fevalues);
       fill_local_data(gradients,
                       split_fevalues);
     //fill_local_data(hessians,
     //                split_fevalues);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, typename FEVALUESBASE>
template<typename DHCellIterator, typename DHFaceIterator>
inline
void
IntegrationInfo<dim, FEVALUESBASE>::reinit(const DHCellIterator& c,
                                           const DHFaceIterator& f,
                                           const unsigned int    fn)
{
       DoFInfo<dim>::reinit(c,
                            f,
                            fn);

       for(unsigned int i = 0; i < fevalv.size(); ++i)
       {
         FEVALUESBASE&      fe_values_base = *fevalv[i];
         FEFaceValues<dim>& fe_face_values = dynamic_cast<FEFaceValues<dim>&>(fe_values_base);
         fe_face_values.reinit(this->cell,
                               fn);
       }

       FEVALUESBASE&      fe_values_base_unsplit = *fe_val_unsplit;
       FEFaceValues<dim>& fe_face_values_unsplit = dynamic_cast<FEFaceValues<dim>&>(fe_values_base_unsplit);
       fe_face_values_unsplit.reinit(this->dof_active_cell,
                                     fn);

       const bool split_fevalues = this->block_info->local_renumbering.size() != 0;
       fill_local_data(values,
                       split_fevalues);
       fill_local_data(gradients,
                       split_fevalues);
     //fill_local_data(hessians,
     //                split_fevalues);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, typename FEVALUESBASE>
template<typename DHCellIterator, typename DHFaceIterator>
inline
void
IntegrationInfo<dim, FEVALUESBASE>::reinit(const DHCellIterator& c,
                                           const DHFaceIterator& f,
                                           const unsigned int    fn,
                                           const unsigned int    sn)
{
       DoFInfo<dim>::reinit(c, f, fn, sn);

       for(unsigned int i = 0; i < fevalv.size(); ++i)
       {
         FEVALUESBASE&         fe_values_base    = *fevalv[i];
         FESubfaceValues<dim>& fe_subface_values = dynamic_cast<FESubfaceValues<dim>&>(fe_values_base);
         fe_subface_values.reinit(this->cell,
                                  fn,
                                  sn);
       }

       FEVALUESBASE&         fe_values_base_unsplit    = *fe_val_unsplit;
       FESubfaceValues<dim>& fe_subface_values_unsplit = dynamic_cast<FESubfaceValues<dim>&>(fe_values_base_unsplit);
       fe_subface_values_unsplit.reinit(this->dof_active_cell,
                                        fn,
                                        sn);

       const bool split_fevalues = this->block_info->local_renumbering.size() != 0;
       fill_local_data(values,
                       split_fevalues);
       fill_local_data(gradients,
                       split_fevalues);
     //fill_local_data(hessians,
     //                split_fevalues);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, typename FEVALUESBASE>
template<typename TYPE>
inline
void
IntegrationInfo<dim, FEVALUESBASE>::fill_local_data(std::vector< std::vector< std::vector<TYPE> > >& data,
                                                    bool                                             split_fevalues) const
{
       if(split_fevalues)
       {
              std::vector< std::vector<TYPE> > local_data;

              unsigned int comp = 0;
              for(unsigned int b = 0; b < this->block_info->base_element.size(); ++b)
              {
                     const unsigned int       no_fe          = this->block_info->base_element[b];
                     const FEValuesBase<dim>& fe_values_base = this->fe(no_fe);

                     const unsigned int n_comp = fe_values_base.get_fe().n_components();
                     local_data.resize(n_comp,
                                       std::vector<TYPE>(fe_values_base.n_quadrature_points));

                     for(unsigned int i = 0; i < this->global_data->n_vectors(); ++i)
                     {
                       const FEVector& src = this->global_data->vector(i);
                       fill_data(fe_values_base,
                                 src,
                                 this->indices,
                                 this->block_info->local.block_start(b),
                                 fe_values_base.dofs_per_cell,
                                 local_data);

                       for(unsigned int c = 0; c < local_data.size(); ++c)
                         for(unsigned int k = 0; k < local_data[c].size(); ++k)
                           data[i][comp+c][k] = local_data[c][k];
                     }

                     comp += n_comp;
              }
       }
       else
       {
              for(unsigned int i = 0; i < this->global_data->n_vectors(); ++i)
              {
                     const FEVector&          src            = this->global_data->vector(i);
                     const FEValuesBase<dim>& fe_values_base = this->fe();

                     fill_data(fe_values_base,
                               src,
                               this->indices,
                               0,
                               fe_values_base.get_fe().dofs_per_cell,
                               data[i]);
              }
       }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, typename FEVALUESBASE>
inline
const FEVALUESBASE&
IntegrationInfo<dim, FEVALUESBASE>::fe() const
{
       AssertDimension(fevalv.size(), 1);
       return *fevalv[0];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, typename FEVALUESBASE>
inline
const FEVALUESBASE&
IntegrationInfo<dim, FEVALUESBASE>::fe(unsigned int i) const
{
       Assert( i < fevalv.size(), ExcIndexRange(i, 0, fevalv.size()) );
       return *fevalv[i];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, typename FEVALUESBASE>
inline
const FEVALUESBASE&
IntegrationInfo<dim, FEVALUESBASE>::get_fe_val_unsplit() const
{
       return *fe_val_unsplit;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim, typename FEVALUESBASE>
inline
void
IntegrationInfo<dim, FEVALUESBASE>::clear()
{
       fevalv.clear();
       fevalv.resize(0);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

} // ApplicationCore

} // FuelCell

#endif