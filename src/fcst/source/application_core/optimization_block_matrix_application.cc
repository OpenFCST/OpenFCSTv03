// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: optimization_block_matrix_application.cc
// - Description: This is a base class for all applications that provide optimization information
// - Developers: Marc Secanell, University of Alberta
// - $Id: optimization_block_matrix_application.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <application_core/optimization_block_matrix_application.h>

using namespace FuelCell::ApplicationCore;

template <int dim>
const FuelCell::ApplicationCore::Event OptimizationBlockMatrixApplication<dim>::sensitivity_analysis = FuelCell::ApplicationCore::Event::assign("Sensitivity Analysis");


//---------------------------------------------------------------------------
//---------------- OPTIMIZATIONBLOCKMATRIXAPPLICATION -----------------------
//---------------------------------------------------------------------------
template <int dim>
OptimizationBlockMatrixApplication<dim>::OptimizationBlockMatrixApplication(boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data)
:
FuelCell::ApplicationCore::BlockMatrixApplication<dim>(data)
{
    FcstUtilities::log << "->OptimizationBlockMatrix"<<std::endl;
    optimization = false;
    boundary_responses = false;
    this->set_all_response_names();
}

//---------------------------------------------------------------------------
template <int dim>
OptimizationBlockMatrixApplication<dim>::OptimizationBlockMatrixApplication(FuelCell::ApplicationCore::DoFApplication<dim>& other,
                                                                            bool triangulation_only)
:
FuelCell::ApplicationCore::BlockMatrixApplication<dim>(other, triangulation_only)
{
    optimization = false;
    boundary_responses = false;
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::declare_parameters(ParameterHandler& param)
{
    FuelCell::ApplicationCore::BlockMatrixApplication<dim>::declare_parameters(param);

    param.enter_subsection("Output Variables");
    {
        param.declare_entry("Compute boundary responses",
                            "false",
                             Patterns::Bool(),
                             "This flag determines if boundary responses will be computed or not");
        param.declare_entry("num_output_vars",
                            "1",
                            Patterns::Integer(),
                            "Number of output variables to be computed");
        param.declare_entry("Output boundary id",
                            "2",
                            Patterns::List( Patterns::Integer(0) ),
                            "Boundary_id(s) on which the boundary response is computed. "
                            "This is only used for boundary responses. "
                            "Provide a comma-separated list of unsigned integers.");
                            
                            

        std::string joinedString = "no_name|"+boost::algorithm::join(all_response_names, "|");
        
        const unsigned int num_output_vars = 40;
        for (unsigned int index=0; index<num_output_vars; ++index)
        {
            std::ostringstream streamOut;
            streamOut << index;
            std::string name = "Output_var_" + streamOut.str();
            param.declare_entry(name.c_str(),
                                "no_name",
                                Patterns::Selection(joinedString));
        }
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::initialize(ParameterHandler& param)
{
    param.enter_subsection("Output Variables");
    boundary_responses = param.get_bool("Compute boundary responses");
    // For solver:
    this->data->enter_flag("boundary_responses", boundary_responses);
    unsigned int num_output_vars = param.get_integer("num_output_vars");
    this->user_input_bdry = FcstUtilities::string_to_number<unsigned int>( Utilities::split_string_list( param.get("Output boundary id") ) );
    name_output_var.clear();
    name_output_var.resize(num_output_vars);

    unsigned int numOfVarWithOutBoundary = 0;
    for (unsigned int index=0; index<num_output_vars; ++index)
    {
        std::ostringstream streamOut;
        streamOut << index;
        std::string name = "Output_var_" + streamOut.str();
        name_output_var[index] = param.get(name.c_str());
        if(name_output_var[index] == "current" || name_output_var[index] == "cathode_current")
            numOfVarWithOutBoundary++;
    }
    param.leave_subsection();

    if(!optimization)
    {
        n_resp = (num_output_vars - numOfVarWithOutBoundary) * user_input_bdry.size() + numOfVarWithOutBoundary;
        name_responses.clear();
        name_responses.resize(n_resp);
        for(unsigned int i = 0, index = 0; i < num_output_vars; ++i) //Loop through name_output_var
        {
//             if(name_output_var[i] == "current" || name_output_var[i] == "cathode_current") //Need else polarization_curve.cc becomes broken and all tests fail
//             if(!boundary_responses)
//             {
//                 name_responses[index] = name_output_var[i];
//                 index++;
//             }
//             else
//             {
//                 for(unsigned int j = 0; j < user_input_bdry.size(); ++j) //Loop through the boundary IDs to calculate response
//                 {
//                     name_responses[index] = name_output_var[i] + " at boundary " + std::to_string(user_input_bdry[j]);
//                     index++;
//                 }
//             }
            
            if(name_output_var[i].compare(0,24,"total_mass_flux_species_") == 0)
            {
                for(unsigned int j = 0; j < user_input_bdry.size(); ++j) //Loop through the boundary IDs to calculate response
                {
                    name_responses[index] = name_output_var[i] + " at boundary " + std::to_string(user_input_bdry[j]);
                    index++;
                }
            }
            else
            {
                name_responses[index] = name_output_var[i];
                index++;
            }
        }
    }

    FuelCell::ApplicationCore::BlockMatrixApplication<dim>::initialize(param);
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::responses (std::vector<double>& dst,
                                                    const FuelCell::ApplicationCore::FEVectors& src)
{
  // -- clear and resize destination vector:
  dst.clear();
  dst.resize(n_resp);

  // --- assertions ---
  AssertThrow(this->quadrature_residual_cell.size() != 0, ExcNotInitialized());
  AssertThrow(this->quadrature_residual_bdry.size() != 0, ExcNotInitialized());

  // --- types of FEVALUES objects we will use further ---
  // --- note: we only need types ---
  FEValues<dim>*     fe_values      = 0;
  FEFaceValues<dim>* fe_face_values = 0;

  // --- info structures construction ---
  typename DoFApplication<dim>::CellInfo cell_info(src,
                                                   this->block_info);
  typename DoFApplication<dim>::FaceInfo bdry_info(src,
                                                   this->block_info);

  // --- info structures initialization ---
  cell_info.initialize(
                       fe_values,
                      *this->element,
                      *this->mapping,
                       this->quadrature_residual_cell,
                       UpdateFlags(update_q_points | update_values | update_gradients | update_JxW_values)
                      );

  bdry_info.initialize(
                       fe_face_values,
                      *this->element,
                      *this->mapping,
                       this->quadrature_residual_bdry,
                       UpdateFlags(update_q_points | update_values | update_gradients | update_normal_vectors | update_JxW_values)
                      );

  // --- loop ---
  this->tr->clear_user_flags();

  typename DoFHandler<dim>::active_cell_iterator
  cell = this->dof->begin_active(),
  endc = this->dof->end();

  for( ; cell != endc; ++cell)
  {
         cell_info.reinit(cell);

         cell_responses(dst,
                        cell_info,
                        src.vector(0));

         if( this->boundary_responses )
         {
                for(unsigned int no_face = 0; no_face < GeometryInfo<dim>::faces_per_cell; ++no_face)
                {
                       typename DoFHandler<dim>::face_iterator face = cell->face(no_face);

                       if( face->at_boundary() )
                       {
                              bdry_info.reinit(cell,
                                               face,
                                               no_face);

                              bdry_responses(dst,
                                             bdry_info,
                                             src.vector(0));
                       }
                }
         }
  }

  // --- calculate all other responses that are not a functional, i.e. they do not require a loop over cells:
  global_responses(dst,
                   src.vector(0));

}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::check_responses()
{
    // Print info:
    const std::type_info& info = typeid(*this);
    FcstUtilities::log << "Pure function " << __FUNCTION__
    << " called in Class "
    << info.name() << std::endl;

    // Throw an exception:
    ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::cell_responses(std::vector<double>& ,
                                                        const typename DoFApplication<dim>::CellInfo& ,
                                                        const FuelCell::ApplicationCore::FEVector& )
{
    // Print info:
    //const std::type_info& info = typeid(*this);
    //FcstUtilities::log << "Pure function " << __FUNCTION__
    //<< " called in Class "
    //<< info.name() << std::endl;

    // Throw an exception:
    ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::bdry_responses(std::vector<double>& ,
                                                        const typename DoFApplication<dim>::FaceInfo& ,
                                                        const FuelCell::ApplicationCore::FEVector& )
{
    // Print info:
    //const std::type_info& info = typeid(*this);
    //FcstUtilities::log << "Pure function " << __FUNCTION__
    //<< " called in Class "
    //<< info.name() << std::endl;

    // Throw an exception:
    ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::global_responses (std::vector<double>& ,
                                                           const FuelCell::ApplicationCore::FEVector& )
{
    // Print info:
    //const std::type_info& info = typeid(*this);
    //FcstUtilities::log << "Pure function " << __FUNCTION__
    //<< " called in Class "
    //<< info.name() << std::endl;

    // Throw an exception:
    ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::print_responses (std::vector<double>& )
{
    // Print info:
    const std::type_info& info = typeid(*this);
    FcstUtilities::log << "Pure function " << __FUNCTION__
    << " called in Class "
    << info.name() << std::endl;

    // Throw an exception:
    ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::dresponses_dl (std::vector<std::vector<double> >& drespon_dl,
                                                        const FuelCell::ApplicationCore::FEVectors& src)
{



 /*







    Assert(drespon_dl.size() == n_resp,
           ExcDimensionMismatch(drespon_dl.size(), n_resp));
    Assert(drespon_dl[0].size() == n_dvar,
           ExcDimensionMismatch(drespon_dl[0].size(), n_dvar));

    // Assert (this->quadrature_residual_cell, ExcNotInitialized());
    const Quadrature<dim>& cell_quadrature = this->quadrature_residual_cell;

    LocalResidual<dim> integrator;
    integrator.boundary_fluxes = this->boundary_fluxes;
    integrator.interior_fluxes = this->interior_fluxes;
    FEVectors aux;
    //aux.add_vector(dst, "Residual");
    integrator.initialize(this, this->block_info, aux);
    integrator.add_update_flags(update_values | update_gradients | update_q_points);
    MeshWorker::VectorSelector tmp;

    for (unsigned int i=0;i<src.n_vectors();++i)
        tmp.add_vector(src.vector_name(i), true, true, false);
    integrator.initialize_selectors(tmp, tmp, tmp);
    integrator.cell_quadrature = this->quadrature_residual_cell;
    integrator.bdry_quadrature = this->quadrature_residual_bdry;
    integrator.face_quadrature = this->quadrature_residual_face;

    MeshWorker::InfoObjects::IntegrationInfoBox<dim> box(this->block_info);
    box.initialize(integrator, *this->element, *this->mapping, src);


    // Initialize response vector:
    for (unsigned int i=0; i<drespon_dl.size(); ++i)
        for (unsigned int j=0; j<drespon_dl[i].size(); ++j)
            drespon_dl[i][j] = 0.0;

        this->tr->clear_user_flags ();

    typename DoFHandler<dim>::active_cell_iterator c;
    for (c = this->dof->begin_active();  c != this->dof->end();    ++c)
    {
        // Initialize local structures
        box.cell_info.reinit(c);

        // Compute dresp_dl
        cell_dresponses_dl(drespon_dl,
                           box.cell_info, //instead of cell_info
                           src.vector(0));
    } // cells

    global_dresponses_dl(drespon_dl,
                         src.vector(0));















                         */

}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::cell_dresponses_dl(std::vector<std::vector<double> >& ,
                                                            const typename DoFApplication<dim>::CellInfo& ,
                                                            const FuelCell::ApplicationCore::FEVector& )
{
    // Print info:
    const std::type_info& info = typeid(*this);
    FcstUtilities::log << "Pure function " << __FUNCTION__
    << " called in Class "
    << info.name() << std::endl;

    // Throw an exception:
    ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::global_dresponses_dl(std::vector<std::vector<double> >& ,
                                                              const FuelCell::ApplicationCore::FEVector& )
{
    // Print info:
    const std::type_info& info = typeid(*this);
    FcstUtilities::log << "Pure function " << __FUNCTION__
    << " called in Class "
    << info.name() << std::endl;

    // Throw an exception:
    ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector >& df_du,
                                                       const FuelCell::ApplicationCore::FEVectors& sol)
{
    bool dresponses_du = true;
    bool dresidual_dlambda = false;
    this->dfunction(df_du, sol, dresponses_du, dresidual_dlambda);

    global_dresponses_du(df_du,
                         sol.vector(0));
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::cell_dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector >& ,
                                                            const typename DoFApplication<dim>::CellInfo& ,
                                                            std::vector<std::vector<double> >& )
{
    // Print info:
    const std::type_info& info = typeid(*this);
    FcstUtilities::log << "Pure function " << __FUNCTION__
    << " called in Class "
    << info.name() << std::endl;

    // Throw an exception:
    ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::global_dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector >& ,
                                                              const FuelCell::ApplicationCore::FEVector& )
{
    // Print info:
    const std::type_info& info = typeid(*this);
    FcstUtilities::log << "Pure function " << __FUNCTION__
    << " called in Class "
    << info.name() << std::endl;

    // Throw an exception:
    ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::dresidual_dlambda(std::vector<FuelCell::ApplicationCore::FEVector >& dst,
                                                           const FuelCell::ApplicationCore::FEVectors& src)
{
    bool dresponses_du = false;
    bool dresidual_dlambda = true;
    this->dfunction(dst, src, dresponses_du, dresidual_dlambda);
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::cell_dresidual_dlambda(std::vector<FuelCell::ApplicationCore::FEVector >& ,
                                                                const typename DoFApplication<dim>::CellInfo& ,
                                                                std::vector<std::vector<double> >& )
{
    // Print info:
    const std::type_info& info = typeid(*this);
    FcstUtilities::log << "Pure function " << __FUNCTION__
    << " called in Class "
    << info.name() << std::endl;

    // Throw an exception:
    ExcNotImplemented ();
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::solve_direct (std::vector<std::vector<double> >& df_dl,
                                                       const FuelCell::ApplicationCore::FEVectors& vectors)
{

    // Assertations
    Assert(df_dl[0].size() == n_dvar,
           ExcDimensionMismatch(df_dl[0].size(), n_dvar));

    // Set flag to assemble dR/du (Done in solve now)
    this->notify (sensitivity_analysis); // this->assemble(vectors);

    // Assemble dR/dl
    std::vector<unsigned int> block_sizes(this->block_info.local.size());
    for (unsigned int i=0; i<block_sizes.size(); ++i)
        block_sizes[i] = this->block_info.local.block_size(i);
    std::vector<FuelCell::ApplicationCore::FEVector >dR_dl(n_dvar, FuelCell::ApplicationCore::FEVector (this->block_info.global));
    this->dresidual_dlambda(dR_dl,
                            vectors);
    // Solve system to obtain -du_dl by solving -dR/du*(-du/dl) = -dR_dl
    std::vector<FuelCell::ApplicationCore::FEVector >du_dl(n_dvar, FuelCell::ApplicationCore::FEVector (this->block_info.global));
    for (unsigned int i=0; i<n_dvar; ++i)
    {
        FuelCell::ApplicationCore::FEVectors aux;
        aux.add_vector(dR_dl[i],"Newton residual");
        unsigned int ind = vectors.find_vector("Solution");
        aux.add_vector(vectors.vector(ind), "Newton iterate");
        this->solve(du_dl[i],
                    aux);//dR_dl[i]);
    }

    // Obtain df_dl
    std::vector<std::vector<double> > pdf_pdl(this->n_resp, std::vector<double> (n_dvar));
    OptimizationBlockMatrixApplication<dim>::dresponses_dl(pdf_pdl,
                                                           vectors);//sol);

    //print_dresponses_dl(pdf_pdl);
    std::vector<FuelCell::ApplicationCore::FEVector > df_du(this->n_resp, FuelCell::ApplicationCore::FEVector (this->block_info.global));
    OptimizationBlockMatrixApplication<dim>::dresponses_du(df_du, vectors);//sol);
    //print_dresponses_du(df_du);

    for (unsigned int i=0; i < df_dl.size(); ++i)
        for (unsigned int j=0; j<df_dl[0].size(); ++j)
        {
            df_dl[i][j] = pdf_pdl[i][j] - df_du[i]*du_dl[j]; //Note the negative sign is due to solving -dR/du*(-du/dl) = -dR_dl
        }

        // Write gradents:
        for (unsigned int i=0; i<n_resp; ++i)
            for (unsigned int j=0; j<n_dvar; ++j)
                FcstUtilities::log<<"Df"<<i<<"/Dl"<<j<<" is equal to "<<df_dl[i][j]<<std::endl;

}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::solve_adjoint (std::vector<std::vector<double> >& /*df_dl*/,
                                                        const FuelCell::ApplicationCore::FEVector& /*sol*/)
{

    //For the adjoint I need the transpose of the system matrix... => Guido added a transponse matrix to LAC

}

//---------------------------------------------------------------------------
template <int dim>
unsigned int
OptimizationBlockMatrixApplication<dim>::get_n_resp() const
{
    return n_resp;
}

//---------------------------------------------------------------------------
template <int dim>
unsigned int
OptimizationBlockMatrixApplication<dim>::get_n_dvar() const
{
    return n_dvar;
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::set_optimization_parameters(unsigned int &n_dvar, unsigned int &n_resp,
                                                                     std::vector<std::string> &name_design_var, std::vector<std::string> &name_responses)
{
    this->n_dvar = n_dvar;
    this->n_resp = n_resp;
    this->name_design_var = name_design_var;
    set_output_variables(name_responses);
    optimization = true;
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::dfunction(std::vector<FuelCell::ApplicationCore::FEVector >& dst,
                                                   const FuelCell::ApplicationCore::FEVectors& src,
                                                   bool dfunctional_du,
                                                   bool dresidual_dlambda)
{
  /*
  FuelCell::ApplicationCore::FEVectors dummy;

  const Quadrature<dim>& cell_quadrature = this->quadrature_residual_cell;
  const Quadrature<dim-1>& face_quadrature = this->quadrature_residual_face;

  // Loop over all objective function and constraints (functionals)
  for (unsigned int f = 0; f < dst.size(); ++f)
    dst[f].reinit(this->block_info.global);

       // Generate FEValues object for all
       // base elements.
  std::vector<FEValues<dim>*> fe(this->element->n_base_elements());
  std::vector<FEFaceValues<dim>*> int_fe(this->element->n_base_elements());
  std::vector<FESubfaceValues<dim>*> int_fe_sub(this->element->n_base_elements());
  std::vector<FEFaceValues<dim>*> ext_fe(this->element->n_base_elements());
       // and initialize them
  for (unsigned int i=0;i<fe.size();++i)
    {
      fe[i] = new FEValues<dim> (
 *this->mapping,
 this->element->base_element(i),
 cell_quadrature,
 UpdateFlags(update_q_points
      | update_values
      | update_gradients
      | update_JxW_values));

      int_fe[i] = new FEFaceValues<dim> (
 *this->mapping,
 this->element->base_element(i),
 face_quadrature,
 UpdateFlags(update_q_points
      | update_values
      | update_gradients
      | update_normal_vectors
      | update_JxW_values));
      int_fe_sub[i] = new FESubfaceValues<dim> (
 *this->mapping,
 this->element->base_element(i),
 face_quadrature,
 UpdateFlags(update_q_points
      | update_values
      | update_gradients
      | update_normal_vectors
      | update_JxW_values));
      ext_fe[i] = new FEFaceValues<dim> (
 *this->mapping,
 this->element->base_element(i),
 face_quadrature,
 UpdateFlags(update_values
      | update_gradients));
    }

  const unsigned int n_dofs = this->element->dofs_per_cell;
  const unsigned int n_blocks = this->element->n_blocks();
  const unsigned int n_components = this->element->n_components();

       // Mapping from local to global
       // indices
  std::vector<unsigned int> indices(n_dofs);
  std::vector<unsigned int> indices_org(n_dofs);
  std::vector<unsigned int> neighbor_indices(n_dofs);
  std::vector<unsigned int> neighbor_indices_org(n_dofs);

                                   // Objects for cell and face data
       // handed down to the local
       // routines.
  typename DoFApplication<dim>::CellInfo cell_info(src, fe, indices);

  // Initialize local data vectors
  cell_info.values.resize(src.n_vectors(), std::vector<std::vector<double> >(
     n_components, std::vector<double>(
    cell_quadrature.size())));

  cell_info.derivatives.resize(src.n_vectors(), std::vector<std::vector<Tensor<1,dim> > >(
     n_components, std::vector<Tensor<1,dim> >(
       cell_quadrature.size())));

       // Result vectors for the cell wise
       // functions
  std::vector<FuelCell::ApplicationCore::FEVector > intvec (dst.size(), FuelCell::ApplicationCore::FEVector (n_blocks));
  std::vector<FuelCell::ApplicationCore::FEVector > extvec (dst.size(), FuelCell::ApplicationCore::FEVector (n_blocks));

       // Cell wise source vector
       // functions
  std::vector<std::vector<double> > intsrc;
  std::vector<std::vector<double> > extsrc;

  for (unsigned int f = 0; f < dst.size(); ++f)
    {
      for (unsigned int i=0;i<intvec[f].n_blocks();++i)
 {
   const unsigned int s = this->element->base_element(this->element->block_to_base_index(i).first).dofs_per_cell;
   intvec[f].block(i).reinit(s);
   extvec[f].block(i).reinit(s);
 }
      intvec[f].collect_sizes();
      extvec[f].collect_sizes();
    }

  this->tr->clear_user_flags ();

  typename DoFHandler<dim>::active_cell_iterator c;
  for (c = this->dof->begin_active();
       c != this->dof->end();
       ++c)
    {
      // Initialize local structures
      cell_info.reinit(c);
      c->get_dof_indices(indices_org);
      for (unsigned int i=0;i<indices.size();++i)
 indices[this->block_info.local_renumbering[i]] = indices_org[i];

      for (unsigned int f = 0; f < dst.size(); ++f)
 {
   intvec[f] = 0.;
   extvec[f] = 0.;
 }

      for (unsigned int i=0;i<fe.size();++i)
 fe[i]->reinit((typename Triangulation<dim>::cell_iterator) c);

      // Fill local data vectors
      this->fill_local_data(cell_info, cell_info.values, src, true);
      this->fill_local_data(cell_info, cell_info.derivatives, src, true);
      // Choose the matrix we want to assemble. I put this if because this way I don't need to
      // code the same thing twice.
      if (dfunctional_du == true)
 cell_dresponses_du(intvec, cell_info,
      intsrc);
      else if (dresidual_dlambda == true)
 cell_dresidual_dlambda(intvec, cell_info,
          intsrc);
      else
 {
   FcstUtilities::log<<"Nothing to be done in optimization_block_matrix_application.cc"<<std::endl;
   abort();
   }
      // Insted of:
      for (unsigned int f = 0; f < dst.size(); ++f)
 for (unsigned int i=0;i<n_dofs;++i)
   dst[f](indices[i]) += intvec[f](i);
      // Use:
      // Can't because intvec is a BlockVector and this routine needs a Vector
      //hanging_node_constraints.distribute_local_to_global(intvec, indices, dst);
    } // cells

  //////////////////////////////////////////////////////////////////////
  // Condense rhs
  //////////////////////////////////////////////////////////////////////
  for (unsigned int f = 0; f < dst.size(); ++f)
    this->hanging_node_constraints.condense(dst[f]);

  //////////////////////////////////////////////////////////////////////
  // Apply Dirichlet boundary conditions
  //////////////////////////////////////////////////////////////////////
  //for (unsigned int f = 0; f < dst.size(); ++f)
  //  this->filtered_matrix.apply_constraints(dst[f],
  //         false);
  */
}

//---------------------------------------------------------------------------
template <int dim>
const std::string
OptimizationBlockMatrixApplication<dim>::extend_filename(const std::string& basename, const int precision) const
{
    /*
    // Create string that will be appended to the .vtk filename. This string will contain the
    // name and value of the design variables.
    std::string filename;
    
    //loop over all design variables
    for (unsigned int i=0; i<name_design_var.size(); i++)
    {
        //create object that can print out numbers
        std::ostringstream streamOut;
        //Used to set how many digits after the decimal point are used. Insert the required number into the brackets.
        streamOut.precision(precision);
        //If there are fewer digits in the design variable than are requested, the remaining digits can be filled with zeroes.
        //Use std::fixed to force this: streamOut <<  std::fixed << design_var_value[design_var.size()-i-1];
        //This can be used to get paraview to group the vtk files correctly.
        
        //Note that the design variables are printed in reverse, so that file managers will group the vtk files by the last design
        // variable first before moving to the second last etc.
        streamOut <<  std::fixed <<  design_var_value[name_design_var.size()-i-1];
        
        //Append design variable name and value to filename1
        filename += "__" + name_design_var[name_design_var.size()-i-1] + "_" + streamOut.str();
    }
    
    const std::string suffix = this->d_out.default_suffix();
    if (suffix == "") return;
    
    //construct the full filename, where basename is given by the adaptive refinement loop and suffix is the file extension (i.e. .vtk)
    const std::string filename = basename + filename + suffix;
    FcstUtilities::log << "Datafile:" << filename << std::endl;
    std::ofstream of (filename.c_str ());  
    */
}

//---------------------------------------------------------------------------
template <int dim>
void
OptimizationBlockMatrixApplication<dim>::set_all_response_names()
{
    all_response_names.push_back("current");
    all_response_names.push_back("cathode_current");
    all_response_names.push_back("anode_current");
    all_response_names.push_back("OH_coverage");
    all_response_names.push_back("O_coverage");
    all_response_names.push_back("water_cathode");
    all_response_names.push_back("water_anode");
    all_response_names.push_back("cathode_reaction_heat");
    all_response_names.push_back("cathode_irrev_heat");
    all_response_names.push_back("cathode_rev_heat");
    all_response_names.push_back("cathode_watervap_heat");
    all_response_names.push_back("anode_reaction_heat");
    all_response_names.push_back("anode_irrev_heat");
    all_response_names.push_back("anode_rev_heat");
    all_response_names.push_back("sorption_heat_cathode");
    all_response_names.push_back("sorption_heat_anode");
    all_response_names.push_back("electron_ohmic_heat");
    all_response_names.push_back("proton_ohmic_heat");
    all_response_names.push_back("PEM_proton_ohmic_heat");    
    all_response_names.push_back("max_temperature");
    all_response_names.push_back("max_RH");    
    all_response_names.push_back("t_m_Pt");        
    all_response_names.push_back("epsilon_V_cat_c:4");
    all_response_names.push_back("epsilon_S_cat_c:4");
    all_response_names.push_back("epsilon_N_cat_c:4");
    all_response_names.push_back("epsilon_V_cat_c:5");
    all_response_names.push_back("epsilon_S_cat_c:5");
    all_response_names.push_back("epsilon_N_cat_c:5");
    all_response_names.push_back("water_prod_CCL");
    all_response_names.push_back("water_prod_ACL");
    all_response_names.push_back("water_evap_CCL");
    all_response_names.push_back("water_evap_ACL");
    all_response_names.push_back("water_evap_CMPL");
    all_response_names.push_back("water_evap_CGDL");
    all_response_names.push_back("water_evap_overall");
    all_response_names.push_back("heat_evap_CCL");
    all_response_names.push_back("heat_evap_ACL");
    all_response_names.push_back("heat_evap_overall");
    all_response_names.push_back("water_cond_overall");
    all_response_names.push_back("heat_cond_overall");
    all_response_names.push_back("oxygen_molar_fraction");
    all_response_names.push_back("electronic_electrical_potential");  
    all_response_names.push_back("total_mass_flux_species_1");
    all_response_names.push_back("total_mass_flux_species_2");
    all_response_names.push_back("total_mass_flux_species_3");
    all_response_names.push_back("total_mass_flux_species_4");
    all_response_names.push_back("total_mass_flux_species_5");
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template class FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<deal_II_dimension>;
