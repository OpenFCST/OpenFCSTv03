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
// - Class: block_matrix_application.cc
// - Description:
// - Developers: Guido Kanschat, Texas A&M University
//               Marc Secanell,  University of Alberta
//
// ----------------------------------------------------------------------------

#include <application_core/block_matrix_application.h>

//------------------------------
template<int dim>
BlockMatrixApplication<dim>::BlockMatrixApplication(
        boost::shared_ptr<ApplicationData> data) :
        DoFApplication<dim>(data)
{
    repair_diagonal = false; //false as standard unless set by child
    FcstUtilities::log << "->BlockMatrix";
}

//---------------------------------------------------------------------------
template<int dim>
BlockMatrixApplication<dim>::BlockMatrixApplication(DoFApplication<dim>& other,
        bool triangulation_only) :
        DoFApplication<dim>(other, triangulation_only)
{
    repair_diagonal = false; //false as standard unless set by child
    FcstUtilities::log << "->BlockMatrix";
}

//---------------------------------------------------------------------------
template<int dim>
void BlockMatrixApplication<dim>::declare_parameters(ParameterHandler& param) 
{
    
    // First, declare all parameters from parent class:
    DoFApplication<dim>::declare_parameters(param);
    
    // Declare application specific parameters:
    param.enter_subsection("Discretization");
    {
        param.enter_subsection("Matrix");
        {
            param.declare_entry("Quadrature cell", "-1", Patterns::Integer(),
                                "Gauss formula used for integrating the "
                                "bilinear form on a cell. A negative value means "
                                "adding its absolute value to polynomial degree of "
                                "the shape functions.");
            param.declare_entry("Quadrature face", "-1", Patterns::Integer(),
                                "Gauss formula used for integrating the bilinear "
                                "form on faces. See above for more details.");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    param.enter_subsection("Linear Solver");
    {
        param.declare_entry("Print matrix and rhs", "false", Patterns::Bool(),
                            "Specify if you would like to print the matrix and rhs (only for debugging purposes)");
        param.declare_entry("Assemble numerically", "false", Patterns::Bool(),
                            "Specify if you would like to assemble the Jacobian analytically or numerically");
        
        param.declare_entry("Symmetric matrix",
                            "false", 
                            Patterns::Bool(),
                            "Flag for specifying if matrix is symmetric to the linear solver");
        
        
        
        param.declare_entry("Allocate additional memory for MUMPS", "false", Patterns::Bool(),
                            "Configure MUMPS solver to use additional memory during solving if needed.");
        
        param.declare_entry("Output system assembling time",
                            "false", 
                            Patterns::Bool(),
                            "Flag for specifying if you want to see how much time the linear system assembling takes.");        
        
        SolverControl::declare_parameters(param);
        
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
template<int dim>
void BlockMatrixApplication<dim>::initialize(ParameterHandler& param) 
{
    // Initialize parent class:
    DoFApplication<dim>::initialize(param);

    
    _initialize(param);
    

}

//---------------------------------------------------------------------------
template<int dim>
void BlockMatrixApplication<dim>::_initialize(ParameterHandler& param) 
{

    // Initialize default values for flux_couplings (Note that this can be redefined in your application)
    this->flux_couplings.reinit(this->element->n_blocks(),
            this->element->n_blocks());
    for (unsigned int i = 0; i < this->element->n_blocks(); ++i) {
        for (unsigned int j = 0; j < this->element->n_blocks(); ++j) {
            this->flux_couplings(i, j) = DoFTools::none;
        }
    }

    // Get relevant parameter values
    // and initialize objects
    // accordingly.
    param.enter_subsection("Discretization");
    {
        param.enter_subsection("Matrix");
        {
            if (true) {
                int n = param.get_integer("Quadrature cell");
                if (n <= 0)
                    n = this->element->degree - n;
                boost::shared_ptr<Quadrature<dim> > newq(new QGauss<dim>(n));
                quadrature_assemble_cell = newq;
            }
            if (true) {
                int n = param.get_integer("Quadrature face");
                if (n <= 0)
                    n = this->element->degree - n;
                boost::shared_ptr<Quadrature<dim - 1> > newq(
                        new QGauss<dim - 1>(n));
                quadrature_assemble_face = newq;
            }
        }
        param.leave_subsection();
    }
    param.leave_subsection();

    // Declare application specific parameters:
    param.enter_subsection("Linear Solver");
    {
        print_debug = param.get_bool("Print matrix and rhs");
        assemble_numerically_flag = param.get_bool("Assemble numerically");
        mumps_additional_mem = param.get_bool("Allocate additional memory for MUMPS");
        symmetric_matrix_flag = param.get_bool("Symmetric matrix");
        output_system_assembling_time = param.get_bool("Output system assembling time");
        solver_control.parse_parameters(param);
        solver_control.log_history(false);
        solver_control.log_result(false);

    }
    param.leave_subsection();
    
}

//---------------------------------------------------------------------------
template<int dim>
void BlockMatrixApplication<dim>::remesh_matrices() 
{
    
    this->block_info.initialize_local(*this->dof);

    // clear matrices
    matrix.clear();
    // Make the list of constraints associated with hanging nodes
    this->hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints(*this->dof, this->hanging_node_constraints);
    this->hanging_node_constraints.close();

#ifdef OPENFCST_WITH_PETSC
    #if deal_II_dimension < 3 //For 2D simulations can use SparsityPattern for computation speed up
        SparsityPattern sparsity(this->dof->n_dofs(), this->dof->n_dofs(),this->dof->max_couplings_between_dofs());
        DoFTools::make_sparsity_pattern(*this->dof, this->cell_couplings, sparsity);
        this->hanging_node_constraints.condense(sparsity);
    #else //Due to memory allocation issues for 3D problems use CompressedSimpleSparsityPattern instead
        CompressedSimpleSparsityPattern compressed_pattern(this->dof->n_dofs());
        DoFTools::make_sparsity_pattern(*this->dof, compressed_pattern);
        this->hanging_node_constraints.condense(compressed_pattern);
        
        SparsityPattern sparsity;
        sparsity.copy_from(compressed_pattern);
    #endif

    std::vector<unsigned int> local_rows_colums_per_process;

    for(unsigned int i = 0; i < this->n_mpi_processes; i++)
        local_rows_colums_per_process.push_back( DoFTools::count_dofs_with_subdomain_association(*this->dof, i));


    matrix.reinit (this->mpi_communicator,
                   sparsity,
                   local_rows_colums_per_process,
                   local_rows_colums_per_process,
                   this->this_mpi_process);


    boundary_values.clear();
    this->dirichlet_bc(boundary_values);

#else

    const unsigned int n_blocks = this->element->n_blocks();

    BlockCompressedSparsityPattern c_sparsity(n_blocks, n_blocks);
    for (unsigned int i = 0; i < n_blocks; ++i)
        for (unsigned int j = 0; j < n_blocks; ++j)
            c_sparsity.block(i, j).reinit(this->block_info.global.block_size(i), this->block_info.global.block_size(j));
    c_sparsity.collect_sizes();

    if (this->interior_fluxes)
        DoFTools::make_flux_sparsity_pattern(*this->dof, c_sparsity, this->cell_couplings, this->flux_couplings);
    else
        DoFTools::make_sparsity_pattern(*this->dof, this->cell_couplings, c_sparsity);

    // Condense sparsity pattern to account for hanging nodes
    this->hanging_node_constraints.condense(c_sparsity);

    sparsities.copy_from(c_sparsity);
    sparsities.compress();

    matrix.reinit(sparsities);

    //////////////////////////////////////////////////////////////////////
    // Calculate Dirichlet boundary conditions
    //////////////////////////////////////////////////////////////////////
    boundary_values.clear();
    this->dirichlet_bc(boundary_values);

#endif

}

//----------------------------
template<int dim>
void BlockMatrixApplication<dim>::remesh() {
    DoFApplication<dim>::remesh();
    remesh_matrices();
}


//---------------------------------------------------------------------------
#ifdef OPENFCST_WITH_PETSC
template<int dim>
void BlockMatrixApplication<dim>::PETSc_solve(FuelCell::ApplicationCore::FEVector system_rhs,
                                              FEVector& solution, 
                                              const FEVectors& src) 
{    
    //Make parallel copies of serial vectors
    PETScWrappers::MPI::Vector del_sol, sys_rhs;
    
    const types::global_dof_index n_local_dofs = DoFTools::count_dofs_with_subdomain_association (*this->dof,this->this_mpi_process);
    del_sol.reinit(this->mpi_communicator, this->dof->n_dofs(), n_local_dofs);
    sys_rhs.reinit(this->mpi_communicator, this->dof->n_dofs(), n_local_dofs);
    
    FEValues<dim>*     fe_values      = 0;
    FuelCell::ApplicationCore::IntegrationInfo< dim, FEValuesBase<dim> >  cell_info(src, this->block_info);
 
    cell_info.initialize(
        fe_values,
        *this->element,
        *this->mapping,
        this->quadrature_residual_cell,
        UpdateFlags(update_q_points | update_values | update_gradients | update_JxW_values)
    );
    
    typename DoFHandler<dim>::active_cell_iterator
    cell = this->dof->begin_active(),
    endc = this->dof->end();
    
    for( ; cell != endc; ++cell)
    {
        cell_info.reinit(cell);      
        
        if(cell->subdomain_id() == this->this_mpi_process){            
            for(unsigned int i = 0; i < this->element->dofs_per_cell; ++i)
            {
                del_sol(cell_info.indices[i]) = solution(cell_info.indices[i]);
                sys_rhs(cell_info.indices[i]) = system_rhs(cell_info.indices[i]);
            }
            
        }
    }
    
    del_sol.compress(VectorOperation::insert);
    sys_rhs.compress(VectorOperation::insert);
        
    MatrixTools::apply_boundary_values (boundary_values, matrix, del_sol, sys_rhs, false);
    
    this->print_matrix_and_rhs(sys_rhs);
    
    if (this->data->get_linear_solver() == FuelCell::ApplicationCore::LinearSolver::CG) {
        FcstUtilities::log << "Solving linear system with CG..." << std::endl;
        
        PETScWrappers::PreconditionJacobi prec(matrix);        
        PETScWrappers::SolverCG solver(solver_control, this->mpi_communicator);        
        solver.solve(this->matrix, del_sol, sys_rhs, prec);
    }
    else if (this->data->get_linear_solver() == FuelCell::ApplicationCore::LinearSolver::BICGSTAB) {
        FcstUtilities::log << "Solving linear system with Bicgstab..." << std::endl;
        
        PETScWrappers::PreconditionJacobi prec(matrix);        
        PETScWrappers::SolverBicgstab solver(solver_control, this->mpi_communicator);        
        solver.solve(this->matrix, del_sol, sys_rhs, prec);
    }
    else if (this->data->get_linear_solver() == FuelCell::ApplicationCore::LinearSolver::ILU_GMRES) {
        FcstUtilities::log << "Solving linear system with ILU-GMRES..." << std::endl;
        
        PETScWrappers::PreconditionNone prec(matrix);
        PETScWrappers::SolverGMRES solver(solver_control, this->mpi_communicator);
        solver.solve(this->matrix,
                     del_sol,
                     sys_rhs,
                     prec);
    }
    else if (this->data->get_linear_solver() == FuelCell::ApplicationCore::LinearSolver::MUMPS) {
        FcstUtilities::log << "Solving linear system with MUMPS..." << std::endl;
        
        PETScWrappers::SparseDirectMUMPS solver(solver_control, this->mpi_communicator);
        solver.set_symmetric_mode(symmetric_matrix_flag);
        
        if (mumps_additional_mem)
            solver.set_prefix(std::string("mat_mumps_icntl_14 100"));
        
        solver.solve(this->matrix,
                     del_sol,
                     sys_rhs);
    }
    else {
        
        const std::type_info& info = typeid (*this);
        FcstUtilities::log << "Solver not implemented in class " << info.name() << " member function solve()" << std::endl;
        abort();
    }   
    
    this->hanging_node_constraints.distribute(del_sol);
    
    //Copy to linear dealii vector
    const PETScWrappers::Vector localized_solution(del_sol);
    solution = localized_solution;
                                                  
}
#else
//---------------------------------------------------------------------------
template<int dim>
void BlockMatrixApplication<dim>::serial_solve(FuelCell::ApplicationCore::FEVector system_rhs,
                                               FEVector& solution) 
{     
    // --- Apply boundary conditions ---
    MatrixTools::apply_boundary_values(boundary_values, matrix, solution, system_rhs, false);

    this->print_matrix_and_rhs(system_rhs);

            
    if (this->data->get_linear_solver() == FuelCell::ApplicationCore::LinearSolver::CG) {
        
        SolverCG<FEVector> solver(solver_control, this->data->block_vector_pool);
        solver.solve(this->matrix, solution, system_rhs, PreconditionIdentity());

    }
    else if (this->data->get_linear_solver() == FuelCell::ApplicationCore::LinearSolver::BICGSTAB) {

        SolverBicgstab<FEVector> solver(solver_control, this->data->block_vector_pool);
        solver.solve(this->matrix, solution, system_rhs, PreconditionIdentity());

    }
    else if (this->data->get_linear_solver() == FuelCell::ApplicationCore::LinearSolver::ILU_GMRES) {

        LinearSolvers::ILUPreconditioner prec(this->matrix);
        LinearSolvers::GMRESSolver solver;
        solver.solve(solver_control, this->matrix, solution, system_rhs, prec.preconditioner);
    }
    else if (this->data->get_linear_solver() == FuelCell::ApplicationCore::LinearSolver::UMFPACK) {
        
        LinearSolvers::SparseDirectUMFPACKSolver solver;
        solver.solve(this->matrix, solution, system_rhs);
    }
    
    else {
        const std::type_info& info = typeid (*this);
        FcstUtilities::log << "Linear Solver" <<"not implemented in class " << info.name() << " member function solve()" << std::endl;
        abort();
    }

    // --- Finally apply hanging node constraints to solution ---
    this->hanging_node_constraints.distribute(solution);
}
#endif

//---------------------------------------------------------------------------
template<int dim>
void BlockMatrixApplication<dim>::solve(FEVector& solution,
        const FEVectors& src) {

    unsigned int index;

    timer.restart();
    
    // --- Matrix assembly ---
    if (this->notifications.any()) 
    {

        // --- Source the solution ---
        FuelCell::ApplicationCore::FEVectors sol;
        std::string solution_vector_name = this->data->get_solution_vector_name(this->data->get_nonlinear_solver());
        
        if (src.count_vector(solution_vector_name))
        {
            index = src.find_vector(solution_vector_name);
            sol.add_vector(src.vector(index), solution_vector_name);
        }

        else
            throw std::runtime_error("BlockMatrixApplication<dim>::solve "
                    "cannot find solution from FEVectors& src.");

        // --- Assemble ---
        if (this->assemble_numerically_flag)
            this->assemble_numerically(sol);
        else
            this->assemble(sol); //Note the second component of

    }

    // --- Source the rhs ---
    FuelCell::ApplicationCore::FEVector system_rhs;
    system_rhs.reinit(this->block_info.global);
    
    std::string residual_vector_name = this->data->get_residual_vector_name(this->data->get_nonlinear_solver());

    if (src.count_vector(residual_vector_name))
    {
        index = src.find_vector(residual_vector_name);
        system_rhs =  src.vector(index);
    }
    else
        throw std::runtime_error("BlockMatrixApplication<dim>::solve "
                "cannot find residual from FEVectors& src.");


    // --- Repair diagonal elements ---
    if (repair_diagonal)
        SolverUtils::repair_diagonal(matrix);
/*
     FcstUtilities::log<<"== PRINTING MATRIX =="<<std::endl;
     std::ofstream file;
     file.open("Matrix.dat");
     this->matrix.print(file,true);
     file.close();
     FcstUtilities::log<<"============================"<<std::endl;
     exit(-1);
     */

    timer.stop();
    
    if (output_system_assembling_time)
        if (this->notifications.any())
        {
            FcstUtilities::log << "The linear system was assembled in " << timer.wall_time() << " seconds." << std::endl;
            FcstUtilities::log << "Jacobians were computed ";
            if (this->assemble_numerically_flag)
            {
                FcstUtilities::log << "numerically because <Assemble numerically> flag was set to <true>." << std::endl;
            }
            else
            {
                FcstUtilities::log << "analytically because <Assemble numerically> flag was set to <false>." << std::endl;
            }
            this->notifications.clear();
        }
        else
            FcstUtilities::log << "The right hand side was assembled in " << timer.wall_time() << " seconds." << std::endl;
    else
        if (this->notifications.any())
            this->notifications.clear();

    // --- solve ---
    #ifdef OPENFCST_WITH_PETSC
        PETSc_solve(system_rhs, solution, src);
    #else
        serial_solve(system_rhs, solution);
    #endif



}

//------------------------------
template<int dim>
void BlockMatrixApplication<dim>::assemble(const FEVectors& src) {
    // Clear old matrices
    matrix = 0;

    #ifdef OPENFCST_WITH_PETSC
        PETSc_assemble(src);
    #else
        serial_assemble(src);
    #endif

}
//------------------------------
//------------------------------

#ifndef OPENFCST_WITH_PETSC
// Define these so the multigrid and single grid codes become more
// similar
#define MATRIX this->matrix.block(block_row, block_col)
#define INDEX(i,l) this->block_info.global.global_to_local(i).second

template<int dim>
void BlockMatrixApplication<dim>::serial_assemble(const FEVectors& src) {

    typedef typename DoFHandler<dim>::active_cell_iterator active_cell_iterator;

    active_cell_iterator begin = this->dof->begin_active();
    active_cell_iterator end = this->dof->end();

    FcstUtilities::log.push("Assembly");
    Assert(quadrature_assemble_cell, ExcNotInitialized());
    Assert(quadrature_assemble_face, ExcNotInitialized());
    const Quadrature<dim>& cell_quadrature = *quadrature_assemble_cell;
    const Quadrature<dim - 1>& face_quadrature = *quadrature_assemble_face;

    // Objects for cell and face data
    // handed down to the local
    // routines.
    FEValues<dim>* fevalues = 0;
    FEFaceValues<dim>* fefacevalues = 0;
    FESubfaceValues<dim>* fesubfacevalues = 0;
    typename DoFApplication<dim>::CellInfo cell_info(src, this->block_info);
    typename DoFApplication<dim>::FaceInfo bdry_info(src, this->block_info);
    typename DoFApplication<dim>::FaceInfo face_info(src, this->block_info);
    typename DoFApplication<dim>::FaceInfo subface_info(src, this->block_info);
    typename DoFApplication<dim>::FaceInfo neighbor_info(src, this->block_info);

    cell_info.initialize(fevalues, *this->element, *this->mapping,
            cell_quadrature,
            UpdateFlags(
                    update_q_points | update_values | update_gradients
                            | update_JxW_values));
    bdry_info.initialize(fefacevalues, *this->element, *this->mapping,
            face_quadrature,
            UpdateFlags(
                    update_q_points | update_values | update_gradients
                            | update_normal_vectors | update_JxW_values));
    face_info.initialize(fefacevalues, *this->element, *this->mapping,
            face_quadrature,
            UpdateFlags(
                    update_q_points | update_values | update_gradients
                            | update_normal_vectors | update_JxW_values));
    subface_info.initialize(fesubfacevalues, *this->element, *this->mapping,
            face_quadrature,
            UpdateFlags(
                    update_q_points | update_values | update_gradients
                            | update_normal_vectors | update_JxW_values));
    neighbor_info.initialize(fefacevalues, *this->element, *this->mapping,
            face_quadrature,
            UpdateFlags(
                    update_values | update_gradients | update_normal_vectors
                            | update_JxW_values));

    // Initialize local data
    MatrixVector intint;
    MatrixVector intext;
    MatrixVector extint;
    MatrixVector extext;

    for (unsigned int i = 0; i < matrix.n_block_rows(); ++i)
        for (unsigned int j = 0; j < matrix.n_block_rows(); ++j) {
            if (this->cell_couplings(i, j) == DoFTools::none
                    && this->flux_couplings(i, j) == DoFTools::none)
                continue;

            const unsigned int rows = this->block_info.local.block_size(i);
            const unsigned int cols = this->block_info.local.block_size(j);
            MatrixBlock<FullMatrix<double> > block(i, j);
            intint.push_back(block);
            intext.push_back(block);
            extint.push_back(block);
            extext.push_back(block);
            intint[intint.size() - 1].matrix.reinit(rows, cols);
            intext[intint.size() - 1].matrix.reinit(rows, cols);
            extint[intint.size() - 1].matrix.reinit(rows, cols);
            extext[intint.size() - 1].matrix.reinit(rows, cols);
        }

    this->tr->clear_user_flags();

    typename DoFHandler<dim>::active_cell_iterator c;
    for (c = begin; c != end; ++c) {
        for (unsigned int i = 0; i < intint.size(); ++i)
            intint[i].matrix = 0.;
        // Initialize local structures
        cell_info.reinit(c);

        // Fill local data vectors
        cell_info.fill_local_data(cell_info.values, true);
        cell_info.fill_local_data(cell_info.derivatives, true);
        cell_matrix(intint, cell_info);


        for (unsigned int i = 0; i < intint.size(); ++i) {
            for (unsigned int j = 0; j < intint[i].matrix.n_rows(); ++j)
                for (unsigned int k = 0; k < intint[i].matrix.n_cols(); ++k)
                    if (fabs(intint[i].matrix(j, k)) > 1.e-15) {
                        const unsigned int block_row = intint[i].row;
                        const unsigned int block_col = intint[i].column;
                        const unsigned int jcell =
                                this->block_info.local.local_to_global(
                                        block_row, j);
                        const unsigned int kcell =
                                this->block_info.local.local_to_global(
                                        block_col, k);

                        MATRIX.add(INDEX(cell_info.indices[jcell], level),
                                INDEX(cell_info.indices[kcell], level),
                                intint[i].matrix(j, k));
                    }
        }
    }
    this->post_cell_assemble();


    #if deal_II_dimension > 1

    // //////////////////////////////////////////////////////////////////////
    // // Integration of fluxes
    // //////////////////////////////////////////////////////////////////////
    if (this->interior_fluxes || this->boundary_fluxes)

    for (c = begin; c != end; ++c)
    {
        for (unsigned int face_nr=0;
                face_nr < GeometryInfo<dim>::faces_per_cell;
                ++ face_nr)
        {
            typename DoFHandler<dim>::face_iterator face = c->face(face_nr);

            // Avoid computing this face
            // a second time.
            if (face->user_flag_set ())
            continue;

            if (c->at_boundary(face_nr))
            {
                for (unsigned int i=0;i<intint.size();++i)
                intint[i].matrix = 0.;

                bdry_info.reinit(c, face, face_nr);
                // Boundary values
                // Fill local data vectors
                bdry_info.fill_local_data(bdry_info.values, true);
                bdry_info.fill_local_data(bdry_info.derivatives, true);
                this->bdry_matrix(intint, bdry_info);

                for (unsigned int i=0;i<intint.size();++i)
                {
                    for (unsigned int j=0;j<intint[i].matrix.n_rows();++j)
                    for (unsigned int k=0;k<intint[i].matrix.n_cols();++k)
                    if (fabs(intint[i].matrix(j,k)) > 1.e-12)
                    {
                        const unsigned int block_row = intint[i].row;
                        const unsigned int block_col = intint[i].column;
                        const unsigned int
                        jcell = this->block_info.local.local_to_global(block_row, j);
                        const unsigned int
                        kcell = this->block_info.local.local_to_global(block_col, k);

                        MATRIX.add(INDEX(bdry_info.indices[jcell], level),
                                INDEX(bdry_info.indices[kcell], level),
                                intint[i].matrix(j,k));
                    }
                }
            }
            else if (this->interior_fluxes)
            {
                // Interior face
                typename DoFHandler<dim>::cell_iterator
                neighbor = c->neighbor(face_nr);
                // Do refinement face
                // from the coarse side
                if (neighbor->level() < c->level())
                continue;

                unsigned int neighbor_face_nr = c->neighbor_of_neighbor(face_nr);
                // Make sure this
                // face is not worked
                // on from the other
                // side
                if (face->user_flag_set ())
                    continue;

                face->set_user_flag ();
                if (!neighbor->has_children())
                neighbor->face(neighbor_face_nr)->set_user_flag();

                if (neighbor->has_children())
                {
                    // Refinement face
                    // needs additional
                    // treatment
                    for (unsigned int sub_nr = 0; sub_nr != face->number_of_children(); ++sub_nr)
                    {
                        typename DoFHandler<dim>::cell_iterator sub_neighbor
                        = c->neighbor_child_on_subface(face_nr, sub_nr);

                        for (unsigned int i=0;i<intint.size();++i)
                        {
                            intint[i].matrix = 0.;
                            intext[i].matrix = 0.;
                            extint[i].matrix = 0.;
                            extext[i].matrix = 0.;
                        }

                        subface_info.reinit(c, face, face_nr, sub_nr);
                        neighbor_info.reinit(sub_neighbor,
                                sub_neighbor->face(neighbor_face_nr),
                                neighbor_face_nr);

                        // Fill local data vectors
                        subface_info.fill_local_data(subface_info.values, true);
                        subface_info.fill_local_data(subface_info.derivatives, true);
                        neighbor_info.fill_local_data(neighbor_info.values, true);
                        neighbor_info.fill_local_data(neighbor_info.derivatives, true);
                        this->face_matrix(intint, intext, extint, extext, subface_info, neighbor_info);

                        sub_neighbor->face(neighbor_face_nr)->set_user_flag ();
                        for (unsigned int i=0;i<intint.size();++i)
                        {
                            for (unsigned int j=0;j<intint[i].matrix.n_rows();++j)
                            for (unsigned int k=0;k<intint[i].matrix.n_cols();++k)
                            {
                                const unsigned int block_row = intint[i].row;
                                const unsigned int block_col = intint[i].column;
                                const unsigned int
                                jcell = this->block_info.local.local_to_global(block_row, j);
                                const unsigned int
                                kcell = this->block_info.local.local_to_global(block_col, k);

                                if (fabs(intint[i].matrix(j,k)) > 1.e-12)
                                MATRIX.add(INDEX(subface_info.indices[jcell], level),
                                        INDEX(subface_info.indices[kcell], level),
                                        intint[i].matrix(j,k));
                                if (fabs(intext[i].matrix(j,k)) > 1.e-12)
                                MATRIX.add(INDEX(subface_info.indices[jcell], level),
                                        INDEX(neighbor_info.indices[kcell], level),
                                        intext[i].matrix(j,k));
                                if (fabs(extint[i].matrix(j,k)) > 1.e-12)
                                MATRIX.add(INDEX(neighbor_info.indices[jcell], level),
                                        INDEX(subface_info.indices[kcell], level),
                                        extint[i].matrix(j,k));
                                if (fabs(extext[i].matrix(j,k)) > 1.e-12)
                                MATRIX.add(INDEX(neighbor_info.indices[jcell], level),
                                        INDEX(neighbor_info.indices[kcell], level),
                                        extext[i].matrix(j,k));
                            }
                        }

                    }
                } else {
                    // Regular interior face
                    face_info.reinit(c, face, face_nr);
                    neighbor_info.reinit(neighbor,
                            neighbor->face(neighbor_face_nr),
                            neighbor_face_nr);

                    for (unsigned int i=0;i<intint.size();++i)
                    {
                        intint[i].matrix = 0.;
                        intext[i].matrix = 0.;
                        extint[i].matrix = 0.;
                        extext[i].matrix = 0.;
                    }

                    // Fill local data vectors
                    face_info.fill_local_data(face_info.values, true);
                    face_info.fill_local_data(face_info.derivatives, true);
                    neighbor_info.fill_local_data(neighbor_info.values, true);
                    neighbor_info.fill_local_data(neighbor_info.derivatives, true);
                    this->face_matrix(intint, intext, extint, extext, face_info, neighbor_info);


                    for (unsigned int i=0;i<intint.size();++i)
                    {
                        for (unsigned int j=0;j<intint[i].matrix.n_rows();++j)
                        for (unsigned int k=0;k<intint[i].matrix.n_cols();++k)
                        {
                            const unsigned int block_row = intint[i].row;
                            const unsigned int block_col = intint[i].column;
                            const unsigned int
                            jcell = this->block_info.local.local_to_global(block_row, j);
                            const unsigned int
                            kcell = this->block_info.local.local_to_global(block_col, k);

                            if (fabs(intint[i].matrix(j,k)) > 1.e-12)
                            MATRIX.add(INDEX(face_info.indices[jcell], level),
                                    INDEX(face_info.indices[kcell], level),
                                    intint[i].matrix(j,k));
                            if (fabs(intext[i].matrix(j,k)) > 1.e-12)
                            MATRIX.add(INDEX(face_info.indices[jcell], level),
                                    INDEX(neighbor_info.indices[kcell], level),
                                    intext[i].matrix(j,k));
                            if (fabs(extint[i].matrix(j,k)) > 1.e-12)
                            MATRIX.add(INDEX(neighbor_info.indices[jcell], level),
                                    INDEX(face_info.indices[kcell], level),
                                    extint[i].matrix(j,k));
                            if (fabs(extext[i].matrix(j,k)) > 1.e-12)
                            MATRIX.add(INDEX(neighbor_info.indices[jcell], level),
                                    INDEX(neighbor_info.indices[kcell], level),
                                    extext[i].matrix(j,k));
                        }
                    }
                }
            }

        } // faces
    } // cells


    #endif

    //////////////////////////////////////////////////////////////////////
    // Condense matrix
    //////////////////////////////////////////////////////////////////////
    this->hanging_node_constraints.condense(this->matrix);

    FcstUtilities::log.pop();

}

#undef MATRIX
#undef INDEX
#else
template<int dim>
void BlockMatrixApplication<dim>::PETSc_assemble(const FEVectors& src) 
{
    typedef typename DoFHandler<dim>::active_cell_iterator active_cell_iterator;

    active_cell_iterator begin = this->dof->begin_active();
    active_cell_iterator end = this->dof->end();

    FcstUtilities::log.push("Assembly");
    Assert(quadrature_assemble_cell, ExcNotInitialized());
    const Quadrature<dim>& cell_quadrature = *quadrature_assemble_cell;
    const Quadrature<dim - 1>& face_quadrature = *quadrature_assemble_face;

    // Objects for cell and face data handed down to the local routines.
    FEValues<dim>* fevalues = 0;
    FEFaceValues<dim>* fefacevalues = 0;
    
    typename DoFApplication<dim>::CellInfo cell_info(src, this->block_info);
    typename DoFApplication<dim>::FaceInfo bdry_info(src, this->block_info);

    bdry_info.initialize(fefacevalues, *this->element, *this->mapping, face_quadrature,
                         UpdateFlags(update_q_points | update_values | update_gradients
                                     | update_normal_vectors | update_JxW_values));

    cell_info.initialize(fevalues, *this->element, *this->mapping, cell_quadrature,
                         UpdateFlags(update_q_points | update_values | update_gradients
                                     | update_JxW_values));

    // Initialize local data
    MatrixVector intint;

    for (unsigned int i = 0; i < this->block_info.global.size(); ++i){
        for (unsigned int j = 0; j < this->block_info.global.size(); ++j) {

            if (this->cell_couplings(i, j) == DoFTools::none
                    && this->flux_couplings(i, j) == DoFTools::none)
                continue;

            MatrixBlock<FullMatrix<double> > block(i, j);
            intint.push_back(block);
            intint[intint.size()-1].matrix.reinit (this->block_info.local.block_size(i),
                                                   this->block_info.local.block_size(j));
        }
    }

    // Auxiliary objects used in loop over cells to store row and column indices:
    std::vector<unsigned int> row_indices;
    std::vector<unsigned int> col_indices;

    this->tr->clear_user_flags();

    // Start loops over cells
    typename DoFHandler<dim>::active_cell_iterator c;
    
    // First loop over cells:
    for (c = begin ; c != end ; ++c)
    {
        if(c->subdomain_id() == this->this_mpi_process)
        {
            for (unsigned int i=0;i<intint.size();++i)
                intint[i].matrix = 0.;

            // Initialize local structures
            cell_info.reinit(c);

            // Fill local data vectors
            cell_info.fill_local_data(cell_info.values, true);
            cell_info.fill_local_data(cell_info.derivatives, true);
            cell_matrix(intint, cell_info);

            //Local matrices are created assembled above and put in block in below loop
            for (unsigned int i=0;i<intint.size();++i)
            {
                const unsigned int n = intint[i].matrix.n_rows();
                const unsigned int m = intint[i].matrix.n_cols();
                row_indices.resize(n);
                col_indices.resize(m);

                //Calculate global indices
                for (unsigned int j=0;j < n;++j){

                    row_indices[j] = cell_info.indices
                            [this->block_info.local.local_to_global(intint[i].row, j)];
                }
                for(unsigned int k=0;k< m ;++k){
                    col_indices[k] = cell_info.indices
                            [this->block_info.local.local_to_global(intint[i].column, k)];
                }

                this->hanging_node_constraints.distribute_local_to_global
                (intint[i].matrix ,row_indices, col_indices, matrix);               
            }
        }
    }

    this->post_cell_assemble();
    
#if deal_II_dimension > 1

    // Now, loop over cells and then, over boundaries:
    if (this->boundary_fluxes)
    {
        for (c = begin; c != end; ++c)
        {
            if(c->subdomain_id() == this->this_mpi_process)
            {
                for (unsigned int face_nr=0; face_nr < GeometryInfo<dim>::faces_per_cell; ++ face_nr)
                {
                    typename DoFHandler<dim>::face_iterator face = c->face(face_nr);
                    
                    // Avoid computing this face
                    // a second time.
                    if (face->user_flag_set ())
                        continue;
                    
                    if (c->at_boundary(face_nr))
                    {
                        for (unsigned int i=0;i<intint.size();++i)
                            intint[i].matrix = 0.;
                        
                        bdry_info.reinit(c, face, face_nr);
                        // Boundary values
                        // Fill local data vectors
                        bdry_info.fill_local_data(bdry_info.values, true);
                        bdry_info.fill_local_data(bdry_info.derivatives, true);
                        this->bdry_matrix(intint, bdry_info);
                        
                        for (unsigned int i=0;i<intint.size();++i)
                        {   
                            const unsigned int n = intint[i].matrix.n_rows();
                            const unsigned int m = intint[i].matrix.n_cols();
                            row_indices.resize(n);
                            col_indices.resize(m);                            
                            
                            //Calculate global indices
                            for (unsigned int j=0;j < n;++j){                                
                                row_indices[j] = bdry_info.indices
                                [this->block_info.local.local_to_global(intint[i].row, j)];
                            }                          
                            
                            for(unsigned int k=0;k< m ;++k){
                                col_indices[k] = bdry_info.indices
                                [this->block_info.local.local_to_global(intint[i].column, k)];
                                
                            }
                            // Distribute (See deal.II ConstraintMatrix info):                            
                            this->hanging_node_constraints.distribute_local_to_global
                            (intint[i].matrix ,row_indices, col_indices, matrix);

                        }
                    }
                }
            }
        }
    }
    else if(this->interior_fluxes)
    {
        throw std::runtime_error("Assembly of interior_fluxes is not implemented for PETSc implementation.");
        
    }
    #endif

    // Finish assembling the global matrix:
    matrix.compress(VectorOperation::add);
    /*
    // Fix eigenvalues of constraints (See deal.II ConstraintMatrix info):
    for (unsigned int r=0;r<matrix.m();++r)
        if (this->hanging_node_constraints.is_constrained(r))
            matrix.set(r,r,1.);
     */   
    FcstUtilities::log.pop();
}

#endif //Not OPENFCST_WITH_PETSC


//------------------------------
template<int dim>
void BlockMatrixApplication<dim>::assemble_numerically(const FEVectors& src, const double delta) 
{
    #ifdef OPENFCST_WITH_PETSC
        throw std::runtime_error("BlockMatrixApplication<dim>::assemble_numerically are not implemented for PETSc implementation.");
    #else
        this->matrix = 0.0;

        double l2_norm_residual;

        // Find the index where the solution vector is stored:
        unsigned int ind = src.find_vector(this->data->get_solution_vector_name(this->data->get_nonlinear_solver()));

        // Vector where the pertrubed solution is stored:
        FEVector solution_delta;
        solution_delta = src.vector(ind);

        // FEVectors object that is passed to residual and that includes the perturbed solution:
        FEVectors solution_copy;
        solution_copy.add_vector(solution_delta, this->data->get_solution_vector_name(this->data->get_nonlinear_solver()));

        // Create a vector where the residual is stored and initialize
        // it with the residual:
        FEVector residual;
        residual.reinit(src.vector(ind));
        l2_norm_residual = this->residual(residual, solution_copy, false);

        // Create and initialize a vector where the residual obtained
        // using the perturbed solution will be stored:
        FEVector residual_copy;
        residual_copy.reinit(src.vector(ind));

        // Loop over DOFs (columns)
        for (unsigned int j = 0; j < residual.size(); j++) {
            // Perturb solution:
            solution_delta[j] += delta;

            // Compute residual
            l2_norm_residual = this->residual(residual_copy, solution_copy, false);

            // Loop over DOFs (rows) and copy entries in the Jacobian as required:
            for (unsigned int i = 0; i < residual.size(); i++) {
                this->matrix.set(i, j, (residual_copy(i) - residual(i)) / delta);
            }
            // Remove perturbation:
            solution_delta[j] -= delta;
        }

        //////////////////////////////////////////////////////////////////////
        // Condense matrix, i.e. correct for hanging nodes
        //////////////////////////////////////////////////////////////////////
        this->hanging_node_constraints.condense(this->matrix);

        //////////////////////////////////////////////////////////////////////
        // Apply Dirichlet boundary conditions
        //////////////////////////////////////////////////////////////////////

        // Note: In solve, if we want to apply BC to matrix, rhs and sol use:
        MatrixTools::apply_boundary_values(boundary_values, this->matrix,
                solution_delta, solution_delta);

        /*
        // ========================================
        // USE FOR DEBUGGING ONLY:
        // ========================================

        FcstUtilities::log<<"== PRINTING NUMERICAL MATRIX =="<<std::endl;
        std::ofstream file;
        file.open("Numerical_matrix.dat");
        this->matrix.print_formatted(file, 5, true, 0, " 0 ");
        file.close();
        FcstUtilities::log<<"============================"<<std::endl;
        // ========================================
        */

        // ========================================
        // USE FOR DEBUGGING ONLY:
        // ========================================
        /*
        //Analytical sensitivities:
        this->matrix=0.0;
        FEVectors sol;
        unsigned int ind = src.find_vector("Newton iterate");
        sol.add_vector(src.vector(ind),"Newton iterate");
        this->assemble(sol); //Note the second component of

        FEVector solution_delta;
        solution_delta = src.vector(ind);

        std::map<unsigned int, double> boundary_values;
        this->dirichlet_bc(boundary_values);

        MatrixTools::apply_boundary_values(boundary_values,
        this->matrix,
        solution_delta,
        solution_delta);

        FcstUtilities::log<<"== PRINTING ANALYTICAL MATRIX =="<<std::endl;

        std::ofstream file;
        file.open("Analytical_matrix.dat");
        this->matrix.print_formatted(file, 5, true, 0, " 0 ");
        file.close();
        FcstUtilities::log<<"============================"<<std::endl;
        */
        // ========================================
    #endif

}

//------------------------------
template<int dim>
void BlockMatrixApplication<dim>::residual_constraints(FEVector& dst) const {

    //////////////////////////////////////////////////////////////////////
    // Condense rhs
    //////////////////////////////////////////////////////////////////////
    this->hanging_node_constraints.condense(dst);

    //////////////////////////////////////////////////////////////////////
    // Apply Dirichlet boundary conditions
    //////////////////////////////////////////////////////////////////////
    for(auto m: boundary_values)
            dst(m.first) = m.second;
}

//------------------------------
#ifdef OPENFCST_WITH_PETSC
template<int dim>
void BlockMatrixApplication<dim>::residual_constraints(PETScWrappers::MPI::Vector& dst, const std::vector<unsigned int>& idx_list) const {

	//Only modify the dofs in local process (idx_list)
    for(auto i: idx_list)
        if (boundary_values.count(i))
            dst(i) = boundary_values.at(i);
}
#endif

//------------------------------
template<int dim>
void BlockMatrixApplication<dim>::post_cell_assemble() {

}

//------------------------------
template<int dim>
void BlockMatrixApplication<dim>::cell_matrix(MatrixVector& matrices,
        const typename DoFApplication<dim>::CellInfo&) {
    FcstUtilities::log << "Void function " << __FUNCTION__ << " called" << std::endl;
    for (unsigned int i = 0; i < matrices.size(); ++i)
        FcstUtilities::log << "Matrix[" << i << "](" << matrices[i].row << ','

                << matrices[i].column << ")\t" << matrices[i].matrix.n_rows()
                << 'x' << matrices[i].matrix.n_cols() << std::endl;
}

//------------------------------
template<int dim>
void BlockMatrixApplication<dim>::bdry_matrix(MatrixVector& matrices,
        const typename DoFApplication<dim>::FaceInfo&) {
    FcstUtilities::log << "Void function " << __FUNCTION__ << " called" << std::endl;
    for (unsigned int i = 0; i < matrices.size(); ++i)
        FcstUtilities::log << "Matrix[" << i << "](" << matrices[i].row << ','
                << matrices[i].column << ")\t" << matrices[i].matrix.n_rows()
                << 'x' << matrices[i].matrix.n_cols() << std::endl;
}

//------------------------------
template<int dim>
void BlockMatrixApplication<dim>::face_matrix(MatrixVector& matrices11,
        MatrixVector& matrices12, MatrixVector& matrices21,
        MatrixVector& matrices22, const typename DoFApplication<dim>::FaceInfo&,
        const typename DoFApplication<dim>::FaceInfo&) {

    FcstUtilities::log << "Void function " << __FUNCTION__ << " called" << std::endl;
    Assert(matrices11.size() == matrices12.size(),
            ExcDimensionMismatch(matrices11.size(), matrices12.size()));
    Assert(matrices11.size() == matrices21.size(),
            ExcDimensionMismatch(matrices11.size(), matrices21.size()));
    Assert(matrices11.size() == matrices22.size(),
            ExcDimensionMismatch(matrices11.size(), matrices22.size()));

    for (unsigned int i = 0; i < matrices11.size(); ++i) {
        Assert(matrices11[i].matrix.n_rows() == matrices12[i].matrix.n_rows(),
                ExcDimensionMismatch(matrices11[i].matrix.n_rows(), matrices12[i].matrix.n_rows()));
        Assert(matrices11[i].matrix.n_cols() == matrices12[i].matrix.n_cols(),
                ExcDimensionMismatch(matrices11[i].matrix.n_cols(), matrices12[i].matrix.n_cols()));
        Assert(matrices11[i].matrix.n_rows() == matrices21[i].matrix.n_rows(),
                ExcDimensionMismatch(matrices11[i].matrix.n_rows(), matrices21[i].matrix.n_rows()));
        Assert(matrices11[i].matrix.n_cols() == matrices21[i].matrix.n_cols(),
                ExcDimensionMismatch(matrices11[i].matrix.n_cols(), matrices21[i].matrix.n_cols()));
        Assert(matrices11[i].matrix.n_rows() == matrices22[i].matrix.n_rows(),
                ExcDimensionMismatch(matrices11[i].matrix.n_rows(), matrices22[i].matrix.n_rows()));
        Assert(matrices11[i].matrix.n_cols() == matrices22[i].matrix.n_cols(),
                ExcDimensionMismatch(matrices11[i].matrix.n_cols(), matrices22[i].matrix.n_cols()));


        FcstUtilities::log << "Matrix[" << i << "](" << matrices11[i].row << ','
                << matrices11[i].column << ")\t"
                << matrices11[i].matrix.n_rows() << 'x'
                << matrices11[i].matrix.n_cols() << std::endl;
    }
}

//--------------------
template<int dim>
void BlockMatrixApplication<dim>::dirichlet_bc(
        std::map<unsigned int, double>& /*boundary_values*/) const {
    const std::type_info& info = typeid(*this);
    FcstUtilities::log << "Pure function " << __FUNCTION__ << " called in Class "
            << info.name() << std::endl;

}

//-------------------------
template<int dim>
template<typename Vector >
void BlockMatrixApplication<dim>::print_matrix_and_rhs(Vector& sys_rhs) const
{
    if (print_debug)
    {
        std::ofstream file;
        unsigned int iterNum = 0;
        std::string fileName = "";
        bool fileExists      = true;
        
        while( fileExists )
        {
            iterNum++;
            fileName = "Solver_globalMatrix_Iter"+std::to_string(iterNum)+".dat";
            if( !boost::filesystem::exists( fileName ) )
                fileExists = false;
        }
        file.open(fileName);    
        
        #ifdef OPENFCST_WITH_PETSC
        for(unsigned int i = 0; i < this->matrix.m(); ++i)
        {
            for(unsigned int j = 0; j < this->matrix.m(); ++j)
            {
                file << std::setprecision(7) << std::scientific << this->matrix(i,j) << "\t";
            }
            file << std::endl;
        }
        #else
        this->matrix.print_formatted(file, 7, true, 14, "0");    
        #endif
        
        file.close();
        
        //RHS
        fileExists = true;
        while( fileExists )
        {
            iterNum++;
            fileName = "Solver_globalRHS_Iter"+std::to_string(iterNum)+".dat";
            if( !boost::filesystem::exists( fileName ) )
                fileExists = false;
        }
        file.open(fileName);
        for(unsigned int i = 0; i < sys_rhs.size(); ++i)
            file << std::setprecision(7) << std::scientific << sys_rhs(i) << std::endl;
        file.close();
    }
}

//------------------------------
//------------------------------
//------------------------------
template class BlockMatrixApplication<deal_II_dimension> ;
