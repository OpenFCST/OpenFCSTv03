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
// - Class: dof_application.cc
// - Description: This class implements:
//                - triangulation
//                - finite element spaces
//                - dof handler
//                - distribution of dofs over triangulation
//                - residual of nonlinear system of equations
//                - estimations of cell-wise errors for adaptive refinement
//                - and more
// - Developers: Guido Kanschat,     Texas A&M University
//               Valentin N. Zingan, University of Alberta
//               Marc Secanell,      University of Alberta
//               Peter Dobson,       University of Alberta
//               Aslan Kosakian,     University of Alberta
//
// ----------------------------------------------------------------------------

#include <application_core/dof_application.h>

template<int dim>
DoFApplication<dim>::DoFApplication(boost::shared_ptr<ApplicationData> data)
:
ApplicationBase(data),
verbosity(0),
tr(new Triangulation<dim>),
n_ref(100),
boundary_fluxes(true),
interior_fluxes(true),
system_management(this->block_info,
                  cell_couplings  ,
                  flux_couplings  ),
output_materials_and_levels(true),
output_actual_degree(true),
print_solution(false),
print_postprocessing(false)
{
    FcstUtilities::log << "->DoF";
    boost::shared_ptr<DoFHandler<dim> >
    newdof(new DoFHandler<dim>(*tr));
    dof = newdof;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
DoFApplication<dim>::DoFApplication()
:
ApplicationBase(boost::shared_ptr<ApplicationData>()),
tr(new Triangulation<dim>),
boundary_fluxes(true),
interior_fluxes(true),
system_management(this->block_info,
                  cell_couplings  ,
                  flux_couplings  ),
output_materials_and_levels(true),
output_actual_degree(true),
print_solution(false),
print_postprocessing(false)
{
    FcstUtilities::log << "->DoF";
    boost::shared_ptr<DoFHandler<dim> >
    newdof(new DoFHandler<dim>(*tr));
    dof = newdof;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
DoFApplication<dim>::DoFApplication(DoFApplication<dim>& other,
                                    bool triangulation_only)
:
ApplicationBase(other),
tr(other.tr),
system_management(this->block_info,
                  cell_couplings  ,
                  flux_couplings  ),
output_actual_degree(true),
print_solution(false),
print_postprocessing(false)
{
    tr = other.tr;
    if (triangulation_only)
    {
        boost::shared_ptr<DoFHandler<dim> >
        newdof(new DoFHandler<dim>(*tr));
        dof = newdof;
    }
    else
    {
        dof = other.dof;
    }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
DoFApplication<dim>::~DoFApplication()
{ 
    mesh_generator.reset();
}
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
void
DoFApplication<dim>::declare_parameters(ParameterHandler& param)
{
    // Declare parameters in mesh generator:
    FuelCellShop::Geometry::GridBase<dim>::declare_GridGenerator_parameters(param);
    
    //==
    param.enter_subsection("Grid generation");
    {
        param.declare_entry("Initial refinement", 
                            "1", 
                            Patterns::Integer(),
                            "Enter the number of times each cell in the original mesh should be divided into four cells.");
        /*
        param.declare_entry("Sort by component", 
                            "true", 
                            Patterns::Bool()
                            "Organize the degree of freedom numbering for the mesh by component");
        */
        param.declare_entry("Sort Cuthill-McKee", 
                            "false", 
                            Patterns::Bool(),
                            "Organize the degree of freedom numbering for the mesh using the Cuthill-McKee algorithm");        
    }
    param.leave_subsection();
    param.enter_subsection("Adaptive refinement");
    {
        /*
         * Here, we allow arbitrarily 100 refinements so cells will practically never be marked as refinement-saturated
         * (see @remesh function of this class). This number will be used only by applications that do not use adaptive
         * refinement class and want refinement saturation to not happen.
         */
        param.declare_entry ("Number of Refinements",
                             "100",
                             Patterns::Integer(),
                             "This parameter is used to define the number of times the mesh will be adaptively refined. \n"
                             "The minimum value is one, i.e., only the original mesh is solved. At each adaptive refinement level, \n"
                             "either all the cell (global) or 30% of the cells with largest error (computed using an error estimator) are split into four. \n"
                             "The process is repeated at each refinement level. ");        
        param.declare_entry("Refinement", 
                            "global",
                            Patterns::Selection("global|adaptive"),
                            "Global refinement will break each cell into four cells, \n"
                            "Adaptive refinement will only refine the percentace of cells in Refinement threshold with largest error.");
        param.declare_entry("Refinement threshold", 
                            "0.3", 
                            Patterns::Double(),
                            "For adaptive refinement, the percentage of cells with largest error that should be refined");
        param.declare_entry("Coarsening threshold", 
                             "0", 
                            Patterns::Double(),
                            "For adaptive refinement, the percentage of cells with smallest error that should be coarsened");
        param.declare_entry("Curved boundary ID", 
                            "0", 
                            Patterns::Integer());
    }
    param.leave_subsection();
    //==
    param.enter_subsection("Discretization");
    {
        param.declare_entry("Element", "FE_Q(1)", Patterns::Anything(),
                            "The finite element used for discretization. "
                            "Any input for FETools::get_element_by_name is possible."
                            "FE_Q(1) is a Lagrange element of order 1");
        param.declare_entry("Mapping degree",
                            "1",
                            Patterns::Integer(1),
                            "Degree used for polynomial mapping of arbitrary order");
        param.declare_entry("Boundary fluxes",
                            "false",
                            Patterns::Bool(),
                            "Do you have any Neumann boundary conditions?");
        param.declare_entry("Interior fluxes",
                            "false",
                            Patterns::Bool(),
                            "Do you have any interior fluxes (usually applies to DG)?");
        
        param.enter_subsection("Residual");
        {
            param.declare_entry("Quadrature cell", "-1", Patterns::Integer(),
                                "Gauss formula for integrating the residual on cells."
                                " See above for more details.");
            param.declare_entry("Quadrature bdry", "-1", Patterns::Integer(),
                                "Gauss formula used for integrating the residual "
                                "on boundary faces. See above for more details.");
            param.declare_entry("Quadrature face", "-1", Patterns::Integer(),
                                "Gauss formula for integrating the residual on faces."
                                " See above for more details.");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    //== Declare parameters in system management:
    this->system_management.declare_parameters(param);
    
    //==
    param.enter_subsection("Initial Solution");
    {
        // Initial solution parameters:
        param.declare_entry("Read in initial solution from file",
                            "false",
                            Patterns::Bool(),
                            "Check whether a stored solution should be read from file and applied to the grid");
        param.declare_entry("Use pre-defined initial solution",
                            "true",
                            Patterns::Bool(),
                            "Use a developer pre-defined routine to setup initial solution instead"
                            "of using the piece-wise function defined using Equations>>Initial Data."
                            "This might be beneficial in some cases as a more appropriate initial solution might be used.");
        param.declare_entry("Output solution for transfer",
                            "false",
                            Patterns::Bool(),
                            "Check whether a solution on a refined grid should be output to file on a coarse mesh");
        param.declare_entry("Output initial solution",
                            "false",
                            Patterns::Bool(),
                            "Set flag to true if you want to output the initial solution to file");
        param.declare_entry("Initial solution output filename",
                            "initial_solution",
                            Patterns::Anything(),
                            "File where the initial solution will be output");
    }
    param.leave_subsection();

    //==
    param.enter_subsection ("Output");
    {
        param.enter_subsection ("Grid");
        {
            GridOut::declare_parameters(param);
        }
        param.leave_subsection();

        //--
        param.enter_subsection("Data");
        {
            param.declare_entry("Output actual degree",
                                "true",
                                Patterns::Bool());
            param.declare_entry("Print solution",
                                "false",
                                Patterns::Bool());
            param.declare_entry("Print postprocessing",
                                "false",
                                Patterns::Bool());
            param.declare_entry("Solution printing indices",
                                """",
                                Patterns::List( Patterns::Integer(0) ),
                                "");
            param.declare_entry("Postprocessing printing indices",
                                """",
                                Patterns::List( Patterns::Integer(0) ),
                                "");
            param.declare_entry("Print blocks instead of indices",
                                "false",
                                Patterns::Bool());
            DataOutInterface<dim>::declare_parameters(param);
            // Default output format for openFCST is vtu not gnuplot
            param.set("Output format","vtu");
            param.declare_entry("Debug output of system matrices and right hand sides",
                                "false",
                                Patterns::Bool(),
                                "WARNING: lots of output files."
                                );            
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
void
DoFApplication<dim>::initialize (ParameterHandler& param)
{
    _initialize(param);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
void
DoFApplication<dim>::_initialize (ParameterHandler& param)
{
    // --- then read additional information:
    param.enter_subsection ("Grid generation");
    {        
        initial_refinement = param.get_integer("Initial refinement");             
        sort_cuthill = param.get_bool("Sort Cuthill-McKee");
    }
    param.leave_subsection();
    param.enter_subsection("Adaptive refinement");
    {
        refinement = param.get("Refinement");
        // For safe initialization of n_ref in case adaptive refinement is not used.
        this->n_ref = param.get_integer("Number of Refinements");
        refinement_threshold = param.get_double("Refinement threshold");
        coarsening_threshold = param.get_double("Coarsening threshold");        
        curved_bdry_id = param.get_integer("Curved boundary ID");
    }
    param.leave_subsection();
    //--
    param.enter_subsection("Discretization");
    {
        boost::shared_ptr<FiniteElement<dim> >
        newelement(FETools::get_fe_from_name<dim>(param.get("Element")));
        element = newelement;
        FcstUtilities::log << element->get_name() << std::endl;

        mapping_degree = param.get_integer("Mapping degree");

        if(mapping == 0)
        {
            if(mapping_degree == 1)
            {
                boost::shared_ptr<Mapping<dim> > newmapping(new MappingQ1<dim>());
                mapping = newmapping;
            }
            else
            {
                boost::shared_ptr<Mapping<dim> > newmapping(new MappingQ<dim>(mapping_degree));
                mapping = newmapping;
            }
        }

        boundary_fluxes = param.get_bool("Boundary fluxes");
        interior_fluxes = param.get_bool("Interior fluxes");

        //--
        param.enter_subsection("Residual");
        {
            if (true)
            {
                int n = param.get_integer("Quadrature cell");
                if (n<=0) n = this->element->degree - n;
                this->quadrature_residual_cell = QGauss<dim>(n);
            }
            if (true)
            {
                int n = param.get_integer("Quadrature bdry");
                if (n<=0) n = this->element->degree - n;
                this->quadrature_residual_bdry = QGauss<dim-1>(n);
            }
            if (true)
            {
                int n = param.get_integer("Quadrature face");
                if (n<=0) n = this->element->degree - n;
                this->quadrature_residual_face = QGauss<dim-1>(n);
            }
        }
        param.leave_subsection();
    }
    param.leave_subsection();

    // -- Initialize parameters in system management:
    this->system_management.initialize(param);
    
    // and make sure that the number of input finite elements matches with the problem:
    Assert (this->element->n_blocks() == this->system_management.get_number_of_solution_names(),
            ExcDimensionMismatch(this->element->n_blocks(), this->system_management.get_number_of_solution_names()));

    //==
    param.enter_subsection("Initial Solution");
    {
        read_in_initial_solution = param.get_bool("Read in initial solution from file");
        use_predefined_solution = param.get_bool("Use pre-defined initial solution");
        filename_initial_sol = param.get("Initial solution output filename");
        output_initial_sol = param.get_bool("Output initial solution");
        output_coarse_solution = param.get_bool("Output solution for transfer");
    }
    param.leave_subsection();

    // -- Parse output parameters
    param.enter_subsection ("Output");
    {
        param.enter_subsection ("Grid");
        {
            g_out.parse_parameters(param);
        }
        param.leave_subsection();

        param.enter_subsection("Data");
        {
            output_actual_degree = param.get_bool("Output actual degree");
            print_solution       = param.get_bool("Print solution");
            print_postprocessing = param.get_bool("Print postprocessing");
            print_blocks_instead_of_indices = param.get_bool("Print blocks instead of indices");

            if( !param.get("Solution printing indices").empty() )
            {
                solution_printing_indices = FcstUtilities::string_to_number<unsigned int>( Utilities::split_string_list( param.get("Solution printing indices") ) );
            }

            if( !param.get("Postprocessing printing indices").empty() )
                {
                    postprocessing_printing_indices = FcstUtilities::string_to_number<unsigned int>( Utilities::split_string_list( param.get("Postprocessing printing indices") ) );
                }
                d_out.parse_parameters(param);
            output_matrices_and_rhs = param.get_bool("Debug output of system matrices and right hand sides");            
        }
        param.leave_subsection();
    }
    param.leave_subsection();

    // --- first read the mesh:
     this->initialize_triangulation(param);

     //--- then distribute dofs:
     this->remesh_dofs();
     
     // This controls data output (search for "already_run" in this file).
     this->data->enter_flag("already_run",false);
}

//===========================================================================================
template <int dim>
void
DoFApplication<dim>::initialize_triangulation(ParameterHandler& param)
{
    // Create a grid object:
    mesh_generator = FuelCellShop::Geometry::GridBase<dim>::create_GridGenerator(param);

    // Check to see whether we are reading the mesh_generator from a file
    // Must read grid if we are reading a stored solution
    bool read_grid_from_file = mesh_generator->read_from_file;

    // Read in new mesh from file:
    if (read_grid_from_file && !this->read_in_initial_solution){
        FcstUtilities::log << "Reading new grid..." << std::endl;
        if(curved_boundary)
           mesh_generator->generate_grid_with_curved_boundaries( *this->tr,
                                                                  curved_bdry_id,
                                                                  curved_boundary );
        else
           mesh_generator->generate_grid(*this->tr);
    }
    // or load mesh associated with previous solution:
    else if (this->read_in_initial_solution) {

        std::ifstream grid_file("._mesh_to_be_transfered.msh");

        if(!grid_file.fail())
        {
            FcstUtilities::log << "Reading grid from previous solution" << std::endl;
            GridIn<dim> grid_in;
            grid_in.attach_triangulation(*this->tr);
            grid_in.read_msh(grid_file);
            grid_file.close();
            read_grid_from_file = true;
        }
        else
        {
            // if grid fails to read, generate standard grid, don't read initial solution
            FcstUtilities::log << "Grid file not available for reading" << std::endl;
            read_grid_from_file = false;
            this->read_in_initial_solution = false;
        }
    }

    // If a grid has not been read from file:
    if(!read_grid_from_file)
    {
        FcstUtilities::log << "Generating new grid..." << std::endl;
        mesh_generator->generate_grid(*this->tr);
    }
    
    if (mesh_generator->grid_in->field_data.size()>0)
        data->field_data=mesh_generator->grid_in->field_data;
    
    // If running in parallel, then subdivide the mesh by DOF.
#ifdef OPENFCST_WITH_PETSC
    GridTools::partition_triangulation (n_mpi_processes, *tr);
#endif

}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
void
DoFApplication<dim>::add_vector_for_transfer(FEVector* v)
{
    transfer_vectors.push_back(v);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
void
DoFApplication<dim>::delete_vector_for_transfer()
{
    transfer_vectors.pop_back();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <int dim>
void
DoFApplication<dim>::init_vector (FEVector& dst) const
{
    dst.reinit(block_info.global);
}

//===========================================================================================
template <int dim>
void
DoFApplication<dim>::initialize_solution (FEVector& initial_guess,
                                          std::shared_ptr<Function<dim> > initial_function)
{

    // Load old solution
    bool solution_read_success = false;

    // Two options
    // a) You need to create a completely new solution either from input file
    //  or from data file
    // b) You already have a solution and you simply want to modify B.C.

    // Option a)
    if( initial_guess.size() == 0 || initial_guess.all_zero() ) {

        // Initialize solution vector size:
        initial_guess.reinit(this->block_info.global);

        // Load old solution if possible
        // Read initial solution from file if read_in_initial_solution set to true:
        if (this->read_in_initial_solution)
        {
            FcstUtilities::log << "Reading initial guess from file" << std::endl;
            this->read_init_solution(initial_guess, solution_read_success);

            // If loading old solution successful, modify initial guess to match B.C.
            if (solution_read_success)
            {
                // --- here we only apply boundary values over the initial guess which has been already prepared ---
                FcstUtilities::log << "Correcting the initial guess" << std::endl;
                FuelCell::InitialAndBoundaryData::apply_piece_wise_constant_DirichletBCs( initial_guess,
                                                                                          *this->mapping,
                                                                                          *this->dof,
                                                                                          this->system_management,
                                                                                          this->component_boundaryID_value_maps );
            }
        }

        // If you do not have an old solution to load or you could not load it, create a new one:
        if (!this->read_in_initial_solution || !solution_read_success)
        {

            FcstUtilities::log << "Generating the initial guess" << std::endl;
            // --- we first form the initial guess ---
            if (initial_function && use_predefined_solution)
            {
                //If I have a function that I can use, then I use the function:
                VectorTools::interpolate (*this->dof, *initial_function.get(), initial_guess);

            }
            else
            {
                //Otherwise, use a piece_wise constant function
                FuelCell::InitialAndBoundaryData::make_piece_wise_constant_initial_data( initial_guess,
                                                                                         *this->mapping,
                                                                                         *this->dof,
                                                                                         this->system_management,
                                                                                         this->component_materialID_value_maps );
            }

            // --- and then apply boundary values over the initial guess we have just formed ---
            FuelCell::InitialAndBoundaryData::apply_piece_wise_constant_DirichletBCs( initial_guess,
                                                                                        *this->mapping,
                                                                                        *this->dof,
                                                                                        this->system_management,
                                                                                        this->component_boundaryID_value_maps );        
                                                                                    
        }
    }
    //
    else{
        // Check that the solution entered is the right size
        AssertThrow(initial_guess.size() == this->block_info.global.total_size(),
                    ExcDimensionMismatch(initial_guess.size(), this->block_info.global.total_size()) );

        // --- here we only apply boundary values over the initial guess which has been already prepared ---
        FcstUtilities::log << "Correcting the initial guess" << std::endl;
        FuelCell::InitialAndBoundaryData::apply_piece_wise_constant_DirichletBCs( initial_guess,
                                                                                  *this->mapping,
                                                                                  *this->dof,
                                                                                  this->system_management,
                                                                                  this->component_boundaryID_value_maps );
    }

    // Output the initial solution if desired:
    if(output_initial_sol)
    {
        FuelCell::ApplicationCore::FEVectors vectors;
        vectors.add_vector(initial_guess,
                           "Solution");
        if (!this->data->flag("already_run"))
        {
            data_out(filename_initial_sol, vectors);
            this->data->enter_flag("already_run",true);
        }
    }
}

//===========================================================================================

template <int dim>
void
DoFApplication<dim>::remesh()
{
    if (refinement == "global")
    {
        for (typename Triangulation<dim>::active_cell_iterator
            cell = tr->begin_active();
        cell != tr->end(); ++cell)
            cell->set_refine_flag();
    }
    else if (refinement == "adaptive")
    {
        GridRefinement::refine_and_coarsen_fixed_number(
            *tr, cell_errors, refinement_threshold, coarsening_threshold);            
    }
    else
    {
        Assert (false, ExcNotImplemented());
    }

    if (transfer_vectors.size() > 0)
    {
        // Initialize SolutionTransfer
        SolutionTransfer<dim, BlockVector<double> > transfer(*dof);
        
        std::vector<BlockVector<double>> aux1(transfer_vectors.size());
        
        // Initialize aux with FEVector-type objects
        for (int i=0; i<transfer_vectors.size(); i++)
        {
            this->init_vector(aux1[i]);
            aux1[i]=*transfer_vectors[i];
            FcstUtilities::log << "LOOP i=" << i << std::endl;
        }
        
        // Prepare mesh (in case some cells need to be refined due to cells having double refinement)
        tr->prepare_coarsening_and_refinement ();
        
        // Prepare solution
        transfer.prepare_for_coarsening_and_refinement(aux1);
        
        // Exectue refinement
        tr->execute_coarsening_and_refinement();
        
        // Apply to vectors
        remesh_dofs();
        
        // Transfer to new mesh (next few lines)
        std::vector<BlockVector<double>> aux2(transfer_vectors.size());
        
        for (int i=0; i<aux1.size(); i++)
        {
            this->init_vector(aux2[i]);
        }
        // aux1 in, aux2 out
        transfer.interpolate(aux1,aux2);
        for (int i=0; i<aux1.size(); i++)
        {
            this->init_vector(*transfer_vectors[i]);
            *transfer_vectors[i] = aux2[i];
        }            
    }
    else
    {
        tr->prepare_coarsening_and_refinement ();
        tr->execute_coarsening_and_refinement();
        remesh_dofs();
    }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
void
DoFApplication<dim>::grid_out(const std::string& name)
{

    const std::string suffix = g_out.default_suffix();
    const std::string filename = name + suffix;
    if (suffix == std::string(""))
    {
        FcstUtilities::log << "Not writing " << filename << std::endl;
        return;
    }

    FcstUtilities::log << "Gridfile:" << filename << std::endl;
    std::ofstream out(filename.c_str());

    g_out.write(*tr, out);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
void
DoFApplication<dim>::data_out(const std::string& basename,
                              const FEVectors& src)
{

    unsigned int index = src.find_vector("Solution");

    if (index != static_cast<unsigned int>(-1))
    {
        FuelCell::ApplicationCore::FEVector solution = src.vector(index);

        // Fill solution names with system_management!!!
        std::vector<std::string> solution_names;
        for (unsigned int i=0; i<solution.n_blocks(); ++i)
            solution_names.push_back("not provided");

        FuelCell::ApplicationCore::FEVector postprocessing;
        std::vector<std::string> postprocessing_names;

        data_out(basename,
                 solution,
                 solution_names,
                 postprocessing,
                 postprocessing_names);
    }
    else
        FcstUtilities::log<<"Warning: Solution vector could not be found in DoFApplication<dim>::data_out"<<std::endl;


}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
void
DoFApplication<dim>::data_out(const std::string&              basename,

                              const FEVector&                 solution,
                              const std::vector<std::string>& solution_names,

                              const FEVector&                 postprocessing,
                              const std::vector<std::string>& postprocessing_names)
{
    //-- Data_out should only be performed on a single CPU

    //-- Setup a variable to see if are running either in series or as process zero
    // in parallel
    bool zero_mpi_process = true;

    #ifdef OPENFCST_WITH_PETSC
    //Only Process 0 writes data out
    if(this->this_mpi_process != 0){
        zero_mpi_process = false;
    }
    #endif
    //--

    //-- Then, assemble data_out vector only for zero_mpi_process:
    if (zero_mpi_process)
    {
        d_out.attach_dof_handler(*this->dof);
        d_out.add_data_vector(solution,
                              solution_names,
                              DataOut<dim>::type_dof_data,
                              solution_interpretations);

        if( postprocessing.size() != 0 )
        {
            d_out.add_data_vector(postprocessing,
                                  postprocessing_names,
                                  DataOut<dim>::type_dof_data,
                                  postprocessing_interpretations);
        }

        output_materials.reinit(this->tr->n_active_cells());
        output_levels.reinit(this->tr->n_active_cells());

        // -- In parallel, I also want to output subdomain_id
        std::vector<unsigned int> partition_int (this->tr->n_active_cells());
        GridTools::get_subdomain_association (*this->tr, partition_int);
        const Vector<double> output_partitioning(partition_int.begin(), partition_int.end());

        unsigned int i = 0;
        for( typename Triangulation<dim>::active_cell_iterator cell  = this->tr->begin_active();
            cell != this->tr->end();  ++cell,   ++i )
            {
                output_materials(i) = cell->material_id();
                output_levels(i)    = cell->level();
            }

            if( output_materials_and_levels )
            {
                d_out.add_data_vector(output_materials,
                                      "cell_material_ids",
                                      DataOut<dim>::type_cell_data);
                d_out.add_data_vector(output_levels,
                                      "cell_levels",
                                      DataOut<dim>::type_cell_data);
                d_out.add_data_vector (output_partitioning, "subdomain_id");
            }

            if( output_actual_degree )
                d_out.build_patches(this->element->degree);
            else
                d_out.build_patches();

            std::ofstream output( (basename + d_out.default_suffix()).c_str() );
            d_out.write(output);
            d_out.clear();

            if( print_solution )
            {
                print("solution.dat",
                solution,
                solution_printing_indices);
            }

            if( print_postprocessing )
            {
                print("postprocessing.dat",
                postprocessing,
                postprocessing_printing_indices);
            }

    }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
void
DoFApplication<dim>::data_out(const std::string&                            basename,
                              const FEVector&                               solution,
                              const std::vector<std::string>&               solution_names,
                              const std::vector< DataPostprocessor<dim>* >& PostProcessing)
{
    bool zero_mpi_process = true;

    #ifdef OPENFCST_WITH_PETSC
    //Only Process 0 writes data out
    if(this->this_mpi_process != 0){
        zero_mpi_process = false;
    }
    #endif

    if (zero_mpi_process)
    {
        d_out.attach_dof_handler(*this->dof);
        d_out.add_data_vector(solution,
                              solution_names,
                              DataOut<dim>::type_dof_data,
                              solution_interpretations);

        for(unsigned int i = 0; i < PostProcessing.size(); ++i)
            d_out.add_data_vector(solution, *PostProcessing[i]);

        output_materials.reinit(this->tr->n_active_cells());
        output_levels.reinit(this->tr->n_active_cells());

        // -- In parallel, I also want to output subdomain_id
        std::vector<unsigned int> partition_int (this->tr->n_active_cells());
        GridTools::get_subdomain_association (*this->tr, partition_int);
        const Vector<double> output_partitioning(partition_int.begin(), partition_int.end());

        unsigned int i = 0;
        for( typename Triangulation<dim>::active_cell_iterator cell  = this->tr->begin_active();
            cell != this->tr->end(); ++cell, ++i )
            {
                output_materials(i) = cell->material_id();
                output_levels(i)    = cell->level();
            }

            if( output_materials_and_levels )
            {
                d_out.add_data_vector(output_materials,
                                      "cell_material_ids",
                                      DataOut<dim>::type_cell_data);
                d_out.add_data_vector(output_levels,
                                      "cell_levels",
                                      DataOut<dim>::type_cell_data);
                d_out.add_data_vector (output_partitioning, "subdomain_id");
            }

            if( output_actual_degree )
                d_out.build_patches(this->element->degree);
            else
                d_out.build_patches();

            std::ofstream output( (basename + d_out.default_suffix()).c_str() );
            d_out.write(output);
            d_out.clear();
            
            if( print_solution )
            {
                print("solution.dat",
                solution,
                solution_printing_indices);
            }
    }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
void
DoFApplication<dim>::print(const std::string&               basename,
                           const FEVector&                  src,
                           const std::vector<unsigned int>& src_indices) const
{
    std::ofstream myfile;
    myfile.open(basename, std::ofstream::out | std::ofstream::trunc);

    if( src_indices.size() == 0 )
    {
        for(unsigned int i = 0; i < src.size(); ++i)
            myfile << src(i) << "\n";
    }
    else
    {
        if( !print_blocks_instead_of_indices )
        {
            for(unsigned int i = 0; i < src_indices.size(); ++i)
            {
                AssertThrow( src_indices[i] < src.size() , ExcIndexRange( src_indices[i] , 0 , src.size() ) );
                myfile << src(src_indices[i]) << "\n";
            }
        }
        else
        {
            for(unsigned int i = 0; i < src_indices.size(); ++i)
            {
                AssertThrow( src_indices[i] < src.n_blocks() , ExcIndexRange( src_indices[i] , 0 , src.n_blocks() ) );
                for(unsigned int j = 0; j < src.block(src_indices[i]).size(); ++j)
                    myfile << src.block(src_indices[i])(j) << "\n";
            }
        }
    }

    myfile.close();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
void
DoFApplication<dim>::remesh_dofs()
{
    FcstUtilities::log << "Levels:" << dof->get_tria().n_levels() << std::endl
    << "Cells:" << dof->get_tria().n_active_cells() << std::endl;
    
    dof->clear();
    dof->distribute_dofs (*element);
    DoFRenumbering::subdomain_wise (*dof);
    
    block_info.initialize(*dof);
    
    FcstUtilities::log << "Dofs:" << dof->n_dofs();
    for (unsigned int i=0;i<block_info.global.size();++i)
        FcstUtilities::log << ((i==0) ? '(' : ',')
        << block_info.global.block_size(i);

    FcstUtilities::log << ')';

    sort_dofs(dof.get());

    // Make the list of constraints associated with hanging nodes
    hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints (*this->dof,
                                             hanging_node_constraints);
    hanging_node_constraints.close();


    // Provide indices for cells and faces to be used in error
    // estimation
    tr->clear_user_data();
    unsigned int k=0;
    for (typename Triangulation<dim>::active_cell_iterator i = tr->begin_active();
            i != tr->end(); ++i)
    {
        i->set_user_index(k++);
    }

    k=0;
    for (typename Triangulation<dim>::face_iterator i = tr->begin_face();
            i != tr->end_face(); ++i)
    {
        i->set_user_index(k++);
    }

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
void
DoFApplication<dim>::sort_dofs(DoFHandler<dim>* dof_handler) const
{
    if (sort_direction.norm_square() != 0.)
    {
        FcstUtilities::log << " sorting downstream";
        DoFRenumbering::downstream(*dof_handler, sort_direction);
    }

    if (sort_cuthill)
    {
        FcstUtilities::log << " sorting Cuthill-McKee";
        DoFRenumbering::Cuthill_McKee(*dof_handler);
    }

    FcstUtilities::log << " sorting by component";
    DoFRenumbering::component_wise(*dof_handler);
    FcstUtilities::log << std::endl;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
void
DoFApplication<dim>::transfer_solution_to_coarse_mesh(Triangulation<dim>& tr_coarse,
                                                      FEVector& coarse_solution,
                                                      FEVector& refined_solution)
{

    // Initialize a dof handler object with the triangulation
    DoFHandler<dim> dof_handler_coarse;
    const FiniteElement<dim>& fe_coarse(*this->element);
    dof_handler_coarse.initialize(tr_coarse,fe_coarse);
    dof_handler_coarse.distribute_dofs(fe_coarse);
    DoFRenumbering::subdomain_wise (dof_handler_coarse);

    sort_dofs(&dof_handler_coarse);

    // Obtain the common cells between the two triangulations
    typedef std::list<std::pair<typename DoFHandler<dim>::cell_iterator,
    typename DoFHandler<dim>::cell_iterator> > CellList;
    CellList finest_common_cells = GridTools::get_finest_common_cells (dof_handler_coarse, *this->dof);

    // Assign a vector for the dof values on the cell to transfer between solutions
    Vector<double> local_dof_values(fe_coarse.dofs_per_cell);

    // Loop over the common cells and transfer the solution values
    for (typename CellList::const_iterator  cell_pair = finest_common_cells.begin(); cell_pair != finest_common_cells.end(); ++cell_pair)
    {
        cell_pair->second->get_interpolated_dof_values(refined_solution, local_dof_values);
        cell_pair->first->set_dof_values_by_interpolation(local_dof_values, coarse_solution);
    }

    // Solution already created in coarse_solution, if we want to print it to a file, do so now:
    if (output_coarse_solution)
    {
        // Write coarse mesh to a file:
        std::ofstream grid_file("._mesh_to_be_transfered.msh");
        GridOut grid_out;
        GridOutFlags::Msh msh_flags;
        msh_flags.write_faces = true;
        grid_out.set_flags(msh_flags);
        grid_out.write_msh(tr_coarse,grid_file);
        grid_file.close();

        // Writ solution to a file:
        std::ofstream vector_out(".transfer_solution.FEVector");
        coarse_solution.block_write(vector_out);
        vector_out.close();

        FcstUtilities::log << "Coarse mesh and solution written to file" << std::endl;
    }

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
void
DoFApplication<dim>::read_init_solution(FEVector& dst, bool &good_solution) const
{
    std::ifstream vector_in(".transfer_solution.FEVector");
    if(!vector_in.fail())
    {
        FcstUtilities::log << "Initializing solution with available data" << std::endl;
        dst.block_read(vector_in);
        vector_in.close();
        good_solution = true;
    }
    else
    {
        FcstUtilities::log << "Cannot read file containing initial solution" << std::endl;
        good_solution = false;
    }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
void
DoFApplication<dim>::constrain_boundary(FEVector& v, bool homogeneous) const
{
    for (typename std::map<unsigned int, double>::const_iterator
        i = boundary_constraints.begin();
    i != boundary_constraints.end();
    ++i)
        v(i->first) = (homogeneous) ? 0. : (i->second);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
double
DoFApplication<dim>::evaluate(const FEVectors&)
{
    FcstUtilities::log.push("evaluate");
    unsigned int n = block_info.global.size();
    unsigned int mem_tr = tr->memory_consumption();
    unsigned int mem_dof = dof->memory_consumption();
    FcstUtilities::log << "DoFApplication "
    << memory_consumption() << " bytes" << std::endl
    << "Triangulation  " << mem_tr << " bytes" << std::endl
    << "DoFHandler     " << mem_dof << " bytes" << std::endl
    << "Dofs  \t" << dof->n_dofs() << std::endl
    << "Blocks\t" << block_info.global.size();
    for (unsigned int i=0;i<n;++i)
        FcstUtilities::log << '\t' << block_info.global.block_size(i);
    FcstUtilities::log << std::endl
    << "Start \t";
    for (unsigned int i=0;i<n;++i)
        FcstUtilities::log << '\t' << block_info.global.local_to_global(i,0);
    FcstUtilities::log << std::endl;
    FcstUtilities::log.pop();
    return 0.;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
unsigned int
DoFApplication<dim>::memory_consumption() const
{
    unsigned int result = sizeof (*this)
    + tr->memory_consumption()
    + dof->memory_consumption()
    + block_info.global.memory_consumption()
    + block_info.local.memory_consumption()
    + g_out.memory_consumption()
    + d_out.memory_consumption()
    - sizeof(block_info)
    - sizeof (d_out)
    - sizeof (g_out);
    return result;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
double
DoFApplication<dim>::estimate(const FEVectors& src)
{
    unsigned int index;

    // --- solution ---

    FuelCell::ApplicationCore::FEVector solution;
    solution.reinit( this->block_info.global );
    index = src.find_vector("Solution");
    solution = src.vector(index);

    // --- components to be used ---

    std::vector<bool> component_mask( this->element->n_blocks(), true );

    // --- cell errors ---

    this->cell_errors.reinit( this->tr->n_active_cells() );
    KellyErrorEstimator<dim>::estimate(*this->mapping,
                                       *this->dof,
                                       this->quadrature_residual_face,
                                       typename FunctionMap<dim>::type(),
                                       solution,
                                       this->cell_errors,
                                       component_mask);

    return 0.0;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
double
DoFApplication<dim>::residual(FEVector&        dst,
                              const FEVectors& src,
                              bool             apply_boundaries)
{
  // --- assertions ---
  AssertThrow(this->quadrature_residual_cell.size() != 0, ExcNotInitialized());
  AssertThrow(this->quadrature_residual_bdry.size() != 0, ExcNotInitialized());

  // --- global residual ---
  dst.reinit(this->block_info.global);

  // --- local residual ---
  FEVector local_residual(this->block_info.local);

  // --- types of FEVALUES objects we will use further ---
  // --- note: we only need types ---
  FEValues<dim>*     fe_values      = 0;
  FEFaceValues<dim>* fe_face_values = 0;

  // --- info structures construction ---
  CellInfo cell_info(src,
                     this->block_info);
  FaceInfo bdry_info(src,
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


#ifdef OPENFCST_WITH_PETSC
  // -- PETSc parallel global residual --

  std::vector<unsigned int> idx_list;

  PETScWrappers::MPI::Vector DST;
  const types::global_dof_index n_local_dofs
  = DoFTools::count_dofs_with_subdomain_association (*this->dof,this_mpi_process);


  DST.reinit (mpi_communicator, this->dof->n_dofs(), n_local_dofs);


  for( ; cell != endc; ++cell)
  {

      if(cell->subdomain_id() == this->this_mpi_process){
          local_residual = 0;
          local_residual.reinit(this->block_info.local);

          cell_info.reinit(cell);

          cell_residual(local_residual, cell_info);

          if( this->boundary_fluxes )
          {
              for(unsigned int no_face = 0; no_face < GeometryInfo<dim>::faces_per_cell; ++no_face)
              {

                  typename DoFHandler<dim>::face_iterator face = cell->face(no_face);
                  if( face->at_boundary() )
                  {
                      bdry_info.reinit(cell,
                              face,
                              no_face);

                      bdry_residual(local_residual,
                              bdry_info);
                  }
              }
          }

          hanging_node_constraints.distribute_local_to_global(local_residual,cell_info.indices, DST);

          for (auto i: cell_info.indices)
              idx_list.push_back(i);
      }

  }
  DST.compress(VectorOperation::add);

  if( apply_boundaries == true )
      residual_constraints(DST, idx_list);

  DST.compress(VectorOperation::insert);

  //Copy to linear dealii vector
  dst = DST;

#else
  for( ; cell != endc; ++cell)
  {
      local_residual = 0;
      local_residual.reinit(this->block_info.local);

      cell_info.reinit(cell);

      cell_residual(local_residual,
              cell_info);

      if( this->boundary_fluxes )
      {
          for(unsigned int no_face = 0; no_face < GeometryInfo<dim>::faces_per_cell; ++no_face)
          {

              typename DoFHandler<dim>::face_iterator face = cell->face(no_face);
              if( face->at_boundary() )
              {
                  bdry_info.reinit(cell,
                          face,
                          no_face);

                  bdry_residual(local_residual,
                          bdry_info);
              }
          }
      }

      for(unsigned int i = 0; i < this->element->dofs_per_cell; ++i)
          dst(cell_info.indices[i]) += local_residual(i);
  }

  if( apply_boundaries == true )
      residual_constraints(dst);

#endif

  return dst.l2_norm();

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
void
DoFApplication<dim>::residual_constraints(FEVector&) const
{}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef OPENFCST_WITH_PETSC
template <int dim>
void
DoFApplication<dim>::residual_constraints(PETScWrappers::MPI::Vector&, const std::vector<unsigned int>&) const
{}
#endif


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
void
DoFApplication<dim>::cell_residual(
    FEVector& vectors,
    const CellInfo&)
{
    FcstUtilities::log << "Void function " << __FUNCTION__ << " called" << std::endl;
    for (unsigned int i=0;i<vectors.n_blocks();++i)
        FcstUtilities::log << "Vector[" << i << "]\t"
        << vectors.block(i).size() << std::endl;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
void
DoFApplication<dim>::bdry_residual(
    FEVector& vectors,
    const FaceInfo&)
{
    FcstUtilities::log << "Void function " << __FUNCTION__ << " called" << std::endl;
    for (unsigned int i=0;i<vectors.n_blocks();++i)
        FcstUtilities::log << "Vector[" << i << "]\t"
        << vectors.block(i).size() << std::endl;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
void
DoFApplication<dim>::face_residual(
    FEVector& face_vector1,
    FEVector& face_vector2,
    const FaceInfo&,
    const FaceInfo&)
{
    Assert (face_vector1.n_blocks() == face_vector2.n_blocks(),
            ExcDimensionMismatch(face_vector1.n_blocks(), face_vector2.n_blocks()));

    FcstUtilities::log << "Void function " << __FUNCTION__ << " called" << std::endl;
    for (unsigned int i=0;i<face_vector1.n_blocks();++i)
    {
        Assert (face_vector1.block(i).size() == face_vector2.block(i).size(),
                ExcDimensionMismatch(face_vector1.block(i).size(), face_vector2.block(i).size()));
        FcstUtilities::log << "Vector[" << i << "]\t"
        << face_vector1.block(i).size() << std::endl;
    }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
double
DoFApplication<dim>::cell_estimate(const CellInfo&)
{
    return 0.;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
double
DoFApplication<dim>::bdry_estimate(const FaceInfo&)
{
    return 0.;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
double
DoFApplication<dim>::face_estimate(const FaceInfo& info1,
                                   const FaceInfo& info2)
{
    const FEFaceValuesBase<dim>& fe1 = info1.fe(0);
    //  const FEFaceValuesBase<dim>& fe2 = info2.fe(0);
    const std::vector<std::vector<Tensor<1,dim> > >& Du1 = info1.derivatives[0];
    const std::vector<std::vector<Tensor<1,dim> > >& Du2 = info2.derivatives[0];

    double result = 0.;


    double h = info1.face->diameter();
    for (unsigned k=0;k<fe1.n_quadrature_points;++k)
    {
        const double dx = fe1.JxW(k);
        //const Point<dim>& n = fe1.normal_vector(k);
        const dealii::Tensor<1, dim>& n = fe1.normal_vector(k);

        // Loop over all components
        for (unsigned int d=0;d<Du1.size();++d)
        {
            const double diff = (n * Du1[d][k]) - (n * Du2[d][k]);
            result += dx * h * diff * diff;
        }
    }
    return result;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <int dim>
void
DoFApplication<dim>::distribute_face_to_cell_errors()
{
    typename DoFHandler<dim>::active_cell_iterator c;
    for (c = this->dof->begin_active();
         c != this->dof->end();
    ++c)
         {
             const unsigned int cell_index = c->user_index();
             for (unsigned int f = 0; f< GeometryInfo<dim>::faces_per_cell;++f)
             {
                 const unsigned int face_index = c->face(f)->user_index();
                 if (c->at_boundary(f))
                     cell_errors(cell_index) += face_errors(face_index);
                 else
                     cell_errors(cell_index) += 0.5 * face_errors(face_index);
             }
         }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
double
DoFApplication<dim>::global_from_local_errors() const
{
    double sum = 0.;
    for (unsigned int i=0;i<cell_errors.size();++i)
        sum += cell_errors(i);
    return std::sqrt(sum);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<int dim>
void
DoFApplication<dim>::store_triangulation(Triangulation<dim>& new_tr)
{
    new_tr.copy_triangulation(*this->tr);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template class DoFApplication<deal_II_dimension>;