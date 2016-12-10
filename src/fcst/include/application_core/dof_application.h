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
// - Class: dof_application.h
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
//
// ----------------------------------------------------------------------------

#ifndef _FUEL_CELL_APPLICATION_CORE_DOF_APPLICATION_H_
#define _FUEL_CELL_APPLICATION_CORE_DOF_APPLICATION_H_
//-- dealII
#include <deal.II/base/data_out_base.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/distributed/solution_transfer.h>

//-- OpenFCST
#include <grid/geometry.h>
#include <grid/geometries.h>
#include <application_core/application_wrapper.h>
#include <application_core/mesh_loop_info_objects.h>
#include <utils/fcst_utilities.h>
#include <application_core/initial_and_boundary_data.h>

//--C++ Standard Libraries
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <typeinfo>

using namespace dealii;
using namespace FuelCell::ApplicationCore;

namespace dealii
{
    template<int, int> class Triangulation;
    template<int, int> class Mapping;
    template<int, int> class FiniteElement;
    template<int, int> class DoFHandler;
}

namespace FuelCell
{
    namespace ApplicationCore
    {

        /**
         * Base class for all linear applications, i.e., all applications requiring Triangulation and
         * DoFHandler classes.
         *
         * The mesh as well as the dof handler may be created by this class,
         * which is the default, or they may be provided by another object, in
         * which case they must be specified in the constructor.
         *
         * Note that in this class the received or created dof handler is associated
         * to the finite element given by the argument "Element" on the parameter
         * file. Therefore, this class in not responsible to generate the system of equations to
         * be solved, only to initialize the dof handler.
         *
         * Element can either be a single element or a FESystem. In the latter case, the
         * nomenclature used in the paramter file is:
         *
         * set Element = FESystem[element1_type(element1_degree)^number_of_elements1-...-elementN_type(elementN_degree)^number_of_elementsN]
         * Example: set Element = FESystem[FE_DGQ(0)-FE_Q(1)^2]
         *
         * TODO: Improve Parallelization of FEVectors. Currently all FEVectors are stored in each MPI process thereby making the
         * calculations expensive. All FEVectors should be PETScWrappers::BlockVector so that only a part of the solution is stored. Once
         * we do this data_out and transfer_solution_to_coarse_mesh will have to change as solution will not store everything. (Priority: low)
         *
         * @author Guido Kanschat
         * @author Valentin N. Zingan
         * @author Marc Secanell
         * @author Peter Dobson
         */
        template<int dim>
        class DoFApplication : public ApplicationBase
        {
        public:

            /**
             * Shortcut.
             */
            typedef typename FuelCell::ApplicationCore::IntegrationInfo< dim, FEValuesBase<dim> >     CellInfo;

            /**
             * Shortcut.
             */
            typedef typename FuelCell::ApplicationCore::IntegrationInfo< dim, FEFaceValuesBase<dim> > FaceInfo;

            ///@name Constructor, destructor and initialization:
            ///@{
            /**
             * Constructor for an object
             * owning its own mesh and dof handler.
             *
             * If <tt>data</tt> is null, a new dof handler for
             * application data is generated.
             */
            DoFApplication(boost::shared_ptr<ApplicationData> data      = boost::shared_ptr<ApplicationData>());

            /**
             * Constructor for an object owning its own mesh and dof
             * handler and creating new ApplicationData.
             */
            DoFApplication();

            /**
             * Constructor for an object, borrowing mesh and dof
             * handler from another object. If <tt>triangulation_only</tt>
             * is true, only the triangulation is borrowed.
             *
             * @note A pointer to <tt>DH2</tt> must be
             * convertible to a pointer to <tt>DH</tt>.
             */
            DoFApplication(DoFApplication<dim>& dof_app,
                           bool                 triangulation_only);

            /**
             * Destructor which deletes
             * owned objects.
             */
            ~DoFApplication();

            /**
             * Declare parameters related to mesh generation and
             * finite elements. Most importantly, here we declare "Element".
             * "Element" is declared under the subsection "Discretization".
             * "Element" is the element that will be associated
             * to the dof handler (by default FE_Q(1)). This can either be a single element of a FESystem.
             * In the latter case, the
             * nomenclature used in the paramter file is:
             * set Element = FESystem[element1_type(element1_degree)^number_of_elements1-
             * ...-elementN_type(elementN_degree)^number_of_elementsN]
             * Example: set Element = FESystem[FE_DGQ(0)-FE_Q(1)^2]
             *
             * @note The function
             * declare_parameters of
             * derived applications should
             * always call the one of its
             * base class.
             */
            virtual void declare_parameters(ParameterHandler& param);

            /**
             * Initialize application data based on #param object. Once all paramters have been
             * initialized, call #initialize_grid in order to generate the mesh for your domain.
             */
            virtual void initialize(ParameterHandler& param);

            /**
             * Initialize dof handler, count the dofs in each block
             * and renumber the dofs. Different numbering schemes can be easily
             * implemented by reimplementing this function
             * in a derived class and do the sorting after having
             * called this function. Remember though to
             * always sort by component in the end.
             */
            virtual void remesh_dofs();

            /**
             * Refine grid accordingly to the Refinement options under "Grid Generation"
             * in the parameter file.
             */
            virtual void remesh();

            /**
             * Reinitialize the BlockVector dst such that it contains
             * block_size.size() blocks. Each block is
             * reinitialized to dimensions block_sizes[i]. Note that
             * block_sizes contains the number of degrees of freedom
             * per block. So, init_vector could be used to
             * reinitialize the solution vector.
             *
             * @note: init_vector can only be called after
             * remesh_dofs has been called otherwise block_sizes will
             * not have been initialized.
             */
            virtual void init_vector(FEVector& dst) const;

            /**
             * For nonlinear applications, an intial solution is needed in order to
             * assemble the residual and the global matrix. This member function is used to either
             * generate a solution based on initial data or to modify an existing solution in order
             * to meet the boundary conditions imposed in the initial data file.
             *
             * There are three possibilities,
             * - The solution is read from a file and boundary conditions are then applied to it.
             * - If solution is not read from a file and  \p initial_function is specified,
             * the function that is passed to initialize_solution  is used to generate the initial solution
             * and then the boundaries are applied.
             * - If \p initial_function is not specified, then the information in the input file
             * is used to generate a piecewise initial solution, with a constant value per material_ID.
             *
             * In order to correct for boundary conditions and compute the initial solution,
             * DoFApplication will need component_boundaryID_value_maps and
             * component_boundaryID_value_maps respectively
             * to be initialized in the initialization function of the application
             * using the equation classes used in the application, for example
             * @code
             * // Now, initialize object that are used to setup initial solution and boundary conditions:
             * component_materialID_value_maps.push_back( ficks_transport_equation.get_component_materialID_value()    );
             * component_materialID_value_maps.push_back( electron_transport_equation.get_component_materialID_value() );
             * component_materialID_value_maps.push_back( proton_transport_equation.get_component_materialID_value()   );
             *
             * component_boundaryID_value_maps.push_back( ficks_transport_equation.get_component_boundaryID_value() );
             * component_boundaryID_value_maps.push_back( electron_transport_equation.get_component_boundaryID_value() );
             * component_boundaryID_value_maps.push_back( proton_transport_equation.get_component_boundaryID_value()   );
             * @endcode
             *
             * The value of the initial solution in each material id and each boundary id is given in the parameter file in
             * sections like the following for, in this case, equation #FicksTransportEquation:
             * @code
             * subsection Equations
             *   subsection Ficks Transport Equation - oxygen
             *     subsection Initial data
             *       set oxygen_molar_fraction = 2: 1.0, 4: 1.0
             *     end
             *     subsection Boundary data
             *       set oxygen_molar_fraction = 3: 1.0
             *     end
             *   end
             * end
             *
             * @endcode
             *
             * The initial solution could also be read from a previous simulation. In order to do so, the
             * previous simulation would have been called with the following flag in the input file:             *
             * @code
             * subsection Adaptive refinement
             *   set Output solution for transfer = true
             * end
             * @endcode
             * Then, a hidden file .transfer_solution.FEVector is created and loaded afterward by setting
             * the flag
             * @code
             * subsection Adaptive refinement
             *   set Read in initial solution from file = true
             * end
             * @endcode
             *
             * @note: initialize_solution can only be called after
             * remesh_dofs has been called otherwise block_sizes will
             * not have been initialized.
             *
             * <h3>Usage </h3>
             * This class can be used directly or, if you would like to specify an initial solution function,
             * you would need to reimplement the function in the child class and call it as follows:
             * @code
             * template <int dim>
             * void NAME::AppPemfc<dim>::init_solution(FuelCell::ApplicationCore::FEVector& initial_guess) {
             *    std::shared_ptr< Function<dim> > initial_solution (new FuelCell::InitialSolution::AppPemfcIC<dim> (&OC, grid));
             *    DoFApplication<dim>::initialize_solution(initial_guess, initial_solution);
             * }
             * @endcode
             * For an example see @AppPemfc.
             */
            virtual void initialize_solution (FEVector& initial_guess,
                                              std::shared_ptr<Function<dim> > initial_function = std::shared_ptr<Function<dim> >());

            ///@}

            ///@name Other
            ///@{
            /**
             * Estimate the error. By default, the KellyErrorEstimator is used in order to estimate the error
             * in every cell. The error for all components of the solution is added.
             *
             * In general, here a loop over all cells and faces, using the virtual local functions
             * - cell_residual(),
             * - face_residual() and
             * - bdry_residual().
             * could be implemented.
             */
            virtual double estimate(const FEVectors& src);

            /**
             * Evaluate a functional during postprocessing such as drag or lift on an aerodynamics problem
             */
            virtual double evaluate(const FEVectors& src);

            /**
             * Loop over all cells to compute the residual. Uses the virtual local functions
             * - cell_residual(),
             * - face_residual() and
             * - bdry_residual()
             * to be implemented by derived classes.
             */
            virtual double residual(FEVector&        dst,
                                    const FEVectors& src,
                                    bool             apply_boundaries = true);

            /**
             * Apply boundary conditions and hanging node constraints
             * to a residual vector after it has been computed by
             * residual(). This function is used in residual() and its
             * default implementation is void.
             *
             * @note Reimplement this function if you want to
             * constrain the residual. A typical implementation would be
             * <pre>
             * hanging_node_constraints.condense(dst);
             * constrain_boundary(v, true);
             * </pre>
             */
            virtual void residual_constraints(FEVector& dst) const;

            #ifdef OPENFCST_WITH_PETSC
            virtual void residual_constraints(PETScWrappers::MPI::Vector& dst, const std::vector<unsigned int>&) const;
            #endif

            /**
             * Function to copy a triangulation object for use after refinement.
             */
            void store_triangulation(Triangulation<dim>& new_tr);

            /**
             * Add the vector to be transfered from one mesh to the next.
             */
            void add_vector_for_transfer(FEVector* src);

            /**
             * Delete the vector to be transfered from one mesh to the next.
             */
            void delete_vector_for_transfer();

            /**
             * Function to perform the transfer of a solution on a refined grid to the initial coarse grid.
             */
            void transfer_solution_to_coarse_mesh(Triangulation<dim>& tr_coarse,
                                                  FEVector&           coarse_solution,
                                                  FEVector&           refined_solution);

            /**
             * Compute the amount of memory needed by this object.
             */
            unsigned int memory_consumption() const;
            ///@}

            ///@name Initial solution management
            ///@{
            /**
             * Filename where to output the initial grid
             *
             * This value is initialized via the input file using:
             * @code
             * subsection Adaptive refinement
             *   set Initial solution filename = initial_sol             # Set flag to true if you want to output the initial solution to file
             * end
             * @endcode
             */
            std::string filename_initial_sol;

            /** Flag to output the initial solution used to start the solving process. This flag
             * is set to true to make sure the initial solution satisfies initial and boundary conditions.
             * For nonlinear problems, the final solution might also be dependent on the initial guess, so
             * this flag can be used to asses if that is the case.
             *
             * This flag is initialized via the input file using:
             * @code
             * subsection Adaptive refinement
             *   set Output initial solution = false             # Set flag to true if you want to output the initial solution to file
             * end
             * @endcode
             */
            bool output_initial_sol;
            /**
             * Bool flag used to specify if the initial solution to the problem, specially important for non-linear problems,
             * should be read from file. This flag is specified in the input file using:
             * @code
             * subsection Adaptive refinement
             *   set Read in initial solution from file = false  # Check whether a stored solution should be read from file and applied to the grid
             * end
             * @endcode
             */
            bool read_in_initial_solution;

            /**
             * Use user pre-defined initial solution
             */
            bool use_predefined_solution;

            /**
             * Bool flag used to specify if the final solution should be stored in the coarse mesh
             * in order to be used later as an initial solution to solve another problem using the
             * flag #read_in_initial_solution
             *
             * @code
             * subsection Adaptive refinement
             *   set Output solution for transfer = false        # Check whether a solution on a refined grid should be output to file on a coarse mesh
             * end
             * @endcode
             */
            bool output_coarse_solution;

            ///@}

            ///@name Output member functions:
            ///@{
            /**
             * Output the grid used to solve the problem.
             */
            virtual void grid_out(const std::string& basename);

            /**
             * Write data in the format specified by the ParameterHandler.
             */
            virtual void data_out(const std::string& basename,
                                  const FEVectors&   src);

            /**
             * This function prints
             * @p FEVector @p src to a text file @p basename.
             *
             * @note If @p src_indices is empty, the whole
             * @p FEVector @p src will be printed to a text file @p basename.
             */
            void print(const std::string&               basename,
                       const FEVector&                  src,
                       const std::vector<unsigned int>& src_indices = std::vector<unsigned int>()) const;
            ///@}

            ///@name Output data:
            ///@{
            /**
             * \p solution_interpretations identifies whether some \p solution_names are scalars or parts of a vector.
             * Note that if \p solution_interpretations is empty, all \p solution_names are treated as scalars.
             */
            std::vector< DataComponentInterpretation::DataComponentInterpretation > solution_interpretations;

            /**
             * \p postprocessing_interpretations identifies whether some \p postprocessing_names are scalars or parts of a vector.
             * Note that if \p postprocessing_interpretations is empty, all \p postprocessing_names are treated as scalars.
             */
            std::vector< DataComponentInterpretation::DataComponentInterpretation > postprocessing_interpretations;

            /**
             * @deprecated Vector for adjusting output to vector data.
             */
            std::vector<DataComponentInterpretation::DataComponentInterpretation> data_interpretation;

            /**
             * \p output_materials_and_levels if \p true then visualized, otherwise suppressed.
             * Initialized as \p true in the constructors.
             */
            bool output_materials_and_levels;

            /**
             * Vector that will be used by \p data_out to store the
             * material ids.
             */
            Vector<double> output_materials;

            /**
             * Vector that will be used by \p data_out to store the
             * refinement levels at each cell.
             */
            Vector<double> output_levels;

            /**
             * @p true, if you want to output
             * solution and postprocessing data
             * using actual finite element fields Q_n
             * with n >= 1.
             *
             * @p false, if only Q_1 is used
             * even for higher degree data outputs.
             */
            bool output_actual_degree;

            /**
             * @p true, if you want to print
             * @p FEVector @p solution to a text file.
             *
             * @p false, otherwise.
             */
            bool print_solution;

            /**
             * @p true, if you want to print
             * @p FEVector @p postprocessing to a text file.
             *
             * @p false, otherwise.
             */
            bool print_postprocessing;

            /**
             * The indices of the @p FEVector @p solution
             * to be printed to a text file.
             *
             * @note If @p print_solution=true and @p solution_printing_indices is not
             * specified in the parameter file, the whole @p FEVector @p solution will be printed to a text file.
             */
            std::vector<unsigned int> solution_printing_indices;

            /**
             * The indices of the @p FEVector @p postprocessing
             * to be printed to a text file.
             *
             * @note If @p print_postprocessing=true and @p postprocessing_printing_indices is not
             * specified in the parameter file, the whole @p FEVector @p postprocessing will be printed to a text file.
             */
            std::vector<unsigned int> postprocessing_printing_indices;

            /**
             * @p true, if the whole blocks of @p FEVector @p solution or @p FEVector @p postprocessing
             * to be printed to a text file instead of separate indices.
             *
             * @p false, otherwise.
             */
            bool print_blocks_instead_of_indices;
            
            /**
             * If true, all cell matrices and right hand sides will be output.
             * 
             * @warning
             * Might lead to huge output depending on number of degrees of freedom.
             * @endwarning
             */
            bool output_matrices_and_rhs;
            ///@}

            /**
             * Curved boundary.
             */
            boost::shared_ptr< Boundary<dim> > curved_boundary;

            /**
             * Curved boundary ID.
             */
            types::boundary_id curved_bdry_id;

        protected:
            ///@name Initialization of application and residual
            ///@{
            /**
             * Initialize from parameter values. Called by initialize.
             *
             * @deprecated
             */
            void _initialize(ParameterHandler& param);

            /**
             * Function used to read in a mesh and hand it over to the boost::shared_ptr<Triangulation<dim> > tr
             * object. This object stores the mesh in the application and it is used to initialize
             * the residual and global matrix sizes.
             *
             * @note This routine is usually called by initialize or _initialize
             *
             * @note under development
             */
            virtual void initialize_triangulation(ParameterHandler& param);

            /**
             * Create a mesh and assign it to object #tr. This member function is usually called by #initialize
             */
            //virtual void initialize_grid(ParameterHandler& param){};

            /**
             * Function to transfer a solution on a refined grid to the initial coarse grid
             *
             * Set \p read to true if you would like the application to read the initial solution from
             * a previous simulation. The initial solution from a previous iteration is stored in the hiddern
             * file .transfer_solution.FEVector. In order to generate this file, you need to setup the
             * following parameter in the input file to true:
             * @code
             * subsection Adaptive refinement
             *   set Output solution for transfer = true
             * end
             * @endcode
             *
             * @note This function is called by initialize_solution
             * @note This function should be private
             */
            void read_init_solution(FEVector& dst,
                                    bool&     good_solution) const;
            ///@}

            ///@name Local assembly routines
            ///@{
            /**
             * Local integration.
             */
            virtual void cell_residual(FEVector&       cell_vector,
                                       const CellInfo& cell);

            /**
             * Local integration.
             */
            virtual void bdry_residual(FEVector&       face_vector,
                                       const FaceInfo& face);

            /**
             * Local integration.
             */
            virtual void face_residual(FEVector&       face_vector1,
                                       FEVector&       face_vector2,
                                       const FaceInfo& face1,
                                       const FaceInfo& face2);

            /**
             * Local estimation.
             */
            virtual double cell_estimate(const CellInfo& src);

            /**
             * Local estimation.
             */
            virtual double bdry_estimate(const FaceInfo& src);

            /**
             * Local estimation.
             */
            virtual double face_estimate(const FaceInfo& src1,
                                         const FaceInfo& src2);
            ///@}

            ///@name Initial and boundary data information:
            ///@{
            /**
             * Each entry of this std::vector reflects the following
             * structure (see FuelCell::InitialAndBoundaryData namespace docs):
             *
             * - \p first  \p argument : name of the solution component,
             * - \p second \p argument : material id,
             * - \p third  \p argument : value of the solution component.
             */
            std::vector< component_materialID_value_map > component_materialID_value_maps;

            /**
             * Each entry of this std::vector reflects the following
             * structure (see FuelCell::InitialAndBoundaryData namespace docs):
             *
             * - \p first  \p argument : name of the solution component,
             * - \p second \p argument : boundary id,
             * - \p third  \p argument : value of the solution component.
             */
            std::vector< component_boundaryID_value_map > component_boundaryID_value_maps;

            ///@}
            ///@name Pre-processor objects
            ///@{
            /**
             * Grid.
             */
            boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > mesh_generator;

            /**
             * Pointer to the Triangulation
             * object.
             */
            boost::shared_ptr<Triangulation<dim> > tr;

            ///@}
            /**
             * The object for writing grids.
             */
            GridOut g_out;

            /**
             * The object for writing data.
             */
            DataOut<dim, DoFHandler<dim> > d_out;

            /**
             * This routine is used to write data in the format specified by the ParameterHandler.
             * This routine is usually called inside child applications wheh data_out is
             * specified.
             *
             * Usually, the solution if first processed and this routine is called last (see Usage below).
             *
             * Parameters:
             *
             * - \p basename is a base of an output file name,
             * - \p solution       is the global solution computed in DoFs,
             * - \p solution_names is the component names of the \p solution,
             * - \p PostProcessing is a vector of DataPostprocessor objects used to evaluate other
             * values using the solution vector. Many DataPostprocessor objects have already been developed.
             *
             * <h3>Usage </h3>
             *
             * The code below, from AppCathode, highlights how this routine is usually used inside an application.
             * First, the filename, solution vector and solution_name vectors are obtained. Then,
             * a vector of DataPostprocessor objects is used to compute any necessary additional
             * post-processing results such as current density, relative humidity or others.
             * Finally, this data_out routine is called:
             *
             * @code
             * template<int dim>
             * void NAME::AppCathode<dim>::data_out(const std::string& filename, const FuelCell::ApplicationCore::FEVectors& src){
             *  // --- Find solution ---
             *  FuelCell::ApplicationCore::FEVector solution = src.vector( src.find_vector("Solution") );
             *  // --- Assign solution names ---
             *  std::vector<std::string> solution_names;
             *  solution_names.push_back("oxygen_molar_fraction");
             *  solution_names.push_back("protonic_electrical_potential");
             *  solution_names.push_back("electronic_electrical_potential");
             *  // --- Assign solution interpretations ---
             *  this->solution_interpretations.clear();
             *  this->solution_interpretations.resize(this->element->n_blocks(), DataComponentInterpretation::component_is_scalar);
             *  // --- Create vector of PostProcessing objects ---
             *  std::vector< DataPostprocessor<dim>* > PostProcessing;
             *  // --- current ---
             *  FuelCellShop::PostProcessing::ORRCurrentDensityDataOut<dim> current(&this->system_management, CCL, &OC);
             *  PostProcessing.push_back(&current);
             *  // --- output ---
             *  DoFApplication<dim>::data_out( filename,
             *                                 solution,
             *                                 solution_names,
             *                                 PostProcessing);
             * }
             * @endcode
             */
            virtual void data_out(const std::string&                            basename,
                                  const FEVector&                               solution,
                                  const std::vector<std::string>&               solution_names,
                                  const std::vector< DataPostprocessor<dim>* >& PostProcessing);

            /**
             * This function outputs
             * the results of a computation.
             *
             * Parameters:
             *
             * - \p basename is a base of an output file name,
             *
             * - \p solution       is the global solution computed in DoFs,
             * - \p solution_names is the component names of the \p solution,
             *
             * - \p postprocessing       is what you compute in DoFs (or in cells in the case of piece-wise constant functionals) after the \p solution has been computed,
             * - \p postprocessing_names is the component names of the \p postprocessing.
             */
            virtual void data_out(const std::string&              basename,

                                  const FEVector&                 solution,
                                  const std::vector<std::string>& solution_names,

                                  const FEVector&                 postprocessing       = FEVector(),
                                  const std::vector<std::string>& postprocessing_names = std::vector<std::string>());

            /**
             * Apply either homogeneous or inhomogeneous
             * #boundary_constraints.
             */
            void constrain_boundary(FEVector& v,
                                    bool      homogeneous) const;

            /**
             * This object knows everything about FCST equations, variables,
             * couplings, etc.
             */
            FuelCell::SystemManagement system_management;

            /**
             * List of all nodes constrained by a strong
             * boundary condition, together with a value to be assigned.
             *
             * This map must be filled by the application and is used
             * in constrain_boundary().
             */
            std::map<unsigned int, double> boundary_constraints;

            /**
             * Constraint Matrix object. This object contains a list of the
             * constraints from the hanging nodes.
             * It needs to be used to remove hanging nodes
             * from the system matrix and rhs before the
             * system of equations is solved using
             * @verbatim
             * hanging_node_constraints.condense(system_matrix);
             * hanging_node_constraints.condense(system_rhs);
             * @endverbatim
             * and it needs to be used to add the hanging nodes to the solution
             * once the system is solved using
             * @verbatim
             * hanging_node_constraints.distribute(solution);
             * @endverbatim
             */
            ConstraintMatrix hanging_node_constraints;

            /**
             * The mapping used for the
             * transformation of mesh
             * cells.
             *
             * @note Should be set by
             * initialize() in a derived
             * class before _initialize()
             * is called.
             */
            boost::shared_ptr<Mapping<dim> > mapping;

            /**
             * Degree used for polynomial mapping of arbitrary order.
             *
             * By mapping we mean the transformation between the unit cell (i.e. the unit line, square, or cube) to the cells in real space.
             * Before we have implicitly used linear or d-linear mappings and have not noticed this at all, since this is what happens if we do not do anything special.
             * However, if the domain has curved boundaries, there are cases where the piecewise linear approximation of the boundary (i.e. by straight line segments) is not sufficient,
             * and we may want that the computational domain is an approximation to the real domain using curved boundaries as well.
             *
             * If the boundary approximation uses piecewise quadratic parabolas to approximate the true boundary,
             * then we say that this is a quadratic or \f$ Q_2 \f$ approximation.
             * If we use piecewise graphs of cubic polynomials, then this is a \f$ Q_3 \f$ approximation, and so on.
             *
             * For some differential equations, it is known that piecewise linear approximations of the boundary, i.e. \f$ Q_1 \f$ mappings, are not sufficient
             * if the boundary of the domain is curved.
             * Examples are the biharmonic equation using \f$ C_1 \f$ elements,
             * or the Euler equation on domains with curved reflective boundaries.
             * In these cases, it is necessary to compute the integrals using a higher order mapping.
             * The reason, of course, is that if we do not use a higher order mapping,
             * the order of approximation of the boundary dominates the order of convergence of the entire numerical scheme,
             * irrespective of the order of convergence of the discretization in the interior of the domain.
             */
            unsigned int mapping_degree;

            /**
             * The finite element used in #dof. This can be a single element or a FESystem
             */
            boost::shared_ptr<FiniteElement<dim> > element;

            /**
             * Pointer to the DoFHandler object.
             */
            boost::shared_ptr<DoFHandler<dim> > dof;

            BlockInfo block_info;

            /**
             * The result of error estimation by cell. This variable will be used for adaptive mesh refinement.
             */
            Vector<float> cell_errors;

            /**
             * The result of error estimation by face.
             */
            Vector<float> face_errors;

            /**
             * Refinement parameter from parameter file.
             */
            std::string refinement;

            /**
             * Initial refinement from parameter file
             */
            unsigned int initial_refinement;

            /**
             * Refinement threshold for adaptive method from parameter file.
             */
            double refinement_threshold;

            /**
             * Coarsening threshold for adaptive method from parameter file.
             */
            double coarsening_threshold;

            /**
             * Flag for sorting with Cuthill McKee algorithm.
             *
             * @note Read from parameter file
             */
            bool sort_cuthill;

            /**
             * Direction for downstream sorting. No downstream
             * sorting if this vector is zero.
             *
             * @note Read from parameter file
             */
            Point<dim> sort_direction;

            /**
             * List of vector names to be transfered from one grid to the next.
             */
            std::vector<FEVector*> transfer_vectors;

            /**
             * Quadrature rule for residual computation on cells.
             */
            Quadrature<dim> quadrature_residual_cell;

            /**
             * Quadrature rule for residual computation on boundary faces.
             */
            Quadrature<dim-1> quadrature_residual_bdry;

            /**
             * Quadrature rule for residual
             * computation on faces.
             */
            Quadrature<dim-1> quadrature_residual_face;

            /**
             * Couplings through cell bilinear forms.
             * cell_coupling allows to specify which variables
             * couple with which equations. This is used
             * by DoFTools to generate a sparsity pattern
             * that does not contain elements unless specified in the
             * cell_coupings, therefore, reducing the amount of memory
             * needed.
             *
             * @note cell_coupings is a two dimensional table (i.e. a matrix) of
             * DoFTools::Coupling objects. DoFTools::Coulings takes the values:
             * always (components couple), zero (no coupling) and nonzero
             *
             * Example: For the 2D Stokes equations
             * \f{eqnarray*}
             * \frac{d^2u}{dx^2} +\frac{d^2u}{dy^2} + \frac{dp}{dx} = 0 \\
             * \frac{d^2v}{dx^2} +\frac{d^2v}{dy^2} + \frac{dp}{dy} = 0 \\
             * \frac{du}{dx} +\frac{dv}{dy} = 0
             *     \f}
             * we need to initialize cell_coupling to:
             * @verbatim
             * cell_coupling.reinit(3,3);
             * cell_coupling(0,0) = DoFTools::always;
             * cell_coupling(0,1) = DoFTools::zero;
             * cell_couplings(0,2) = DoFTools::always;
             * cell_couplings(1,0) = DoFTools::zero;
             * ...
             * @endverbatim
             */
            Table<2, DoFTools::Coupling> cell_couplings;

            /**
             * Couplings through flux
             * bilinear forms
             */
            Table<2, DoFTools::Coupling> flux_couplings;

            /**
             * Extend the integration loops in assemble() and residual()
             * also to boundary faces.
             */
            bool boundary_fluxes;

            /**
             * Extend the integration loops in assemble() and residual()
             * also to interior faces.
             */
            bool interior_fluxes;

            /**
             * Controls verbosity of certain functions. Zero is
             * quiet and default is one.
             */
            unsigned int verbosity;

        private:

            /**
             * Sort the degrees of freedom. In a derived class,
             * sorting order and schemes can be changed by
             * overloading this function.
             */
            virtual void sort_dofs(DoFHandler<dim>* dof_handler) const;

            /**
             * After computing the error contributions on faces and
             * cells, this function adds half of the contribution
             * computed by face_estimate() and the contribution
             * computed by bdry_estimate() to the adjacent cell. It
             * deletes #face_errors in the end.
             *
             * This function can be overridden by derived
             * classes to handle both estimates separately.
             */
            virtual void distribute_face_to_cell_errors();

            /**
             * Compute a global estimation value from #cell_errors. In
             * order for this to be useful, #face_errors should be
             * transfered to cells first, for instance by
             * distribute_face_to_cell_errors().
             *
             * This function computes the sum of the #cell_errors and
             * takes the square root in the end.  Reimplementation in
             * derived classes allows for different norms.
             *
             * The function is called by estimate() to compute its
             * return value.
             */
            virtual double global_from_local_errors() const;
            
            /**
             * Number of refinements.
             * 
             * @note Normally, <tt>n_ref</tt> is stored in the data object. But if
             * adaptive refinement class is not used, we need to initialize this
             * variable in this class so that the mesh refinement saturation does
             * not happen (see the source file of this class, search for <tt>n_ref</tt>).
             */
            unsigned int n_ref;
        };

    } // ApplicationCore

} // FuelCell

#endif