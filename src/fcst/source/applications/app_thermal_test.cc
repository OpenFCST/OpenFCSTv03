//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: app_thermal_test.cc 21-05-2013
//    - Description: Class designed to test thermal transport equation.
//    - Developers: Madhur Bhaiya
//    - Id: $Id: app_thermal_test.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <applications/app_thermal_test.h>

namespace NAME = FuelCell::Application;

using namespace FuelCell::ApplicationCore;

//---------------------------------------------------------------------------
template <int dim>
NAME::AppThermalTest<dim>::AppThermalTest(boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data)
:
OptimizationBlockMatrixApplication<dim>(data),
thermal_transport(this->system_management, data)
{
    this->repair_diagonal = true;
    FcstUtilities::log << "FuelCell::Application::AppThermal_test-" << dim<<"d"<<std::endl;
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::AppThermalTest<dim>::declare_parameters(ParameterHandler& param)
{
    OptimizationBlockMatrixApplication<dim>::declare_parameters(param);
    
    // Declaring case type parameter entry
    param.enter_subsection("Thermal test");
    {
        param.declare_entry ("Case type",
                             "Case1",
                             Patterns::Selection("Case1|Case2"),
                             "Case type used in test the thermal transport equation class.");
    }
    param.leave_subsection();
    
    // Declare parameters in system management:
    this->system_management.declare_parameters(param);
    
    // Declare parameters for layer classes:
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Test layer", param);
    
    // Declare parameters for physics classes:
    thermal_transport.declare_parameters(param);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppThermalTest<dim>::initialize(ParameterHandler& param)
{
    OptimizationBlockMatrixApplication<dim>::initialize(param);
    
    // Determining the case type to solve:
    param.enter_subsection("Thermal test");
    {
        case_type = param.get("Case type");
    }
    param.leave_subsection();
    
    //Initialize the problem data:
    this->system_management.initialize(param);
    
    // Make sure that the number of input finite elements matches with the problem:
    AssertThrow (this->element->n_blocks() == this->system_management.get_equation_names().size(),
                 ExcDimensionMismatch(this->element->n_blocks(), this->system_management.get_equation_names().size()));
    AssertThrow (this->element->n_blocks() == this->system_management.get_number_of_solution_names(),
                 ExcDimensionMismatch(this->element->n_blocks(), this->system_management.get_number_of_solution_names()));
    
    // Initialize layer classes:    
    test_layer = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Test layer",param);

    // Initializing Equation classes
    thermal_transport.initialize(param);
    
    std::vector<couplings_map> tmp;
    tmp.push_back( thermal_transport.get_internal_cell_couplings() );    
    this->system_management.make_cell_couplings(tmp);

    // Initialize matrices and spartisity pattern for the whole system
    this->remesh_matrices();
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppThermalTest<dim>::initialize_triangulation(ParameterHandler& param)
{
    Point<dim> p1(0.0, 0.0);
    Point<dim> p2(5.0, 1.0);
    
    GridGenerator::hyper_rectangle(*this->tr, p1, p2, true);
    
    typename Triangulation<dim>::active_cell_iterator cell = (*this->tr).begin_active(), endc = (*this->tr).end();
    
    for (; cell!=endc; ++cell)
    {
        cell->set_material_id(1);
    }
    (*this->tr).refine_global(4);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppThermalTest<dim>::initialize_solution(FuelCell::ApplicationCore::FEVector& dst,
                                               std::shared_ptr<Function<dim> > initial_function)
{
    FuelCell::InitialSolution::AppThermalTestIC<dim> initial_solution(case_type);
        
    // resize vector:
    FcstUtilities::log << "Generating default initial solution" << std::endl;
    dst.reinit(this->block_info.global);
    VectorTools::interpolate (*this->dof, initial_solution, dst);
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::AppThermalTest<dim>::cell_matrix(FuelCell::ApplicationCore::MatrixVector& cell_matrices,
                                       const typename DoFApplication<dim>::CellInfo& info)
{
    thermal_transport.assemble_cell_matrix(cell_matrices, info, test_layer.get());
}
     
//---------------------------------------------------------------------------
template <int dim>
void 
NAME::AppThermalTest<dim>::cell_residual(FuelCell::ApplicationCore::FEVector& cell_vector,
                                         const typename DoFApplication<dim>::CellInfo& info)
{ 
    // -- Assertion before starting routine:
    // Make sure vectors are the right size
    Assert (cell_vector.n_blocks() == this->element->n_blocks(),
            ExcDimensionMismatch (cell_vector.n_blocks(), this->element->n_blocks()));
    
    thermal_transport.assemble_cell_residual(cell_vector, info, test_layer.get());
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppThermalTest<dim>::bdry_matrix(FuelCell::ApplicationCore::MatrixVector& bdry_matrices,
                                       const typename DoFApplication<dim>::FaceInfo& bdry_info)
{
    thermal_transport.assemble_bdry_matrix(bdry_matrices, bdry_info, test_layer.get());
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppThermalTest<dim>::bdry_residual(FuelCell::ApplicationCore::FEVector& bdry_vector,
                                         const typename DoFApplication<dim>::FaceInfo& bdry_info)
{
    thermal_transport.assemble_bdry_residual(bdry_vector, bdry_info, test_layer.get());
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::AppThermalTest<dim>::dirichlet_bc(std::map<unsigned int, double>& boundary_values) const
{
    std::vector<bool> comp_mask(this->element->n_blocks(), true);
    
    if (case_type == "Case1")
    {
        for (unsigned int i=0 ; i<4; ++i)
        {
            VectorTools::interpolate_boundary_values(*this->mapping,
                                                     *this->dof,
                                                     i,
                                                     ZeroFunction<dim>(this->element->n_blocks()),
                                                     boundary_values,
                                                     comp_mask);
        }
    }
    
    else if (case_type == "Case2")
    {
        VectorTools::interpolate_boundary_values(*this->mapping,
                                                 *this->dof,
                                                 1,
                                                 ZeroFunction<dim>(this->element->n_blocks()),
                                                 boundary_values,
                                                 comp_mask);
    }
    
    else
        AssertThrow( false, ExcNotImplemented() );
}

//---------------------------------------------------------------------------
template <int dim>
double 
NAME::AppThermalTest<dim>::estimate(const FuelCell::ApplicationCore::FEVectors& vectors)
{
    // --- solution ---
    FuelCell::ApplicationCore::FEVector solution;
    solution.reinit( this->block_info.global );
    unsigned int index = vectors.find_vector("Solution");
    solution = vectors.vector(index);
    
    // --- components to be used ---
    std::vector<bool> comp_mask (this->element->n_blocks(), true);
    
    // --- cell errors ---
    this->cell_errors.reinit(this->tr->n_active_cells());
    KellyErrorEstimator<dim>::estimate(*this->mapping,
                                       *this->dof,
                                       this->quadrature_residual_face,
                                       typename FunctionMap<dim>::type(),
                                       solution,
                                       this->cell_errors,
                                       comp_mask);
    
    return 0.0; 
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::AppThermalTest<dim>::data_out(const std::string& basename,
                                    const FuelCell::ApplicationCore::FEVectors& vectors)
{
    this->solution_interpretations.clear();
    this->solution_interpretations.resize(this->element->n_blocks(), DataComponentInterpretation::component_is_scalar);

    FcstUtilities::log << "Datafile:" << basename + this->d_out.default_suffix() << std::endl;
    
    // --- data out ---
    DoFApplication<dim>::data_out(basename,
                                  vectors.vector( vectors.find_vector("Solution") ),
                                  this->system_management.get_solution_names());
}
    

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
FuelCell::InitialSolution::AppThermalTestIC<dim>::AppThermalTestIC(const std::string& type)
:
Function<dim> (1),
case_type(type)
{
}

//---------------------------------------------------------------------------
template <int dim>
double
FuelCell::InitialSolution::AppThermalTestIC<dim>::value (const Point<dim> &p,
                                                         const unsigned int component) const
{
    double x = p(0);
    double y = p(1);
    
    if (case_type == "Case1")
    {
        if (y == 1.0)
            return 50.0 + std::sin((Constants::Pi()*x)/5.0)*30.0;
        else
            return 50.0;
    }
    
    else if (case_type == "Case2")
        return 50.0;
    
    else
        AssertThrow( false, ExcNotImplemented() );
}
                                                                                    
                                                                                    
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations. 
template class NAME::AppThermalTest<deal_II_dimension>;
template class FuelCell::InitialSolution::AppThermalTestIC<deal_II_dimension>;
                                                                                        
