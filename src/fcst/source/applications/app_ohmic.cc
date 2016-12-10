// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: app_diffusion.cc
// - Description: This class describes diffusion in fuel cell cathodes
//                Ficks, one gas
// - Developers: Mayank Sabharwal,University of Alberta
//               Marc Secanell, University of Alberta
//
// ----------------------------------------------------------------------------

#include <applications/app_ohmic.h>

namespace NAME  = FuelCell::Application;

// ---              ---
// ---              ---
// --- AppOhmic ---
// ---              ---
// ---              ---

       //////////////////////////////////////////////////
       //////////////////////////////////////////////////
       // CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
       //////////////////////////////////////////////////
       //////////////////////////////////////////////////

// ---             ---
// --- constructor ---
// ---             ---

template<int dim>
NAME::AppOhmic<dim>::AppOhmic( boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data )
:
FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>(data),
electron_transport_equation(this->system_management, data)
{
  this->repair_diagonal = true;
  FcstUtilities::log <<  "->FuelCell::Application::AppOhmic-" << dim << "D" << std::endl;

}

// ---            ---
// --- destructor ---
// ---            ---

template<int dim>
NAME::AppOhmic<dim>::~AppOhmic()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::AppOhmic<dim>::declare_parameters(ParameterHandler& param)
{

    OptimizationBlockMatrixApplication<dim>::declare_parameters(param);

    // Declare parameters in system management:
    this->system_management.declare_parameters(param);

    // Declare parameters in operating conditions:
    OC.declare_parameters(param);
    
    
   
    // Declare layer class:
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Cathode gas diffusion layer", param);
    
    // Declare equation class:
    electron_transport_equation.declare_parameters(param);
    
    
}

// ---            ---
// --- initialize ---
// ---            ---
template<int dim>
void
NAME::AppOhmic<dim>::initialize(ParameterHandler& param)
{
    OptimizationBlockMatrixApplication<dim>::initialize(param);

    // Initialize parameters in system management:FuelCellShop::Material::PureGas
    this->system_management.initialize(param);
    
    // Initialize parameters in operating conditions:
    OC.initialize(param);
            
    // Initialize layer classes:
    CGDL = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Cathode gas diffusion layer", param);
    std::vector<FuelCellShop::Material::PureGas*> gases;
    gases.push_back(&oxygen);
    gases.push_back(&nitrogen);
    CGDL->set_gases_and_compute(gases, OC.get_pc_atm(), OC.get_T());
    
    // Initialize parameters for physics classes:
    electron_transport_equation.initialize(param);
    
    // --- we make cell couplings for this problem ---
    std::vector<couplings_map> tmp;
    tmp.push_back( electron_transport_equation.get_internal_cell_couplings()    );
    this->system_management.make_cell_couplings(tmp);
    
    // Now, initialize object that are used to setup initial solution and boundary conditions:    
    this->component_materialID_value_maps.push_back( electron_transport_equation.get_component_materialID_value()    );
    OC.adjust_initial_solution(this->component_materialID_value_maps, this->mesh_generator);
    
    this->component_boundaryID_value_maps.push_back( electron_transport_equation.get_component_boundaryID_value() );
    OC.adjust_boundary_conditions(this->component_boundaryID_value_maps, this->mesh_generator);
    
    // --- and then allocate memory for matrices ---
    this->remesh_matrices();
    
    // Output options:
    // - system info:
    //this->system_management.print_system_info();
    
    // - layers info:
    CGDL->print_layer_properties();
        
    // - equations info:
    electron_transport_equation.print_equation_info();
    
}

// ---               ---
// --- init_solution ---
// ---               ---

template<int dim>
void
NAME::AppOhmic<dim>::initialize_solution(FuelCell::ApplicationCore::FEVector& initial_guess,
                                             std::shared_ptr<Function<dim> > initial_function)
{
    DoFApplication<dim>::initialize_solution(initial_guess); 
}

       ///////////////////////////////////
       ///////////////////////////////////
       // LOCAL CG FEM BASED ASSEMBLERS //
       ///////////////////////////////////
       ///////////////////////////////////

// ---             ---
// --- cell_matrix ---
// ---             ---

template<int dim>
void
NAME::AppOhmic<dim>::cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                     const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{

    if(      CGDL->belongs_to_material(cell_info.cell->material_id())   )
    {
        electron_transport_equation.assemble_cell_matrix(cell_matrices, cell_info, CGDL.get());
        

    }
    
    else
    {
      
      const int index = this->system_management.matrix_block_index(electron_transport_equation.get_equation_name(), "electronic_electrical_potential");
      cell_matrices[index].matrix.all_zero();
      
    }

}

// ---               ---
// --- cell_residual ---
// ---               ---

template<int dim>
void
NAME::AppOhmic<dim>::cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_res,
                                       const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
    if(      CGDL->belongs_to_material(cell_info.cell->material_id())   )
    {
        cell_res = 0;
    
    }
    else
    {
        cell_res = 0;
        
    }
}


       /////////////////////
       /////////////////////
       // OTHER FUNCTIONS //
       /////////////////////
       /////////////////////

// ---              ---
// --- dirichlet_bc ---
// ---              ---

template<int dim>
void
NAME::AppOhmic<dim>::dirichlet_bc(std::map<unsigned int, double>& boundary_values) const
{
  FuelCell::InitialAndBoundaryData::make_constant_DirichletBC_values( boundary_values,
                                                                      *this->mapping,
                                                                      *this->dof,
                                                                      this->system_management,
                                                                      this->component_boundaryID_value_maps );
} 

///////////////////
/////////////////////
// POST-PROCESSING //
/////////////////////
/////////////////////

// ---                ---
// --- bdry_responses ---
// ---                ---
template<int dim>
void
NAME::AppOhmic<dim>::bdry_responses(std::vector<double>&                                                     dst,
                                      const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                      const FuelCell::ApplicationCore::FEVector& src)
{
    int user_input_bdry = 2;
    
    // -- Find out what material is the cell made of, i.e. MEA layer)
    const unsigned int bdry_id = bdry_info.dof_face->boundary_indicator();
    
    unsigned solution_index = bdry_info.global_data->find_vector("Solution");
    Tensor<2,dim> sigmaSeff_cell;
    
    FuelCellShop::Equation::VariableInfo xi;
    
    if ( this->system_management.solution_in_userlist("electronic_electrical_potential") )
    {
        xi.solution_index = this->system_management.solution_name_to_index("electronic_electrical_potential"); 
        xi.fetype_index = this->system_management.block_info->base_element[xi.solution_index];
        xi.indices_exist = true;
    }
     
    CGDL->effective_electron_conductivity(sigmaSeff_cell);
    
    int n_q_points_bdry = (bdry_info.fe(xi.fetype_index)).n_quadrature_points;
    
    //-------- Looping over Quadrature points ----------------------------
    std::vector<double> JxW_bdry;
    JxW_bdry.resize(n_q_points_bdry);
    
    if (bdry_id==this->user_input_bdry[0])
    {
        for (unsigned int q = 0; q < n_q_points_bdry; ++q)
        {
            JxW_bdry[q] = (bdry_info.fe(xi.fetype_index)).JxW(q);
            
            dst[0]+= -sigmaSeff_cell*bdry_info.gradients[solution_index][xi.solution_index][q] * bdry_info.fe(xi.fetype_index).normal_vector(q)* JxW_bdry[q];
        }   
    }
}

// ---          ---
// --- evaluate ---
// ---          ---
// This function 
template <int dim>
double
NAME::AppOhmic<dim>::evaluate (const FuelCell::ApplicationCore::FEVectors& src)
{
    return 0.0;
}

// ---          ---
// --- data_out ---
// ---          ---

template<int dim>
void
NAME::AppOhmic<dim>::data_out(const std::string& filename,
                                  const FuelCell::ApplicationCore::FEVectors& src)
{
  //////////////
  // SOLUTION //
  //////////////

  // --- Find solution ---
  FuelCell::ApplicationCore::FEVector solution = src.vector( src.find_vector("Solution") );

  // --- Assign solution names ---
  std::vector<std::string> solution_names;

  solution_names.push_back("electronic_electrical_potential");

  // --- Assign solution interpretations ---
  this->solution_interpretations.clear();
  this->solution_interpretations.resize(this->element->n_blocks(),
                                        DataComponentInterpretation::component_is_scalar);

  ///////////////////////////////////
  // Do further POST-PROCESSING    //
  ///////////////////////////////////

  // --- Create vector of PostProcessing objects ---
  std::vector< DataPostprocessor<dim>* > PostProcessing;
  
  // --- output ---
  DoFApplication<dim>::data_out( filename,
                                 solution,
                                 solution_names);
}

/////////////////////////////
// EXPLICIT INSTANTIATIONS //
/////////////////////////////
template class NAME::AppOhmic<deal_II_dimension>;
