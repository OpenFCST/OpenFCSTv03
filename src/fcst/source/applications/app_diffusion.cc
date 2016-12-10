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
//               Aslan Kosakian, University of Alberta
//
// ----------------------------------------------------------------------------

#include <applications/app_diffusion.h>

namespace NAME  = FuelCell::Application;

// ---              ---
// ---              ---
// --- AppDiffusion ---
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
std::string name_section = "oxygen";
template<int dim>
NAME::AppDiffusion<dim>::AppDiffusion( boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data )
:
FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>(data),
ficks_transport_equation(this->system_management,name_section, data),
equation_debug_output(this->system_management, data)
{
  this->repair_diagonal = true;
  FcstUtilities::log <<  "->FuelCell::Application::AppDiffusion-" << dim << "D" << std::endl;

}

// ---            ---
// --- destructor ---
// ---            ---

template<int dim>
NAME::AppDiffusion<dim>::~AppDiffusion()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::AppDiffusion<dim>::declare_parameters(ParameterHandler& param)
{

    OptimizationBlockMatrixApplication<dim>::declare_parameters(param);

    // Declare parameters in system management:
    this->system_management.declare_parameters(param);

    // Declare parameters in operating conditions:
    OC.declare_parameters(param);
    
    // Setup System management parameters that we would like to use if were are using this application:
    param.enter_subsection("System management");
    {
        param.declare_entry("Number of solution variables","1");
        param.enter_subsection("Solution variables");
        {   
            param.declare_entry("Solution variable 1","oxygen_molar_fraction");
        }
        param.leave_subsection();
        param.enter_subsection("Equations");
        {
            param.declare_entry("Equation 1","Ficks Transport Equation - oxygen");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Materials");
        {   param.declare_entry("solute",
                                "oxygen",
                                Patterns::Selection("oxygen | nitrogen | hydrogen | water | air | helium"),
                                "Gas species to be treated as the solute for solving the Fick's Law ");
            param.declare_entry("solvent",
                                "nitrogen",
                                Patterns::Selection("oxygen | nitrogen | hydrogen | water | air | helium"),
                                "Gas species to be treated as the solvent for solving the Fick's Law ");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
   
    // Declare layer class:
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Cathode gas diffusion layer", param);
    
    // Declare equation class:
    ficks_transport_equation.declare_parameters(param);
    
    
}

// ---            ---
// --- initialize ---
// ---            ---
template<int dim>
void
NAME::AppDiffusion<dim>::initialize(ParameterHandler& param)
{
    OptimizationBlockMatrixApplication<dim>::initialize(param);
    
    std::string solutename;
    std::string solventname;
    
    // Initialize parameters in system management:FuelCellShop::Material::PureGas
    this->system_management.initialize(param);
    
    // Initialize parameters in operating conditions:
    OC.initialize(param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Materials");
        {   
            solutename = param.get("solute");
            solventname = param.get("solvent");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    if (solutename=="oxygen")
        solute = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Oxygen);
    else if (solutename=="nitrogen")
        solute = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Nitrogen);
    else if (solutename=="helium")
        solute = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Helium);
    else if (solutename=="hydrogen")
        solute = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Hydrogen);
    else if (solutename=="water")
        solute = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::WaterVapor);
    else
    {
        FcstUtilities::log << "Solute specified does not match any of components defined in the database" << std::endl;
        AssertThrow(false, ExcInternalError());
    }

    if (solventname=="oxygen")
        solvent = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Oxygen);
    else if (solventname=="nitrogen")
        solvent = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Nitrogen);
    else if (solventname=="helium")
        solvent = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Helium);
    else if (solventname=="hydrogen")
        solvent = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Hydrogen);
    else if (solventname=="water")
        solvent = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::WaterVapor);
    else
    {
        FcstUtilities::log << "Solvent specified does not match any of components defined in the database" << std::endl;
        AssertThrow(false, ExcInternalError());
    }
    // Initialize materials and layers:
    solute->declare_parameters(param);
    solvent->declare_parameters(param);
    solute->initialize(param);
    solvent->initialize(param);
    
    // Initialize gases and material classes:  
    std::vector< FuelCellShop::Material::PureGas* > gases;
    gases.push_back(solute.get());
    gases.push_back(solvent.get());
        
    // Initialize layer classes:
    CGDL = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Cathode gas diffusion layer", param);
    CGDL->set_gases_and_compute(gases, OC.get_pc_atm(), OC.get_T());
    
    // Initialize parameters for physics classes:
    ficks_transport_equation.set_solute_and_solvent (solute.get(), solvent.get(), param);
    ficks_transport_equation.initialize(param);
    
    // --- we make cell couplings for this problem ---
    std::vector<couplings_map> tmp;
    tmp.push_back( ficks_transport_equation.get_internal_cell_couplings()    );
    this->system_management.make_cell_couplings(tmp);
    
    // Now, initialize object that are used to setup initial solution and boundary conditions:    
    this->component_materialID_value_maps.push_back( ficks_transport_equation.get_component_materialID_value()    );
    OC.adjust_initial_solution(this->component_materialID_value_maps, this->mesh_generator);
    
    this->component_boundaryID_value_maps.push_back( ficks_transport_equation.get_component_boundaryID_value() );
    OC.adjust_boundary_conditions(this->component_boundaryID_value_maps, this->mesh_generator);
    
    // --- and then allocate memory for matrices ---
    this->remesh_matrices();
    
    // Output options:
    // - system info:
    //this->system_management.print_system_info();
    
    // - layers info:
    //CGDL->print_layer_properties();
        
    // - equations info:
    //ficks_transport_equation.print_equation_info();
}

// ---               ---
// --- init_solution ---
// ---               ---

template<int dim>
void
NAME::AppDiffusion<dim>::initialize_solution(FuelCell::ApplicationCore::FEVector& initial_guess,
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
NAME::AppDiffusion<dim>::cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                     const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{     
    unsigned int index=get_solution_index();
    if(      CGDL->belongs_to_material(cell_info.cell->material_id())   )
    {
        ficks_transport_equation.assemble_cell_matrix(cell_matrices, cell_info, CGDL.get());   
    }
    FuelCellShop::Equation::VariableInfo xi = get_xi();
    if (this->output_matrices_and_rhs)
        equation_debug_output.output_matrix(cell_matrices, cell_info, xi);
}

// ---               ---
// --- cell_residual ---
// ---               ---

template<int dim>
void
NAME::AppDiffusion<dim>::cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_res,
                                       const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
    if ((this->data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NEWTONBASIC) || 
        (this->data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NEWTON3PP) || 
        (this->data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NEWTONLINESEARCH))
    {
        ficks_transport_equation.assemble_cell_residual(cell_res,cell_info,CGDL.get());
    }
    else if ((this->data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::PICARD) || 
        (this->data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NONE))
    {
        if(      CGDL->belongs_to_material(cell_info.cell->material_id())   )
        {       
            cell_res = 0;
            if (this->output_matrices_and_rhs)
                equation_debug_output.output_vector(cell_res);   
        }
    
        else
        {
            cell_res = 0;
            
        }
    }
}

// ---               ---
// --- bdry_matrix   ---
// ---               ---

template<int dim>
void
NAME::AppDiffusion<dim>::bdry_matrix(FuelCell::ApplicationCore::MatrixVector&                               bdry_matrices,
                                       const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info)
{    }


template<int dim>
void
NAME::AppDiffusion<dim>::bdry_residual(FuelCell::ApplicationCore::FEVector&                                     bdry_res,
                                       const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info)
{
    const unsigned int material_id = bdry_info.dof_active_cell->material_id();
    if(CGDL->belongs_to_material(material_id))
        ficks_transport_equation.assemble_bdry_residual(bdry_res,bdry_info,CGDL.get());
    
}

       /////////////////////
       /////////////////////
       // OTHER FUNCTIONS //
       /////////////////////
       /////////////////////
     
// ---                    ---
// --- get_solution_index ---
// ---                    ---     
       
template<int dim>
unsigned int
NAME::AppDiffusion<dim>::get_solution_index()
{
    std::stringstream ss;
    ss <<solute->name_material()<<"_molar_fraction";
    std::string speciesname = ss.str();        
    return this->system_management.matrix_block_index(ficks_transport_equation.get_equation_name(), speciesname);  
}

// ---              ---
// --- dirichlet_bc ---
// ---              ---

template<int dim>
void
NAME::AppDiffusion<dim>::dirichlet_bc(std::map<unsigned int, double>& boundary_values) const
{
  FuelCell::InitialAndBoundaryData::make_constant_DirichletBC_values( boundary_values,
                                                                      *this->mapping,
                                                                      *this->dof,
                                                                      this->system_management,
                                                                      this->component_boundaryID_value_maps );
} 



template<int dim>
FuelCellShop::Equation::VariableInfo
NAME::AppDiffusion<dim>::get_xi()
{
    FuelCellShop::Equation::VariableInfo xi;
    
    if ( this->system_management.solution_in_userlist("oxygen_molar_fraction") )
    {
        xi.solution_index = this->system_management.solution_name_to_index("oxygen_molar_fraction"); 
        xi.fetype_index = this->system_management.block_info->base_element[xi.solution_index];
        xi.indices_exist = true;
    }     
    return xi;  
}

////////////////////
/////////////////////
// POST-PROCESSING //
/////////////////////
/////////////////////

// ---                ---
// --- bdry_responses ---
// ---                ---


template<int dim>
void
NAME::AppDiffusion<dim>::bdry_responses(std::vector<double>&                                                     dst,
                                      const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                      const FuelCell::ApplicationCore::FEVector& src)
{
    const unsigned int material_id = bdry_info.dof_active_cell->material_id();
    if(CGDL->belongs_to_material(material_id))
    {
        
        int index_gas, index_solvent;
        CGDL->get_gas_index(this->solute.get(), index_gas);
        CGDL->get_gas_index(this->solvent.get(), index_solvent);
        // -- Find out what material is the cell made of, i.e. MEA layer)
        const unsigned int bdry_id = bdry_info.dof_face->boundary_indicator();
        
        unsigned solution_index = bdry_info.global_data->find_vector("Solution");
        Table< 2, Tensor<2,dim> > Deff_iso;
        Tensor<2,dim> Deff;
        double T,p; //Temperature[K] and pressure [atm]
        
        FuelCellShop::Equation::VariableInfo xi = get_xi();
        
        CGDL->effective_gas_diffusivity(Deff_iso);
        CGDL->get_T_and_p(T,p);
        double concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T))*Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3);  //Units mol/cm3
        int n_q_points_bdry = (bdry_info.fe(xi.fetype_index)).n_quadrature_points;
        
        if (this->data->flag_exists("knudsen") && this->data->flag("knudsen"))
        {
            double radius,DKnud =0;
            std::string text = boost::lexical_cast<std::string>(bdry_info.cell->id());
            int cell_id =  FcstUtilities::cellId_to_index(text);
            DKnud = ficks_transport_equation.compute_Knudsen_diffusivity(cell_id)*Units::convert(1., Units::UNIT2,Units::C_UNIT2);  //converting form cm2/s to m2/s because D_bulk is in m2/s 
            Deff = ficks_transport_equation.effective_diffusion_coefficient(Deff_iso(index_gas,index_solvent),DKnud)*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);  //converting back to cm2/s
        }
        else
            Deff = Deff_iso(index_gas,index_solvent)*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);
        //-------- Looping over Quadrature points ----------------------------
        std::vector<double> JxW_bdry;
        JxW_bdry.resize(n_q_points_bdry);
        
        if (bdry_id==this->user_input_bdry[0])
        {
            for (unsigned int q = 0; q < n_q_points_bdry; ++q)
            {
                JxW_bdry[q] = (bdry_info.fe(xi.fetype_index)).JxW(q);
                dst[0]+= -Deff*concentration*bdry_info.gradients[solution_index][xi.solution_index][q] * bdry_info.fe(xi.fetype_index).normal_vector(q)* JxW_bdry[q];
            }   
        }
    }

}

// ---          ---
// --- evaluate ---
// ---          ---
// This function 
template <int dim>
double
NAME::AppDiffusion<dim>::evaluate (const FuelCell::ApplicationCore::FEVectors& src)
{
    return 0.0;
}

// ---          ---
// --- data_out ---
// ---          ---

template<int dim>
void
NAME::AppDiffusion<dim>::data_out(const std::string& filename,
                                  const FuelCell::ApplicationCore::FEVectors& src)
{
  //////////////
  // SOLUTION //
  //////////////

  // --- Find solution ---
  FuelCell::ApplicationCore::FEVector solution = src.vector( src.find_vector("Solution") );

  // --- Assign solution names ---
  std::vector<std::string> solution_names;

  solution_names.push_back("oxygen_molar_fraction");

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
/////////////////////////////
// EXPLICIT INSTANTIATIONS //
/////////////////////////////
/////////////////////////////

// ---              ---
// --- AppDiffusion ---
// ---              ---
template class NAME::AppDiffusion<deal_II_dimension>;
