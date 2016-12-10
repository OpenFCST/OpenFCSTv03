// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: app_cathode.cc
// - Description: This class describes diffusion in fuel cell cathodes
//                Ficks, one gas
// - Developers: Marc Secanell and Valentin N. Zingan, University of Alberta
//
// ----------------------------------------------------------------------------

#include <applications/app_cathode.h>

namespace NAME  = FuelCell::Application;

//////////////////////////////////////////////////
//////////////////////////////////////////////////
// CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
//////////////////////////////////////////////////
//////////////////////////////////////////////////

// --- constructor ---

template<int dim>
NAME::AppCathode<dim>::AppCathode( boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data )
:
FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>(data),

ficks_transport_equation_cathode(this->system_management, &solute, &solvent, data),
ficks_transport_equation_anode(this->system_management, &solute_anode, &solvent, data),
electron_transport_equation(this->system_management, data),
proton_transport_equation(this->system_management, data),
reaction_source_terms(this->system_management, data),
ORRCurrent(this->system_management),
HORCurrent(this->system_management),
anode(false)
{
  this->repair_diagonal = true;
  FcstUtilities::log << "->FuelCell::Application::AppCathode-" << dim << "D" << std::endl;
}

// ---            ---
// --- destructor ---
// ---            ---

template<int dim>
NAME::AppCathode<dim>::~AppCathode()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::AppCathode<dim>::declare_parameters(ParameterHandler& param)
{
    OptimizationBlockMatrixApplication<dim>::declare_parameters(param);
    
    // Declare application specific parameters:
    param.enter_subsection("Application");
    {
        param.enter_subsection("Electrode");
        {
            param.declare_entry("Anode simulation",
                                "false",
                                Patterns::Bool(),
                                "Set to true if you want to simulate an anode. Otherwise, it is assumed that the user would like to run a cathode case.");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    // Declare parameters in operating conditions:
    OC.declare_parameters(param);
      
    // Declare layer and material classes:
    solute.declare_parameters(param);
    solute_anode.declare_parameters(param);
    solvent.declare_parameters(param);
    
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Gas diffusion layer", param);
    FuelCellShop::Layer::MicroPorousLayer<dim>::declare_MicroPorousLayer_parameters("Microporous layer", param);
    FuelCellShop::Layer::CatalystLayer<dim>::declare_CatalystLayer_parameters("Catalyst layer", param);
    
    // Declare equation classes:
    ficks_transport_equation_anode.declare_parameters(param);
    ficks_transport_equation_cathode.declare_parameters(param);
    electron_transport_equation.declare_parameters(param);
    proton_transport_equation.declare_parameters(param);
    reaction_source_terms.declare_parameters(param);
    
    // Declare post-processing routines
    ORRCurrent.declare_parameters(param);
    HORCurrent.declare_parameters(param);    
    
    // Set new default variables for variables and equations in SystemManagement
    this->set_default_parameters_for_application(param);
}

// ---            ---
// --- initialize ---
// ---            ---
template<int dim>
void
NAME::AppCathode<dim>::initialize(ParameterHandler& param)
{   
    OptimizationBlockMatrixApplication<dim>::initialize(param);
  
    param.enter_subsection("Application");
    {
        param.enter_subsection("Electrode");
        {
            anode = param.get_bool("Anode simulation");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
      
    // Initialize parameters in operating conditions:
    OC.initialize(param);
    
    // Initialize gas classes:      
    std::vector< FuelCellShop::Material::PureGas* > gases;  
    solute.initialize(param);
    solute_anode.initialize(param);
    solvent.initialize(param);
    if (anode)
        gases.push_back(&solute_anode);
    else
        gases.push_back(&solute);
    gases.push_back(&solvent);
    
    // Initialize parameters for physics classes:
    if (anode)
        ficks_transport_equation = &ficks_transport_equation_anode;
    else
        ficks_transport_equation = &ficks_transport_equation_cathode;        
    ficks_transport_equation->initialize(param);
    electron_transport_equation.initialize(param);
    proton_transport_equation.initialize(param);
    
    // Initialize layer classes:
    double pressure = OC.get_pc_atm();
    if (anode)
        pressure = OC.get_pa_atm();
    CGDL = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Gas diffusion layer", param);
    CGDL->set_gases_and_compute(gases, pressure, OC.get_T());
    
    CMPL = FuelCellShop::Layer::MicroPorousLayer<dim>::create_MicroPorousLayer("Microporous layer", param);
    CMPL->set_gases_and_compute(gases, pressure, OC.get_T());
    
    CCL  = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("Catalyst layer", param);
    CCL->set_gases_and_compute(gases, pressure, OC.get_T());
    
    // Initialize the necessary kinetics parameters in CCL.
    ReactionNames name = ORR;
    if (anode)
        name = HOR;
    CCL->set_reaction_kinetics(name);
    CCL->set_constant_solution(Units::convert(pressure, Units::ATM_to_PA), VariableNames::total_pressure);
    CCL->set_constant_solution(OC.get_T(), VariableNames::temperature_of_REV);
    
    // -- Setting kinetics in the reaction source terms object.    
    if (anode)
        reaction_source_terms.set_anode_kinetics(CCL->get_kinetics());
    else
        reaction_source_terms.set_cathode_kinetics(CCL->get_kinetics());
    
    reaction_source_terms.initialize(param);
        
    // --- second of all we make cell couplings for this problem ---
    std::vector<couplings_map> tmp;
    tmp.push_back( ficks_transport_equation->get_internal_cell_couplings()    );
    tmp.push_back( electron_transport_equation.get_internal_cell_couplings() );
    tmp.push_back( proton_transport_equation.get_internal_cell_couplings()   );
    reaction_source_terms.adjust_internal_cell_couplings(tmp);
    this->system_management.make_cell_couplings(tmp);
    
    // Now, initialize object that are used to setup initial solution and boundary conditions:    
    this->component_materialID_value_maps.push_back( ficks_transport_equation->get_component_materialID_value()    );
    this->component_materialID_value_maps.push_back( electron_transport_equation.get_component_materialID_value() );
    this->component_materialID_value_maps.push_back( proton_transport_equation.get_component_materialID_value()   );  
    OC.adjust_initial_solution(this->component_materialID_value_maps, this->mesh_generator);
    //-- OC.adjust_initial_solution is designed for a cathode, so here we make sure anode is ala
    if (anode) 
    {
        for(unsigned int i = 0; i < this->component_materialID_value_maps.size(); ++i)
        {
            for(component_materialID_value_map::iterator iter  = this->component_materialID_value_maps[i].begin(); iter != this->component_materialID_value_maps[i].end(); ++iter)
            {      
                if ((iter->first.compare("protonic_electrical_potential") == 0))
                {
                    std::vector<unsigned int> material_ID = this->mesh_generator->get_material_id("Cathode CL");
                    for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                        iter->second[material_ID[ind]] = 0.0;
                }
            }
        }
    }
    //--
    this->component_boundaryID_value_maps.push_back( ficks_transport_equation->get_component_boundaryID_value() );
    this->component_boundaryID_value_maps.push_back( electron_transport_equation.get_component_boundaryID_value() );
    this->component_boundaryID_value_maps.push_back( proton_transport_equation.get_component_boundaryID_value()   );
    OC.adjust_boundary_conditions(this->component_boundaryID_value_maps, this->mesh_generator);
    
    // --- and then allocate memory for vectors and matrices ---
    this->remesh_matrices();

    // Initialize post-processing routines:
    ORRCurrent.initialize(param);
    HORCurrent.initialize(param);
    
    // Output options:
    // - system info:
    //this->system_management.print_system_info();
    
    // - layers info:
    CGDL->print_layer_properties();
    CMPL->print_layer_properties();
    CCL->print_layer_properties();
    
    // - equations info:
    /*
    ficks_transport_equation.print_equation_info();
    electron_transport_equation.print_equation_info();
    proton_transport_equation.print_equation_info();
    reaction_source_terms.print_equation_info();
    */
    
    FcstUtilities::log << "Theoretical cell voltage: "<<CCL->get_kinetics()->get_cat ()->voltage_cell_th(OC.get_T())<<" V"<<std::endl;

}

// ---               ---
// --- init_solution ---
// ---               ---

template<int dim>
void
NAME::AppCathode<dim>::initialize_solution(FuelCell::ApplicationCore::FEVector& initial_guess,
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
NAME::AppCathode<dim>::cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                     const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
    if( CGDL->belongs_to_material(cell_info.cell->material_id())   )
    {
        ficks_transport_equation->assemble_cell_matrix(cell_matrices, cell_info, CGDL.get());
        electron_transport_equation.assemble_cell_matrix(cell_matrices, cell_info, CGDL.get());
    }
    else if( CMPL->belongs_to_material(cell_info.cell->material_id())   )
    {
        ficks_transport_equation->assemble_cell_matrix(cell_matrices, cell_info, CMPL.get());
        electron_transport_equation.assemble_cell_matrix(cell_matrices, cell_info, CMPL.get());
    }
    else if( CCL->belongs_to_material(cell_info.cell->material_id())    )
    {
        ficks_transport_equation->assemble_cell_matrix(cell_matrices, cell_info, CCL.get());
        electron_transport_equation.assemble_cell_matrix(cell_matrices, cell_info, CCL.get());
        proton_transport_equation.assemble_cell_matrix(cell_matrices, cell_info, CCL.get());
        reaction_source_terms.assemble_cell_matrix(cell_matrices, cell_info, CCL.get());
    }
    else
    {
        FcstUtilities::log<<"Material id: "    <<cell_info.cell->material_id()<<" does not correspond to any layer"<<std::endl;
        Assert( false , ExcNotImplemented() );
    }
}

// ---               ---
// --- cell_residual ---
// ---               ---

template<int dim>
void
NAME::AppCathode<dim>::cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_res,
                                       const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
    if(      CGDL->belongs_to_material(cell_info.cell->material_id())   )
    {
        ficks_transport_equation->assemble_cell_residual(cell_res, cell_info, CGDL.get());
        electron_transport_equation.assemble_cell_residual(cell_res, cell_info, CGDL.get());
    }
    else if( CMPL->belongs_to_material(cell_info.cell->material_id())   )
    {
        ficks_transport_equation->assemble_cell_residual(cell_res, cell_info, CMPL.get());
        electron_transport_equation.assemble_cell_residual(cell_res, cell_info, CMPL.get());
    }
    else if( CCL->belongs_to_material(cell_info.cell->material_id())    )
    {
        ficks_transport_equation->assemble_cell_residual(cell_res, cell_info, CCL.get());
        electron_transport_equation.assemble_cell_residual(cell_res, cell_info, CCL.get());
        proton_transport_equation.assemble_cell_residual(cell_res, cell_info, CCL.get());
        reaction_source_terms.assemble_cell_residual(cell_res, cell_info, CCL.get());
    }
    else
    {
        FcstUtilities::log<<"Material id: "    <<cell_info.cell->material_id()<<" does not correspond to any layer"<<std::endl;
        Assert( false , ExcNotImplemented() );
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
NAME::AppCathode<dim>::dirichlet_bc(std::map<unsigned int, double>& boundary_values) const
{
  FuelCell::InitialAndBoundaryData::make_zero_boundary_values( boundary_values,
                                                              *this->mapping,
                                                              *this->dof,
                                                               this->system_management,
                                                               this->component_boundaryID_value_maps );
}

/////////////////////
/////////////////////
// POST-PROCESSING //
/////////////////////
/////////////////////

// ---                ---
// --- cell_responses ---
// ---                ---

template<int dim>
void
NAME::AppCathode<dim>::cell_responses(std::vector<double>&                                                     dst,
                                        const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                        const FuelCell::ApplicationCore::FEVector& src)
{
  // -- Find out what material is the cell made of, i.e. MEA layer)
    const unsigned int material_id = cell_info.dof_active_cell->material_id();
    
    // Compute ORR responses in the CL
    std::map<FuelCellShop::PostProcessing::ResponsesNames, double> ORR_responses;
    
    if (CCL->belongs_to_material(material_id)) //the material is the catalyst layer
    {
        if (anode)
            HORCurrent.compute_responses(cell_info, CCL.get(), ORR_responses);        
        else
            ORRCurrent.compute_responses(cell_info, CCL.get(), ORR_responses);        
    }
    
    // Organize responses:
    // -- Used to normalize current density (note necessary since it can be done via input file)
    std::vector<double> L_cat_vec = this->mesh_generator->L_cat_c();
    double L_cat = std::accumulate(L_cat_vec.begin(),L_cat_vec.end(),0.0);
    double area_CL = this->mesh_generator->L_channel_c()/2.0 + this->mesh_generator->L_land_c()/2.0;
    double volume_CL = (this->mesh_generator->L_channel_c()/2.0 + this->mesh_generator->L_land_c()/2.0)*L_cat;
    
    for (unsigned int r = 0; r < this->n_resp; ++r)
    {
        if ( (this->name_responses[r] == "current" || this->name_responses[r] == "cathode_current") && CCL->belongs_to_material(material_id))
        {
            if (anode)
                dst[r] += (1/area_CL)*ORR_responses[FuelCellShop::PostProcessing::ResponsesNames::HOR_current];
            else
                dst[r] += (1/area_CL)*ORR_responses[FuelCellShop::PostProcessing::ResponsesNames::ORR_current];
        }
            //(If using response normalization) resp[r] += ORR_responses[FuelCellShop::PostProcessing::ResponsesNames::current];               
        else if (this->name_responses[r] == "OH_coverage" && CCL->belongs_to_material(material_id))
            dst[r] += (1/volume_CL)*ORR_responses[FuelCellShop::PostProcessing::ResponsesNames::OH_coverage];
            //(If using response normalization) resp[r] += ORR_responses[FuelCellShop::PostProcessing::ResponsesNames::OH_coverage];
        else if (this->name_responses[r] == "O_coverage" && CCL->belongs_to_material(material_id))
            dst[r] += (1/volume_CL)*ORR_responses[FuelCellShop::PostProcessing::ResponsesNames::O_coverage];
    }
    
    cell_responses_aux(dst, cell_info, src);
    
}

//-------------------------------------------------------------------
template<int dim>
void
NAME::AppCathode<dim>::cell_responses_aux(std::vector<double>&                          dst,
                                            const typename DoFApplication<dim>::CellInfo& cell_info,
                                            const FEVector& src)
{
    // -- Find out what material is the cell made of, i.e. MEA layer)
    const unsigned int material_id = cell_info.dof_active_cell->material_id();
    
    // Orgaize responses:
    // -- Used to normalize current density (note necessary since it can be done via input file)
    std::vector<double> L_cat_vec = this->mesh_generator->L_cat_c();
    double L_cat = std::accumulate(L_cat_vec.begin(),L_cat_vec.end(),0.0);
    double volume_CL = (this->mesh_generator->L_channel_c()/2.0 + this->mesh_generator->L_land_c()/2.0)*L_cat;
    
    for (unsigned int r = 0; r < this->n_resp; ++r)
    {
        
        if (this->name_responses[r] == "high_total_coverage" && CCL->belongs_to_material(material_id))
        {    
            // Creating solution variable vector to passed to catalyst layer classes
            // Storing solution indices
            unsigned int solution_cell = cell_info.global_data->find_vector("Solution");
            const unsigned int xO2_index = this->system_management.solution_name_to_index("oxygen_molar_fraction");
            const unsigned int phiM_index = this->system_management.solution_name_to_index("protonic_electrical_potential");
            const unsigned int phiS_index = this->system_management.solution_name_to_index("electronic_electrical_potential");
            unsigned int fetype_index = this->system_management.block_info->base_element[xO2_index];
            unsigned int n_q_points_cell = (cell_info.fe(fetype_index)).n_quadrature_points;
            
            std::vector< FuelCellShop::SolutionVariable > solution_variables;
            solution_variables.push_back( FuelCellShop::SolutionVariable(&cell_info.values[solution_cell][xO2_index], oxygen_molar_fraction) );
            std::vector<double> temp_phiM(n_q_points_cell, 0.0);
            solution_variables.push_back( FuelCellShop::SolutionVariable(temp_phiM, protonic_electrical_potential) );
            std::vector<double> temp_phiS(n_q_points_cell, 0.85);
            solution_variables.push_back( FuelCellShop::SolutionVariable(temp_phiS, electronic_electrical_potential) );
            
            std::map<FuelCellShop::PostProcessing::ResponsesNames, double> ORR_responses_temp;
            ORRCurrent.compute_responses(solution_variables, cell_info, CCL.get(), ORR_responses_temp);  
            dst[r] += (1/volume_CL)*(ORR_responses_temp[FuelCellShop::PostProcessing::ResponsesNames::OH_coverage]+ORR_responses_temp[FuelCellShop::PostProcessing::ResponsesNames::O_coverage]);
            
        }
        else if (this->name_responses[r] == "low_total_coverage" && CCL->belongs_to_material(material_id))
        {      
            // Creating solution variable vector to passed to catalyst layer classes
            // Storing solution indices
            unsigned int solution_cell = cell_info.global_data->find_vector("Solution");
            const unsigned int xO2_index = this->system_management.solution_name_to_index("oxygen_molar_fraction");
            const unsigned int phiM_index = this->system_management.solution_name_to_index("protonic_electrical_potential");
            const unsigned int phiS_index = this->system_management.solution_name_to_index("electronic_electrical_potential");
            unsigned int fetype_index = this->system_management.block_info->base_element[xO2_index];
            unsigned int n_q_points_cell = (cell_info.fe(fetype_index)).n_quadrature_points;
            std::vector< FuelCellShop::SolutionVariable > solution_variables;
            solution_variables.push_back( FuelCellShop::SolutionVariable(&cell_info.values[solution_cell][xO2_index], oxygen_molar_fraction) );
            std::vector<double> temp_phiM(n_q_points_cell, 0.0);
            solution_variables.push_back( FuelCellShop::SolutionVariable(temp_phiM, protonic_electrical_potential) );
            std::vector<double> temp_phiS(n_q_points_cell, 0.4);
            solution_variables.push_back( FuelCellShop::SolutionVariable(temp_phiS, electronic_electrical_potential) );
            
            std::map<FuelCellShop::PostProcessing::ResponsesNames, double> ORR_responses_temp;
            ORRCurrent.compute_responses(solution_variables, cell_info, CCL.get(), ORR_responses_temp);  
            dst[r] += (1/volume_CL)*(ORR_responses_temp[FuelCellShop::PostProcessing::ResponsesNames::OH_coverage]+ORR_responses_temp[FuelCellShop::PostProcessing::ResponsesNames::O_coverage]);
            
        }
    }
    
}

// ---          ---
// --- evaluate ---
// ---          ---
template <int dim>
void
NAME::AppCathode<dim>::global_responses(std::vector<double>& resp,
                                          const FuelCell::ApplicationCore::FEVector& /*src*/)
{
    
    std::vector<unsigned int> cathode_CL_material_ids;
    std::map<std::string, double> volume_fractions;
    cathode_CL_material_ids = CCL->get_material_ids();
    
    for (unsigned int i = 0; i < cathode_CL_material_ids.size(); ++i)
    {
        
        CCL->set_local_material_id(cathode_CL_material_ids[i]);
        CCL->get_volume_fractions(volume_fractions);
        std::string epsilonc = "epsilon_V_cat_c:";
        std::stringstream epsilon_c;
        epsilon_c << epsilonc << cathode_CL_material_ids[i];
        
        std::string epsilons = "epsilon_S_cat_c:";
        std::stringstream epsilon_s;
        epsilon_s << epsilons << cathode_CL_material_ids[i];
        
        std::string epsilonn = "epsilon_N_cat_c:";
        std::stringstream epsilon_n;
        epsilon_n << epsilonn << cathode_CL_material_ids[i];
        
        for (unsigned int r = 0; r < this->n_resp; ++r)
        {
            if (this->name_responses[r].compare(epsilon_c.str()) == 0) {
                resp[r] = volume_fractions.find("Void")->second;
            }
            else if (this->name_responses[r].compare(epsilon_s.str()) == 0) {
                resp[r] = volume_fractions.find("Solid")->second;
            }
            else if (this->name_responses[r].compare(epsilon_n.str()) == 0) {
                resp[r] = volume_fractions.find("Ionomer")->second;
            }
        }
    }
}

// ---          ---
// --- evaluate ---
// ---          ---
template <int dim>
double
NAME::AppCathode<dim>::evaluate (const FuelCell::ApplicationCore::FEVectors& src)
{
    std::vector<double> test(this->n_resp, 0.0);
    this->responses(test,
                    src);  
    
    return -test[0];
}

// ---          ---
// --- data_out ---
// ---          ---

template<int dim>
void
NAME::AppCathode<dim>::data_out(const std::string& filename,
                                  const FuelCell::ApplicationCore::FEVectors& src)
{
  //////////////
  // SOLUTION //
  //////////////

  // --- Find solution ---
  FuelCell::ApplicationCore::FEVector solution = src.vector( src.find_vector("Solution") );

  // --- Assign solution names ---
  std::vector<std::string> solution_names;

  if (anode)
      solution_names.push_back("hydrogen_molar_fraction");
  else
      solution_names.push_back("oxygen_molar_fraction");
  solution_names.push_back("protonic_electrical_potential");
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

  // --- current ---
  FuelCellShop::PostProcessing::ORRCurrentDensityDataOut<dim> currentORR(&this->system_management, CCL, &OC);
  FuelCellShop::PostProcessing::HORCurrentDensityDataOut<dim> currentHOR(&this->system_management, CCL, &OC);
  if (anode)
      PostProcessing.push_back(&currentHOR);
  else
      PostProcessing.push_back(&currentORR);

  // --- output ---
  DoFApplication<dim>::data_out( filename,
                                 solution,
                                 solution_names,
                                 PostProcessing);
}

/////////////////////////////
/////////////////////////////
// EXPLICIT INSTANTIATIONS //
/////////////////////////////
/////////////////////////////

// ---              ---
// --- AppCathode ---
// ---              ---
template class NAME::AppCathode<deal_II_dimension>;
