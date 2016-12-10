//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: app_pemfc_multiscale.cc
//    - Description:
//    - Developers: M. Secanell, P. Dobson, M. Bhaiya
//
//---------------------------------------------------------------------------

#include <applications/app_pemfc.h>

namespace NAME = FuelCell::Application;

//---------------------------------------------------------------------------
template <int dim>
NAME::AppPemfc<dim>::AppPemfc(boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data)
  :
  OptimizationBlockMatrixApplication<dim>(data),
  proton_transport(this->system_management, data),
  lambda_transport(this->system_management, data),
  electron_transport(this->system_management, data),
  reaction_source_terms(this->system_management, data),
  sorption_source_terms(this->system_management, data),
  ficks_oxygen_nitrogen(this->system_management, &oxygen, &nitrogen, data),
  ficks_water_nitrogen(this->system_management, &water, &nitrogen, data),
  ficks_water_hydrogen(this->system_management, &water, &hydrogen, data),
  ORRCurrent(this->system_management),
  HORCurrent(this->system_management),
  WaterSorption(this->system_management, &sorption_source_terms)
{
    this->repair_diagonal = true;
    FcstUtilities::log << "FuelCell::Application::AppPemfc-" << dim<<"d"<<std::endl;
}

//---------------------------------------------------------------------------
template <int dim>
NAME::AppPemfc<dim>::~AppPemfc()
{

}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfc<dim>::declare_parameters(ParameterHandler& param)
{
    OptimizationBlockMatrixApplication<dim>::declare_parameters(param);

    // Declare parameters in system management:
    this->system_management.declare_parameters(param);

    // Declare parameters for operating conditions
    OC.declare_parameters(param);

    // Declare parameters for layer classes:
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Anode gas diffusion layer", param);
    FuelCellShop::Layer::MicroPorousLayer<dim>::declare_MicroPorousLayer_parameters("Anode microporous layer", param);
    FuelCellShop::Layer::CatalystLayer<dim>::declare_CatalystLayer_parameters("Anode catalyst layer", param);
    FuelCellShop::Layer::MembraneLayer<dim>::declare_MembraneLayer_parameters("Membrane layer", param);
    FuelCellShop::Layer::CatalystLayer<dim>::declare_CatalystLayer_parameters("Cathode catalyst layer", param);
    FuelCellShop::Layer::MicroPorousLayer<dim>::declare_MicroPorousLayer_parameters("Cathode microporous layer", param);
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Cathode gas diffusion layer", param);

    // Declare parameters for physics classes:
    proton_transport.declare_parameters(param);
    lambda_transport.declare_parameters(param);
    electron_transport.declare_parameters(param);
    reaction_source_terms.declare_parameters(param);
    sorption_source_terms.declare_parameters(param);
    ficks_oxygen_nitrogen.declare_parameters(param);
    ficks_water_hydrogen.declare_parameters(param);
    ficks_water_nitrogen.declare_parameters(param);
    
    // Declare post-processing routines
    ORRCurrent.declare_parameters(param);
    HORCurrent.declare_parameters(param);
    
    // Set new default variables for variables and equations in SystemManagement
    this->set_default_parameters_for_application(param);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfc<dim>::initialize(ParameterHandler& param)
{
    OptimizationBlockMatrixApplication<dim>::initialize(param);
    
    // Initialize coefficients (NOTE: grid already initialized)
    OC.initialize(param);
    
    // Initialize gases and material classes:
    std::vector< FuelCellShop::Material::PureGas* > cathode_gases;
    cathode_gases.push_back(&oxygen);
    cathode_gases.push_back(&water);
    cathode_gases.push_back(&nitrogen);
    std::vector< FuelCellShop::Material::PureGas* > anode_gases;
    anode_gases.push_back(&water);
    anode_gases.push_back(&hydrogen);
    
    // Initialize layer classes:
    AGDL = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Anode gas diffusion layer",param);
    AGDL->set_gases_and_compute(anode_gases, OC.get_pa_atm (), OC.get_T());
    
    AMPL = FuelCellShop::Layer::MicroPorousLayer<dim>::create_MicroPorousLayer("Anode microporous layer",param);
    AMPL->set_gases_and_compute(anode_gases, OC.get_pa_atm (), OC.get_T());
    
    ACL = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("Anode catalyst layer", param);
    ACL->set_gases_and_compute (anode_gases, OC.get_pa_atm (), OC.get_T());
    
    ML = FuelCellShop::Layer::MembraneLayer<dim>::create_MembraneLayer("Membrane layer", param);
    
    CCL = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("Cathode catalyst layer", param);
    CCL->set_gases_and_compute (cathode_gases, OC.get_pc_atm (), OC.get_T());
    
    CMPL = FuelCellShop::Layer::MicroPorousLayer<dim>::create_MicroPorousLayer("Cathode microporous layer",param);
    CMPL->set_gases_and_compute(cathode_gases, OC.get_pc_atm (), OC.get_T());
    
    CGDL = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Cathode gas diffusion layer",param);
    CGDL->set_gases_and_compute(cathode_gases, OC.get_pc_atm (), OC.get_T());
    
    // Initialize necessary kinetics parameters in catalyst layers.
    ReactionNames anode_rxn_name = HOR;
    ACL->set_reaction_kinetics(anode_rxn_name);
    ReactionNames cathode_rxn_name = ORR;
    CCL->set_reaction_kinetics(cathode_rxn_name);
    
    // Setting required constant solutions, viz., temperature and pressure in the catalyst layer and membrane layer objects.
    ACL->set_constant_solution(OC.get_pa_Pa(), total_pressure);   //Note: Pressure and temp. can only be set after setting kinetics.
    ACL->set_constant_solution(OC.get_T(), temperature_of_REV);
    CCL->set_constant_solution(OC.get_pc_Pa(), total_pressure);   //Note: Pressure and temp. can only be set after setting kinetics.
    CCL->set_constant_solution(OC.get_T(), temperature_of_REV);
    ML->set_constant_solution(OC.get_T(), temperature_of_REV);
    
    // Setting kinetics in the reaction source terms object.
    reaction_source_terms.set_cathode_kinetics(CCL->get_kinetics());
    reaction_source_terms.set_anode_kinetics(ACL->get_kinetics());
    
    // Initialize parameters for physics classes:
    ficks_oxygen_nitrogen.initialize(param);
    ficks_water_hydrogen.initialize(param);
    ficks_water_nitrogen.initialize(param);
    proton_transport.initialize(param);
    lambda_transport.initialize(param);
    electron_transport.initialize(param);
    reaction_source_terms.initialize(param);
    sorption_source_terms.initialize(param);
    
    // Initialize system of equations based on the physical phenomena being solved for:
    std::vector<couplings_map> tmp;
    tmp.push_back(ficks_water_nitrogen.get_internal_cell_couplings() );
    tmp.push_back(ficks_water_hydrogen.get_internal_cell_couplings() );
    tmp.push_back(ficks_oxygen_nitrogen.get_internal_cell_couplings() );
    tmp.push_back(proton_transport.get_internal_cell_couplings() );
    tmp.push_back(electron_transport.get_internal_cell_couplings() );
    tmp.push_back(lambda_transport.get_internal_cell_couplings() );
    reaction_source_terms.adjust_internal_cell_couplings(tmp);
    sorption_source_terms.adjust_internal_cell_couplings(tmp);
    this->system_management.make_cell_couplings(tmp);
    
    // Now, initialize object that are used to setup initial solution and boundary conditions:    
    this->component_materialID_value_maps.push_back( ficks_oxygen_nitrogen.get_component_materialID_value()    );
    this->component_materialID_value_maps.push_back( ficks_water_hydrogen.get_component_materialID_value() );
    this->component_materialID_value_maps.push_back( ficks_water_nitrogen.get_component_materialID_value()   );
    this->component_materialID_value_maps.push_back( proton_transport.get_component_materialID_value()   );
    this->component_materialID_value_maps.push_back( electron_transport.get_component_materialID_value()   );
    this->component_materialID_value_maps.push_back( lambda_transport.get_component_materialID_value()   );
    OC.adjust_initial_solution(this->component_materialID_value_maps, this->mesh_generator);
    
    this->component_boundaryID_value_maps.push_back( ficks_oxygen_nitrogen.get_component_boundaryID_value()    );
    this->component_boundaryID_value_maps.push_back( ficks_water_hydrogen.get_component_boundaryID_value() );
    this->component_boundaryID_value_maps.push_back( ficks_water_nitrogen.get_component_boundaryID_value()   );
    this->component_boundaryID_value_maps.push_back( proton_transport.get_component_boundaryID_value()   );
    this->component_boundaryID_value_maps.push_back( electron_transport.get_component_boundaryID_value()   );
    this->component_boundaryID_value_maps.push_back( lambda_transport.get_component_boundaryID_value()   );
    OC.adjust_boundary_conditions(this->component_boundaryID_value_maps, this->mesh_generator);
    
    // Initialize matrices and spartisity pattern for the whole system
    this->remesh_matrices();
    
    // Initialize post-processing routines:
    ORRCurrent.initialize(param);
    HORCurrent.initialize(param);
    WaterSorption.initialize(param);
    
    // Output CL properties
    OC.print_operating_conditions();
    ACL->print_layer_properties();
    CCL->print_layer_properties();
    
    // FOR DEBUGGING PURPOSES ONLY:
    //   std::ofstream file;
    //   file.open("Modified_parameter_file.prm");
    //   param.print_parameters (file, ParameterHandler::XML);
    //   file.close();
    
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfc<dim>::initialize_solution(FuelCell::ApplicationCore::FEVector& initial_guess,
                                         std::shared_ptr<Function<dim> > initial_function)
{
    std::shared_ptr< Function<dim> > initial_solution (new FuelCell::InitialSolution::AppPemfcIC<dim> (&OC, this->mesh_generator));
    
    DoFApplication<dim>::initialize_solution(initial_guess, initial_solution);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfc<dim>::cell_matrix(MatrixVector& cell_matrices,
                                      const typename DoFApplication<dim>::CellInfo& info)
{

    const unsigned int material_id = info.cell->material_id();

    if ( CGDL->belongs_to_material(material_id) )
    {
        ficks_oxygen_nitrogen.assemble_cell_matrix(cell_matrices, info, CGDL.get());
        ficks_water_nitrogen.assemble_cell_matrix(cell_matrices, info, CGDL.get());
        electron_transport.assemble_cell_matrix(cell_matrices, info, CGDL.get());
    }
    else if ( CMPL->belongs_to_material(material_id) )
    {
        ficks_oxygen_nitrogen.assemble_cell_matrix(cell_matrices, info, CMPL.get());
        ficks_water_nitrogen.assemble_cell_matrix(cell_matrices, info, CMPL.get());
        electron_transport.assemble_cell_matrix(cell_matrices, info, CMPL.get());
    }
    else if ( CCL->belongs_to_material(material_id) )
    {
        ficks_oxygen_nitrogen.assemble_cell_matrix(cell_matrices, info, CCL.get());
        ficks_water_nitrogen.assemble_cell_matrix(cell_matrices, info, CCL.get());
        proton_transport.assemble_cell_matrix(cell_matrices, info, CCL.get());
        lambda_transport.assemble_cell_matrix(cell_matrices, info, CCL.get());
        electron_transport.assemble_cell_matrix(cell_matrices, info, CCL.get());
        reaction_source_terms.assemble_cell_matrix(cell_matrices, info, CCL.get());
        sorption_source_terms.assemble_cell_matrix(cell_matrices, info, CCL.get());
    }
    else if ( ML->belongs_to_material(material_id) )
    {
        proton_transport.assemble_cell_matrix(cell_matrices, info, ML.get());
        lambda_transport.assemble_cell_matrix(cell_matrices, info, ML.get());
    }
    else if ( ACL->belongs_to_material(material_id) )
    {
        ficks_water_hydrogen.assemble_cell_matrix(cell_matrices, info, ACL.get());
        proton_transport.assemble_cell_matrix(cell_matrices, info, ACL.get());
        lambda_transport.assemble_cell_matrix(cell_matrices, info, ACL.get());
        electron_transport.assemble_cell_matrix(cell_matrices, info, ACL.get());
        reaction_source_terms.assemble_cell_matrix(cell_matrices, info, ACL.get());
        sorption_source_terms.assemble_cell_matrix(cell_matrices, info, ACL.get());
    }
    else if ( AMPL->belongs_to_material(material_id) )
    {
        ficks_water_hydrogen.assemble_cell_matrix(cell_matrices, info, AMPL.get());
        electron_transport.assemble_cell_matrix(cell_matrices, info, AMPL.get());
    }
    else if ( AGDL->belongs_to_material(material_id) )
    {
        ficks_water_hydrogen.assemble_cell_matrix(cell_matrices, info, AGDL.get());
        electron_transport.assemble_cell_matrix(cell_matrices, info, AGDL.get());
    }
    else
        FcstUtilities::log<<"Material id: "    <<info.cell->material_id()<<" does not correspond to any layer"<<std::endl;
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfc<dim>::cell_residual(FuelCell::ApplicationCore::FEVector& cell_vector,
                                        const typename DoFApplication<dim>::CellInfo& info)
{

    // -- Assertion before starting routine:
    // Make sure vectors are the right size
    Assert (cell_vector.n_blocks() == this->element->n_blocks(),
                        ExcDimensionMismatch (cell_vector.n_blocks(), this->element->n_blocks()));

    // -- Find out what material is the cell made of, i.e. MEA layer)
    const unsigned int material_id = info.cell->material_id();

    if ( CGDL->belongs_to_material(material_id) )
    {
        ficks_oxygen_nitrogen.assemble_cell_residual(cell_vector, info, CGDL.get());
        ficks_water_nitrogen.assemble_cell_residual(cell_vector, info, CGDL.get());
        electron_transport.assemble_cell_residual(cell_vector, info, CGDL.get());
        
    }
    else if ( CMPL->belongs_to_material(material_id) )
    {
        ficks_oxygen_nitrogen.assemble_cell_residual(cell_vector, info, CMPL.get());
        ficks_water_nitrogen.assemble_cell_residual(cell_vector, info, CMPL.get());
        electron_transport.assemble_cell_residual(cell_vector, info, CMPL.get());
    }
    else if ( CCL->belongs_to_material(material_id) )
    {
        ficks_oxygen_nitrogen.assemble_cell_residual(cell_vector, info, CCL.get());
        ficks_water_nitrogen.assemble_cell_residual(cell_vector, info, CCL.get());
        proton_transport.assemble_cell_residual(cell_vector, info, CCL.get());
        lambda_transport.assemble_cell_residual(cell_vector, info, CCL.get());
        electron_transport.assemble_cell_residual(cell_vector, info, CCL.get());
        reaction_source_terms.assemble_cell_residual(cell_vector, info, CCL.get());
        sorption_source_terms.assemble_cell_residual(cell_vector, info, CCL.get());
    }
    else if ( ML->belongs_to_material(material_id) )
    {
        proton_transport.assemble_cell_residual(cell_vector, info, ML.get());
        lambda_transport.assemble_cell_residual(cell_vector, info, ML.get());
    }
    else if ( ACL->belongs_to_material(material_id) )
    {
        ficks_water_hydrogen.assemble_cell_residual(cell_vector, info, ACL.get());
        proton_transport.assemble_cell_residual(cell_vector, info, ACL.get());
        lambda_transport.assemble_cell_residual(cell_vector, info, ACL.get());
        electron_transport.assemble_cell_residual(cell_vector, info, ACL.get());
        reaction_source_terms.assemble_cell_residual(cell_vector, info, ACL.get());
        sorption_source_terms.assemble_cell_residual(cell_vector, info, ACL.get());
    }
    else if ( AMPL->belongs_to_material(material_id) )
    {
        ficks_water_hydrogen.assemble_cell_residual(cell_vector, info, AMPL.get());
        electron_transport.assemble_cell_residual(cell_vector, info, AMPL.get());
    }
    else if ( AGDL->belongs_to_material(material_id) )
    {
        ficks_water_hydrogen.assemble_cell_residual(cell_vector, info, AGDL.get());
        electron_transport.assemble_cell_residual(cell_vector, info, AGDL.get());
    }
    else
        AssertThrow(false, ExcNotImplemented());

}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfc<dim>::dirichlet_bc(std::map<unsigned int, double>& boundary_values) const
{
      FuelCell::InitialAndBoundaryData::make_zero_boundary_values( boundary_values,
                                                                   *this->mapping,
                                                                   *this->dof,
                                                                   this->system_management,
                                                                   this->component_boundaryID_value_maps );
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//-----------    OPTIMIZATION ROUTINES     ----------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfc<dim>::check_responses()
{
//   for (unsigned int r = 0; r < this->n_resp; ++r)
//     {
//       if (this->name_responses[r] != "current" &&
//    this->name_responses[r] != "m_Pt_c" &&
//    this->name_responses[r] != "epsilon_V_cat_c" &&
//    this->name_responses[r] != "epsilon_N_cat_c" &&
//    this->name_responses[r] != "epsilon_S_cat_c" &&
//    this->name_responses[r] != "m_Pt_a" &&
//    this->name_responses[r] != "epsilon_V_cat_a" &&
//    this->name_responses[r] != "epsilon_N_cat_a" &&
//    this->name_responses[r] != "epsilon_S_cat_a" &&
//    this->name_responses[r] != "t_m_Pt" &&
//    this->name_responses[r] != "l_cat_c" &&
//    this->name_responses[r] != "l_cat_a")
//  {
//    FcstUtilities::log <<"Response "<<this->name_responses[r]<<" is not implemented"
//       <<"Error in file "<<__FILE__<<"line "<<__LINE__<< std::endl;
//    // Throw an exception:
//    ExcNotImplemented ();
//  }
//     }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfc<dim>::cell_responses (std::vector<double>& resp,
                                          const typename DoFApplication<dim>::CellInfo& info,
                                          const FuelCell::ApplicationCore::FEVector& /*src*/)
{
    // TODO: ISSUE HERE
    // Info needed during integration. Note this might be a problem when the geometry is supplied from input file
    l_channel = this->mesh_generator->L_channel_c();
    l_land = this->mesh_generator->L_land_c();
    
    // -- Find out what material is the cell made of, i.e. MEA layer)
    const unsigned int material_id = info.dof_active_cell->material_id();

    // -- Object used to store responses:
    std::map<FuelCellShop::PostProcessing::ResponsesNames, double> responses;
                
    for (unsigned int r = 0; r < this->n_resp; ++r)
    {
        if ( (this->name_responses[r] == "cathode_current" || this->name_responses[r] == "current") && (CCL->belongs_to_material(material_id)) )
        {
            // Compute ORR responses in the CL
            ORRCurrent.compute_responses(info, CCL.get(), responses);   
            resp[r] += responses[FuelCellShop::PostProcessing::ResponsesNames::ORR_current]/ (l_channel/2.0 + l_land/2.0);
        }

        else if ( (this->name_responses[r] == "anode_current") && (ACL->belongs_to_material(material_id)) )
        {
            
            HORCurrent.compute_responses(info, ACL.get(), responses);   
            resp[r] += responses[FuelCellShop::PostProcessing::ResponsesNames::HOR_current]/ (l_channel/2.0 + l_land/2.0);
        }

        else if ( (this->name_responses[r] == "water_cathode") && (CCL->belongs_to_material(material_id)) )
        {
            
            WaterSorption.compute_responses(info, CCL.get(), responses);
            resp[r] += responses[FuelCellShop::PostProcessing::ResponsesNames::sorbed_water]/ (l_channel/2.0 + l_land/2.0);
        }

        else if ( (this->name_responses[r] == "water_anode") && (ACL->belongs_to_material(material_id)) )
        {
            WaterSorption.compute_responses(info, ACL.get(), responses);
            resp[r] += responses[FuelCellShop::PostProcessing::ResponsesNames::sorbed_water]/ (l_channel/2.0 + l_land/2.0);
        }
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfc<dim>::global_responses(std::vector<double>& resp,
        const FuelCell::ApplicationCore::FEVector& /*src*/)
{
/*
    std::map<std::string, double> volume_fractions;
    std::map<std::string, double> loadings;
    const std::vector<unsigned int> cathode_CL_material_ids (CCL->get_material_ids());
    
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
        
        std::string loading_name = "l_cat_c:";
        std::stringstream loading_n;
        loading_n << loading_name << cathode_CL_material_ids[i];
        
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
            else if (this->name_responses[r].compare(loading_n.str()) == 0) {
                resp[r] = loadings.find("V_Pt")->second/this->mesh_generator->L_cat_c();
            }
        }
    }
    
    
    // Anode
    const std::vector<unsigned int> anode_CL_material_ids (ACL->get_material_ids());
    
    for (unsigned int i = 0; i < anode_CL_material_ids.size(); ++i)
    {
        
        ACL->set_local_material_id(anode_CL_material_ids[i]);
        ACL->get_volume_fractions(volume_fractions);
        ACL->get_loadings(loadings);
        std::string epsilonc = "epsilon_V_cat_a:";
        std::stringstream epsilon_c;
        epsilon_c << epsilonc << anode_CL_material_ids[i];
        
        std::string epsilons = "epsilon_S_cat_a:";
        std::stringstream epsilon_s;
        epsilon_s << epsilons << anode_CL_material_ids[i];
        
        std::string epsilonn = "epsilon_N_cat_a:";
        std::stringstream epsilon_n;
        epsilon_n << epsilonn << anode_CL_material_ids[i];
        
        std::string loading_name = "l_cat_a:";
        std::stringstream loading_n;
        loading_n << loading_name << anode_CL_material_ids[i];
        
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
            else if (this->name_responses[r].compare(loading_n.str()) == 0) {
                resp[r] = loadings.find("V_Pt")->second/this->mesh_generator->L_cat_a();
            }
        }
    }
    
    // Total loading:
    */
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfc<dim>::global_dresponses_dl(std::vector<std::vector<double> >& df_dl,
           const FuelCell::ApplicationCore::FEVector& )
{}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfc<dim>::global_dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector >& df_du,
                                               const FuelCell::ApplicationCore::FEVector& /*src*/)
{}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//-----------    POSTPROCESSING ROUTINES     --------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
template <int dim>
double
NAME::AppPemfc<dim>::evaluate (const FuelCell::ApplicationCore::FEVectors& vectors)
{
    std::vector<double> test(this->n_resp, 0.0);
    this->responses(test, vectors);

    for (unsigned int r = 0; r < this->n_resp; ++r)
        FcstUtilities::log << "Response " << this->name_responses[r] << " is: :" << test[r] << std::endl;

    return 0.0;
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfc<dim>::data_out(const std::string &basename,
                              const FuelCell::ApplicationCore::FEVectors &vectors)
{
    // --- output file name:
    // Create string that will be appended to the .vtk filename. This string will contain the name and value of the design variables.
    std::string filename1;
    
    //loop over all design variables
    for (unsigned int i=0; i<design_var.size(); i++)
    {
        //create object that can print out numbers
        std::ostringstream streamOut;
        //Used to set how many digits after the decimal point are used. Insert the required number into the brackets.
        streamOut.precision();
        //If there are fewer digits in the design variable than are requested, the remaining digits can be filled with zeroes.
        //Use std::fixed to force this: streamOut <<  std::fixed << design_var_value[design_var.size()-i-1];
        //This can be used to get paraview to group the vtk files correctly.
        
        //Note that the design variables are printed in reverse, so that file managers will group the vtk files by the last design
        // variable first before moving to the second last etc.
        streamOut << std::fixed <<  design_var_value[design_var.size()-i-1];
        
        //Append design variable name and value to filename1
        filename1 += "__" + design_var[design_var.size()-i-1] + "_" + streamOut.str();
    }
    //construct the full filename, where basename is given by the adaptive refinement loop and suffix is the file extension (i.e. .vtk)
    const std::string fileName = basename + filename1;
    FcstUtilities::log << "Datafile:" << fileName + this->d_out.default_suffix() << std::endl;
    
    // --- Find solution and do post-processing, otherwise skip:
    FuelCell::ApplicationCore::FEVector solution;
    unsigned int index = vectors.find_vector("Solution");
    
    if (index != static_cast<unsigned int>(-1))
    {
        solution = vectors.vector(index);
        // --- Assign solution interpretations
        this->solution_interpretations.clear();
        this->solution_interpretations.resize(this->element->n_blocks(),
                                              DataComponentInterpretation::component_is_scalar);
        
        // --- Create vector of PostProcessing objects ---
        std::vector< DataPostprocessor<dim>* > PostProcessing;
        
        // --- current ---
        FuelCellShop::PostProcessing::ORRCurrentDensityDataOut<dim> currentORR(&this->system_management, CCL, &OC);
        PostProcessing.push_back(&currentORR);
        
        FuelCellShop::PostProcessing::HORCurrentDensityDataOut<dim> currentHOR(&this->system_management, ACL, &OC);
        PostProcessing.push_back(&currentHOR);
        
        // --- relative humidity
        std::vector< boost::shared_ptr< FuelCellShop::Layer::PorousLayer<dim> > > porous_layers;
        porous_layers.push_back(CGDL);
        porous_layers.push_back(CMPL);
        porous_layers.push_back(CCL);
        porous_layers.push_back(ACL);
        porous_layers.push_back(AMPL);
        porous_layers.push_back(AGDL);
        FuelCellShop::PostProcessing::RelativeHumidityDataOut<dim> relative_humidity (&this->system_management, porous_layers, &OC);   
        PostProcessing.push_back(&relative_humidity);
        
        // --- data out ---
        DoFApplication<dim>::data_out( fileName,
                                       solution,
                                       this->system_management.get_solution_names(),
                                       PostProcessing);
    }
    else
        FcstUtilities::log<<"Warning: Solution vector could not be found in AppPemfc<dim>::data_out"<<std::endl;
    
    //-- If solution could not be found try update:
    index = vectors.find_vector("update");
    FuelCell::ApplicationCore::FEVector update;
    
    if (index != static_cast<unsigned int>(-1) && solution.size() == 0)
    {
        update = vectors.vector(index);
        // --- data out ---
        DoFApplication<dim>::data_out( fileName,
                                       update,
                                       this->system_management.get_solution_names());
    }
    //-- If solution and update could not be found try residual:
    index = vectors.find_vector("residual");
    FuelCell::ApplicationCore::FEVector residual;
    
    if (index != static_cast<unsigned int>(-1) && solution.size() == 0 && update.size() == 0)
    {
        residual = vectors.vector(index);
        // --- data out ---
        DoFApplication<dim>::data_out( fileName,
                                       residual,
                                       this->system_management.get_solution_names());
    }
    
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
FuelCell::InitialSolution::AppPemfcIC<dim>::AppPemfcIC(FuelCell::OperatingConditions* OC,
                                                                 boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid)
  :
  Function<dim> (5),
  OC(OC),
  grid(grid)
{

}

//---------------------------------------------------------------------------
template <int dim>
FuelCell::InitialSolution::AppPemfcIC<dim>::~AppPemfcIC()
{}

//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::InitialSolution::AppPemfcIC<dim>::vector_value (const Point<dim> &p,
                                                               Vector<double> &v) const
{
    // size checking for vectors
    Assert(v.size() == this->n_components,
           ExcDimensionMismatch (v.size(), this->n_components));
//
    double RH_c = OC->get_RH_c();
    double RH_a = OC->get_RH_a();
    double V_cell = std::fabs(OC->get_V());
    double x_wv = OC->get_x_wv();
    double x_o2 = OC->get_x_o2();
    double x_h2 = OC->get_x_h2();
    double OCV = OC->get_OCV();

    double l_mem = grid->L_mem();
    double l_cat_a = grid->L_cat_a();
    double l_mpl_a = grid->L_mpl_a();
    double l_gdl_a = grid->L_gdl_a();

    double x = p(0);
    double x_norm = (x - l_gdl_a - l_cat_a - l_mpl_a)/l_mem;
    double delta = -0.0001;

    if (x <= l_gdl_a + l_cat_a + l_mpl_a) //Anode
    {
        v(0) = 0.0;//0.5*x_o2;
        v(1) = 1 - x_h2;  // x_wv_a - Dirichlet B.C.
        v(2) = delta; // phi_m
        v(3) = 0;        // phi_s - Dirichlet B.C.
        if(RH_a <= 0.5)
            v(4) = 4.0;
        else
            v(4) = 7.0; //lamda
    }
    else if (x > l_gdl_a + l_cat_a + l_mpl_a && x < l_gdl_a + l_cat_a + l_mpl_a + l_mem) //membrane
    {
        v(0) = 0.0;//.5*x_o2;
        v(1) = x_wv;
        v(2) = ((V_cell-OCV)/2.0 - (delta))*x_norm + delta;
        v(3) = V_cell * x_norm;
        if(RH_a <= 0.5 && RH_c <= 0.5)
            v(4) = 4.0;
        else if(RH_a <= 0.5 && RH_c > 0.5)
            v(4) = 3.0*x_norm + 4.0;
        else if(RH_a > 0.5 && RH_c <= 0.5)
            v(4) = -3.0*x_norm + 7.0;
        else
            v(4) = 7.0;
    }
    else if (x >= l_gdl_a + l_cat_a + l_mpl_a + l_mem) //Cathode
    {
        v(0) = x_o2;  //Dirichlet B.C.
        v(1) = x_wv;  //Dirichlet B.C.
        v(2) = (V_cell-OCV)/2.0;  // phi_m - i.e. 1/2 of the losses through the CCL
        v(3) = V_cell;   // phi_s - Dirichlet B.C.
        // lamda
        if (RH_c <= 0.5)
            v(4) = 4.0;
        else
            v(4) = 7.0;
    }
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::AppPemfc<deal_II_dimension>;
template class FuelCell::InitialSolution::AppPemfcIC<deal_II_dimension>;
