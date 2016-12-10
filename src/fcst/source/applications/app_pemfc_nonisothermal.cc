//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: app_pemfc_nonisothermal.cc 21-05-2013
//    - Description: Class designed to solve a non-isothermal PEMFC model.
//    - Developers: Madhur Bhaiya
//    - Id: $Id: app_pemfc_nonisothermal.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <applications/app_pemfc_nonisothermal.h>

namespace NAME = FuelCell::Application;

using namespace FuelCell::ApplicationCore;

//---------------------------------------------------------------------------
template <int dim>
NAME::AppPemfcNIThermal<dim>::AppPemfcNIThermal(boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data)
:
OptimizationBlockMatrixApplication<dim>(data),
thermal_transport(this->system_management, data),
proton_transport(this->system_management, data),
lambda_transport(this->system_management, data),
electron_transport(this->system_management, data),
reaction_source_terms(this->system_management, data),
sorption_source_terms(this->system_management, data),
ficks_oxygen_nitrogen(this->system_management, &oxygen, &nitrogen, data),
ficks_water_nitrogen(this->system_management,  &water,  &nitrogen, data),
ficks_water_hydrogen(this->system_management,  &water,  &hydrogen, data),
ORRCurrent(this->system_management),
HORCurrent(this->system_management),
electronOhmicHeat(this->system_management, &thermal_transport),
protonOhmicHeat(this->system_management, &thermal_transport),
sorptionHeat(this->system_management, &sorption_source_terms),
catReactionHeat(this->system_management, &reaction_source_terms),
anReactionHeat(this->system_management, &reaction_source_terms),
waterSorption(this->system_management, &sorption_source_terms)
{
    this->repair_diagonal = true;
    FcstUtilities::log << "FuelCell::Application::AppPemfc_NonIsothermal-" << dim<<"d"<<std::endl;
}

//---------------------------------------------------------------------------
template <int dim>
NAME::AppPemfcNIThermal<dim>::~AppPemfcNIThermal()
{ 
}
     
//---------------------------------------------------------------------------
template <int dim>
void 
NAME::AppPemfcNIThermal<dim>::declare_parameters(ParameterHandler& param)
{
    OptimizationBlockMatrixApplication<dim>::declare_parameters(param);
    
    // Declare parameters in system management:
    this->system_management.declare_parameters(param);
    
    // Declare parameters for operating conditions
    OC.declare_parameters(param);
    
    // Declare parameters for layer classes:
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Cathode gas diffusion layer", param);
    FuelCellShop::Layer::MicroPorousLayer<dim>::declare_MicroPorousLayer_parameters("Cathode microporous layer", param);
    FuelCellShop::Layer::CatalystLayer<dim>::declare_CatalystLayer_parameters("Cathode catalyst layer", param);
    FuelCellShop::Layer::MembraneLayer<dim>::declare_MembraneLayer_parameters("Membrane layer", param);
    FuelCellShop::Layer::CatalystLayer<dim>::declare_CatalystLayer_parameters("Anode catalyst layer", param);
    FuelCellShop::Layer::MicroPorousLayer<dim>::declare_MicroPorousLayer_parameters("Anode microporous layer", param);
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Anode gas diffusion layer", param);
    
    // Declare parameters for physics classes:
    thermal_transport.declare_parameters(param);
    proton_transport.declare_parameters(param);
    lambda_transport.declare_parameters(param);
    electron_transport.declare_parameters(param);
    reaction_source_terms.declare_parameters(param);
    sorption_source_terms.declare_parameters(param);
    ficks_oxygen_nitrogen.declare_parameters(param);
    ficks_water_nitrogen.declare_parameters(param);
    ficks_water_hydrogen.declare_parameters(param);
    
    // Declare current post processing responses classes:
    ORRCurrent.declare_parameters(param);
    HORCurrent.declare_parameters(param);
    
    // Set new default variables for variables and equations in SystemManagement
    this->set_default_parameters_for_application(param);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfcNIThermal<dim>::_initialize(ParameterHandler& param)
{
    //Initialize the problem data:
    this->system_management.initialize(param);
    
    // Make sure that the number of input finite elements matches with the problem:
    AssertThrow (this->element->n_blocks() == this->system_management.get_equation_names().size(),
                 ExcDimensionMismatch(this->element->n_blocks(), this->system_management.get_equation_names().size()));
    AssertThrow (this->element->n_blocks() == this->system_management.get_number_of_solution_names(),
                 ExcDimensionMismatch(this->element->n_blocks(), this->system_management.get_number_of_solution_names()));
    
    // Initialize coefficients
    OC.initialize(param);
    
    // Initialize gases and material classes
    std::vector< FuelCellShop::Material::PureGas* > cathode_gases;
    cathode_gases.push_back(&oxygen);
    cathode_gases.push_back(&water);
    cathode_gases.push_back(&nitrogen);
    std::vector< FuelCellShop::Material::PureGas* > anode_gases;
    anode_gases.push_back(&water);
    anode_gases.push_back(&hydrogen);
    
    // Initialize layer classes:    
    CGDL = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Cathode gas diffusion layer",param);
    CGDL->set_gases (cathode_gases, OC.get_pc_atm());
    
    CMPL = FuelCellShop::Layer::MicroPorousLayer<dim>::create_MicroPorousLayer("Cathode microporous layer",param);
    CMPL->set_gases (cathode_gases, OC.get_pc_atm());

    CCL = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("Cathode catalyst layer", param);
    CCL->set_gases (cathode_gases, OC.get_pc_atm());
    
    ML = FuelCellShop::Layer::MembraneLayer<dim>::create_MembraneLayer("Membrane layer", param);
    
    ACL = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("Anode catalyst layer", param);
    ACL->set_gases (anode_gases, OC.get_pa_atm());

    AMPL = FuelCellShop::Layer::MicroPorousLayer<dim>::create_MicroPorousLayer("Anode microporous layer",param);
    AMPL->set_gases (anode_gases, OC.get_pa_atm());
    
    AGDL = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Anode gas diffusion layer",param);
    AGDL->set_gases (anode_gases, OC.get_pa_atm());
    
    // Initialise the necessary kinetics parameters.
    ReactionNames cathode_rxn_name = ORR;
    CCL->set_reaction_kinetics(cathode_rxn_name);
    ReactionNames  anode_rxn_name = HOR;
    ACL->set_reaction_kinetics(anode_rxn_name);
    // Set the kinetics inside the catalyst layers
    CCL->set_constant_solution(OC.get_pc_Pa(), total_pressure);  // Note: pressure can only be set after setting the kinetics  
    ACL->set_constant_solution(OC.get_pa_Pa(), total_pressure);
    
    // Setting kinetics in the reaction source terms object.
    reaction_source_terms.set_cathode_kinetics(CCL->get_kinetics());
    reaction_source_terms.set_anode_kinetics(ACL->get_kinetics());

    // Initializing Equation classes
    thermal_transport.initialize(param);
    proton_transport.initialize(param);
    lambda_transport.initialize(param);
    electron_transport.initialize(param);
    reaction_source_terms.initialize(param);
    sorption_source_terms.initialize(param);
    ficks_oxygen_nitrogen.initialize(param);
    ficks_water_nitrogen.initialize(param);
    ficks_water_hydrogen.initialize(param);
    
    time_k = sorption_source_terms.get_time_constant();
    l_channel = this->mesh_generator->L_channel_c();
    l_land = this->mesh_generator->L_land_c();
    
    std::vector<couplings_map> tmp;
    tmp.push_back( thermal_transport.get_internal_cell_couplings() );
    tmp.push_back( proton_transport.get_internal_cell_couplings() );
    tmp.push_back( electron_transport.get_internal_cell_couplings() );
    tmp.push_back( lambda_transport.get_internal_cell_couplings() );
    tmp.push_back( ficks_water_nitrogen.get_internal_cell_couplings() );
    tmp.push_back( ficks_water_hydrogen.get_internal_cell_couplings() );
    tmp.push_back( ficks_oxygen_nitrogen.get_internal_cell_couplings() );
    reaction_source_terms.adjust_internal_cell_couplings(tmp);
    sorption_source_terms.adjust_internal_cell_couplings(tmp);
    
    //
    this->system_management.make_cell_couplings(tmp);

    // Now, initialize object that are used to setup initial solution and boundary conditions:    
    this->component_materialID_value_maps.push_back( ficks_oxygen_nitrogen.get_component_materialID_value()    );
    this->component_materialID_value_maps.push_back( ficks_water_hydrogen.get_component_materialID_value() );
    this->component_materialID_value_maps.push_back( ficks_water_nitrogen.get_component_materialID_value()   );
    this->component_materialID_value_maps.push_back( proton_transport.get_component_materialID_value()   );
    this->component_materialID_value_maps.push_back( electron_transport.get_component_materialID_value()   );
    this->component_materialID_value_maps.push_back( lambda_transport.get_component_materialID_value()   );
    this->component_materialID_value_maps.push_back( thermal_transport.get_component_materialID_value()   );
    OC.adjust_initial_solution(this->component_materialID_value_maps, this->mesh_generator);
        
    this->component_boundaryID_value_maps.push_back( ficks_oxygen_nitrogen.get_component_boundaryID_value()    );
    this->component_boundaryID_value_maps.push_back( ficks_water_hydrogen.get_component_boundaryID_value() );
    this->component_boundaryID_value_maps.push_back( ficks_water_nitrogen.get_component_boundaryID_value()   );
    this->component_boundaryID_value_maps.push_back( proton_transport.get_component_boundaryID_value()   );
    this->component_boundaryID_value_maps.push_back( electron_transport.get_component_boundaryID_value()   );
    this->component_boundaryID_value_maps.push_back( lambda_transport.get_component_boundaryID_value()   );
    this->component_boundaryID_value_maps.push_back( thermal_transport.get_component_boundaryID_value()   );
    OC.adjust_boundary_conditions(this->component_boundaryID_value_maps, this->mesh_generator);
    
    OC.print_operating_conditions();
    ACL->print_layer_properties();
    CCL->print_layer_properties();
    
    // Initialize matrices and sparticity pattern for the whole system
    this->remesh_matrices();
    
    // Initialize post-processing routines:
    ORRCurrent.initialize(param);
    HORCurrent.initialize(param);
    electronOhmicHeat.initialize(param);
    protonOhmicHeat.initialize(param);
    sorptionHeat.initialize(param);
    catReactionHeat.initialize(param);
    anReactionHeat.initialize(param);
    waterSorption.initialize(param);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfcNIThermal<dim>::initialize(ParameterHandler& param)
{
    OptimizationBlockMatrixApplication<dim>::initialize(param);
    _initialize(param);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfcNIThermal<dim>::initialize_solution(FuelCell::ApplicationCore::FEVector& initial_guess,
                                                  std::shared_ptr<Function<dim> > initial_function)
{
    std::shared_ptr< Function<dim> > initial_solution (new FuelCell::InitialSolution::AppPemfcNIThermalIC<dim> (&OC, this->mesh_generator, &this->system_management));
    
    DoFApplication<dim>::initialize_solution(initial_guess, initial_solution);
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::AppPemfcNIThermal<dim>::cell_matrix(FuelCell::ApplicationCore::MatrixVector& cell_matrices,
                                                        const typename DoFApplication<dim>::CellInfo& info)
{
    // -- Find out what material is the cell made of, i.e. MEA layer)
    const unsigned int material_id = info.dof_active_cell->material_id();

    //---- Equation Classes -----------------------------------------------------------
    if ( CGDL->belongs_to_material(material_id) )
    {
        ficks_oxygen_nitrogen.assemble_cell_matrix(cell_matrices, info, CGDL.get());
        ficks_water_nitrogen.assemble_cell_matrix(cell_matrices, info, CGDL.get());
        thermal_transport.assemble_cell_matrix(cell_matrices, info, CGDL.get());
        electron_transport.assemble_cell_matrix(cell_matrices, info, CGDL.get());
    }
    else if ( CMPL->belongs_to_material(material_id) )
    {
        ficks_oxygen_nitrogen.assemble_cell_matrix(cell_matrices, info, CMPL.get());
        ficks_water_nitrogen.assemble_cell_matrix(cell_matrices, info, CMPL.get());
        thermal_transport.assemble_cell_matrix(cell_matrices, info, CMPL.get());
        electron_transport.assemble_cell_matrix(cell_matrices, info, CMPL.get());
    }
    else if ( CCL->belongs_to_material(material_id) )
    {
        ficks_oxygen_nitrogen.assemble_cell_matrix(cell_matrices, info, CCL.get());
        ficks_water_nitrogen.assemble_cell_matrix(cell_matrices, info, CCL.get());
        thermal_transport.assemble_cell_matrix(cell_matrices, info, CCL.get());
        proton_transport.assemble_cell_matrix(cell_matrices, info, CCL.get());
        lambda_transport.assemble_cell_matrix(cell_matrices, info, CCL.get());
        electron_transport.assemble_cell_matrix(cell_matrices, info, CCL.get());
        reaction_source_terms.assemble_cell_matrix(cell_matrices, info, CCL.get());
        sorption_source_terms.assemble_cell_matrix(cell_matrices, info, CCL.get());
    }
    else if ( ML->belongs_to_material(material_id) )
    {
        thermal_transport.assemble_cell_matrix(cell_matrices, info, ML.get());
        proton_transport.assemble_cell_matrix(cell_matrices, info, ML.get());
        lambda_transport.assemble_cell_matrix(cell_matrices, info, ML.get());
    }
    else if ( ACL->belongs_to_material(material_id) )
    {
        ficks_water_hydrogen.assemble_cell_matrix(cell_matrices, info, ACL.get());
        thermal_transport.assemble_cell_matrix(cell_matrices, info, ACL.get());
        proton_transport.assemble_cell_matrix(cell_matrices, info, ACL.get());
        lambda_transport.assemble_cell_matrix(cell_matrices, info, ACL.get());
        electron_transport.assemble_cell_matrix(cell_matrices, info, ACL.get());
        reaction_source_terms.assemble_cell_matrix(cell_matrices, info, ACL.get());
        sorption_source_terms.assemble_cell_matrix(cell_matrices, info, ACL.get());
    }
    else if ( AMPL->belongs_to_material(material_id) )
    {
        ficks_water_hydrogen.assemble_cell_matrix(cell_matrices, info, AMPL.get());
        thermal_transport.assemble_cell_matrix(cell_matrices, info, AMPL.get());
        electron_transport.assemble_cell_matrix(cell_matrices, info, AMPL.get());
    }
    else if ( AGDL->belongs_to_material(material_id) )
    {
        ficks_water_hydrogen.assemble_cell_matrix(cell_matrices, info, AGDL.get());
        thermal_transport.assemble_cell_matrix(cell_matrices, info, AGDL.get());
        electron_transport.assemble_cell_matrix(cell_matrices, info, AGDL.get());
    }
    else
        AssertThrow(false, ExcNotImplemented());
    }
     
//---------------------------------------------------------------------------
template <int dim>
void 
NAME::AppPemfcNIThermal<dim>::cell_residual(FuelCell::ApplicationCore::FEVector& cell_vector,
                                                            const typename DoFApplication<dim>::CellInfo& info)
{ 
    // -- Assertion before starting routine:
    // Make sure vectors are the right size
    Assert (cell_vector.n_blocks() == this->element->n_blocks(),
            ExcDimensionMismatch (cell_vector.n_blocks(), this->element->n_blocks()));

    // -- Find out what material is the cell made of, i.e. MEA layer)
    const unsigned int material_id = info.dof_active_cell->material_id();

    //---- Equation Classes -----------------------------------------------------------
    if ( CGDL->belongs_to_material(material_id) )
    {
        ficks_oxygen_nitrogen.assemble_cell_residual(cell_vector, info, CGDL.get());
        ficks_water_nitrogen.assemble_cell_residual(cell_vector, info, CGDL.get());
        thermal_transport.assemble_cell_residual(cell_vector, info, CGDL.get());
        electron_transport.assemble_cell_residual(cell_vector, info, CGDL.get());
    }
    else if ( CMPL->belongs_to_material(material_id) )
    {
        ficks_oxygen_nitrogen.assemble_cell_residual(cell_vector, info, CMPL.get());
        ficks_water_nitrogen.assemble_cell_residual(cell_vector, info, CMPL.get());
        thermal_transport.assemble_cell_residual(cell_vector, info, CMPL.get());
        electron_transport.assemble_cell_residual(cell_vector, info, CMPL.get());
    }
    else if ( CCL->belongs_to_material(material_id) )
    {
        ficks_oxygen_nitrogen.assemble_cell_residual(cell_vector, info, CCL.get());
        ficks_water_nitrogen.assemble_cell_residual(cell_vector, info, CCL.get());
        thermal_transport.assemble_cell_residual(cell_vector, info, CCL.get());
        proton_transport.assemble_cell_residual(cell_vector, info, CCL.get());
        lambda_transport.assemble_cell_residual(cell_vector, info, CCL.get());
        electron_transport.assemble_cell_residual(cell_vector, info, CCL.get());
        reaction_source_terms.assemble_cell_residual(cell_vector, info, CCL.get());
        sorption_source_terms.assemble_cell_residual(cell_vector, info, CCL.get());
    }
    else if ( ML->belongs_to_material(material_id) )
    {
        thermal_transport.assemble_cell_residual(cell_vector, info, ML.get());
        proton_transport.assemble_cell_residual(cell_vector, info, ML.get());
        lambda_transport.assemble_cell_residual(cell_vector, info, ML.get());
    }
    else if ( ACL->belongs_to_material(material_id) )
    {
        ficks_water_hydrogen.assemble_cell_residual(cell_vector, info, ACL.get());
        thermal_transport.assemble_cell_residual(cell_vector, info, ACL.get());
        proton_transport.assemble_cell_residual(cell_vector, info, ACL.get());
        lambda_transport.assemble_cell_residual(cell_vector, info, ACL.get());
        electron_transport.assemble_cell_residual(cell_vector, info, ACL.get());
        reaction_source_terms.assemble_cell_residual(cell_vector, info, ACL.get());
        sorption_source_terms.assemble_cell_residual(cell_vector, info, ACL.get());
    }
    else if ( AMPL->belongs_to_material(material_id) )
    {
        ficks_water_hydrogen.assemble_cell_residual(cell_vector, info, AMPL.get());
        thermal_transport.assemble_cell_residual(cell_vector, info, AMPL.get());
        electron_transport.assemble_cell_residual(cell_vector, info, AMPL.get());
    }
    else if ( AGDL->belongs_to_material(material_id) )
    {
        ficks_water_hydrogen.assemble_cell_residual(cell_vector, info, AGDL.get());
        thermal_transport.assemble_cell_residual(cell_vector, info, AGDL.get());
        electron_transport.assemble_cell_residual(cell_vector, info, AGDL.get());
    }
    else
        AssertThrow(false, ExcNotImplemented());
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfcNIThermal<dim>::bdry_matrix(FuelCell::ApplicationCore::MatrixVector& bdry_matrices,
                                                        const typename DoFApplication<dim>::FaceInfo& bdry_info)
{
    const unsigned int material_id = bdry_info.dof_active_cell->material_id();
    
    if ( CGDL->belongs_to_material(material_id) )
    {
        electron_transport.assemble_bdry_matrix(bdry_matrices, bdry_info, CGDL.get());
        thermal_transport.assemble_bdry_matrix(bdry_matrices, bdry_info, CGDL.get());
    }
    else if ( CMPL->belongs_to_material(material_id) )
    {
        electron_transport.assemble_bdry_matrix(bdry_matrices, bdry_info, CMPL.get());
        thermal_transport.assemble_bdry_matrix(bdry_matrices, bdry_info, CMPL.get());
    }
    else if ( CCL->belongs_to_material(material_id) )
    {
        electron_transport.assemble_bdry_matrix(bdry_matrices, bdry_info, CCL.get());
        thermal_transport.assemble_bdry_matrix(bdry_matrices, bdry_info, CCL.get());
    }
    else if ( ML->belongs_to_material(material_id) )
    {
        thermal_transport.assemble_bdry_matrix(bdry_matrices, bdry_info, ML.get());
    }
    else if ( ACL->belongs_to_material(material_id) )
    {
        electron_transport.assemble_bdry_matrix(bdry_matrices, bdry_info, ACL.get());
        thermal_transport.assemble_bdry_matrix(bdry_matrices, bdry_info, ACL.get());
    }
    else if ( AMPL->belongs_to_material(material_id) )
    {
        electron_transport.assemble_bdry_matrix(bdry_matrices, bdry_info, AMPL.get());
        thermal_transport.assemble_bdry_matrix(bdry_matrices, bdry_info, AMPL.get());
    }
    else if ( AGDL->belongs_to_material(material_id) )
    {
        electron_transport.assemble_bdry_matrix(bdry_matrices, bdry_info, AGDL.get());
        thermal_transport.assemble_bdry_matrix(bdry_matrices, bdry_info, AGDL.get());
    }
    else
        AssertThrow(false, ExcNotImplemented());
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppPemfcNIThermal<dim>::bdry_residual(FuelCell::ApplicationCore::FEVector& bdry_vector,
                                                            const typename DoFApplication<dim>::FaceInfo& bdry_info)
{
    const unsigned int material_id = bdry_info.dof_active_cell->material_id();
    
    if ( CGDL->belongs_to_material(material_id) )
    {
        electron_transport.assemble_bdry_residual(bdry_vector, bdry_info, CGDL.get());
        thermal_transport.assemble_bdry_residual(bdry_vector, bdry_info, CGDL.get());
    }
    else if ( CMPL->belongs_to_material(material_id) )
    {
        electron_transport.assemble_bdry_residual(bdry_vector, bdry_info, CMPL.get());
        thermal_transport.assemble_bdry_residual(bdry_vector, bdry_info, CMPL.get());
    }
    else if ( CCL->belongs_to_material(material_id) )
    {
        electron_transport.assemble_bdry_residual(bdry_vector, bdry_info, CCL.get());
        thermal_transport.assemble_bdry_residual(bdry_vector, bdry_info, CCL.get());
    }
    else if ( ML->belongs_to_material(material_id) )
    {
        thermal_transport.assemble_bdry_residual(bdry_vector, bdry_info, ML.get());
    }
    else if ( ACL->belongs_to_material(material_id) )
    {
        electron_transport.assemble_bdry_residual(bdry_vector, bdry_info, ACL.get());
        thermal_transport.assemble_bdry_residual(bdry_vector, bdry_info, ACL.get());
    }
    else if ( AMPL->belongs_to_material(material_id) )
    {
        electron_transport.assemble_bdry_residual(bdry_vector, bdry_info, AMPL.get());
        thermal_transport.assemble_bdry_residual(bdry_vector, bdry_info, AMPL.get());
    }
    else if ( AGDL->belongs_to_material(material_id) )
    {
        electron_transport.assemble_bdry_residual(bdry_vector, bdry_info, AGDL.get());
        thermal_transport.assemble_bdry_residual(bdry_vector, bdry_info, AGDL.get());
    }
    else
        AssertThrow(false, ExcNotImplemented());
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::AppPemfcNIThermal<dim>::dirichlet_bc(std::map<unsigned int, double>& boundary_values) const
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
//-----------    OPTIMIZATION ROUTINES     ---------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::AppPemfcNIThermal<dim>::check_responses()
{
    // This has not yet been implemented. See app_cathode in the archive folder for
    // the code that used the old coding system.
    const std::type_info& info = typeid(*this);
    FcstUtilities::log << "Pure function " << __FUNCTION__
    << " called in Class "
    << info.name()  << std::endl;
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::AppPemfcNIThermal<dim>::cell_responses (std::vector<double>& resp,
                                              const typename DoFApplication<dim>::CellInfo& info,
                                              const FuelCell::ApplicationCore::FEVector& /*src*/)
{    
    // Find out what material is the cell made of, i.e. MEA layer)
    const unsigned int material_id = info.dof_active_cell->material_id();
    
    // Create a response Map
    std::map<FuelCellShop::PostProcessing::ResponsesNames, double> responseMap;
    
    // Finally add/organize the cell responses into the input response vector to determine overall response value.
    // All the computed total cell response values are normalized against the surface area of the layer.
    for (unsigned int r = 0; r < this->n_resp; ++r)
    {
        if ( (this->name_responses[r] == "cathode_current" || this->name_responses[r] == "current") and CCL->belongs_to_material(material_id) )
        {
            ORRCurrent.compute_responses(info, CCL.get(), responseMap);
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::ORR_current] / (l_channel/2.0 + l_land/2.0);
        }
        
        else if ( this->name_responses[r] == "anode_current" and ACL->belongs_to_material(material_id) )
        {            
            HORCurrent.compute_responses(info, ACL.get(), responseMap);
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::HOR_current] / (l_channel/2.0 + l_land/2.0);
        }
        
        else if ( this->name_responses[r] == "water_cathode" and CCL->belongs_to_material(material_id) )
        {
            waterSorption.compute_responses(info, CCL.get(), responseMap);
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::sorbed_water] / (l_channel/2.0 + l_land/2.0);
        }
        
        else if ( this->name_responses[r] == "water_anode" and ACL->belongs_to_material(material_id) )
        {
            waterSorption.compute_responses(info, ACL.get(), responseMap);
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::sorbed_water] / (l_channel/2.0 + l_land/2.0);
        }
        
        else if ( this->name_responses[r] == "cathode_reaction_heat" and CCL->belongs_to_material(material_id) )
        {
            catReactionHeat.compute_responses(info, CCL.get(), responseMap);
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::ORR_reaction_heat] / (l_channel/2.0 + l_land/2.0);
        }
            
        else if ( this->name_responses[r] == "cathode_irrev_heat" and CCL->belongs_to_material(material_id) )
        {
            catReactionHeat.compute_responses(info, CCL.get(), responseMap);
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::ORR_irrev_heat] / (l_channel/2.0 + l_land/2.0);
        }
        
        else if ( this->name_responses[r] == "cathode_rev_heat" and CCL->belongs_to_material(material_id) )
        {
            catReactionHeat.compute_responses(info, CCL.get(), responseMap);
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::ORR_rev_heat] / (l_channel/2.0 + l_land/2.0);
        }
        
        else if ( this->name_responses[r] == "cathode_watervap_heat" and CCL->belongs_to_material(material_id) )
        {
            catReactionHeat.compute_responses(info, CCL.get(), responseMap);
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::ORR_watervap_heat] / (l_channel/2.0 + l_land/2.0);
        }
        
        else if ( this->name_responses[r] == "anode_reaction_heat" and ACL->belongs_to_material(material_id) )
        {
            anReactionHeat.compute_responses(info, ACL.get(), responseMap);
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::HOR_reaction_heat] / (l_channel/2.0 + l_land/2.0);
        }
        
        else if ( this->name_responses[r] == "anode_irrev_heat" and ACL->belongs_to_material(material_id) )
        {
            anReactionHeat.compute_responses(info, ACL.get(), responseMap);
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::HOR_irrev_heat] / (l_channel/2.0 + l_land/2.0);
        }
        
        else if ( this->name_responses[r] == "anode_rev_heat" and ACL->belongs_to_material(material_id) )
        {
            anReactionHeat.compute_responses(info, ACL.get(), responseMap);
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::HOR_rev_heat] / (l_channel/2.0 + l_land/2.0);
        }
        
        else if ( this->name_responses[r] == "sorption_heat_cathode" and CCL->belongs_to_material(material_id) )
        {
            sorptionHeat.compute_responses(info, CCL.get(), responseMap);
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::sorption_heat] / (l_channel/2.0 + l_land/2.0);
        }
        
        else if ( this->name_responses[r] == "sorption_heat_anode" and ACL->belongs_to_material(material_id) )
        {
            sorptionHeat.compute_responses(info, ACL.get(), responseMap);
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::sorption_heat] / (l_channel/2.0 + l_land/2.0);
        }
        
        else if (this->name_responses[r] == "electron_ohmic_heat")
        {
            // Note that responseMap will be filled with electron_ohmic_heat functional only when we are in GDL/MPL/CL
            // Note that by picking a select or combination of layers below, one can compute partial contribution
            // from a layer or combination of layers to overall ohmic heat generation.
            // In this case, we are computing total ohmic heat generated in all the layers.
            
            if ( CGDL->belongs_to_material(material_id) )
                electronOhmicHeat.compute_responses(info, CGDL.get(), responseMap);
            else if ( CMPL->belongs_to_material(material_id) )
                electronOhmicHeat.compute_responses(info, CMPL.get(), responseMap);
            else if ( CCL->belongs_to_material(material_id) )
                electronOhmicHeat.compute_responses(info, CCL.get(), responseMap);
            else if ( ACL->belongs_to_material(material_id) )
                electronOhmicHeat.compute_responses(info, ACL.get(), responseMap);
            else if ( AMPL->belongs_to_material(material_id) )
                electronOhmicHeat.compute_responses(info, AMPL.get(), responseMap);
            else if ( AGDL->belongs_to_material(material_id) )
                electronOhmicHeat.compute_responses(info, AGDL.get(), responseMap);
            else // If membrane layer, then move on to next iteration.
                continue;
                
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::electron_ohmic_heat] / (l_channel/2.0 + l_land/2.0);
        }
        
        else if (this->name_responses[r] == "proton_ohmic_heat")
        {
            // Note that responseMap will be filled with proton_ohmic_heat functional only when we are in CL/Membrane
            if ( CCL->belongs_to_material(material_id) )
                protonOhmicHeat.compute_responses(info, CCL.get(), responseMap);
            else if ( ML->belongs_to_material(material_id) )
                protonOhmicHeat.compute_responses(info, ML.get(), responseMap);
            else if ( ACL->belongs_to_material(material_id) )
                protonOhmicHeat.compute_responses(info, ACL.get(), responseMap);
            else // For other layers, move on to next iteration
                continue;
                
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::proton_ohmic_heat] / (l_channel/2.0 + l_land/2.0);
        }
        else if (this->name_responses[r] == "PEM_proton_ohmic_heat")
        {
            // Note that responseMap will be filled with proton_ohmic_heat functional only when we are in CL/Membrane
            if ( ML->belongs_to_material(material_id) )
                protonOhmicHeat.compute_responses(info, ML.get(), responseMap);
            else // For other layers, move on to next iteration
                continue;
                
            resp[r] += responseMap[FuelCellShop::PostProcessing::ResponsesNames::proton_ohmic_heat] / (l_channel/2.0 + l_land/2.0);
        }
    }
}
     
//---------------------------------------------------------------------------
template <int dim>
void 
NAME::AppPemfcNIThermal<dim>::global_responses(std::vector<double>& resp,
                                               const FuelCell::ApplicationCore::FEVector& sol)
{
    const unsigned int t_index = this->system_management.solution_name_to_index("temperature_of_REV");
    const unsigned int xH2O_index = this->system_management.solution_name_to_index("water_molar_fraction");
    
    for (unsigned int r=0; r<this->n_resp; ++r)
    {
        if (this->name_responses[r] == "max_temperature")
        {
            double maxTemp(0.0);
            for (unsigned int i=0; i<sol.block(t_index).size(); ++i)
            {
                if (sol.block(t_index)(i) > maxTemp)
                    maxTemp = sol.block(t_index)(i);
            }
            resp[r] = maxTemp;
        }
        
        else if (this->name_responses[r] == "max_RH")
        {
            AssertThrow( OC.get_pc_Pa() == OC.get_pa_Pa(), ExcMessage("max_RH response current works only in the case of equal anode and cathode pressures.") );
            AssertThrow( sol.block(xH2O_index).size() == sol.block(t_index).size(), ExcInternalError() );
            
            double maxRH(0.0);
            for (unsigned int i=0; i<sol.block(t_index).size(); ++i)
            {
                double localRH = ( OC.get_pc_Pa()*sol.block(xH2O_index)(i) ) / ( water.get_water_vapor_saturation_pressure(sol.block(t_index)(i)) );
                if ( localRH > maxRH )
                    maxRH = localRH;
            }
            resp[r] = maxRH;
        }
    }
} 

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::AppPemfcNIThermal<dim>::global_dresponses_dl(std::vector<std::vector<double> >& df_dl,
                                                                const FuelCell::ApplicationCore::FEVector& )
{
    // This has not yet been implemented.
    const std::type_info& info = typeid(*this);
    FcstUtilities::log << "Pure function " << __FUNCTION__
    << " called in Class "
    << info.name()  << std::endl;
}

//---------------------------------------------------------------------------
template <int dim>
void 
NAME::AppPemfcNIThermal<dim>::global_dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector >& df_du,
                                                                const FuelCell::ApplicationCore::FEVector& /*src*/)
{
    // This has not yet been implemented. See app_cathode in the archive folder for
    // the code that used the old coding system.
    const std::type_info& info = typeid(*this);
    FcstUtilities::log << "Pure function " << __FUNCTION__
    << " called in Class "
    << info.name()  << std::endl;
}

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
NAME::AppPemfcNIThermal<dim>::evaluate (const FuelCell::ApplicationCore::FEVectors& vectors)
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
NAME::AppPemfcNIThermal<dim>::data_out(const std::string& filename,
                                       const FuelCell::ApplicationCore::FEVectors& src)
{
    // --- Find solution ---
    FuelCell::ApplicationCore::FEVector solution = src.vector( src.find_vector("Solution") );

    // --- Assign solution interpretations ---
    this->solution_interpretations.clear();
    this->solution_interpretations.resize(this->element->n_blocks(),
                                          DataComponentInterpretation::component_is_scalar);
    
    // DO FURTHER POST-PROCESSING
    
    // --- Create vector of PostProcessing objects ---
    std::vector< DataPostprocessor<dim>* > PostProcessing;
    
    // --- cathode_current ---
    FuelCellShop::PostProcessing::ORRCurrentDensityDataOut<dim> cathode_current(&this->system_management, CCL, &OC);
    PostProcessing.push_back(&cathode_current);

    // --- anode_current ---
    FuelCellShop::PostProcessing::HORCurrentDensityDataOut<dim> anode_current(&this->system_management, ACL, &OC);
    PostProcessing.push_back(&anode_current);
    
    // --- relative_humidity ---
    std::vector< boost::shared_ptr< FuelCellShop::Layer::PorousLayer<dim> > > porous_layers;
    porous_layers.push_back(CGDL);
    porous_layers.push_back(CMPL);
    porous_layers.push_back(CCL);
    porous_layers.push_back(ACL);
    porous_layers.push_back(AMPL);
    porous_layers.push_back(AGDL);
    FuelCellShop::PostProcessing::RelativeHumidityDataOut<dim> relative_humidity(&this->system_management, porous_layers, &OC);
    PostProcessing.push_back(&relative_humidity);

    // --- output ---
    DoFApplication<dim>::data_out(filename,
                                  solution,
                                  this->system_management.get_solution_names(),
                                  PostProcessing);
}
    

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
FuelCell::InitialSolution::AppPemfcNIThermalIC<dim>::AppPemfcNIThermalIC(FuelCell::OperatingConditions* OC, 
                                                                         boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid,
                                                                         FuelCell::SystemManagement* system_mgmt)
:
Function<dim> (system_mgmt->get_number_of_solution_names()),
OC(OC),
grid(grid),
system(system_mgmt)
{}

//---------------------------------------------------------------------------
template <int dim>
FuelCell::InitialSolution::AppPemfcNIThermalIC<dim>::~AppPemfcNIThermalIC()
{}

//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::InitialSolution::AppPemfcNIThermalIC<dim>::vector_value (const Point<dim> &p, 
                                                                   Vector<double> &v) const
{
    // size checking for vectors
    Assert(v.size() == this->n_components,
            ExcDimensionMismatch (v.size(), this->n_components));
    
    double RH_c = OC->get_RH_c();
    double RH_a = OC->get_RH_a();
    double V_cell = std::fabs(OC->get_V());
    double x_wv = OC->get_x_wv();
    double x_o2 = OC->get_x_o2();
    double x_h2 = OC->get_x_h2();
    double OCV = OC->get_OCV();
    double tempCathodeRib = OC->get_T();
    double tempAnodeRib = OC->get_T();
    
    double l_mem = grid->L_mem();
    double l_cat_a = grid->L_cat_a();
    double l_mpl_a = grid->L_mpl_a();
    double l_gdl_a = grid->L_gdl_a();
    
    double x = p(0);
    double x_norm = (x - l_gdl_a - l_cat_a - l_mpl_a)/l_mem;  
    double delta = -0.0001;
    
    if (x <= l_gdl_a + l_cat_a + l_mpl_a)   //Anode
    {
        v( system->solution_name_to_index("oxygen_molar_fraction") ) = 0.0;              // 0.5*x_o2; 
        v( system->solution_name_to_index("water_molar_fraction") ) = 1 - x_h2;          // x_wv_a - Dirichlet B.C. 
        v( system->solution_name_to_index("protonic_electrical_potential") ) = delta;    // phi_m
        v( system->solution_name_to_index("electronic_electrical_potential") ) = 0;      // phi_s - Dirichlet B.C.
        if(RH_a <= 0.5)
            v( system->solution_name_to_index("membrane_water_content") ) = 4.0;
        else
            v( system->solution_name_to_index("membrane_water_content") ) = 7.0;         // lambda
            v( system->solution_name_to_index("temperature_of_REV") ) = tempAnodeRib;
    }
    else if (x > l_gdl_a + l_cat_a + l_mpl_a && x < l_gdl_a + l_cat_a + l_mpl_a + l_mem) //membrane
    {
        v( system->solution_name_to_index("oxygen_molar_fraction") ) = 0.0;
        v( system->solution_name_to_index("water_molar_fraction") ) = x_wv;
        v( system->solution_name_to_index("protonic_electrical_potential") ) = ((V_cell-OCV)/2.0 - (delta))*x_norm + delta; 
        v( system->solution_name_to_index("electronic_electrical_potential") ) = V_cell * x_norm;
        if(RH_a <= 0.5 && RH_c <= 0.5)
            v( system->solution_name_to_index("membrane_water_content") ) = 4.0;
        else if(RH_a <= 0.5 && RH_c > 0.5)
            v( system->solution_name_to_index("membrane_water_content") ) = 3.0*x_norm + 4.0;
        else if(RH_a > 0.5 && RH_c <= 0.5)
            v( system->solution_name_to_index("membrane_water_content") ) = -3.0*x_norm + 7.0;
        else
            v( system->solution_name_to_index("membrane_water_content") ) = 7.0;
        
        v( system->solution_name_to_index("temperature_of_REV") ) = (tempAnodeRib + tempCathodeRib)/2.0;
    }
    else if (x >= l_gdl_a + l_cat_a + l_mpl_a + l_mem)  //Cathode
    {
        v( system->solution_name_to_index("oxygen_molar_fraction") ) = x_o2;
        v( system->solution_name_to_index("water_molar_fraction") ) = x_wv;
        v( system->solution_name_to_index("protonic_electrical_potential") ) = (V_cell-OCV)/2.0;
        v( system->solution_name_to_index("electronic_electrical_potential") ) = V_cell;
        // lambda
        if (RH_c <= 0.5)
            v( system->solution_name_to_index("membrane_water_content") ) = 4.0;
        else
            v( system->solution_name_to_index("membrane_water_content") ) = 7.0;
        
        v( system->solution_name_to_index("temperature_of_REV") ) = tempCathodeRib;
    }
}
                                                                                    
                                                                                    
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations. 
template class NAME::AppPemfcNIThermal<deal_II_dimension>;
template class FuelCell::InitialSolution::AppPemfcNIThermalIC<deal_II_dimension>;
                                                                                        
