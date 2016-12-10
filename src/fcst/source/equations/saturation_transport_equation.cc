//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: saturation_transport_equation.cc
//    - Description: Implementation of saturation_transport_equation.h. This
//      class is used to setup the matrix and rhs for liquid water saturation 
//      transport due to capillary diffusion.
//    - Developers: Madhur Bhaiya
//
//---------------------------------------------------------------------------

#include "equations/saturation_transport_equation.h"

namespace NAME = FuelCellShop::Equation;

// ---             ---
// --- Constructors ---
// ---             ---
template<int dim>
NAME::SaturationTransportEquation<dim>::SaturationTransportEquation(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data)
:
EquationBase<dim>(system_management,data)
{
    this->name_base_variable = "liquid_water_saturation";
    this->equation_name = "Liquid Water Saturation Transport Equation";
    FcstUtilities::log << "->FuelCellShop::Equation::SaturationTransportEquation" << std::endl;
    
    //----Initializing VariableInfo Structs------------------------------------------
    //----Setting indices_exist to false --------------------------------------------
    s_liquid_water.indices_exist = false;
    x_water.indices_exist = false;
    t_rev.indices_exist = false;
    
    //-- Counter to store if generic_data and constant_data need to be initialized.
    this->counter.resize(2, true);
}
// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::SaturationTransportEquation<dim>::~SaturationTransportEquation()
{}

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::SaturationTransportEquation<dim>::declare_parameters(ParameterHandler& param) const
{
    NAME::EquationBase<dim>::declare_parameters(param);
    
    param.enter_subsection("Equations");
    {
        param.enter_subsection("Liquid Water Saturation Transport Equation");
        {             
            
            param.declare_entry("Evaporation rate constant",
                                "1.0e-10",
                                Patterns::Double(),
                                "Evaporation rate constant, [mol/(Pa-cm^2-s)].");
            
            param.declare_entry("Condensation rate constant",
                                "1.0e-7",
                                Patterns::Double(),
                                "Condensation rate constant, [mol/(Pa-cm^2-s)].");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

// ---            ---
// --- initialize ---
// ---            ---

template<int dim>
void
NAME::SaturationTransportEquation<dim>::initialize(ParameterHandler& param)
{
    
    NAME::EquationBase<dim>::initialize(param);
    
    param.enter_subsection("Equations");
    {
        param.enter_subsection("Liquid Water Saturation Transport Equation");
        {  
            k_e = param.get_double("Evaporation rate constant");
            k_c = param.get_double("Condensation rate constant");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    this->make_internal_cell_couplings();
}

// ---                      ---
// --- assemble_cell_matrix ---
// ---                      ---

template<int dim>
void
NAME::SaturationTransportEquation<dim>::assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                                             const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                             FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
    FuelCellShop::Material::WaterVapor water;
	
    if ( this->counter[0] )
    {
        this->make_assemblers_generic_constant_data();
        this->counter[0] = false;
    }
    
    if ( this->counter[1] )
    {
        this->make_assemblers_cell_constant_data(cell_info);
        this->counter[1] = false;
    }
    
    cell_residual_counter = false;
    this->make_assemblers_cell_variable_data(cell_info, layer);
    
    double k_phase(0.);    
    
    //-------- Looping over Quadrature points ----------------------------
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
        // Calculating certain values:
        double p_vap = p_cell * (cell_info.values[last_iter_cell][x_water.solution_index][q]); // water vapour pressure [Pa]
        double psat = water.get_water_vapor_saturation_pressure(cell_info.values[last_iter_cell][t_rev.solution_index][q]);  // saturation pressure [Pa]
        double dpsat_dT = water.get_Dwater_vapor_saturation_pressure_Dtemperature(cell_info.values[last_iter_cell][t_rev.solution_index][q]); // derivative of saturation pressure w.r.t temperature
        double h_lv = FuelCellShop::Material::LiquidWater::latentVap_heat(cell_info.values[last_iter_cell][t_rev.solution_index][q]); // latent heat of vaporization [J/mol]
        double dhlv_dT = FuelCellShop::Material::LiquidWater::deriv_latentVap_heat(cell_info.values[last_iter_cell][t_rev.solution_index][q]); // derivative of latent heat w.r.t temperature
        
        if (psat > p_vap)
            k_phase = k_e;
        else if (psat < p_vap)
            k_phase = k_c;
        else
            k_phase = 0.0;
        
        //---------------Liquid Water Saturation Transport Equation -----------------------------------
        for (unsigned int i=0; i < (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell; ++i)
        {            
            //-----------Assembling Matrix for terms corresponding to "s" BLOCK------------------------
            for (unsigned int j=0; j < (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell; ++j)
            {
                cell_matrices[s_liquid_water.block_index].matrix(i,j) += ( this->JxW_cell[q] * (grad_phi_s_cell[q][j] * rhok_mu_dpcds_cell[q] * grad_phi_s_cell[q][i]) 
                                                                                                +
                                                                           this->JxW_cell[q] * phi_s_cell[q][j] * (grad_phi_s_cell[q][i] * 
                                                                           ds_rhok_mu_dpcds_cell[q] * cell_info.gradients[last_iter_cell][s_liquid_water.solution_index][q]) 
                                                                                                + 
                                                                           this->JxW_cell[q] * (-1.0) * phi_s_cell[q][j] * phi_s_cell[q][i] * k_phase * M_water * darea_lv_ds_cell[q] * (p_vap-psat) );
            }

            //-----Assembling Matrix for terms corresponding to "T" BLOCK------------------------------
            for (unsigned int j=0; j < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++j)
            {
                cell_matrices[t_rev.block_index].matrix(i,j) += ( this->JxW_cell[q] * (-1.0) * phi_s_cell[q][j] * (grad_phi_T_cell[q][i] * 
                                                                  dT_rhok_mu_dpcds_cell[q] * cell_info.gradients[last_iter_cell][s_liquid_water.solution_index][q])
                                                                                                + 
                                                                  this->JxW_cell[q] * phi_s_cell[q][j] * phi_T_cell[q][i] * k_phase * M_water * area_lv_cell[q] * dpsat_dT );
            }
            
            //-----Assembling Matrix for terms corresponding to "x_H2O" BLOCK---------------------------
            for (unsigned int j=0; j < (cell_info.fe(x_water.fetype_index)).dofs_per_cell; ++j)
            {
                cell_matrices[x_water.block_index].matrix(i,j) += ( this->JxW_cell[q] * (-1.0) * phi_s_cell[q][j] * phi_xwater_cell[q][i] * k_phase * M_water * area_lv_cell[q] * p_cell );
            }
            
        }
        
        
        //---------------Thermal Transport Equation ---------------------------------------------------
        for (unsigned int i=0; i < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++i)
        {
            //-----------Assembling Matrix for terms corresponding to "s" BLOCK------------------------
            for (unsigned int j=0; j < (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell; ++j)
            {
                cell_matrices[thermal_blockindex_sliquid].matrix(i,j) += ( this->JxW_cell[q] * (-1.0) * phi_T_cell[q][j] * phi_s_cell[q][i] * k_phase * darea_lv_ds_cell[q] * h_lv * (p_vap-psat) );
            }
            
            //-----Assembling Matrix for terms corresponding to "T" BLOCK------------------------------
            for (unsigned int j=0; j < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++j)
            {
                cell_matrices[thermal_blockindex_trev].matrix(i,j) += ( this->JxW_cell[q] * phi_T_cell[q][j] * phi_T_cell[q][i] * k_phase * area_lv_cell[q] * 
                                                                        ( (h_lv*dpsat_dT) - ((p_vap-psat)*dhlv_dT) ) );
            }
            
            //-----Assembling Matrix for terms corresponding to "x_H2O" BLOCK---------------------------
            for (unsigned int j=0; j < (cell_info.fe(x_water.fetype_index)).dofs_per_cell; ++j)
            {
                cell_matrices[thermal_blockindex_xwater].matrix(i,j) += ( this->JxW_cell[q] * phi_T_cell[q][j] * phi_xwater_cell[q][i] * k_phase * area_lv_cell[q] * p_cell * h_lv );
            }
        }
        
        
        // --- Ficks Transport Equation - water -------------------------------------------------------
        for (unsigned int i=0; i < (cell_info.fe(x_water.fetype_index)).dofs_per_cell; ++i)
        {
            //-----------Assembling Matrix for terms corresponding to "s" BLOCK------------------------
            for (unsigned int j=0; j < (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell; ++j)
            {
                cell_matrices[fickswater_blockindex_sliquid].matrix(i,j) += ( this->JxW_cell[q] * (-1.0) * phi_xwater_cell[q][j] * phi_s_cell[q][i] * k_phase * darea_lv_ds_cell[q] * (psat-p_vap) );
            }
            
            //-----Assembling Matrix for terms corresponding to "T" BLOCK------------------------------
            for (unsigned int j=0; j < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++j)
            {
                cell_matrices[fickswater_blockindex_trev].matrix(i,j) += ( this->JxW_cell[q] * (-1.0) * phi_xwater_cell[q][j] * phi_T_cell[q][i] * k_phase * area_lv_cell[q] * dpsat_dT );
            }
            
            //-----Assembling Matrix for terms corresponding to "x_H2O" BLOCK---------------------------
            for (unsigned int j=0; j < (cell_info.fe(x_water.fetype_index)).dofs_per_cell; ++j)
            {
                cell_matrices[fickswater_blockindex_xwater].matrix(i,j) += ( this->JxW_cell[q] * phi_xwater_cell[q][j] * phi_xwater_cell[q][i] * k_phase * area_lv_cell[q] * p_cell );
            }
        }
    }
}

// ---                        ---
// --- assemble_cell_residual ---
// ---                        ---
template<int dim>
void
NAME::SaturationTransportEquation<dim>::assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
                                                               const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                               FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{   
    FuelCellShop::Material::WaterVapor water;

    if ( this->counter[0] )
    {
        this->make_assemblers_generic_constant_data();
        this->counter[0] = false;
    }
    
    if ( this->counter[1] )
    {
        this->make_assemblers_cell_constant_data(cell_info);
        this->counter[1] = false;
    }
    
    cell_residual_counter = true;
    this->make_assemblers_cell_variable_data(cell_info, layer);
    
    double k_phase(0.);
    
    for (unsigned int q=0; q < this->n_q_points_cell; ++q)
    {
        // Calculating certain values:
        double p_vap = p_cell * (cell_info.values[last_iter_cell][x_water.solution_index][q]); // water vapour pressure [Pa]
        double psat = water.get_water_vapor_saturation_pressure(cell_info.values[last_iter_cell][t_rev.solution_index][q]);  // saturation pressure [Pa]
        double h_lv = FuelCellShop::Material::LiquidWater::latentVap_heat(cell_info.values[last_iter_cell][t_rev.solution_index][q]); // latent heat of vaporization [J/mol]
        
        if (psat > p_vap)
            k_phase = k_e;
        else if (psat < p_vap)
            k_phase = k_c;
        else
            k_phase = 0.0;
        
        //---------------Liquid Water Saturation Transport Equation -----------------------------------
        for (unsigned int i=0; i < (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell; ++i)
        {
            cell_residual.block(s_liquid_water.solution_index)(i) += ( this->JxW_cell[q] * (grad_phi_s_cell[q][i] * rhok_mu_dpcds_cell[q] * 
                                                                       cell_info.gradients[last_iter_cell][s_liquid_water.solution_index][q])
                                                                                     +
                                                                       this->JxW_cell[q] * (-1.0) * phi_s_cell[q][i] * k_phase * M_water * area_lv_cell[q] * (p_vap-psat) );
        }
        
        // --- Ficks Transport Equation - water -------------------------------------------------------
        for (unsigned int i=0; i < (cell_info.fe(x_water.fetype_index)).dofs_per_cell; ++i)
        {
            cell_residual.block(x_water.solution_index)(i) += ( this->JxW_cell[q] * (-1.0) * phi_xwater_cell[q][i] * k_phase * area_lv_cell[q] * (psat-p_vap) );
        }
        
        // --- Thermal Transport Equation -------------------------------------------------------------
        for (unsigned int i=0; i < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++i)
        {
            cell_residual.block(t_rev.solution_index)(i) += ( this->JxW_cell[q] * (-1.0) * phi_T_cell[q][i] * k_phase * area_lv_cell[q] * h_lv * (p_vap-psat) );
        }
    }
}

// ---                              ---
// --- make_internal_cell_couplings ---
// ---                              ---
template<int dim>
void
NAME::SaturationTransportEquation<dim>::make_internal_cell_couplings()
{
    AssertThrow(this->system_management->solution_in_userlist(this->name_base_variable), VariableShouldExistForEquation(this->name_base_variable, this->equation_name) );
    AssertThrow(this->system_management->equation_name_to_index(this->equation_name) == this->system_management->solution_name_to_index(this->name_base_variable), 
              IndexDoNotMatch(this->name_base_variable, this->equation_name) );
    AssertThrow(this->system_management->solution_in_userlist("temperature_of_REV"), ExcMessage("temperature_of_REV should also be solved for, when considering liquid water saturation transport."));
    AssertThrow(this->system_management->solution_in_userlist("water_molar_fraction"), ExcMessage("water_molar_fraction should also be solved for, when considering liquid water saturation transport."));
    
    std::map< std::string, DoFTools::Coupling > tmp;
    
    std::vector< std::string> sol_names = this->system_management->get_solution_names(); 
    
    for (unsigned int i = 0; i < sol_names.size(); ++i)
    {
        if (sol_names[i] == this->name_base_variable)
            tmp[this->name_base_variable] = DoFTools::always;
        
        else if (sol_names[i] == "temperature_of_REV")
            tmp["temperature_of_REV"] = DoFTools::always;
        
        else if (sol_names[i] == "water_molar_fraction")
            tmp["water_molar_fraction"] = DoFTools::always;
        
        else
            tmp[sol_names[i]] = DoFTools::none;
    }
    
    this->internal_cell_couplings[this->equation_name] = tmp;
}

// ---                     ---
// --- print_equation_info ---
// ---                     ---

template<int dim>
void
NAME::SaturationTransportEquation<dim>::print_equation_info() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "-------------------------------------------------------------------------------" << std::endl;
    FcstUtilities::log << std::endl;
    
    FcstUtilities::log << "PARAMETERS FOR \" Liquid Water Saturation Transport Equation \":" << std::endl;
    FcstUtilities::log << std::endl;

    FcstUtilities::log << "Evaporation rate constant:  " << k_e << std::endl;
    FcstUtilities::log << "Condensation rate constant:  " << k_c << std::endl;
    FcstUtilities::log << std::endl;
        
    FcstUtilities::log << "INTERNAL CELL COUPLINGS FOR \" Liquid Water Saturation Transport Equation \":" << std::endl;
    FcstUtilities::log << std::endl;
        
    couplings_map::const_iterator iter;
    
    for( iter = this->internal_cell_couplings.begin(); iter != this->internal_cell_couplings.end(); ++iter )
    {
        FcstUtilities::log << "\"" << iter->first << "\"" << ":";
        
        std::map<std::string, DoFTools::Coupling> int_map = iter->second;
        std::map<std::string, DoFTools::Coupling>::const_iterator int_iter;
        
        for( int_iter = int_map.begin(); int_iter != int_map.end(); ++int_iter )
            FcstUtilities::log << "\"" << int_iter->first << "\"" << " is coupled as " << int_iter->second << std::endl;
        
        FcstUtilities::log << std::endl;
    }
    
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "-------------------------------------------------------------------------------" << std::endl;
    FcstUtilities::log << std::endl;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
// LOCAL CG FEM BASED ASSEMBLERS - make_ FUNCTIONS //
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
       
// ---                                       ---
// --- make_assemblers_generic_constant_data ---
// ---                                       ---

template<int dim>
void
NAME::SaturationTransportEquation<dim>::make_assemblers_generic_constant_data()
{
    FuelCellShop::Material::WaterVapor water;    

    M_water = water.get_molar_mass() * 1.0e3;        // gm/mol
    rho_l = FuelCellShop::Material::LiquidWater::get_density();    // gm/cm^3
    
    //-----------Filling VariableInfo structures------------------------------------------
    //----------liquid_water_saturation---------------------------------------------------
    s_liquid_water.solution_index = this->system_management->solution_name_to_index(this->name_base_variable);
    s_liquid_water.block_index = this->system_management->matrix_block_index(this->equation_name, this->name_base_variable);
    s_liquid_water.fetype_index = this->system_management->block_info->base_element[s_liquid_water.solution_index];
    s_liquid_water.indices_exist = true;
        
    //-----------water_molar_fraction-----------------------------------------------------
    x_water.solution_index = this->system_management->solution_name_to_index("water_molar_fraction");
    x_water.block_index = this->system_management->matrix_block_index(this->equation_name, "water_molar_fraction");
    x_water.fetype_index = this->system_management->block_info->base_element[x_water.solution_index];
    x_water.indices_exist = true;
        
    //----------temperature_of_solid_phase----------------------------------------------------
    t_rev.solution_index = this->system_management->solution_name_to_index("temperature_of_REV");
    t_rev.block_index = this->system_management->matrix_block_index(this->equation_name, "temperature_of_REV");
    t_rev.fetype_index = this->system_management->block_info->base_element[t_rev.solution_index];
    t_rev.indices_exist = true;
    
    //---Filling block indices corresponding to "Thermal Transport Equation" and "Ficks Transport Equation - water"-----------------------------
    thermal_blockindex_xwater = this->system_management->matrix_block_index("Thermal Transport Equation", "water_molar_fraction");
    thermal_blockindex_sliquid = this->system_management->matrix_block_index("Thermal Transport Equation", "liquid_water_saturation");
    thermal_blockindex_trev = this->system_management->matrix_block_index("Thermal Transport Equation", "temperature_of_REV");
    
    fickswater_blockindex_xwater = this->system_management->matrix_block_index("Ficks Transport Equation - water", "water_molar_fraction");
    fickswater_blockindex_sliquid = this->system_management->matrix_block_index("Ficks Transport Equation - water", "liquid_water_saturation");
    fickswater_blockindex_trev = this->system_management->matrix_block_index("Ficks Transport Equation - water", "temperature_of_REV");
}

// ---                                    ---
// --- make_assemblers_cell_constant_data ---
// ---                                    ---

template<int dim>
void
NAME::SaturationTransportEquation<dim>::make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
    Assert( ((s_liquid_water.indices_exist) && (x_water.indices_exist) && (t_rev.indices_exist)), 
            ExcMessage("make_assemblers_generic_constant_data function not called before.") );
    
    this->n_q_points_cell = (cell_info.fe(s_liquid_water.fetype_index)).n_quadrature_points;
    last_iter_cell = cell_info.global_data->find_vector(this->solution_vector_name);
    
    //-------------Allocation------------------------------------------
    phi_s_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell ) );
    phi_xwater_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(x_water.fetype_index)).dofs_per_cell ) );
    phi_T_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(t_rev.fetype_index)).dofs_per_cell ) );
    grad_phi_s_cell.resize( this->n_q_points_cell, std::vector< Tensor<1,dim> >( (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell ) );
    grad_phi_T_cell.resize( this->n_q_points_cell, std::vector< Tensor<1,dim> >( (cell_info.fe(t_rev.fetype_index)).dofs_per_cell ) );
    
    //-----------------------------------------------------------------
    this->JxW_cell.resize(this->n_q_points_cell);

    rhok_mu_dpcds_cell.resize( this->n_q_points_cell, Tensor<2,dim>() );
    ds_rhok_mu_dpcds_cell.resize( this->n_q_points_cell, Tensor<2,dim>() );
    dT_rhok_mu_dpcds_cell.resize( this->n_q_points_cell, Tensor<2,dim>() );
    area_lv_cell.resize( this->n_q_points_cell, 0.0);
}

// ---                                    ---
// --- make_assemblers_cell_variable_data ---
// ---                                    ---

template<int dim>
void
NAME::SaturationTransportEquation<dim>::make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                         FuelCellShop::Layer::BaseLayer<dim>* const layer)
{
    // ----- type infos -------------
    const std::type_info& GasDiffusionLayer = typeid(FuelCellShop::Layer::GasDiffusionLayer<dim>);
    const std::type_info& MicroPorousLayer = typeid(FuelCellShop::Layer::MicroPorousLayer<dim>);
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    
    const std::type_info& base_layer = layer->get_base_type();
    
    std::vector< Tensor<2,dim> > k_l(this->n_q_points_cell, Tensor<2,dim>() );
    std::vector<double> dpc_ds(this->n_q_points_cell, 0.0);
    
    std::map< VariableNames, std::vector< Tensor<2,dim> > > deriv_kl;
    std::map< VariableNames, std::vector<double> > deriv_dpc_ds;
    
    // ----- dynamic cast and filling the containers -----------------
    try
    {
        if (base_layer == GasDiffusionLayer)
        {
            FuelCellShop::Layer::GasDiffusionLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::GasDiffusionLayer<dim>* >(layer);
            //-------------------
            ptr->set_saturation( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][s_liquid_water.solution_index],liquid_water_saturation) );
            ptr->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index],temperature_of_REV) );
            //-------------------
            ptr->get_p(p_cell);
            p_cell *= Units::convert(1.,Units::ATM_to_PA); // converting into pascals
            ptr->liquid_permeablity(k_l);                  // cm^2
            ptr->dpcapillary_dsat(dpc_ds);                 // dyne/cm^2
            ptr->interfacial_surface_area(area_lv_cell);   // cm^2/cm^3
            
            // --- Derivative effective property for cell_matrix --------------------
            if ( !cell_residual_counter )
            {
                std::vector<VariableNames> deriv_flags;
                deriv_flags.push_back(liquid_water_saturation);    
                deriv_flags.push_back(temperature_of_REV);
                ptr->set_derivative_flags(deriv_flags);
                //-------------------
                ptr->derivative_liquid_permeablity(deriv_kl);
                ptr->derivative_dpcapillary_dsat(deriv_dpc_ds);
                std::map< VariableNames, std::vector<double> > deriv_alv;
                ptr->derivative_interfacial_surface_area(deriv_alv);
                darea_lv_ds_cell = deriv_alv[liquid_water_saturation];
            }
        }
        
        else if (base_layer == MicroPorousLayer)
        {
            FuelCellShop::Layer::MicroPorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::MicroPorousLayer<dim>* >(layer);
            //-------------------
            ptr->set_saturation( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][s_liquid_water.solution_index],liquid_water_saturation) );
            ptr->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index],temperature_of_REV) );
            //-------------------
            ptr->get_p(p_cell);
            p_cell *= Units::convert(1.,Units::ATM_to_PA); // converting into pascals
            ptr->liquid_permeablity(k_l);                  // cm^2
            ptr->dpcapillary_dsat(dpc_ds);                 // dyne/cm^2
            ptr->interfacial_surface_area(area_lv_cell);   // cm^2/cm^3
            
            // --- Derivative effective property for cell_matrix --------------------
            if ( !cell_residual_counter )
            {
                std::vector<VariableNames> deriv_flags;
                deriv_flags.push_back(liquid_water_saturation);    
                deriv_flags.push_back(temperature_of_REV);
                ptr->set_derivative_flags(deriv_flags);
                //-------------------
                ptr->derivative_liquid_permeablity(deriv_kl);
                ptr->derivative_dpcapillary_dsat(deriv_dpc_ds);
                std::map< VariableNames, std::vector<double> > deriv_alv;
                ptr->derivative_interfacial_surface_area(deriv_alv);
                darea_lv_ds_cell = deriv_alv[liquid_water_saturation];
            }
        }
        
        else if (base_layer == CatalystLayer)
        {
            FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);
            //-------------------
            ptr->set_saturation( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][s_liquid_water.solution_index],liquid_water_saturation) );
            ptr->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index],temperature_of_REV) );
            //-------------------
            ptr->get_p(p_cell);
            p_cell *= Units::convert(1.,Units::ATM_to_PA); // converting into pascals
            ptr->liquid_permeablity(k_l);                  // cm^2
            ptr->dpcapillary_dsat(dpc_ds);                 // dyne/cm^2
            ptr->interfacial_surface_area(area_lv_cell);   // cm^2/cm^3
            
            // --- Derivative effective property for cell_matrix --------------------
            if ( !cell_residual_counter )
            {
                std::vector<VariableNames> deriv_flags;
                deriv_flags.push_back(liquid_water_saturation);    
                deriv_flags.push_back(temperature_of_REV);
                ptr->set_derivative_flags(deriv_flags);
                //-------------------
                ptr->derivative_liquid_permeablity(deriv_kl);
                ptr->derivative_dpcapillary_dsat(deriv_dpc_ds);
                std::map< VariableNames, std::vector<double> > deriv_alv;
                ptr->derivative_interfacial_surface_area(deriv_alv);
                darea_lv_ds_cell = deriv_alv[liquid_water_saturation];
            }
        }
        
        else
            AssertThrow( false, ExcNotImplemented() );
    }
    catch(const std::bad_cast& e)
    {
        const std::type_info& info = typeid(*layer);  
        FcstUtilities::log << "Object of type "<<info.name()<<" not implemented"<< std::endl;
        FcstUtilities::log << e.what() << std::endl;
    }
    
    //---------------------------------------------------------------------------------------------------------------
    //------------Looping over quadrature points in the cell --------------------------------------------------------
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
        double mu_l = FuelCellShop::Material::LiquidWater::viscosity( cell_info.values[last_iter_cell][t_rev.solution_index][q] );  // gm/(cm-s)
        rhok_mu_dpcds_cell[q] = ((rho_l * dpc_ds[q])/ mu_l) * k_l[q];
        
        if (!cell_residual_counter)
        {          
            ds_rhok_mu_dpcds_cell[q] = (rho_l/mu_l)*( (k_l[q]*deriv_dpc_ds[liquid_water_saturation][q])+(dpc_ds[q]*deriv_kl[liquid_water_saturation][q]) );
            
            double dmu_dT = FuelCellShop::Material::LiquidWater::deriv_viscosity( cell_info.values[last_iter_cell][t_rev.solution_index][q] ); // gm/(cm-s-K)
            
            dT_rhok_mu_dpcds_cell[q] = ((rho_l*((mu_l*deriv_dpc_ds[temperature_of_REV][q])-(dpc_ds[q]*dmu_dT)))/(mu_l*mu_l)) * k_l[q];
        }
        
        //---------------------------------------------------------------------------------------------------------------
        //-------JxW----------
        this->JxW_cell[q] = (cell_info.fe(s_liquid_water.fetype_index)).JxW(q);
        
        //------ Filling shape functions etc ----------------------------------------------------------------------
        //------ This avoids recalculating shape functions etc for efficiency -------------------------------------
        for (unsigned int k=0; k < (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell; ++k)
        {
            phi_s_cell[q][k] = (cell_info.fe(s_liquid_water.fetype_index)).shape_value(k,q);
            grad_phi_s_cell[q][k] = (cell_info.fe(s_liquid_water.fetype_index)).shape_grad(k,q);
        }
        
        for (unsigned int k=0; k < (cell_info.fe(x_water.fetype_index)).dofs_per_cell; ++k)
            phi_xwater_cell[q][k] = (cell_info.fe(x_water.fetype_index)).shape_value(k,q);
        
        for (unsigned int k=0; k < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++k)
        {
            phi_T_cell[q][k] = (cell_info.fe(t_rev.fetype_index)).shape_value(k,q);
            
            if (!cell_residual_counter)
                grad_phi_T_cell[q][k] = (cell_info.fe(t_rev.fetype_index)).shape_grad(k,q);
        }
    }
}

// ---                         ---
// --- explicit instantiations ---
// ---                         ---

template class NAME::SaturationTransportEquation<deal_II_dimension>;
