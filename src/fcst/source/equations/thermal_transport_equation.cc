//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: thermal_transport_equation.cc 18-04-2013
//    - Description: Equation class for Thermal transport
//    - Developers: Madhur Bhaiya
//
//---------------------------------------------------------------------------

#include "equations/thermal_transport_equation.h"

namespace NAME = FuelCellShop::Equation;

// ---             ---
// --- Constructor ---
// ---             ---

template <int dim>
NAME::ThermalTransportEquation<dim>::ThermalTransportEquation(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data)
:
NAME::EquationBase<dim>(system_management,data)
{
    this->name_base_variable = "temperature_of_REV";
    this->equation_name = "Thermal Transport Equation";
    FcstUtilities::log << "->FuelCellShop::Equation::ThermalTransportEquation" << std::endl;
    
    //----Initializing VariableInfo Structs------------------------------------------
    //----Setting indices_exist to false --------------------------------------------
    t_rev.indices_exist = false;
    phi_s.indices_exist = false;
    phi_m.indices_exist = false;
    lambda.indices_exist = false;
    s_liquid_water.indices_exist = false;
    p_liquid_water.indices_exist = false;
    
    this->counter.resize(3, true);
}

// ---            ---
// --- Destructor ---
// ---            ---

template <int dim>
NAME::ThermalTransportEquation<dim>::~ThermalTransportEquation()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template <int dim>
void
NAME::ThermalTransportEquation<dim>::declare_parameters(ParameterHandler& param) const
{
    
    NAME::EquationBase<dim>::declare_parameters(param);
    
    param.enter_subsection("Equations");
    {
        param.enter_subsection(this->equation_name);
        {
            param.enter_subsection("Boolean flags");
            {
                param.declare_entry("Electronic ohmic heat in GDL",
                                    "true",
                                    Patterns::Bool(),
                                    "Flag to include electronic ohmic heating source term in Gas diffusion layer(s).");
                param.declare_entry("Electronic ohmic heat in MPL",
                                    "true",
                                    Patterns::Bool(),
                                    "Flag to include electronic ohmic heating source term in Microporous layer(s).");
                param.declare_entry("Electronic ohmic heat in CL",
                                    "true",
                                    Patterns::Bool(),
                                    "Flag to include electronic ohmic heating source term in Catalyst layer(s).");
                param.declare_entry("Protonic ohmic heat in CL",
                                    "true",
                                    Patterns::Bool(),
                                    "Flag to include protonic ohmic heating source term in Catalyst layer(s).");
                param.declare_entry("Protonic ohmic heat in ML",
                                    "true",
                                    Patterns::Bool(),
                                    "Flag to include protonic ohmic heating source term in Membrane layer(s).");
                
                param.declare_entry("Enthalpy transport due to fickian diffusion of gases",
                                    "false",
                                    Patterns::Bool(),
                                    "Flag to include term(s) corresponding to enthalpy transport due to fickian diffusion of gases in the thermal transport equation.");
                param.declare_entry("Enthalpy transport associated with lambda transport",
                                    "false",
                                    Patterns::Bool(),
                                    "Flag to include term(s) enthalpy transport due to lambda (sorbed water) transport inside the polymer electrolyte due to electroosmosis, diffusion and thermoosmosis.");
            }
            param.leave_subsection();
            
            param.enter_subsection("Boundary conditions");
            {
                param.declare_entry("Constant Heat Flux Boundary Conditions",
                                    "",
                                    Patterns::Map( Patterns::Integer(0), Patterns::Double() ),
                                    "A comma-separated list of Constant heat flux boundaries, with prescribed values of heat fluxes [W/cm^2]."
                                    "\n"
                                    "Correct format is of a map, given as \"id1: value1, id2: value2, id3: value3, ...\" where "
                                    "each 'id' corresponds to boundary_id, and each 'value' refers to constant heat flux [W/cm^2] (positive for heat leaving out; negative for coming in).");
                
                param.declare_entry("Convective Heat Flux Boundary Conditions",
                                    "",
                                    Patterns::Map( Patterns::Integer(0), Patterns::Anything()),
                                    "A comma-separated list of Convective heat flux boundary conditions."
                                    "\n"
                                    "Correct format is of a map, given as \"id1: htc1;amb1, id2: htc2;amb2, id3: htc3;amb3, ...\" where "
                                    "each 'id' corresponds to boundary_id, and 'htc' refers to heat transfer coeffcient [W/(cm^2-K)] and  'amb' refers to ambient temperature [K].");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

// ---            ---
// --- initialize ---
// ---            ---

template <int dim>
void
NAME::ThermalTransportEquation<dim>::initialize(ParameterHandler& param)
{
    
    NAME::EquationBase<dim>::initialize(param);        
        
    param.enter_subsection("Equations");
    {
        param.enter_subsection(this->equation_name);
        {
            param.enter_subsection("Boolean flags");
            {
                electron_ohmic_heat_gdl = param.get_bool("Electronic ohmic heat in GDL");
                electron_ohmic_heat_mpl = param.get_bool("Electronic ohmic heat in MPL");
                electron_ohmic_heat_cl = param.get_bool("Electronic ohmic heat in CL");
                proton_ohmic_heat_cl = param.get_bool("Protonic ohmic heat in CL");
                proton_ohmic_heat_ml = param.get_bool("Protonic ohmic heat in ML");
                
                enthalpy_fickian_transport = param.get_bool("Enthalpy transport due to fickian diffusion of gases");
                enthalpy_lambda_transport = param.get_bool("Enthalpy transport associated with lambda transport");
            }
            param.leave_subsection();
            
            param.enter_subsection("Boundary conditions");
            {
                const_heat_flux_map = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Constant Heat Flux Boundary Conditions")) );
                conv_heat_flux_map = FcstUtilities::split_mapvalue_list<unsigned int, double>( FcstUtilities::string_to_map<unsigned int, std::string>( param.get("Convective Heat Flux Boundary Conditions") ) , ';');
            }
            param.leave_subsection();
        }
        param.leave_subsection();
        
        if (enthalpy_lambda_transport)
        {
            AssertThrow(this->system_management->solution_in_userlist("membrane_water_content"), ExcMessage("membrane_water_content should be one of the solution variables, if enthalpy transport associated with lambda transport if ON."));
            
            // Reading whether Thermal Osmosis is ON or OFF in the Lambda Transport Equation class
            param.enter_subsection("Membrane Water Content Transport Equation");
            {
                param.enter_subsection("Boolean flags for lambda transport modes");
                {
                    flag_thermoosmosis = param.get_bool("Thermo-osmosis");
                }
                param.leave_subsection();
            }
            param.leave_subsection();
        }
    }
    param.leave_subsection();

    // Checking that convective heat flux boundary conditions are passed correctly in the parameter file.
    for (std::map< unsigned int, std::vector<double> >::const_iterator iter = conv_heat_flux_map.begin(); iter != conv_heat_flux_map.end(); ++iter)
    {
        AssertThrow(iter->second.size() == 2, ExcMessage("For convective heat flux boundary conditions, exactly two parameters viz. heat transfer coeff and ambient temperature should be specified."));
        AssertThrow(((iter->second[0]>=0.)&&(iter->second[1]>=0.)), ExcMessage("In convective heat flux boundary conditions, heat transfer coeff/ambient temperature should be positive."));
    }
    
    this->make_internal_cell_couplings();
}

// ---                      ---
// --- assemble_cell_matrix ---
// ---                      ---

template <int dim>
void
NAME::ThermalTransportEquation<dim>::assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                                           const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                           FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
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
    
    
    //-------- Looping over Quadrature points ----------------------------
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
        //---------------LOOP over i -----------------------------------------------------------------
        for (unsigned int i=0; i < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++i)
        {
            //--------------LOOP(s) over j-------------------------------------------------------------
            
            //-----------Assembling Matrix for terms corresponding to "T" BLOCK------------------------
            //----------- If some term is not be included for a particular layer, it will automatically calculate out to zero--------------------
            for (unsigned int j=0; j < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++j)
            {
                cell_matrices[t_rev.block_index].matrix(i,j) += ( this->JxW_cell[q] * (grad_phi_T_cell[q][i] * keff_cell[q] * grad_phi_T_cell[q][j]) 
                                                                                        + 
                                                                  this->JxW_cell[q] * phi_T_cell[q][j] * 
                                                                  (grad_phi_T_cell[q][i] * dkeff_dT_cell[q] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q]) 
                                                                                        + 
                                                                  this->JxW_cell[q] * (-1.0) * dsigmaMeff_dT_cell[q] * phi_T_cell[q][i] * phi_T_cell[q][j] * 
                                                                  (grad_phiM_cell_old[q] * grad_phiM_cell_old[q]) 
                                                                                        + 
                                                                  this->JxW_cell[q] * (-1.0) * delectroosmotic_dhdT_dT[q] * phi_T_cell[q][i] * phi_T_cell[q][j] * 
                                                                  (grad_phiM_cell_old[q] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q]) 
                                                                                        + 
                                                                  this->JxW_cell[q] * (-1.0) * dbackdiff_dhdT_dT[q] * phi_T_cell[q][i] * phi_T_cell[q][j] * 
                                                                  (grad_lambda_cell_old[q] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q]) 
                                                                                        + 
                                                                  this->JxW_cell[q] * (-1.0) * dthermoosmotic_dhdT_dT[q] * phi_T_cell[q][i] * phi_T_cell[q][j] * 
                                                                  (cell_info.gradients[last_iter_cell][t_rev.solution_index][q] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q]) 
                                                                                        + 
                                                                  this->JxW_cell[q] * (-1.0) * electroosmotic_dhdT[q] * phi_T_cell[q][i] * 
                                                                  (grad_phi_T_cell[q][j] * grad_phiM_cell_old[q]) 
                                                                                        + 
                                                                  this->JxW_cell[q] * (-1.0) * backdiff_dhdT[q] * phi_T_cell[q][i] * 
                                                                  (grad_phi_T_cell[q][j] * grad_lambda_cell_old[q]) 
                                                                                        + 
                                                                  this->JxW_cell[q] * (-2.0) * thermoosmotic_dhdT[q] * phi_T_cell[q][i] * 
                                                                  (grad_phi_T_cell[q][j] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q]) );
                
                if (enthalpy_fickian_transport)
                {
                    for (typename std::map< std::string, std::vector< Tensor<2,dim> > >::const_iterator iter = conc_Deff_dHdT_map.begin(); iter!=conc_Deff_dHdT_map.end(); ++iter)
                    {
                        Assert( xi_map.find(iter->first)!=xi_map.end() && (xi_map.at(iter->first)).indices_exist, ExcInternalError() );
                        Assert( dT_concDeffdHdT_map.find(iter->first) != dT_concDeffdHdT_map.end(), ExcInternalError() );
                        
                        cell_matrices[t_rev.block_index].matrix(i,j) += ( this->JxW_cell[q] * (-1.0) * phi_T_cell[q][i] * phi_T_cell[q][j] * 
                                                                          (cell_info.gradients[last_iter_cell][(xi_map.at(iter->first)).solution_index][q] * 
                                                                          dT_concDeffdHdT_map.at(iter->first)[q] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q])
                                                                                                +
                                                                          this->JxW_cell[q] * (-1.0) * phi_T_cell[q][i] * 
                                                                          (cell_info.gradients[last_iter_cell][(xi_map.at(iter->first)).solution_index][q] * 
                                                                          conc_Deff_dHdT_map.at(iter->first)[q] * grad_phi_T_cell[q][j]) );
                    }
                }
            }
            
            //----------ELECTRONIC OHMIC HEAT ---------------------------------------------------------
            //----------TERM CORRESPONDING TO "PHI_S" BLOCK--------------------------------------------
            if (phi_s.indices_exist)
            {
                //-------Assembling Matrix-------------------------------------------------------------
                for (unsigned int j=0; j < (cell_info.fe(phi_s.fetype_index)).dofs_per_cell; ++j)
                {
                    cell_matrices[phi_s.block_index].matrix(i,j) += ( this->JxW_cell[q] * (-2.0) * phi_T_cell[q][i] * 
                                                                      (grad_phi_phiS_cell[q][j] * sigmaSeff_cell * grad_phiS_cell_old[q]) );
                }
            }
            
            //----------PROTONIC OHMIC HEAT -----------------------------------------------------------
            //----------TERM CORRESPONDING TO "PHI_M" BLOCK--------------------------------------------
            if (phi_m.indices_exist)
            {
                //-------Assembling Matrix-------------------------------------------------------------
                for (unsigned int j=0; j < (cell_info.fe(phi_m.fetype_index)).dofs_per_cell; ++j)
                {
                    cell_matrices[phi_m.block_index].matrix(i,j) += ( this->JxW_cell[q] * (-2.0) * phi_T_cell[q][i] * sigmaMeff_cell[q] * 
                                                                      (grad_phi_phiM_cell[q][j] * grad_phiM_cell_old[q]) 
                                                                                            + 
                                                                      this->JxW_cell[q] * (-1.0) * electroosmotic_dhdT[q] * phi_T_cell[q][i] * 
                                                                      (grad_phi_phiM_cell[q][j] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q]) );
                }
            }
            
            //----------TERM CORRESPONDING TO "LAMBDA" BLOCK-------------------------------------------
            if (lambda.indices_exist)
            {
                //-----Assembling Matrix---------------------------------------------------------------
                for (unsigned int j=0; j < (cell_info.fe(lambda.fetype_index)).dofs_per_cell; ++j)
                {
                    cell_matrices[lambda.block_index].matrix(i,j) += ( this->JxW_cell[q] * (-1.0) * dsigmaMeff_dlambda_cell[q] * phi_T_cell[q][i] * phi_lambda_cell[q][j] * 
                                                                       (grad_phiM_cell_old[q] * grad_phiM_cell_old[q]) 
                                                                                            + 
                                                                       this->JxW_cell[q] * (-1.0) * delectroosmotic_dhdT_dlambda[q] * phi_T_cell[q][i] * phi_lambda_cell[q][j] * 
                                                                       (grad_phiM_cell_old[q] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q]) 
                                                                                            + 
                                                                       this->JxW_cell[q] * (-1.0) * dbackdiff_dhdT_dlambda[q] * phi_T_cell[q][i] * phi_lambda_cell[q][j] * 
                                                                       (grad_lambda_cell_old[q] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q]) 
                                                                                            + 
                                                                       this->JxW_cell[q] * (-1.0) * backdiff_dhdT[q] * phi_T_cell[q][i] * 
                                                                       (grad_phi_lambda_cell[q][j] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q]) );
                }
            }
            
            //---------Enthalpy transport via fickian diffusion----------------------------------------
            if (enthalpy_fickian_transport)
            {
                //-------TERM CORRESPONDING TO "XI" BLOCK ---------------------------------------------
                for (typename std::map< std::string, std::vector< Tensor<2,dim> > >::const_iterator iter = conc_Deff_dHdT_map.begin(); iter!=conc_Deff_dHdT_map.end(); ++iter)
                {
                    Assert( xi_map.find(iter->first)!=xi_map.end() && (xi_map.at(iter->first)).indices_exist, ExcInternalError() );
                    Assert( grad_phi_xi_map.find(iter->first)!=grad_phi_xi_map.end(), ExcInternalError() );
                    
                    //-----Assembling Matrix-----------------------------------------------------------
                    for (unsigned int j=0; j < (cell_info.fe((xi_map.at(iter->first)).fetype_index)).dofs_per_cell; ++j)
                    {
                        cell_matrices[(xi_map.at(iter->first)).block_index].matrix(i,j) += ( this->JxW_cell[q] * (-1.0) * phi_T_cell[q][i] * 
                                                                                             (cell_info.gradients[last_iter_cell][t_rev.solution_index][q] * conc_Deff_dHdT_map.at(iter->first)[q] * 
                                                                                             grad_phi_xi_map.at(iter->first)[q][j]) );
                    }
                    
                    //--------TERM CORRESPONDING TO "s" BLOCK-----------------------------------------------
                    if (s_liquid_water.indices_exist)
                    {
                        //-----Assembling Matrix---------------------------------------------------------------
                        for (unsigned int j=0; j < (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell; ++j)
                        {
                            cell_matrices[s_liquid_water.block_index].matrix(i,j) += ( this->JxW_cell[q] * (-1.0) * phi_T_cell[q][i] * phi_s_cell[q][j] * 
                                                                                       (cell_info.gradients[last_iter_cell][(xi_map.at(iter->first)).solution_index][q] * 
                                                                                       ds_concDeffdHdT_map.at(iter->first)[q] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q]) );
                        }
                    }
                    
                    
                    //--------TERM CORRESPONDING TO "p" BLOCK-----------------------------------------------
                    if (p_liquid_water.indices_exist)
                    {
                        //-----Assembling Matrix---------------------------------------------------------------
                        for (unsigned int j=0; j < (cell_info.fe(p_liquid_water.fetype_index)).dofs_per_cell; ++j)
                        {
                            cell_matrices[p_liquid_water.block_index].matrix(i,j) += 0.0 *( this->JxW_cell[q] * (-1.0) * phi_T_cell[q][i] * phi_p_cell[q][j] * 
                            (cell_info.gradients[last_iter_cell][(xi_map.at(iter->first)).solution_index][q] * 
                            dp_concDeffdHdT_map.at(iter->first)[q] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q])  );
                        }
                    }
                }
                
                
                
            }
            
        } // End Loop over "i"
    } // End Loop Over Quadrature Points
}

// ---                        ---
// --- assemble_cell_residual ---
// ---                        ---

template <int dim>
void
NAME::ThermalTransportEquation<dim>::assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
                                                             const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                             FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
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
    
    
    for (unsigned int q=0; q < this->n_q_points_cell; ++q)
    {
        for (unsigned int i=0; i < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++i)
        {
            cell_residual.block(t_rev.solution_index)(i) += ( this->JxW_cell[q] * 
                                                              (grad_phi_T_cell[q][i] * keff_cell[q] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q])
                                                                                + 
                                                              this->JxW_cell[q] * (-1.0) * phi_T_cell[q][i] * 
                                                              (grad_phiS_cell_old[q] * sigmaSeff_cell * grad_phiS_cell_old[q]) 
                                                                                + 
                                                              this->JxW_cell[q] * (-1.0) * phi_T_cell[q][i] * sigmaMeff_cell[q] * 
                                                              (grad_phiM_cell_old[q] * grad_phiM_cell_old[q]) 
                                                                                + 
                                                              this->JxW_cell[q] * (-1.0) * phi_T_cell[q][i] * electroosmotic_dhdT[q] * 
                                                              (grad_phiM_cell_old[q] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q]) 
                                                                                + 
                                                              this->JxW_cell[q] * (-1.0) * phi_T_cell[q][i] * backdiff_dhdT[q] * 
                                                              (grad_lambda_cell_old[q] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q]) 
                                                                                + 
                                                              this->JxW_cell[q] * (-1.0) * phi_T_cell[q][i] * thermoosmotic_dhdT[q] * 
                                                              (cell_info.gradients[last_iter_cell][t_rev.solution_index][q] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q]) );
            
            if (enthalpy_fickian_transport)
            {
                for (typename std::map< std::string, std::vector< Tensor<2,dim> > >::const_iterator iter = conc_Deff_dHdT_map.begin(); iter!=conc_Deff_dHdT_map.end(); ++iter)
                {
                    Assert( xi_map.find(iter->first)!=xi_map.end() && (xi_map.at(iter->first)).indices_exist, ExcInternalError() );
                    cell_residual.block(t_rev.solution_index)(i) +=   ( this->JxW_cell[q] * (-1.0) * phi_T_cell[q][i] * 
                                                                      (cell_info.gradients[last_iter_cell][(xi_map.at(iter->first)).solution_index][q] * 
                                                                      conc_Deff_dHdT_map.at(iter->first)[q] * cell_info.gradients[last_iter_cell][t_rev.solution_index][q]) );
                }
            }
        }
    }
}

// ---                      ---
// --- assemble_bdry_matrix ---
// ---                      ---

template <int dim>
void
NAME::ThermalTransportEquation<dim>::assemble_bdry_matrix(FuelCell::ApplicationCore::MatrixVector&                                 bdry_matrices,
                                                           const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                           FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
    if ( this->counter[0] )
    {
        this->make_assemblers_generic_constant_data();
        this->counter[0] = false;
    }
    
    if ( this->counter[2] )
    {
        this->make_assemblers_bdry_constant_data(bdry_info);
        this->counter[2] = false;
    }
    
    //------ Convective Heat Flux boundaries --------------------------------------------------------------
    std::map<unsigned int, std::vector<double>>::const_iterator iter_conv_heat = conv_heat_flux_map.find( bdry_info.dof_face->boundary_indicator() );
    
    if ( iter_conv_heat != conv_heat_flux_map.end() )
    {        
        bdry_residual_counter = false;
        this->make_assemblers_bdry_variable_data(bdry_info, layer);
        
        //-------- Looping over Quadrature points ----------------------------
        for (unsigned int q = 0; q < this->n_q_points_bdry; ++q)
        {
            //---------------LOOP over i -----------------------------------------------------------------
            for (unsigned int i = 0; i < (bdry_info.fe(t_rev.fetype_index)).dofs_per_cell; ++i)
            {
                //--------------LOOP over j-------------------------------------------------------------
                for (unsigned int j = 0; j < (bdry_info.fe(t_rev.fetype_index)).dofs_per_cell; ++j)
                {
                    bdry_matrices[t_rev.block_index].matrix(i,j) += ( this->JxW_bdry[q] * phi_T_bdry[q][i] * phi_T_bdry[q][j] * (iter_conv_heat->second)[0] );
                }
            }
        } 
    }
}

// ---                        ---
// --- assemble_bdry_residual ---
// ---                        ---

template <int dim>
void
NAME::ThermalTransportEquation<dim>::assemble_bdry_residual(FuelCell::ApplicationCore::FEVector&                                     bdry_residual,
                                                             const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                             FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
    if ( this->counter[0] )
    {
        this->make_assemblers_generic_constant_data();
        this->counter[0] = false;
    }
    
    if ( this->counter[2] )
    {
        this->make_assemblers_bdry_constant_data(bdry_info);
        this->counter[2] = false;
    }
    
    //------ Constant Heat Flux boundaries --------------------------------------------------------------
    std::map<unsigned int, double>::const_iterator iter_const_heat = const_heat_flux_map.find( bdry_info.dof_face->boundary_indicator() );
    //------ Convective Heat Flux boundaries --------------------------------------------------------------
    std::map<unsigned int, std::vector<double>>::const_iterator iter_conv_heat = conv_heat_flux_map.find( bdry_info.dof_face->boundary_indicator() );
    
    if ( (iter_const_heat != const_heat_flux_map.end()) || (iter_conv_heat != conv_heat_flux_map.end()) )
    {
        AssertThrow( !((iter_const_heat!=const_heat_flux_map.end())&&(iter_conv_heat!=conv_heat_flux_map.end())), 
                     ExcMessage("Boundary_id(s) should not match for constant heat flux and convective heat flux boundary conditions.") );
        
        bdry_residual_counter = true;
        this->make_assemblers_bdry_variable_data(bdry_info, layer);
        
        //-------- Looping over Quadrature points ----------------------------
        for (unsigned int q = 0; q < this->n_q_points_bdry; ++q)
        {
            for (unsigned int i = 0; i < (bdry_info.fe(t_rev.fetype_index)).dofs_per_cell; ++i)
            {
                if ( iter_const_heat != const_heat_flux_map.end() )
                {
                    bdry_residual.block(t_rev.solution_index)(i) += ( this->JxW_bdry[q] * phi_T_bdry[q][i] * iter_const_heat->second );
                }
                else if ( iter_conv_heat != conv_heat_flux_map.end() )
                {
                    bdry_residual.block(t_rev.solution_index)(i) += ( this->JxW_bdry[q] * phi_T_bdry[q][i] * (iter_conv_heat->second)[0] * 
                                                                        (bdry_info.values[last_iter_bdry][t_rev.solution_index][q] - (iter_conv_heat->second)[1]) );
                }
            }
        }
    }
}

// ---                              ---
// --- make_internal_cell_couplings ---
// ---                              ---

template <int dim>
void
NAME::ThermalTransportEquation<dim>::make_internal_cell_couplings()
{
    AssertThrow(this->system_management->solution_in_userlist(this->name_base_variable), VariableShouldExistForEquation(this->name_base_variable, this->equation_name) );
    AssertThrow(this->system_management->equation_name_to_index(this->equation_name) == this->system_management->solution_name_to_index(this->name_base_variable), 
                IndexDoNotMatch(this->name_base_variable, this->equation_name) );
    if (electron_ohmic_heat_cl || electron_ohmic_heat_gdl || electron_ohmic_heat_mpl)
        AssertThrow(this->system_management->solution_in_userlist("electronic_electrical_potential"), VariableNotFoundForOhmicHeat("Electronic", "electronic_electrical_potential") );
    
    if (proton_ohmic_heat_cl || proton_ohmic_heat_ml)
        AssertThrow(this->system_management->solution_in_userlist("protonic_electrical_potential"), VariableNotFoundForOhmicHeat("Protonic", "protonic_electrical_potential") );

    std::map< std::string, DoFTools::Coupling > tmp;
    std::vector< std::string> sol_names = this->system_management->get_solution_names();
    
    for (unsigned int i = 0; i < sol_names.size(); ++i)
    {
        if (sol_names[i] == this->name_base_variable)
            tmp[this->name_base_variable] = DoFTools::always;
        
        else if (sol_names[i] == "electronic_electrical_potential")
        {
            if (electron_ohmic_heat_cl || electron_ohmic_heat_gdl || electron_ohmic_heat_mpl)
                tmp["electronic_electrical_potential"] = DoFTools::always;
            else
                tmp["electronic_electrical_potential"] = DoFTools::none;
        }
        
        else if (sol_names[i] == "protonic_electrical_potential")
        {
            if (proton_ohmic_heat_cl || proton_ohmic_heat_ml || enthalpy_lambda_transport)
                tmp["protonic_electrical_potential"] = DoFTools::always;
            else
                tmp["protonic_electrical_potential"] = DoFTools::none;
        }
        
        else if (sol_names[i] == "membrane_water_content")
        {
            if (proton_ohmic_heat_cl || proton_ohmic_heat_ml || enthalpy_lambda_transport)
                tmp["membrane_water_content"] = DoFTools::always;
            else
                tmp["membrane_water_content"] = DoFTools::none;
        }
        
        else if (sol_names[i] == "liquid_water_saturation")
            tmp["liquid_water_saturation"] = DoFTools::always;
        
        else if (sol_names[i] == "capillary_pressure")
            tmp["capillary_pressure"] = DoFTools::always;
        
        else if ( (sol_names[i].find("_molar_fraction") != std::string::npos) && (enthalpy_fickian_transport) )
            tmp[sol_names[i]] = DoFTools::always;
        
        else
            tmp[sol_names[i]] = DoFTools::none;
    }
    
    this->internal_cell_couplings[this->equation_name] = tmp;
}

// ---                     ---
// --- print_equation_info ---
// ---                     ---

template <int dim>
void
NAME::ThermalTransportEquation<dim>::print_equation_info() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "-------------------------------------------------------------------------------" << std::endl;
    FcstUtilities::log << std::endl;
    
    FcstUtilities::log << "PARAMETERS FOR \" Thermal Transport Equation \":" << std::endl;
    FcstUtilities::log << std::endl;
    
    FcstUtilities::log << "Boolean flags for heat source terms:" << std::endl;
    FcstUtilities::log << std::endl;
    
    FcstUtilities::log << "Electronic ohmic heat in GDL: " << electron_ohmic_heat_gdl << std::endl;
    FcstUtilities::log << "Electronic ohmic heat in MPL: " << electron_ohmic_heat_mpl << std::endl;
    FcstUtilities::log << "Electronic ohmic heat in CL:  " << electron_ohmic_heat_cl << std::endl;
    FcstUtilities::log << "Protonic ohmic heat in CL:    " << proton_ohmic_heat_cl << std::endl;
    FcstUtilities::log << "Protonic ohmic heat in ML:    " << proton_ohmic_heat_ml << std::endl;
    
    FcstUtilities::log << "Enthalpy transport via fickian gas diffusion: " << enthalpy_fickian_transport << std::endl;
    FcstUtilities::log << "Enthalpy transport associated with lambda transport: " << enthalpy_lambda_transport << std::endl;
    
    FcstUtilities::log << std::endl;
    
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "INTERNAL CELL COUPLINGS FOR \" Thermal Transport Equation \":" << std::endl;
    FcstUtilities::log << std::endl;
    
    couplings_map::const_iterator iter;
    
    for( iter = this->internal_cell_couplings.begin(); iter != this->internal_cell_couplings.end(); ++iter )
    {
        FcstUtilities::log << "\"" << iter->first << "\"" << ":" << std::endl;
        
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

template <int dim>
void
NAME::ThermalTransportEquation<dim>::make_assemblers_generic_constant_data()
{
    //-----------Filling VariableInfo structures------------------------------------------
    //-----------temperature_of_solid_phase-----------------------------------------------
    t_rev.solution_index = this->system_management->solution_name_to_index(this->name_base_variable);
    t_rev.block_index = this->system_management->matrix_block_index(this->equation_name, this->name_base_variable);
    t_rev.fetype_index = this->system_management->block_info->base_element[t_rev.solution_index];
    t_rev.indices_exist = true;
    
    //-----------electronic_electrical_potential------------------------------------------
    if ( electron_ohmic_heat_gdl || electron_ohmic_heat_mpl || electron_ohmic_heat_cl )
    {
        phi_s.solution_index = this->system_management->solution_name_to_index("electronic_electrical_potential");
        phi_s.block_index = this->system_management->matrix_block_index(this->equation_name, "electronic_electrical_potential");
        phi_s.fetype_index = this->system_management->block_info->base_element[phi_s.solution_index];
        phi_s.indices_exist = true;
    }
    
    //----------protonic_electrical_potential---------------------------------------------
    if ( proton_ohmic_heat_cl || proton_ohmic_heat_ml || (this->system_management->solution_in_userlist("protonic_electrical_potential") && enthalpy_lambda_transport) )
    {
        phi_m.solution_index = this->system_management->solution_name_to_index("protonic_electrical_potential");
        phi_m.block_index = this->system_management->matrix_block_index(this->equation_name, "protonic_electrical_potential");
        phi_m.fetype_index = this->system_management->block_info->base_element[phi_m.solution_index];
        phi_m.indices_exist = true;
    }
        
    //----------membrane_water_content-----------------------------------------------
    if ( this->system_management->solution_in_userlist("membrane_water_content") && (proton_ohmic_heat_cl || proton_ohmic_heat_ml || enthalpy_lambda_transport) )
    {
        lambda.solution_index = this->system_management->solution_name_to_index("membrane_water_content");
        lambda.block_index = this->system_management->matrix_block_index(this->equation_name, "membrane_water_content");
        lambda.fetype_index = this->system_management->block_info->base_element[lambda.solution_index];
        lambda.indices_exist = true;
    }
    
    //----------liquid_water_saturation-----------------------------------------------
    if ( this->system_management->solution_in_userlist("liquid_water_saturation") )
    {
        s_liquid_water.solution_index = this->system_management->solution_name_to_index("liquid_water_saturation");
        s_liquid_water.block_index = this->system_management->matrix_block_index(this->equation_name, "liquid_water_saturation");
        s_liquid_water.fetype_index = this->system_management->block_info->base_element[s_liquid_water.solution_index];
        s_liquid_water.indices_exist = true;
    }
        //----------capillary_pressure-----------------------------------------------
    if ( this->system_management->solution_in_userlist("capillary_pressure") )
    {
        p_liquid_water.solution_index = this->system_management->solution_name_to_index("capillary_pressure");
        p_liquid_water.block_index = this->system_management->matrix_block_index(this->equation_name, "capillary_pressure");
        p_liquid_water.fetype_index = this->system_management->block_info->base_element[p_liquid_water.solution_index];
        p_liquid_water.indices_exist = true;
    }
    
    //--- Filling VariableInfo structures corresponding to gases whose transport via Fickian diffusion, hence also enthalpy transport is being considered
    if (enthalpy_fickian_transport)
    {
        std::vector<std::string> sol_names = this->system_management->get_solution_names();
        for (unsigned int i=0; i<sol_names.size(); ++i)
        {
            if (sol_names[i].find("_molar_fraction") != std::string::npos)
            {
                VariableInfo temp;
                temp.solution_index = this->system_management->solution_name_to_index(sol_names[i]);
                temp.block_index = this->system_management->matrix_block_index(this->equation_name, sol_names[i]);
                temp.fetype_index = this->system_management->block_info->base_element[temp.solution_index];
                temp.indices_exist = true;
                
                xi_map[ sol_names[i] ] = temp;
            }
        }
    }
}

// ---                                    ---
// --- make_assemblers_cell_constant_data ---
// ---                                    ---

template <int dim>
void
NAME::ThermalTransportEquation<dim>::make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
    Assert( t_rev.indices_exist, ExcMessage("make_assemblers_generic_constant_data function not called before.") );

    this->n_q_points_cell = (cell_info.fe(t_rev.fetype_index)).n_quadrature_points;
    last_iter_cell = cell_info.global_data->find_vector(this->solution_vector_name);

    //-------------Allocation------------------------------------------
    phi_T_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(t_rev.fetype_index)).dofs_per_cell ) );
    grad_phi_T_cell.resize( this->n_q_points_cell, std::vector< Tensor<1,dim> >( (cell_info.fe(t_rev.fetype_index)).dofs_per_cell ) );

    if (phi_s.indices_exist)
        grad_phi_phiS_cell.resize( this->n_q_points_cell, std::vector< Tensor<1,dim> >( (cell_info.fe(phi_s.fetype_index)).dofs_per_cell ) );

    if (phi_m.indices_exist)
        grad_phi_phiM_cell.resize( this->n_q_points_cell, std::vector< Tensor<1,dim> >( (cell_info.fe(phi_m.fetype_index)).dofs_per_cell ) );

    if (lambda.indices_exist)
    {
        phi_lambda_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(lambda.fetype_index)).dofs_per_cell ) );
        grad_phi_lambda_cell.resize( this->n_q_points_cell, std::vector< Tensor<1,dim> >( (cell_info.fe(lambda.fetype_index)).dofs_per_cell ) );
    }
    
    if (s_liquid_water.indices_exist)
        phi_s_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell ) );
    
    if (p_liquid_water.indices_exist)
        phi_p_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(p_liquid_water.fetype_index)).dofs_per_cell ) );

    //-----------------------------------------------------------------
    this->JxW_cell.resize(this->n_q_points_cell);

    keff_cell.resize(this->n_q_points_cell);    // Only keff_cell can be allocated here as rest needs to be cleared everytime (since they depend on flags)
}

// ---                                    ---
// --- make_assemblers_bdry_constant_data ---
// ---                                    ---

template <int dim>
void
NAME::ThermalTransportEquation<dim>::make_assemblers_bdry_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info)
{
    Assert( t_rev.indices_exist, ExcMessage("make_assemblers_generic_constant_data function not called before.") );
    
    this->n_q_points_bdry = (bdry_info.fe(t_rev.fetype_index)).n_quadrature_points;
    last_iter_bdry = bdry_info.global_data->find_vector(this->solution_vector_name);
    
    //-------------Allocation------------------------------------------
    phi_T_bdry.resize( this->n_q_points_bdry, std::vector<double>( (bdry_info.fe(t_rev.fetype_index)).dofs_per_cell ) );
    
    //-----------------------------------------------------------------
    this->JxW_bdry.resize(this->n_q_points_bdry);
    this->normal_vectors.resize(this->n_q_points_bdry);
}

// ---                                    ---
// --- make_assemblers_cell_variable_data ---
// ---                                    ---

template <int dim>
void
NAME::ThermalTransportEquation<dim>::make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                        FuelCellShop::Layer::BaseLayer<dim>* const layer)
{
    Assert( this->n_q_points_cell != 0, ExcMessage("make_assemblers_cell_constant_data function not called before.") );
    
    //---------------Effective Transport Properties---------------------------------------------------------------
    // ----- type infos -------------
    const std::type_info& GasDiffusionLayer = typeid(FuelCellShop::Layer::GasDiffusionLayer<dim>);
    const std::type_info& MicroPorousLayer  = typeid(FuelCellShop::Layer::MicroPorousLayer<dim>);
    const std::type_info& CatalystLayer     = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    const std::type_info& MembraneLayer     = typeid(FuelCellShop::Layer::MembraneLayer<dim>);
    
    const std::type_info& base_layer = layer->get_base_type();
    
    // ----- Assigning the containers -----------------------------------------
    // ----- All containers intialized to zero by default -------------------------------------
    
    sigmaSeff_cell.clear();     // Initializes it to zero tensor
    dkeff_dT_cell.assign(this->n_q_points_cell, Tensor<2,dim>());
    sigmaMeff_cell.assign(this->n_q_points_cell,0.0);
    dsigmaMeff_dT_cell.assign(this->n_q_points_cell,0.0);
    dsigmaMeff_dlambda_cell.assign(this->n_q_points_cell,0.0);
    
    conc_Deff_dHdT_map.clear();
    dT_concDeffdHdT_map.clear();
    ds_concDeffdHdT_map.clear();
    dp_concDeffdHdT_map.clear();
    grad_phi_xi_map.clear();

    electroosmotic_dhdT.assign(this->n_q_points_cell,0.0);
    delectroosmotic_dhdT_dlambda.assign(this->n_q_points_cell,0.0);
    delectroosmotic_dhdT_dT.assign(this->n_q_points_cell,0.0);
    backdiff_dhdT.assign(this->n_q_points_cell,0.0);
    dbackdiff_dhdT_dlambda.assign(this->n_q_points_cell,0.0);
    dbackdiff_dhdT_dT.assign(this->n_q_points_cell,0.0);
    thermoosmotic_dhdT.assign(this->n_q_points_cell,0.0);
    dthermoosmotic_dhdT_dT.assign(this->n_q_points_cell,0.0);
    
    grad_phiM_cell_old.assign(this->n_q_points_cell, Tensor<1,dim>());
    grad_phiS_cell_old.assign(this->n_q_points_cell, Tensor<1,dim>());
    grad_lambda_cell_old.assign(this->n_q_points_cell, Tensor<1,dim>());
    
    // ----- dynamic cast and filling the containers -----------------
    try
    {
        if (base_layer == GasDiffusionLayer)
        {
            FuelCellShop::Layer::GasDiffusionLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::GasDiffusionLayer<dim>* >(layer);

            ptr->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index], temperature_of_REV) );
            ptr->effective_thermal_conductivity(keff_cell);
                        
            if (electron_ohmic_heat_gdl)
            {
                grad_phiS_cell_old = cell_info.gradients[last_iter_cell][phi_s.solution_index];
                ptr->effective_electron_conductivity(sigmaSeff_cell);
            }
            
            if (!cell_residual_counter)
                ptr->derivative_effective_thermal_conductivity(dkeff_dT_cell);
            
            if (enthalpy_fickian_transport)
                gas_enthalpy_transport_assemblers_cell_data(cell_info, ptr);
        }
        
        else if (base_layer == MicroPorousLayer)
        {
            FuelCellShop::Layer::MicroPorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::MicroPorousLayer<dim>* >(layer);
            
            ptr->effective_thermal_conductivity(keff_cell);
            
            if (electron_ohmic_heat_mpl)
            {
                grad_phiS_cell_old = cell_info.gradients[last_iter_cell][phi_s.solution_index];
                ptr->effective_electron_conductivity(sigmaSeff_cell);
            }
            
            if (enthalpy_fickian_transport)
                gas_enthalpy_transport_assemblers_cell_data(cell_info, ptr);
        }
        
        else if (base_layer == CatalystLayer)
        {
            FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);

            ptr->effective_thermal_conductivity(keff_cell);
                
            if (electron_ohmic_heat_cl)
            {
                grad_phiS_cell_old = cell_info.gradients[last_iter_cell][phi_s.solution_index];
                ptr->effective_electron_conductivity(sigmaSeff_cell);
            }
            
            if (proton_ohmic_heat_cl || enthalpy_lambda_transport)
            {
                if (phi_m.indices_exist)
                    grad_phiM_cell_old = cell_info.gradients[last_iter_cell][phi_m.solution_index];
                
                std::vector<VariableNames> deriv_flags;
                ptr->get_electrolyte()->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index],
                                                         temperature_of_REV) );
                deriv_flags.push_back(temperature_of_REV);
                
                if (lambda.indices_exist)
                {
                    ptr->get_electrolyte()->set_membrane_water_content( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][lambda.solution_index],
                                                                        membrane_water_content) );
                    deriv_flags.push_back(membrane_water_content);
                    grad_lambda_cell_old = cell_info.gradients[last_iter_cell][lambda.solution_index];
                }
                
                if (!cell_residual_counter)
                    ptr->set_derivative_flags(deriv_flags);
                
                if (proton_ohmic_heat_cl)
                {
                    ptr->effective_proton_conductivity(sigmaMeff_cell);
                
                    if (!cell_residual_counter)
                    {
                        std::map<VariableNames, std::vector<double> > dsigmaMeff_cell;
                        ptr->derivative_effective_proton_conductivity(dsigmaMeff_cell);
                
                        dsigmaMeff_dT_cell = dsigmaMeff_cell[temperature_of_REV];
                        if (lambda.indices_exist)
                            dsigmaMeff_dlambda_cell = dsigmaMeff_cell[membrane_water_content];
                    }
                }
                
                if (enthalpy_lambda_transport)
                    lambda_enthalpy_transport_assemblers_cell_data(cell_info, ptr);
            }
            
            if (enthalpy_fickian_transport)
                gas_enthalpy_transport_assemblers_cell_data(cell_info, ptr);
        }
        
        else if (base_layer == MembraneLayer)
        {
            FuelCellShop::Layer::MembraneLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::MembraneLayer<dim>* >(layer);
            
            ptr->effective_thermal_conductivity(keff_cell);
                
            if (proton_ohmic_heat_ml || enthalpy_lambda_transport)
            {
                if (phi_m.indices_exist)
                    grad_phiM_cell_old = cell_info.gradients[last_iter_cell][phi_m.solution_index];
                
                std::vector<VariableNames> deriv_flags;
                ptr->get_electrolyte()->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index],
                                                         temperature_of_REV) );
                deriv_flags.push_back(temperature_of_REV);
                
                if (lambda.indices_exist)
                {
                    ptr->get_electrolyte()->set_membrane_water_content( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][lambda.solution_index],
                                                                        membrane_water_content) );
                    deriv_flags.push_back(membrane_water_content);
                    grad_lambda_cell_old = cell_info.gradients[last_iter_cell][lambda.solution_index];
                }
                
                if (!cell_residual_counter)
                    ptr->set_derivative_flags(deriv_flags);
                
                if (proton_ohmic_heat_ml)
                {
                    ptr->effective_proton_conductivity(sigmaMeff_cell);
                
                    if (!cell_residual_counter)
                    {
                        std::map<VariableNames, std::vector<double> > dsigmaMeff_cell;
                        ptr->derivative_effective_proton_conductivity(dsigmaMeff_cell);
                
                        dsigmaMeff_dT_cell = dsigmaMeff_cell[temperature_of_REV];
                        if (lambda.indices_exist)
                            dsigmaMeff_dlambda_cell = dsigmaMeff_cell[membrane_water_content];
                    }
                }
                
                if (enthalpy_lambda_transport)
                    lambda_enthalpy_transport_assemblers_cell_data(cell_info, ptr);
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
        //-------JxW----------
        this->JxW_cell[q] = (cell_info.fe(t_rev.fetype_index)).JxW(q);
        
        //------ Filling shape functions etc ----------------------------------------------------------------------
        //------ This avoids recalculating shape functions etc for efficiency -------------------------------------
        
        for (unsigned int k=0; k < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++k)
        {
            phi_T_cell[q][k] = (cell_info.fe(t_rev.fetype_index)).shape_value(k,q);
            grad_phi_T_cell[q][k] = (cell_info.fe(t_rev.fetype_index)).shape_grad(k,q);
        }
        
        if (!cell_residual_counter)
        {
            //------- Checking based on boolean flags for non-base fe elements------- ---------------------------------
            if ( phi_s.indices_exist )
                for (unsigned int k=0; k < (cell_info.fe(phi_s.fetype_index)).dofs_per_cell; ++k)
                    grad_phi_phiS_cell[q][k] = (cell_info.fe(phi_s.fetype_index)).shape_grad(k,q);
                
            if ( phi_m.indices_exist )
                for (unsigned int k=0; k < (cell_info.fe(phi_m.fetype_index)).dofs_per_cell; ++k)
                    grad_phi_phiM_cell[q][k] = (cell_info.fe(phi_m.fetype_index)).shape_grad(k,q);
                    
            if ( lambda.indices_exist )
            {
                for (unsigned int k=0; k < (cell_info.fe(lambda.fetype_index)).dofs_per_cell; ++k)
                {
                    phi_lambda_cell[q][k] = (cell_info.fe(lambda.fetype_index)).shape_value(k,q);
                    grad_phi_lambda_cell[q][k] = (cell_info.fe(lambda.fetype_index)).shape_grad(k,q);
                }
            }
        }
    }
    
}

// ---                                             ---
// --- gas_enthalpy_transport_assemblers_cell_data ---
// ---                                             ---
template <int dim>
template <typename porelayer>
void
NAME::ThermalTransportEquation<dim>::gas_enthalpy_transport_assemblers_cell_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                                 porelayer* const layer)
{   
    std::vector<FuelCellShop::Material::PureGas*> gases_cell;
    gases_cell = layer->get_gases();
    
    double pressure;
    std::vector< Tensor<2,dim> > Deff;
    std::vector< Tensor<2,dim> > dDeff_dT;
    std::vector< Tensor<2,dim> > dDeff_ds;
    std::vector< Tensor<2,dim> > dDeff_dp;
    
    std::vector<VariableNames> deriv_flags;
    std::map< VariableNames, std::vector< Tensor<2,dim> > > dDeff_du;
    layer->get_p(pressure);        // in atmospheres
    
    for (unsigned int i=0; i < gases_cell.size()-1; ++i)
    {
        std::stringstream gas_var_name;
        gas_var_name << gases_cell[i]->name_material()<<"_molar_fraction";

        layer->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index], temperature_of_REV) );
        deriv_flags.push_back(temperature_of_REV);
        if (s_liquid_water.indices_exist)
        {
            layer->set_saturation( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][s_liquid_water.solution_index], liquid_water_saturation) );
            deriv_flags.push_back(liquid_water_saturation);
        }
        
        if (p_liquid_water.indices_exist)
        {
            layer->set_capillary_pressure( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][p_liquid_water.solution_index], capillary_pressure) );
            deriv_flags.push_back(capillary_pressure);
        }
        
        layer->compute_gas_diffusion(gases_cell[i], gases_cell.back());
        layer->effective_gas_diffusivity(Deff);                          // Filling the vector first with effective diffusivity
        
        layer->set_derivative_flags(deriv_flags);
        if (!cell_residual_counter)
            layer->derivative_effective_gas_diffusivity(dDeff_du);       // Filling the vector first with derivative of effective diffusivity

        dDeff_dT = dDeff_du[temperature_of_REV];
        
        if (s_liquid_water.indices_exist)
            dDeff_ds = dDeff_du[liquid_water_saturation];
        
        if (p_liquid_water.indices_exist)
            dDeff_dp = dDeff_du[capillary_pressure];
        
            
        for ( unsigned int q=0; q < this->n_q_points_cell; ++q)
        {
            const double& T = cell_info.values[last_iter_cell][t_rev.solution_index][q];
            const double concentration = ((pressure*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T))*Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3);  // mol/cm^3
            
            Deff[q] *= Units::convert(1.,Units::C_UNIT2, Units::UNIT2); // converting effective diffusivity from m^2/s to cm^2/s
            
            if (!cell_residual_counter)
            {
                const double dc_dT = (-1.)*((pressure*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T*T))*Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3); // mol/(cm^3-K)
                
                // Modifying the vector directly, avoiding allocating new memory
                dDeff_dT[q] *= ( Units::convert(1.,Units::C_UNIT2, Units::UNIT2) * concentration * (gases_cell[i]->get_Dmolar_enthalpy_Dtemperature(T) - gases_cell.back()->get_Dmolar_enthalpy_Dtemperature(T)) );
                dDeff_dT[q] += ( ( Deff[q] * dc_dT * (gases_cell[i]->get_Dmolar_enthalpy_Dtemperature(T) - gases_cell.back()->get_Dmolar_enthalpy_Dtemperature(T)) ) 
                                                            + 
                                   ( Deff[q] * concentration * (gases_cell[i]->get_D2molar_enthalpy_Dtemperature2(T) - gases_cell.back()->get_D2molar_enthalpy_Dtemperature2(T)) ) );
                
                // liquid_water_saturation term
                if (s_liquid_water.indices_exist)
                {
                    dDeff_ds[q] *= ( Units::convert(1.,Units::C_UNIT2, Units::UNIT2) * concentration * (gases_cell[i]->get_Dmolar_enthalpy_Dtemperature(T) - gases_cell.back()->get_Dmolar_enthalpy_Dtemperature(T)) );
                    
                    for (unsigned int k=0; k < (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell; ++k)
                        phi_s_cell[q][k] = (cell_info.fe(s_liquid_water.fetype_index)).shape_value(k,q);
                }
                
                
                // capillary_pressure term
                if (p_liquid_water.indices_exist)
                {
                    dDeff_dp[q] *= ( Units::convert(1.,Units::C_UNIT2, Units::UNIT2) * concentration * (gases_cell[i]->get_Dmolar_enthalpy_Dtemperature(T) - gases_cell.back()->get_Dmolar_enthalpy_Dtemperature(T)) );
                    
                    for (unsigned int k=0; k < (cell_info.fe(p_liquid_water.fetype_index)).dofs_per_cell; ++k)
                        phi_p_cell[q][k] = (cell_info.fe(p_liquid_water.fetype_index)).shape_value(k,q);
                }
                
                
                
                //Filling grad_phi_xi structure
                Assert( (xi_map.at(gas_var_name.str())).indices_exist, ExcInternalError() );
                (grad_phi_xi_map[gas_var_name.str()]).resize( this->n_q_points_cell, std::vector< Tensor<1,dim> >( (cell_info.fe((xi_map.at(gas_var_name.str())).fetype_index)).dofs_per_cell ) );
                for (unsigned int k=0; k < (cell_info.fe((xi_map.at(gas_var_name.str())).fetype_index)).dofs_per_cell; ++k)
                    grad_phi_xi_map[gas_var_name.str()][q][k] = (cell_info.fe((xi_map.at(gas_var_name.str())).fetype_index)).shape_grad(k,q);
            }
            
            // Modifying the vector directly, avoiding allocating new memory
            Deff[q] *= ( concentration * (gases_cell[i]->get_Dmolar_enthalpy_Dtemperature(T) - gases_cell.back()->get_Dmolar_enthalpy_Dtemperature(T)) );
        }
        
        conc_Deff_dHdT_map[ gas_var_name.str() ] = Deff;
        dT_concDeffdHdT_map[ gas_var_name.str() ] = dDeff_dT;
        ds_concDeffdHdT_map[ gas_var_name.str() ] = dDeff_ds;
        dp_concDeffdHdT_map[ gas_var_name.str() ] = dDeff_dp;
    }
}

// ---                                                ---
// --- lambda_enthalpy_transport_assemblers_cell_data ---
// ---                                                ---
template <int dim>
template <typename polymer>
void
NAME::ThermalTransportEquation<dim>::lambda_enthalpy_transport_assemblers_cell_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                                    const polymer* const layer)
{
    FuelCellShop::Material::WaterVapor water;

    std::vector<double> D_lambda;       // water diffusion coefficient
    std::vector<double> n_d;            // electroosmotic drag coefficient
    std::vector<double> sigma_m;        // protonic conductivity
    std::vector<double> D_T;            // thermo-osmosis coefficient

    std::map< VariableNames, std::vector<double> > dD_lambda;   // derivative of water diffusion coefficient
    std::map< VariableNames, std::vector<double> > dn_d;        // derivative of electroosmotic drag coefficient
    std::map< VariableNames, std::vector<double> > dsigma_m;    // derivative of protonic conductivity
    std::map< VariableNames, std::vector<double> > dD_T;        // derivative of thermo-osmosis coefficient
    
    const double rho_dry = layer->get_electrolyte()->get_density();
    const double EW = layer->get_electrolyte()->get_EW();
    layer->effective_water_diffusivity(D_lambda);
    
    if (phi_m.indices_exist)
    {
        layer->effective_proton_conductivity(sigma_m);
        layer->get_electrolyte()->electroosmotic_drag(n_d);
    }
    
    if (flag_thermoosmosis)
        layer->effective_thermoosmotic_diffusivity(D_T);
    
    if (!cell_residual_counter)
    {
        layer->derivative_effective_water_diffusivity(dD_lambda);
        
        if (phi_m.indices_exist)
        {
            layer->get_electrolyte()->electroosmotic_drag_derivative(dn_d);
            layer->derivative_effective_proton_conductivity(dsigma_m);
        }
        
        if (flag_thermoosmosis)
            layer->derivative_effective_thermoosmotic_diffusivity(dD_T);
    }

    // Looping over quadrature points to fill the vector corresponding to complex terms formed out of various effective transport properties:
    for (unsigned int q = 0; q<this->n_q_points_cell; ++q)
    {
        const double& T = cell_info.values[last_iter_cell][t_rev.solution_index][q];
        
        backdiff_dhdT[q] = (rho_dry/EW) * D_lambda.at(q) * layer->get_electrolyte()->get_dHlambda_dT(T);
        
        if (phi_m.indices_exist)
            electroosmotic_dhdT[q] = ( n_d.at(q)/Constants::F() ) * sigma_m.at(q) * layer->get_electrolyte()->get_dHlambda_dT(T);
        
        if (flag_thermoosmosis)
            thermoosmotic_dhdT[q] = ( D_T.at(q)/(water.get_molar_mass()*1.0e3) ) * layer->get_electrolyte()->get_dHlambda_dT(T);
        
        // -- terms for cell matrix
        if (!cell_residual_counter)
        {
            dbackdiff_dhdT_dlambda[q] = (rho_dry/EW) * layer->get_electrolyte()->get_dHlambda_dT(T) * dD_lambda.at(membrane_water_content)[q];
            dbackdiff_dhdT_dT[q] = (rho_dry/EW) * ( layer->get_electrolyte()->get_dHlambda_dT(T) * dD_lambda.at(temperature_of_REV)[q] 
                                                                            + 
                                                    layer->get_electrolyte()->get_d2Hlambda_dT2(T) * D_lambda.at(q) );
            
            if (phi_m.indices_exist)
            {
                delectroosmotic_dhdT_dlambda[q] = ( layer->get_electrolyte()->get_dHlambda_dT(T)/Constants::F() ) * ( n_d.at(q) * dsigma_m.at(membrane_water_content)[q] 
                                                                                                                                        + 
                                                                                                                      dn_d.at(membrane_water_content)[q] * sigma_m.at(q) );
                delectroosmotic_dhdT_dT[q] = ( n_d.at(q)/Constants::F() ) * ( layer->get_electrolyte()->get_dHlambda_dT(T) * dsigma_m.at(temperature_of_REV)[q] 
                                                                                                                + 
                                                                              layer->get_electrolyte()->get_d2Hlambda_dT2(T) * sigma_m.at(q) );
            }
            
            if (flag_thermoosmosis)
                dthermoosmotic_dhdT_dT[q] = ( 1.0/(water.get_molar_mass()*1.0e3) ) * ( layer->get_electrolyte()->get_dHlambda_dT(T) * dD_T.at(temperature_of_REV)[q] 
                                                                                                                                    + 
                                                                                                     layer->get_electrolyte()->get_d2Hlambda_dT2(T) * D_T.at(q) );
        }
    }
}

// ---                                    ---
// --- make_assemblers_bdry_variable_data ---
// ---                                    ---

template <int dim>
void
NAME::ThermalTransportEquation<dim>::make_assemblers_bdry_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                                        FuelCellShop::Layer::BaseLayer<dim>* const layer)
{
    Assert( this->n_q_points_bdry != 0, ExcMessage("make_assemblers_bdry_constant_data function not called before.") );

    const std::type_info& GasDiffusionLayer = typeid(FuelCellShop::Layer::GasDiffusionLayer<dim>);
    
    const std::type_info& base_layer = layer->get_base_type();
    
    // ----- Assigning the containers -----------------------------------------
    // ----- All containers intialized to zero by default -------------------------------------
    dkeff_dT_bdry.assign(this->n_q_points_bdry, Tensor<2,dim>());
    
    // ----- dynamic cast and filling the containers -----------------
    if (!bdry_residual_counter)
    {
        try
        {
            if ( base_layer == GasDiffusionLayer)
            {
                FuelCellShop::Layer::GasDiffusionLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::GasDiffusionLayer<dim>* >(layer);

                ptr->set_temperature( FuelCellShop::SolutionVariable(&bdry_info.values[last_iter_bdry][t_rev.solution_index], temperature_of_REV) );
                ptr->derivative_effective_thermal_conductivity(dkeff_dT_bdry);
            }
            else
            { /* Do nothing for other layers, as dKeff_dT = 0.0 for other layers */ }
        }
        catch(const std::bad_cast& e)
        {
            const std::type_info& info = typeid(*layer);  
            FcstUtilities::log << "Object of type "<<info.name()<<" not implemented"<< std::endl;
            FcstUtilities::log << e.what() << std::endl;
        }
    }
    
    //---------------------------------------------------------------------------------------------------------------
    //------------Looping over quadrature points in the cell --------------------------------------------------------
    for (unsigned int q=0; q < this->n_q_points_bdry; ++q)
    {
        //-------JxW & normal_vectors----------------------------------------------------------------------------
        this->JxW_bdry[q] = (bdry_info.fe(t_rev.fetype_index)).JxW(q);
        //------ Filling shape functions etc --------------------------------------------------------------------
        //------ This avoids recalculating shape functions etc for efficiency -----------------------------------
        for (unsigned int k=0; k < (bdry_info.fe(t_rev.fetype_index)).dofs_per_cell; ++k)
            phi_T_bdry[q][k] = (bdry_info.fe(t_rev.fetype_index)).shape_value(k,q);
        
        if (!bdry_residual_counter)
            this->normal_vectors[q] = Point<dim>((bdry_info.fe(t_rev.fetype_index)).normal_vector(q));
    }
}
       
// ---                           ---
// ---  EXPLICIT INSTANTIATIONS  ---
// ---                           ---

template class NAME::ThermalTransportEquation<deal_II_dimension>;
//----------------------------------------------------------------
//----------------------------------------------------------------
#define LAYER FuelCellShop::Layer::CatalystLayer<deal_II_dimension>
template void NAME::ThermalTransportEquation<deal_II_dimension>::gas_enthalpy_transport_assemblers_cell_data< LAYER >(const typename FuelCell::ApplicationCore::DoFApplication<deal_II_dimension>::CellInfo&,
                                                                                                                      LAYER* const);
#undef LAYER
//----------------------------------------------------------------
#define LAYER FuelCellShop::Layer::MicroPorousLayer<deal_II_dimension>
template void NAME::ThermalTransportEquation<deal_II_dimension>::gas_enthalpy_transport_assemblers_cell_data< LAYER >(const typename FuelCell::ApplicationCore::DoFApplication<deal_II_dimension>::CellInfo&,
                                                                                                                      LAYER* const);
#undef LAYER
//----------------------------------------------------------------
#define LAYER FuelCellShop::Layer::GasDiffusionLayer<deal_II_dimension>
template void NAME::ThermalTransportEquation<deal_II_dimension>::gas_enthalpy_transport_assemblers_cell_data< LAYER >(const typename FuelCell::ApplicationCore::DoFApplication<deal_II_dimension>::CellInfo&,
                                                                                                                      LAYER* const);
#undef LAYER
//----------------------------------------------------------------
//----------------------------------------------------------------
#define LAYER FuelCellShop::Layer::CatalystLayer<deal_II_dimension>
template void NAME::ThermalTransportEquation<deal_II_dimension>::lambda_enthalpy_transport_assemblers_cell_data< LAYER >(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo&,
                                                                                                                         const LAYER* const);
#undef LAYER
//----------------------------------------------------------------
#define LAYER FuelCellShop::Layer::MembraneLayer<deal_II_dimension>
template void NAME::ThermalTransportEquation<deal_II_dimension>::lambda_enthalpy_transport_assemblers_cell_data< LAYER >(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo&,
                                                                                                                         const LAYER* const);
#undef LAYER
