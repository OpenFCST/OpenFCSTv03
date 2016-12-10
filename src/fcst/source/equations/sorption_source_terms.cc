//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: sorption_source_terms.cc
//    - Description: This class is used to assemble cell matrix and cell residual
//      corresponding to sorption/desorption of water inside the catalyst layer.
//    - Developers: Madhur Bhaiya
//
//---------------------------------------------------------------------------

#include "equations/sorption_source_terms.h"

namespace NAME = FuelCellShop::Equation;

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::SorptionSourceTerms<dim>::SorptionSourceTerms(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data)
:
NAME::EquationBase<dim>(system_management,data)
{
    FcstUtilities::log << "->FuelCellShop::Equation::SorptionSourceTerms" << std::endl;
    
    //----Initializing VariableInfo Structs------------------------------------------
    //----Setting indices_exist to false --------------------------------------------
    x_water.indices_exist = false;
    lambda.indices_exist = false;
    t_rev.indices_exist = false;
    
    this->counter.resize(2, true);
}

// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::SorptionSourceTerms<dim>::~SorptionSourceTerms()
{}

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::SorptionSourceTerms<dim>::declare_parameters(ParameterHandler& param) const
{
    param.enter_subsection("Sorption Source Terms");
    {
        param.declare_entry("Water soption time constant [1/s]",
                            "10000.0",
                            Patterns::Double(0.0),
                            "Time constant for sorption isotherm. Units [1/s]");
        
        param.declare_entry("Heat source/sink due to sorption/desorption",
                            "false",
                            Patterns::Bool(),
                            "Flag to include heat release/absorption due to sorption/desorption of water inside the catalyst layer.");
    }
    param.leave_subsection();
}

// ---            ---
// --- initialize ---
// ---            ---

template<int dim>
void
NAME::SorptionSourceTerms<dim>::initialize(ParameterHandler& param)
{
    param.enter_subsection("Sorption Source Terms");
    {
        time_constant = param.get_double("Water soption time constant [1/s]");
        flag_sorp_heat_cl = param.get_bool("Heat source/sink due to sorption/desorption");
    }
    param.leave_subsection();
    
    //-------Assertion check that x_water and lambda should be the solution variables to account for sorption/desorption process --------------------
    AssertThrow( this->system_management->solution_in_userlist("water_molar_fraction"),VariableNotFoundForSorption("water_molar_fraction", "sorption/desorption") );
    AssertThrow( this->system_management->solution_in_userlist("membrane_water_content"), VariableNotFoundForSorption("membrane_water_content", "sorption/desorption") );
    
    if (flag_sorp_heat_cl)
        AssertThrow( this->system_management->solution_in_userlist("temperature_of_REV"), VariableNotFoundForSorption("temperature_of_REV", "heat source/sink due to sorption/desorption") );
}

// ---                      ---
// --- assemble_cell_matrix ---
// ---                      ---

template<int dim>
void
NAME::SorptionSourceTerms<dim>::assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
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
        
    assemble_matrix_for_equation(cell_matrices, cell_info, "Ficks Transport Equation - water", cell_info.fe(x_water.fetype_index), phi_xWater_cell, 1.);
    assemble_matrix_for_equation(cell_matrices, cell_info, "Membrane Water Content Transport Equation", cell_info.fe(lambda.fetype_index), phi_lambda_cell, -1.);
    if (flag_sorp_heat_cl)
        assemble_matrix_for_equation(cell_matrices, cell_info, "Thermal Transport Equation", cell_info.fe(t_rev.fetype_index), phi_T_cell, -1.);
}

// ---                        ---
// --- assemble_cell_residual ---
// ---                        ---

template<int dim>
void
NAME::SorptionSourceTerms<dim>::assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
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
        // ---- Ficks Transport Equation - water ------------------------------
        for (unsigned int i=0; i < (cell_info.fe(x_water.fetype_index)).dofs_per_cell; ++i)
            cell_residual.block(x_water.solution_index)(i) += ( this->JxW_cell[q] * phi_xWater_cell[q][i] * ((time_constant*rho_dry_cell)/EW_cell) * 
                                                                (lambda_eq_cell[q] - (cell_info.values[last_iter_cell][lambda.solution_index][q])) );
        
        // ---- Membrane Water Content Transport Equation ------------------------------
        for (unsigned int i=0; i < (cell_info.fe(lambda.fetype_index)).dofs_per_cell; ++i)
            cell_residual.block(lambda.solution_index)(i) += ( this->JxW_cell[q] * phi_lambda_cell[q][i] * ((time_constant*rho_dry_cell)/EW_cell) * (-1.0) * 
                                                                (lambda_eq_cell[q] - (cell_info.values[last_iter_cell][lambda.solution_index][q])) );
        
        // ---- Thermal Transport Equation ------------------------
        if (flag_sorp_heat_cl)
            for (unsigned int i=0; i < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++i)
                cell_residual.block(t_rev.solution_index)(i) += ( this->JxW_cell[q] * phi_T_cell[q][i] * ((time_constant*rho_dry_cell)/EW_cell) * (-1.0) * h_sorp_cell[q] * 
                                                                    (lambda_eq_cell[q] - (cell_info.values[last_iter_cell][lambda.solution_index][q])) );
    } // End Loop Over Quadrature Points
}

// ---                              ---
// --- adjust_internal_cell_couplings ---
// ---                              ---

template<int dim>
void
NAME::SorptionSourceTerms<dim>::adjust_internal_cell_couplings(std::vector< couplings_map >& equation_map) const
{
    Assert( equation_map.size() != 0, ExcMessage("Vector size should be greater than zero in SorptionSourceTerms::adjust_internal_cell_couplings.") );
    
    for (unsigned int i=0; i<equation_map.size(); ++i)
    {
        for ( couplings_map::iterator iter = equation_map[i].begin(); iter != equation_map[i].end(); ++iter )
        {
            if ( (iter->first == "Ficks Transport Equation - water") || (iter->first == "Membrane Water Content Transport Equation") )
            {
                (iter->second)["water_molar_fraction"] = DoFTools::always;
                (iter->second)["membrane_water_content"] = DoFTools::always;
                
                if (flag_sorp_heat_cl)
                    (iter->second)["temperature_of_REV"] = DoFTools::always;
            }
            
            else if ( iter->first == "Thermal Transport Equation" )
            {
                (iter->second)["temperature_of_REV"] = DoFTools::always;
                
                if (flag_sorp_heat_cl)
                {
                    (iter->second)["water_molar_fraction"] = DoFTools::always;
                    (iter->second)["membrane_water_content"] = DoFTools::always;
                }
            }
        }
    }
    
}

// ---                     ---
// --- print_equation_info ---
// ---                     ---

template<int dim>
void
NAME::SorptionSourceTerms<dim>::print_equation_info() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "-------------------------------------------------------------------------------" << std::endl;
    FcstUtilities::log << std::endl;
    
    FcstUtilities::log << "PARAMETERS for \"Sorption Source Terms\":" << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "Water soption time constant [1/s]: " << time_constant << std::endl;
    FcstUtilities::log << "Heat source/sink due to sorption/desorption: " << flag_sorp_heat_cl << std::endl;
    
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
NAME::SorptionSourceTerms<dim>::make_assemblers_generic_constant_data()
{    
    //--------Block indices can't be filled here, as they depend on what equation we are in -------------------------------------------
    //--------While doing cell_matrix assembly, developers need to be wary of the fact that block_indices are still not filled yet-----
    
    //-----------Filling VariableInfo structures----------------------------------------------------------
    //----------water_molar_fraction--------------------------------------------------------------
    x_water.solution_index = this->system_management->solution_name_to_index("water_molar_fraction");
    x_water.fetype_index = this->system_management->block_info->base_element[x_water.solution_index];
    x_water.indices_exist = true;
    
    //----------membrane_water_content--------------------------------------------------------------
    lambda.solution_index = this->system_management->solution_name_to_index("membrane_water_content");
    lambda.fetype_index = this->system_management->block_info->base_element[lambda.solution_index];
    lambda.indices_exist = true;
    
    if (flag_sorp_heat_cl) // It indirectly checks whether temperature_of_solid_phase is in user-defined list or not.
    {
        //-----------temperature_of_solid_phase-------------------------------------------------------
        t_rev.solution_index = this->system_management->solution_name_to_index("temperature_of_REV");
        t_rev.fetype_index = this->system_management->block_info->base_element[t_rev.solution_index];
        t_rev.indices_exist = true;
    }
}

// ---                                    ---
// --- make_assemblers_cell_constant_data ---
// ---                                    ---

template<int dim>
void
NAME::SorptionSourceTerms<dim>::make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
    Assert( (x_water.indices_exist && lambda.indices_exist), ExcMessage("make_assemblers_generic_constant_data function not called before.") );
    
    this->n_q_points_cell = (cell_info.fe(x_water.fetype_index)).n_quadrature_points;
    last_iter_cell = cell_info.global_data->find_vector(this->solution_vector_name);
    
    //-------Allocation------------------------------------------------------------------------
    // ----- All containers intialized to zero by default -------------------------------------
    phi_xWater_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(x_water.fetype_index)).dofs_per_cell ) );
    phi_lambda_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(lambda.fetype_index)).dofs_per_cell ) );
    if (t_rev.indices_exist)
        phi_T_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(t_rev.fetype_index)).dofs_per_cell ) );
    
    //-----------------------------------------------------------------
    this->JxW_cell.resize(this->n_q_points_cell);
    
    lambda_eq_cell.resize(this->n_q_points_cell);
    dlambdaEq_dxWater_cell.resize(this->n_q_points_cell);
    dlambdaEq_dT_cell.resize(this->n_q_points_cell);
    
    h_sorp_cell.resize(this->n_q_points_cell);
    dhsorp_dT_cell.resize(this->n_q_points_cell);
}

// ---                                    ---
// --- make_assemblers_cell_variable_data ---
// ---                                    ---

template<int dim>
void
NAME::SorptionSourceTerms<dim>::make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                   FuelCellShop::Layer::BaseLayer<dim>* const layer)
{
    Assert( this->n_q_points_cell != 0, ExcMessage("make_assemblers_cell_constant_data function not called before.") );
    
    // ----- type infos -------------
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    
    const std::type_info& base_layer = layer->get_base_type();
    
    // ----- dynamic cast and filling the containers -----------------
    try
    {
        if ( base_layer == CatalystLayer )
        {
            FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);
            
            rho_dry_cell = ptr->get_electrolyte()->get_density();
            EW_cell = ptr->get_electrolyte()->get_EW();
            
            std::vector<VariableNames> deriv_flags;
            
            ptr->get_electrolyte()->set_water_molar_fraction( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][x_water.solution_index],
                                                              water_molar_fraction) );
            deriv_flags.push_back(water_molar_fraction);
            
            if (t_rev.indices_exist)
            {
                ptr->get_electrolyte()->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index],
                                                         temperature_of_REV) );
                deriv_flags.push_back(temperature_of_REV);
            }
            
            ptr->get_electrolyte()->sorption_isotherm(lambda_eq_cell);
            if (flag_sorp_heat_cl)
                ptr->get_electrolyte()->sorption_enthalpy(h_sorp_cell);
            
            if (!cell_residual_counter)
            {
                ptr->set_derivative_flags(deriv_flags);
                std::map< VariableNames, std::vector<double> > dlambdaEq;
                ptr->get_electrolyte()->sorption_isotherm_derivative(dlambdaEq);
                dlambdaEq_dxWater_cell = dlambdaEq[water_molar_fraction];
                if (t_rev.indices_exist)
                    dlambdaEq_dT_cell = dlambdaEq[temperature_of_REV];
                
                if (flag_sorp_heat_cl)
                {
                    std::map< VariableNames, std::vector<double> > dhsorp;
                    ptr->get_electrolyte()->sorption_enthalpy_derivative(dhsorp);
                    dhsorp_dT_cell = dhsorp[temperature_of_REV];
                }
            }
            
        }
        else
            AssertThrow( false, ExcNotImplemented() );
        
    }
    
    catch(const std::bad_cast& e)
    {
        const std::type_info& info = typeid(*layer);  
        FcstUtilities::log << "Object of type " << info.name() << " not implemented" << std::endl;
        FcstUtilities::log << e.what() << std::endl;
    }
    
    //---------------------------------------------------------------------------------------------------------------
    //------------Looping over quadrature points in the cell --------------------------------------------------------
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
        //-------JxW----------
        this->JxW_cell[q] = (cell_info.fe(x_water.fetype_index)).JxW(q);
        
        //------ Filling shape functions etc ----------------------------------------------------------------------
        //------ This avoids recalculating shape functions etc for efficiency -------------------------------------
        
        for (unsigned int k=0; k < (cell_info.fe(x_water.fetype_index)).dofs_per_cell; ++k)
            phi_xWater_cell[q][k] = (cell_info.fe(x_water.fetype_index)).shape_value(k,q);
        
        for (unsigned int k=0; k < (cell_info.fe(lambda.fetype_index)).dofs_per_cell; ++k)
            phi_lambda_cell[q][k] = (cell_info.fe(lambda.fetype_index)).shape_value(k,q);
        
        //------- Checking based on boolean flags for other fe elements--------------------------------------------
        if (t_rev.indices_exist)
            for (unsigned int k=0; k < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++k)
                phi_T_cell[q][k] = (cell_info.fe(t_rev.fetype_index)).shape_value(k,q);
    }
    
}

// ---                              ---
// --- assemble_matrix_for_equation ---
// ---                              ---

template<int dim>
void
NAME::SorptionSourceTerms<dim>::assemble_matrix_for_equation(FuelCell::ApplicationCore::MatrixVector& cell_matrices,
                                                            const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                            const std::string& eq_name,
                                                            const FEValuesBase<dim>& test_fe,
                                                            const std::vector< std::vector<double> >& test_shape_functions,
                                                            const double& sourceterm_factor)
{
    Assert( !cell_residual_counter, ExcInternalError() );
    Assert( this->n_q_points_cell != 0, ExcMessage("make_assemblers_cell_constant_data function not called before.") );
    Assert( ((eq_name == "Ficks Transport Equation - water") || (eq_name == "Membrane Water Content Transport Equation") || (eq_name == "Thermal Transport Equation")), ExcNotImplemented() );
    
    // --- Filling block indices --------------------------------------------
    x_water.block_index = this->system_management->matrix_block_index(eq_name,"water_molar_fraction");
    lambda.block_index = this->system_management->matrix_block_index(eq_name,"membrane_water_content");
    if (t_rev.indices_exist)
        t_rev.block_index = this->system_management->matrix_block_index(eq_name,"temperature_of_REV");
    
    //-------- Looping over Quadrature points ----------------------------
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
        //---------------LOOP over i -----------------------------------------------------------------
        for (unsigned int i=0; i < test_fe.dofs_per_cell; ++i)
        {
            if (eq_name != "Thermal Transport Equation")        // Ficks Transport Equation - water OR Membrane Water Content Transport Equation
            {
                //--------------LOOP(s) over j-------------------------------------------------------------
                for (unsigned int j=0; j < (cell_info.fe(x_water.fetype_index)).dofs_per_cell; ++j)
                    cell_matrices[x_water.block_index].matrix(i,j) += ( this->JxW_cell[q] * sourceterm_factor * ((time_constant*rho_dry_cell)/EW_cell) * dlambdaEq_dxWater_cell[q] * 
                                                                        test_shape_functions[q][i] * phi_xWater_cell[q][j] );
                
                for (unsigned int j=0; j < (cell_info.fe(lambda.fetype_index)).dofs_per_cell; ++j)
                    cell_matrices[lambda.block_index].matrix(i,j) += ( this->JxW_cell[q] * sourceterm_factor * (-1.0) * ((time_constant*rho_dry_cell)/EW_cell) * 
                                                                        test_shape_functions[q][i] * phi_lambda_cell[q][j] );
                
                if (t_rev.indices_exist)
                    for (unsigned int j=0; j < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++j)
                        cell_matrices[t_rev.block_index].matrix(i,j) += ( this->JxW_cell[q] * sourceterm_factor * ((time_constant*rho_dry_cell)/EW_cell) * dlambdaEq_dT_cell[q] * 
                                                                            test_shape_functions[q][i] * phi_T_cell[q][j] );
            }
            
            else // Thermal Transport Equation
            {
                //--------------LOOP(s) over j-------------------------------------------------------------
                for (unsigned int j=0; j < (cell_info.fe(x_water.fetype_index)).dofs_per_cell; ++j)
                    cell_matrices[x_water.block_index].matrix(i,j) += ( this->JxW_cell[q] * sourceterm_factor * ((time_constant*rho_dry_cell)/EW_cell) * dlambdaEq_dxWater_cell[q] * h_sorp_cell[q] * 
                                                                        test_shape_functions[q][i] * phi_xWater_cell[q][j] );
                
                for (unsigned int j=0; j < (cell_info.fe(lambda.fetype_index)).dofs_per_cell; ++j)
                    cell_matrices[lambda.block_index].matrix(i,j) += ( this->JxW_cell[q] * sourceterm_factor * ((time_constant*rho_dry_cell)/EW_cell) * 
                                                                        ( (-1.0) * h_sorp_cell[q]) * 
                                                                        test_shape_functions[q][i] * phi_lambda_cell[q][j] );
                
                for (unsigned int j=0; j < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++j)
                    cell_matrices[t_rev.block_index].matrix(i,j) += ( this->JxW_cell[q] * sourceterm_factor * ((time_constant*rho_dry_cell)/EW_cell) * 
                                                                        ((lambda_eq_cell[q]-cell_info.values[last_iter_cell][lambda.solution_index][q])*dhsorp_dT_cell[q] + dlambdaEq_dT_cell[q]*h_sorp_cell[q]) * 
                                                                        test_shape_functions[q][i] * phi_T_cell[q][j] );  
            }
            
        } // End Loop over "i"
    } // End Loop over "q"
}

// ---                           ---
// ---  EXPLICIT INSTANTIATIONS  ---
// ---                           ---

template class NAME::SorptionSourceTerms<deal_II_dimension>;
