//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: proton_transport_equation.cc
//    - Description: Implementation of proton_transport_equation.h. This
//      class is used to setup the matrix and rhs for proton transport based
//      on Ohms law.
//    - Developers: Madhur Bhaiya, M. Secanell, Valentin N. Zingan
//
//---------------------------------------------------------------------------

#include "equations/proton_transport_equation.h"

namespace NAME = FuelCellShop::Equation;

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::ProtonTransportEquation<dim>::ProtonTransportEquation(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data)
:
EquationBase<dim>(system_management,data)
{
    this->name_base_variable = "protonic_electrical_potential";
    this->equation_name = "Proton Transport Equation";
    FcstUtilities::log << "->FuelCellShop::Equation::ProtonTransportEquation" << std::endl;

    //----Initializing VariableInfo Structs------------------------------------------
    //----Setting indices_exist to false --------------------------------------------
    phi_m.indices_exist = false;
    lambda.indices_exist = false;
    t_rev.indices_exist = false;

    //-- Counter to store if generic_data and constant_data need to be initialized.
    this->counter.resize(3, true);
}

// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::ProtonTransportEquation<dim>::~ProtonTransportEquation()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::ProtonTransportEquation<dim>::declare_parameters(ParameterHandler& param) const
{
    NAME::EquationBase<dim>::declare_parameters(param);
    param.enter_subsection("Equations");
    {
        param.enter_subsection("Proton Transport Equation");
        {
            param.enter_subsection("Boundary data");
            {
                
                param.declare_entry("current_flux",
                                    """",
                                    Patterns::Map(   Patterns::Integer(0) , Patterns::Double()   ),
                                    " ");
            }
            param.leave_subsection();
            
            
            param.enter_subsection("Boundary conditions");
            {
                param.declare_entry("Dirichlet Boundary Conditions",
                                    "",
                                    Patterns::Map( Patterns::Integer(0), Patterns::Double() ),
                                    "A comma-separated list of Dirichlet boundary_id(s) with prescribed values of Potential [V]."
                                    "\n"
                                    "Correct format is of a map, given as ''id1: value1, id2: value2, id3: value3, ...'' where "
                                    "each id must be an unsigned integer, and each value can be a double corresponding to potential value.");
                
                param.declare_entry("Constant Proton Current Flux Boundary Conditions",
                                    "",
                                    Patterns::Map( Patterns::Integer(0), Patterns::Double() ),
                                    "A comma-separated list of boundaries with Galvanostatic b.c., with prescribed values of protonic current fluxes [A/cm^2]."
                                    "\n"
                                    "Correct format is of a map, given as ''id1: value1, id2: value2, id3: value3, ...'' where "
                                    "each id must be an unsigned integer, and each value can be a double (positive for protonic current flux leaving out and negative for coming in).");
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

template<int dim>
void
NAME::ProtonTransportEquation<dim>::initialize(ParameterHandler& param)
{
    NAME::EquationBase<dim>::initialize(param);
    
    param.enter_subsection("Equations");
    {
        param.enter_subsection("Proton Transport Equation");
        {
            param.enter_subsection("Boundary data");
            {
                
                if( !param.get("current_flux").empty() )
                {
                    // proton_current_flux_map = FcstUtilities::string_to_map<unsigned int, double>( param.get("current_flux") ); // will be uncommented later on
                }
            }
            param.leave_subsection();
            
            param.enter_subsection("Boundary conditions");
            {
                dirichlet_bdry_map = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Dirichlet Boundary Conditions") ) );
                proton_current_flux_map = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Constant Proton Current Flux Boundary Conditions") ) );
            }
            param.leave_subsection();
            
        }
        param.leave_subsection();
    }
    param.leave_subsection();

    this->make_internal_cell_couplings();
    this->make_boundary_types();
}

// ---                      ---
// --- assemble_cell_matrix ---
// ---                      ---

template<int dim>
void
NAME::ProtonTransportEquation<dim>::assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
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
        for (unsigned int i=0; i < (cell_info.fe(phi_m.fetype_index)).dofs_per_cell; ++i)
        {
            //--------------LOOP(s) over j-------------------------------------------------------------

            //-----------Assembling Matrix for terms corresponding to "phi_m" BLOCK------------------------
            //----------- If some term is not be included for a particular layer, it will automatically calculate out to zero--------------------
            for (unsigned int j=0; j < (cell_info.fe(phi_m.fetype_index)).dofs_per_cell; ++j)
            {
                cell_matrices[phi_m.block_index].matrix(i,j) += grad_phi_phiM_cell[q][i] * sigmaMeff_cell[q] * grad_phi_phiM_cell[q][j] * this->JxW_cell[q];
            }

            //----------TERM CORRESPONDING TO "LAMBDA" BLOCK-----------------------------------------
            if ( lambda.indices_exist )
            {
                //-----Assembling Matrix------------------------------------------------------------
                for (unsigned int j=0; j < (cell_info.fe(lambda.fetype_index)).dofs_per_cell; ++j)
                {
                    cell_matrices[lambda.block_index].matrix(i,j) += grad_phi_phiM_cell[q][i] * cell_info.gradients[last_iter_cell][phi_m.solution_index][q] * dsigmaMeff_dlambda_cell[q] *
                                                                     phi_lambda_cell[q][j] * this->JxW_cell[q];
                }
            }

            //----------TERM CORRESPONDING TO "T" BLOCK-----------------------------------------
            if ( t_rev.indices_exist )
            {
                //-----Assembling Matrix------------------------------------------------------------
                for (unsigned int j=0; j < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++j)
                {
                    cell_matrices[t_rev.block_index].matrix(i,j) += grad_phi_phiM_cell[q][i] * cell_info.gradients[last_iter_cell][phi_m.solution_index][q] * dsigmaMeff_dT_cell[q] *
                                                                    phi_T_cell[q][j] * this->JxW_cell[q];
                }
            }

        }
    }
}

// ---                        ---
// --- assemble_cell_residual ---
// ---                        ---
template<int dim>
void
NAME::ProtonTransportEquation<dim>::assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
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
        for (unsigned int i=0; i < (cell_info.fe(phi_m.fetype_index)).dofs_per_cell; ++i)
        {
            cell_residual.block(phi_m.solution_index)(i) += grad_phi_phiM_cell[q][i] * sigmaMeff_cell[q] * cell_info.gradients[last_iter_cell][phi_m.solution_index][q] * this->JxW_cell[q];
        }
    }
}

// ---                      ---
// --- assemble_bdry_matrix ---
// ---                      ---

template<int dim>
void
NAME::ProtonTransportEquation<dim>::assemble_bdry_matrix(FuelCell::ApplicationCore::MatrixVector&                                 bdry_matrices,
                                                         const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                         FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
       // The boundary integral " q * [ ... ] " is always zero because of the boundary conditions we use:
       // - if phi_m = known then q = 0
       // - if [-sigma_m * grad_phi_m * n] = known then [ ... ] = 0
}

// ---                        ---
// --- assemble_bdry_residual ---
// ---                        ---

template<int dim>
void
NAME::ProtonTransportEquation<dim>::assemble_bdry_residual(FuelCell::ApplicationCore::FEVector&                                     bdry_residual,
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

    //------ Constant Proton Current Flux boundaries --------------------------------------------------------------
    std::map<unsigned int, double>::const_iterator iter = proton_current_flux_map.find( bdry_info.dof_face->boundary_indicator() );
    if (iter != proton_current_flux_map.end() )
    {
        const double& proton_current_value = iter->second;

        bdry_residual_counter = true;
        this->make_assemblers_bdry_variable_data(bdry_info, layer);

        //-------- Looping over Quadrature points ----------------------------
        for (unsigned int q = 0; q < this->n_q_points_bdry; ++q)
        {
            for (unsigned int i = 0; i < (bdry_info.fe(phi_m.fetype_index)).dofs_per_cell; ++i)
            {
                bdry_residual.block(phi_m.solution_index)(i) += phi_phiM_bdry[q][i] * proton_current_value * this->JxW_bdry[q];
            }
        }
    }
}

// ---                              ---
// --- make_internal_cell_couplings ---
// ---                              ---
template<int dim>
void
NAME::ProtonTransportEquation<dim>::make_internal_cell_couplings()
{
    AssertThrow(this->system_management->solution_in_userlist(this->name_base_variable), VariableShouldExistForEquation(this->name_base_variable, this->equation_name) );
    AssertThrow(this->system_management->equation_name_to_index(this->equation_name) == this->system_management->solution_name_to_index(this->name_base_variable),
              IndexDoNotMatch(this->name_base_variable, this->equation_name) );

    std::map< std::string, DoFTools::Coupling > tmp;

    std::vector< std::string> sol_names = this->system_management->get_solution_names();

    for (unsigned int i = 0; i < sol_names.size(); ++i)
    {
        if (sol_names[i] == this->name_base_variable)
            tmp[this->name_base_variable] = DoFTools::always;

        else if (sol_names[i] == "membrane_water_content")
            tmp["membrane_water_content"] = DoFTools::always;

        else if (sol_names[i] == "temperature_of_REV")
            tmp["temperature_of_REV"] = DoFTools::always;

        else
            tmp[sol_names[i]] = DoFTools::none;
    }

    this->internal_cell_couplings[this->equation_name] = tmp;
}

// ---                     ---
// --- make_boundary_types ---
// ---                     ---

template<int dim>
void
NAME::ProtonTransportEquation<dim>::make_boundary_types()
{
    unsigned int i = 0;
    for ( std::map<unsigned int,double>::const_iterator iter = dirichlet_bdry_map.begin(); iter != dirichlet_bdry_map.end(); ++iter)
    {
        i++;
        BoundaryType temp_dirich;

        std::ostringstream streamOut;
        streamOut << i;
        std::string temp_name = "Dirichlet_" + streamOut.str();

        temp_dirich.boundary_name = temp_name;
        temp_dirich.boundary_id = iter->first;
        temp_dirich.boundary_condition = "Dirichlet";

        this->boundary_types.push_back(temp_dirich);
    }
}

// ---                     ---
// --- print_equation_info ---
// ---                     ---

template<int dim>
void
NAME::ProtonTransportEquation<dim>::print_equation_info() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "-------------------------------------------------------------------------------" << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "INTERNAL CELL COUPLINGS FOR \" Proton Transport Equation \":" << std::endl;
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

  FcstUtilities::log << "Initial data:";
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;

  for(component_materialID_value_map::const_iterator iter  = this->component_materialID_value.begin();
                                                     iter != this->component_materialID_value.end();
                                                   ++iter)
  {
         std::map<types::material_id, double> tmp = iter->second;

         for(std::map<types::material_id, double>::const_iterator iter2  = tmp.begin();
                                                                  iter2 != tmp.end();
                                                                ++iter2)
         {
                FcstUtilities::log << "Name of the solution component  = " << iter->first   << std::endl;
                FcstUtilities::log << "Material id                     = " << iter2->first  << std::endl;
                FcstUtilities::log << "Value of the solution component = " << iter2->second << std::endl;
                FcstUtilities::log << std::endl;
         }
  }

  FcstUtilities::log << "Boundary data:";
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;

  for(component_boundaryID_value_map::const_iterator iter  = this->component_boundaryID_value.begin();
                                                     iter != this->component_boundaryID_value.end();
                                                   ++iter)
  {
         std::map<types::boundary_id, double> tmp = iter->second;

         for(std::map<types::boundary_id, double>::const_iterator iter2  = tmp.begin();
                                                                  iter2 != tmp.end();
                                                                ++iter2)
         {
                FcstUtilities::log << "Name of the solution component  = " << iter->first   << std::endl;
                FcstUtilities::log << "Boundary id                     = " << iter2->first  << std::endl;
                FcstUtilities::log << "Value of the solution component = " << iter2->second << std::endl;
                FcstUtilities::log << std::endl;
         }
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
NAME::ProtonTransportEquation<dim>::make_assemblers_generic_constant_data()
{
    //-----------Filling VariableInfo structures------------------------------------------
    //----------protonic_electrical_potential---------------------------------------------
    phi_m.solution_index = this->system_management->solution_name_to_index(this->name_base_variable);
    phi_m.block_index = this->system_management->matrix_block_index(this->equation_name, this->name_base_variable);
    phi_m.fetype_index = this->system_management->block_info->base_element[phi_m.solution_index];
    phi_m.indices_exist = true;

    //----------membrane_water_content----------------------------------------------------
    if ( this->system_management->solution_in_userlist("membrane_water_content") )
    {
        lambda.solution_index = this->system_management->solution_name_to_index("membrane_water_content");
        lambda.block_index = this->system_management->matrix_block_index(this->equation_name, "membrane_water_content");
        lambda.fetype_index = this->system_management->block_info->base_element[lambda.solution_index];
        lambda.indices_exist = true;
    }

    //-----------temperature_of_solid_phase-----------------------------------------------
    if ( this->system_management->solution_in_userlist("temperature_of_REV") )
    {
        t_rev.solution_index = this->system_management->solution_name_to_index("temperature_of_REV");
        t_rev.block_index = this->system_management->matrix_block_index(this->equation_name, "temperature_of_REV");
        t_rev.fetype_index = this->system_management->block_info->base_element[t_rev.solution_index];
        t_rev.indices_exist = true;
    }
}

// ---                                    ---
// --- make_assemblers_cell_constant_data ---
// ---                                    ---

template<int dim>
void
NAME::ProtonTransportEquation<dim>::make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
    Assert( phi_m.indices_exist, ExcMessage("make_assemblers_generic_constant_data function not called before.") );

    this->n_q_points_cell = (cell_info.fe(phi_m.fetype_index)).n_quadrature_points;
    last_iter_cell = cell_info.global_data->find_vector(this->solution_vector_name);

    //-------------Allocation------------------------------------------

    sigmaMeff_cell.resize(this->n_q_points_cell);
    dsigmaMeff_dlambda_cell.resize(this->n_q_points_cell);
    dsigmaMeff_dT_cell.resize(this->n_q_points_cell);

    grad_phi_phiM_cell.resize( this->n_q_points_cell, std::vector< Tensor<1,dim> >( (cell_info.fe(phi_m.fetype_index)).dofs_per_cell ) );

    if ( lambda.indices_exist )
    {
        phi_lambda_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(lambda.fetype_index)).dofs_per_cell ) );
    }

    if ( t_rev.indices_exist )
    {
        phi_T_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(t_rev.fetype_index)).dofs_per_cell ) );
    }

    this->JxW_cell.resize(this->n_q_points_cell);
}

// ---                                    ---
// --- make_assemblers_bdry_constant_data ---
// ---                                    ---

template<int dim>
void
NAME::ProtonTransportEquation<dim>::make_assemblers_bdry_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info)
{
    Assert( phi_m.indices_exist, ExcMessage("make_assemblers_generic_constant_data function not called before.") );

    this->n_q_points_bdry = (bdry_info.fe(phi_m.fetype_index)).n_quadrature_points;
    last_iter_bdry = bdry_info.global_data->find_vector(this->solution_vector_name);

    //-------------Allocation------------------------------------------
    phi_phiM_bdry.resize( this->n_q_points_bdry, std::vector<double>( (bdry_info.fe(phi_m.fetype_index)).dofs_per_cell ) );

    //-----------------------------------------------------------------
    this->JxW_bdry.resize(this->n_q_points_bdry);
}

// ---                                    ---
// --- make_assemblers_cell_variable_data ---
// ---                                    ---

template<int dim>
void
NAME::ProtonTransportEquation<dim>::make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                         FuelCellShop::Layer::BaseLayer<dim>* const layer)
{
    Assert( this->n_q_points_cell != 0, ExcMessage("make_assemblers_cell_constant_data function not called before.") );

    //---------------Effective Transport Properties---------------------------------------------------------------
    // ----- type infos -------------
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    const std::type_info& MembraneLayer = typeid(FuelCellShop::Layer::MembraneLayer<dim>);

    const std::type_info& base_layer = layer->get_base_type();

    // ----- dynamic cast and filling the containers -----------------
    try
    {
        if (base_layer == CatalystLayer)
        {
            FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);

            if ( lambda.indices_exist )
                ptr->get_electrolyte()->set_membrane_water_content( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][lambda.solution_index], membrane_water_content) );
            if ( t_rev.indices_exist )
                ptr->get_electrolyte()->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index], temperature_of_REV) );

            ptr->effective_proton_conductivity(sigmaMeff_cell);
            //std::cout<<"Effective proton conductivity in CL : "<<sigmaMeff_cell<<std::endl;

            // --- Derivative effective property for cell_matrix --------------------
            if ( !cell_residual_counter && (lambda.indices_exist || t_rev.indices_exist) )
            {
                std::vector<VariableNames> deriv_flags;
                if (lambda.indices_exist)
                    deriv_flags.push_back(membrane_water_content);
                if (t_rev.indices_exist)
                    deriv_flags.push_back(temperature_of_REV);

                std::map< VariableNames, std::vector<double> > DsigmaMeff;
                ptr->set_derivative_flags(deriv_flags);
                ptr->derivative_effective_proton_conductivity(DsigmaMeff);

                if (lambda.indices_exist)
                    dsigmaMeff_dlambda_cell = DsigmaMeff[membrane_water_content];

                if (t_rev.indices_exist)
                    dsigmaMeff_dT_cell = DsigmaMeff[temperature_of_REV];
            }
        }

        else if (base_layer == MembraneLayer)
        {
            FuelCellShop::Layer::MembraneLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::MembraneLayer<dim>* >(layer);

            if ( lambda.indices_exist )
                ptr->get_electrolyte()->set_membrane_water_content( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][lambda.solution_index], membrane_water_content) );
            if ( t_rev.indices_exist )
                ptr->get_electrolyte()->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index], temperature_of_REV) );

            ptr->effective_proton_conductivity(sigmaMeff_cell);

            // --- Derivative effective property for cell_matrix --------------------
            if ( !cell_residual_counter && (lambda.indices_exist || t_rev.indices_exist) )
            {
                std::vector<VariableNames> deriv_flags;
                if (lambda.indices_exist)
                    deriv_flags.push_back(membrane_water_content);
                if (t_rev.indices_exist)
                    deriv_flags.push_back(temperature_of_REV);

                std::map< VariableNames, std::vector<double> > DsigmaMeff;
                ptr->set_derivative_flags(deriv_flags);
                ptr->derivative_effective_proton_conductivity(DsigmaMeff);

                if (lambda.indices_exist)
                    dsigmaMeff_dlambda_cell = DsigmaMeff[membrane_water_content];

                if (t_rev.indices_exist)
                    dsigmaMeff_dT_cell = DsigmaMeff[temperature_of_REV];
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
        this->JxW_cell[q] = (cell_info.fe(phi_m.fetype_index)).JxW(q);

        //------ Filling shape functions etc ----------------------------------------------------------------------
        //------ This avoids recalculating shape functions etc for efficiency -------------------------------------
        for (unsigned int k=0; k < (cell_info.fe(phi_m.fetype_index)).dofs_per_cell; ++k)
        {
            grad_phi_phiM_cell[q][k] = (cell_info.fe(phi_m.fetype_index)).shape_grad(k,q);
        }

        if ( !cell_residual_counter ) // This avoids these loops below if we are not assembling for cell matrix. Same can be done for phi_phiM_cell, but there isn't much of an improvement because of extra "if".
        {
            if ( lambda.indices_exist )
                for (unsigned int k=0; k < (cell_info.fe(lambda.fetype_index)).dofs_per_cell; ++k)
                    phi_lambda_cell[q][k] = (cell_info.fe(lambda.fetype_index)).shape_value(k,q);

            if ( t_rev.indices_exist )
                for (unsigned int k=0; k < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++k)
                    phi_T_cell[q][k] = (cell_info.fe(t_rev.fetype_index)).shape_value(k,q);
        }
    }
}

// ---                                    ---
// --- make_assemblers_bdry_variable_data ---
// ---                                    ---

template<int dim>
void
NAME::ProtonTransportEquation<dim>::make_assemblers_bdry_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                                                        FuelCellShop::Layer::BaseLayer<dim>* const layer)
{
    Assert( this->n_q_points_bdry != 0, ExcMessage("make_assemblers_bdry_constant_data function not called before.") );

    //---------------------------------------------------------------------------------------------------------------
    //------------Looping over quadrature points in the cell --------------------------------------------------------
    for (unsigned int q=0; q < this->n_q_points_bdry; ++q)
    {
        this->JxW_bdry[q] = (bdry_info.fe(phi_m.fetype_index)).JxW(q);

        for (unsigned int k=0; k < (bdry_info.fe(phi_m.fetype_index)).dofs_per_cell; ++k)
            phi_phiM_bdry[q][k] = (bdry_info.fe(phi_m.fetype_index)).shape_value(k,q);
    }
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
// CLASS TEST
// Author: M. Secanell
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////


template<int dim>
void
NAME::ProtonTransportEquation<dim>::class_test()
{

    FuelCell::ApplicationCore::BlockInfo block_info;
    Table< 2, DoFTools::Coupling > cell_couplings;
    Table< 2, DoFTools::Coupling > flux_couplings;
    FuelCell::SystemManagement sys(block_info, cell_couplings, flux_couplings);
    FcstUtilities::log<<"Create object ProtonTransportEquation:"<<std::endl;
    FuelCellShop::Equation::ProtonTransportEquation<dim> test(sys);
    test.print_equation_info();

}

// ---                         ---
// --- explicit instantiations ---
// ---                         ---

template class NAME::ProtonTransportEquation<deal_II_dimension>;
