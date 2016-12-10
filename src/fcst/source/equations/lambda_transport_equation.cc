//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: lambda_transport_equation.cc
//    - Description: Implementation of lambda_transport_equation.h. This
//      class is used to setup the matrix and rhs for lambda transport based
//      on Springer model + Thermal Osmosis.
//    - Developers: Madhur Bhaiya, Valentin N. Zingan
//
//---------------------------------------------------------------------------

#include "equations/lambda_transport_equation.h"

namespace NAME = FuelCellShop::Equation;

// ---             ---
// --- Constructors ---
// ---             ---
template<int dim>
NAME::LambdaTransportEquation<dim>::LambdaTransportEquation(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data)
:
EquationBase<dim>(system_management,data)
{
    this->equation_name = "Membrane Water Content Transport Equation";
    this->name_base_variable = "membrane_water_content";
    FcstUtilities::log << "->FuelCellShop::Equation::LambdaTransportEquation" << std::endl;

    //----Initializing VariableInfo Structs------------------------------------------
    //----Setting indices_exist to false --------------------------------------------
    lambda.indices_exist = false;
    phi_m.indices_exist = false;
    t_rev.indices_exist = false;

    //-- Counter to store if generic_data and constant_data need to be initialized.
    this->counter.resize(2, true);
}
// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::LambdaTransportEquation<dim>::~LambdaTransportEquation()
{}

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::LambdaTransportEquation<dim>::declare_parameters(ParameterHandler& param) const
{
    NAME::EquationBase<dim>::declare_parameters(param);
    
    param.enter_subsection("Equations");
    {
        param.enter_subsection(this->equation_name);
        {
            param.enter_subsection("Boolean flags for lambda transport modes");
            {
                param.declare_entry("Thermo-osmosis",
                                    "false",
                                    Patterns::Bool(),
                                    "Flag to include lambda (sorbed water) transport by Thermo-osmosis.");
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
NAME::LambdaTransportEquation<dim>::initialize(ParameterHandler& param)
{
    NAME::EquationBase<dim>::initialize(param);
    
    param.enter_subsection("Equations");
    {
        param.enter_subsection(this->equation_name);
        {
            param.enter_subsection("Boolean flags for lambda transport modes");
            {
                flag_thermoosmosis = param.get_bool("Thermo-osmosis");
            }
            param.leave_subsection();
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
NAME::LambdaTransportEquation<dim>::assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
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
        for (unsigned int i=0; i < (cell_info.fe(lambda.fetype_index)).dofs_per_cell; ++i)
        {
            //--------------LOOP(s) over j-------------------------------------------------------------

            //-----------Assembling Matrix for terms corresponding to "lambda" BLOCK------------------------
            for (unsigned int j=0; j < (cell_info.fe(lambda.fetype_index)).dofs_per_cell; ++j)
            {
                Tensor<1,dim> dF_L;
                dF_L = ( (nDrag_cell[q]/F)*cell_info.gradients[last_iter_cell][phi_m.solution_index][q]*dsigmaMeff_dlambda_cell[q] + (rho_dry_cell/EW_cell)*cell_info.gradients[last_iter_cell][lambda.solution_index][q]*dDlambdaeff_dlambda_cell[q] )*phi_lambda_cell[q][j]
                       +
                       (rho_dry_cell/EW_cell)*Dlambdaeff_cell[q]*grad_phi_lambda_cell[q][j];

                cell_matrices[lambda.block_index].matrix(i,j) += grad_phi_lambda_cell[q][i] * dF_L * this->JxW_cell[q];
            }

            //----------TERM CORRESPONDING TO "phi_m" BLOCK-----------------------------------------
            if ( phi_m.indices_exist )
            {
                //-----Assembling Matrix------------------------------------------------------------
                for (unsigned int j=0; j < (cell_info.fe(phi_m.fetype_index)).dofs_per_cell; ++j)
                {
                    Tensor<1,dim> dF_Fm;
                    dF_Fm = (nDrag_cell[q]/F) * sigmaMeff_cell[q] * grad_phi_phiM_cell[q][j];

                    cell_matrices[phi_m.block_index].matrix(i,j) += grad_phi_lambda_cell[q][i] * dF_Fm * this->JxW_cell[q];
                }
            }

            //----------TERM CORRESPONDING TO "T" BLOCK-----------------------------------------
            if ( t_rev.indices_exist )
            {
                //-----Assembling Matrix------------------------------------------------------------
                for (unsigned int j=0; j < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++j)
                {
                    Tensor<1,dim> dF_T;
                    dF_T = ( (nDrag_cell[q]/F)*cell_info.gradients[last_iter_cell][phi_m.solution_index][q]*dsigmaMeff_dT_cell[q]+
                             (rho_dry_cell/EW_cell)*cell_info.gradients[last_iter_cell][lambda.solution_index][q]*dDlambdaeff_dT_cell[q]+
                             (1.0/M_water)*cell_info.gradients[last_iter_cell][t_rev.solution_index][q]*dDTeff_dT_cell[q] )*phi_T_cell[q][j]+
                             (1.0/M_water)*DTeff_cell[q]*grad_phi_T_cell[q][j];

                    cell_matrices[t_rev.block_index].matrix(i,j) += grad_phi_lambda_cell[q][i] * dF_T * this->JxW_cell[q];
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
NAME::LambdaTransportEquation<dim>::assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
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
        for (unsigned int i=0; i < (cell_info.fe(lambda.fetype_index)).dofs_per_cell; ++i)
        {
            Tensor<1,dim> F_L_old;
            F_L_old = (rho_dry_cell/EW_cell)*Dlambdaeff_cell[q]*cell_info.gradients[last_iter_cell][lambda.solution_index][q];

            cell_residual.block(lambda.solution_index)(i) += grad_phi_lambda_cell[q][i] * F_L_old * this->JxW_cell[q];
        }

        if ( phi_m.indices_exist )
        {
               for (unsigned int i=0; i < (cell_info.fe(lambda.fetype_index)).dofs_per_cell; ++i)
               {
                   Tensor<1,dim> F_Fm_old;
                   F_Fm_old = (nDrag_cell[q]/F)*sigmaMeff_cell[q]*cell_info.gradients[last_iter_cell][phi_m.solution_index][q];

                   cell_residual.block(lambda.solution_index)(i) += grad_phi_lambda_cell[q][i] * F_Fm_old * this->JxW_cell[q];
               }
        }

        if ( t_rev.indices_exist )
        {
               for (unsigned int i=0; i < (cell_info.fe(lambda.fetype_index)).dofs_per_cell; ++i)
               {
                   Tensor<1,dim> F_T_old;
                   F_T_old = (1.0/M_water)*DTeff_cell[q]*cell_info.gradients[last_iter_cell][t_rev.solution_index][q];

                   cell_residual.block(lambda.solution_index)(i) += grad_phi_lambda_cell[q][i] * F_T_old * this->JxW_cell[q];
               }
        }
    }
}

// ---                              ---
// --- make_internal_cell_couplings ---
// ---                              ---
template<int dim>
void
NAME::LambdaTransportEquation<dim>::make_internal_cell_couplings()
{
    AssertThrow(this->system_management->solution_in_userlist(this->name_base_variable), VariableShouldExistForEquation(this->name_base_variable, this->equation_name) );
    AssertThrow(this->system_management->equation_name_to_index(this->equation_name) == this->system_management->solution_name_to_index(this->name_base_variable),
              IndexDoNotMatch(this->name_base_variable, this->equation_name) );

    if (flag_thermoosmosis)
        AssertThrow(this->system_management->solution_in_userlist("temperature_of_REV"), ExcMessage("temperature_of_solid_phase should be defined in the user list for thermo-osmosis.") );

    std::map< std::string, DoFTools::Coupling > tmp;

    std::vector< std::string> sol_names = this->system_management->get_solution_names();

    for (unsigned int i = 0; i < sol_names.size(); ++i)
    {
        if (sol_names[i] == this->name_base_variable)
            tmp[this->name_base_variable] = DoFTools::always;

        else if (sol_names[i] == "temperature_of_REV")
            tmp["temperature_of_REV"] = DoFTools::always;

        else if (sol_names[i] == "protonic_electrical_potential")
            tmp["protonic_electrical_potential"] = DoFTools::always;

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
NAME::LambdaTransportEquation<dim>::print_equation_info() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "-------------------------------------------------------------------------------" << std::endl;
    FcstUtilities::log << std::endl;

    FcstUtilities::log << "PARAMETERS FOR \" Membrane Water Content Transport Equation \":" << std::endl;
    FcstUtilities::log << std::endl;

    FcstUtilities::log << "Boolean flags for lambda transport modes:" << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "Lambda transport via Thermo-osmosis:  " << flag_thermoosmosis << std::endl;
    FcstUtilities::log << std::endl;

    FcstUtilities::log << "INTERNAL CELL COUPLINGS FOR \" Membrane Water Content Transport Equation \":" << std::endl;
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
NAME::LambdaTransportEquation<dim>::make_assemblers_generic_constant_data()
{
    FuelCellShop::Material::WaterVapor water;

    F = Constants::F();
    M_water = water.get_molar_mass() * 1.0e3;

    //-----------Filling VariableInfo structures------------------------------------------
    //----------membrane_water_content----------------------------------------------------
    lambda.solution_index = this->system_management->solution_name_to_index(this->name_base_variable);
    lambda.block_index = this->system_management->matrix_block_index(this->equation_name, this->name_base_variable);
    lambda.fetype_index = this->system_management->block_info->base_element[lambda.solution_index];
    lambda.indices_exist = true;

    //-----------protonic_electrical_potential------------------------------------------------
    if ( this->system_management->solution_in_userlist("protonic_electrical_potential") )
    {
        phi_m.solution_index = this->system_management->solution_name_to_index("protonic_electrical_potential");
        phi_m.block_index = this->system_management->matrix_block_index(this->equation_name, "protonic_electrical_potential");
        phi_m.fetype_index = this->system_management->block_info->base_element[phi_m.solution_index];
        phi_m.indices_exist = true;
    }

    //----------temperature_of_solid_phase----------------------------------------------------
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
NAME::LambdaTransportEquation<dim>::make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
    Assert( lambda.indices_exist, ExcMessage("make_assemblers_generic_constant_data function not called before.") );

    this->n_q_points_cell = (cell_info.fe(lambda.fetype_index)).n_quadrature_points;
    last_iter_cell = cell_info.global_data->find_vector(this->solution_vector_name);

    //-------------Allocation------------------------------------------
    phi_lambda_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(lambda.fetype_index)).dofs_per_cell ) );
    grad_phi_lambda_cell.resize( this->n_q_points_cell, std::vector< Tensor<1,dim> >( (cell_info.fe(lambda.fetype_index)).dofs_per_cell ) );

    if ( phi_m.indices_exist )
        grad_phi_phiM_cell.resize( this->n_q_points_cell, std::vector< Tensor<1,dim> >( (cell_info.fe(phi_m.fetype_index)).dofs_per_cell ) );

    if ( t_rev.indices_exist )
    {
        phi_T_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(t_rev.fetype_index)).dofs_per_cell ) );
        grad_phi_T_cell.resize( this->n_q_points_cell, std::vector< Tensor<1,dim> >( (cell_info.fe(t_rev.fetype_index)).dofs_per_cell ) );
    }

    //-----------------------------------------------------------------
    this->JxW_cell.resize(this->n_q_points_cell);

    Dlambdaeff_cell.resize(this->n_q_points_cell);
    dDlambdaeff_dlambda_cell.resize(this->n_q_points_cell);
    dDlambdaeff_dT_cell.resize(this->n_q_points_cell);

    sigmaMeff_cell.resize(this->n_q_points_cell);
    dsigmaMeff_dlambda_cell.resize(this->n_q_points_cell);
    dsigmaMeff_dT_cell.resize(this->n_q_points_cell);

    nDrag_cell.resize(this->n_q_points_cell);
    dnDrag_dlambda_cell.resize(this->n_q_points_cell);

    DTeff_cell.resize(this->n_q_points_cell);
    dDTeff_dT_cell.resize(this->n_q_points_cell);
}

// ---                                    ---
// --- make_assemblers_cell_variable_data ---
// ---                                    ---

template<int dim>
void
NAME::LambdaTransportEquation<dim>::make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
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
            //-------------------
            rho_dry_cell = ptr->get_electrolyte()->get_density();
            EW_cell = ptr->get_electrolyte()->get_EW();
            //-------------------
            ptr->get_electrolyte()->set_membrane_water_content( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][lambda.solution_index],
                                                                membrane_water_content) );
            if (t_rev.indices_exist)
            {
                ptr->get_electrolyte()->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index],
                                                         temperature_of_REV) );
            }
            //-------------------
            ptr->effective_water_diffusivity(Dlambdaeff_cell);
            if (phi_m.indices_exist)
            {
                ptr->effective_proton_conductivity(sigmaMeff_cell);
                ptr->get_electrolyte()->electroosmotic_drag(nDrag_cell);
            }
            if (flag_thermoosmosis)
                ptr->effective_thermoosmotic_diffusivity(DTeff_cell);

            // --- Derivative effective property for cell_matrix --------------------
            if ( !cell_residual_counter )
            {
                std::vector<VariableNames> deriv_flags;
                deriv_flags.push_back(membrane_water_content);
                if ( t_rev.indices_exist )
                    deriv_flags.push_back(temperature_of_REV);
                ptr->set_derivative_flags(deriv_flags);
                //-------------------
                std::map< VariableNames, std::vector<double> > dDlambdaeff;
                ptr->derivative_effective_water_diffusivity(dDlambdaeff);
                dDlambdaeff_dlambda_cell = dDlambdaeff[membrane_water_content];
                if (t_rev.indices_exist)
                    dDlambdaeff_dT_cell = dDlambdaeff[temperature_of_REV];
                //-------------------
                if (phi_m.indices_exist)
                {
                    std::map< VariableNames, std::vector<double> > DsigmaMeff;
                    std::map< VariableNames, std::vector<double> > dnDrag;
                    ptr->get_electrolyte()->electroosmotic_drag_derivative(dnDrag);
                    dnDrag_dlambda_cell = dnDrag[membrane_water_content];
                    ptr->derivative_effective_proton_conductivity(DsigmaMeff);
                    dsigmaMeff_dlambda_cell = DsigmaMeff[membrane_water_content];
                    if (t_rev.indices_exist)
                        dsigmaMeff_dT_cell = DsigmaMeff[temperature_of_REV];
                }
                //-------------------
                if (flag_thermoosmosis)
                {
                    std::map< VariableNames, std::vector<double> > dDTeff;
                    ptr->derivative_effective_thermoosmotic_diffusivity(dDTeff);
                    dDTeff_dT_cell = dDTeff[temperature_of_REV];
                }
            }
        }

        else if (base_layer == MembraneLayer)
        {
            FuelCellShop::Layer::MembraneLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::MembraneLayer<dim>* >(layer);
            //-------------------
            rho_dry_cell = ptr->get_electrolyte()->get_density();
            EW_cell = ptr->get_electrolyte()->get_EW();
            //-------------------
            ptr->get_electrolyte()->set_membrane_water_content( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][lambda.solution_index],
                                                                membrane_water_content) );
            if (t_rev.indices_exist)
            {
                ptr->get_electrolyte()->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index],
                                                         temperature_of_REV) );
            }
            //-------------------
            ptr->effective_water_diffusivity(Dlambdaeff_cell);
            if (phi_m.indices_exist)
            {
                ptr->effective_proton_conductivity(sigmaMeff_cell);
                ptr->get_electrolyte()->electroosmotic_drag(nDrag_cell);
            }
            if (flag_thermoosmosis)
                ptr->effective_thermoosmotic_diffusivity(DTeff_cell);

            // --- Derivative effective property for cell_matrix --------------------
            if ( !cell_residual_counter )
            {
                std::vector<VariableNames> deriv_flags;
                deriv_flags.push_back(membrane_water_content);
                if (t_rev.indices_exist)
                    deriv_flags.push_back(temperature_of_REV);
                ptr->set_derivative_flags(deriv_flags);
                //-------------------
                std::map< VariableNames, std::vector<double> > dDlambdaeff;
                ptr->derivative_effective_water_diffusivity(dDlambdaeff);
                dDlambdaeff_dlambda_cell = dDlambdaeff[membrane_water_content];
                if (t_rev.indices_exist)
                    dDlambdaeff_dT_cell = dDlambdaeff[temperature_of_REV];
                //-------------------
                if (phi_m.indices_exist)
                {
                    std::map< VariableNames, std::vector<double> > DsigmaMeff;
                    std::map< VariableNames, std::vector<double> > dnDrag;
                    ptr->get_electrolyte()->electroosmotic_drag_derivative(dnDrag);
                    dnDrag_dlambda_cell = dnDrag[membrane_water_content];
                    ptr->derivative_effective_proton_conductivity(DsigmaMeff);
                    dsigmaMeff_dlambda_cell = DsigmaMeff[membrane_water_content];
                    if (t_rev.indices_exist)
                        dsigmaMeff_dT_cell = DsigmaMeff[temperature_of_REV];
                }
                //-------------------
                if (flag_thermoosmosis)
                {
                    std::map< VariableNames, std::vector<double> > dDTeff;
                    ptr->derivative_effective_thermoosmotic_diffusivity(dDTeff);
                    dDTeff_dT_cell = dDTeff[temperature_of_REV];
                }
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
        this->JxW_cell[q] = (cell_info.fe(lambda.fetype_index)).JxW(q);

        //------ Filling shape functions etc ----------------------------------------------------------------------
        //------ This avoids recalculating shape functions etc for efficiency -------------------------------------
        for (unsigned int k=0; k < (cell_info.fe(lambda.fetype_index)).dofs_per_cell; ++k)
        {
            phi_lambda_cell[q][k] = (cell_info.fe(lambda.fetype_index)).shape_value(k,q);
            grad_phi_lambda_cell[q][k] = (cell_info.fe(lambda.fetype_index)).shape_grad(k,q);
        }

        if ( phi_m.indices_exist )
            for (unsigned int k=0; k < (cell_info.fe(phi_m.fetype_index)).dofs_per_cell; ++k)
                grad_phi_phiM_cell[q][k] = (cell_info.fe(phi_m.fetype_index)).shape_grad(k,q);

        if ( t_rev.indices_exist)
            for (unsigned int k=0; k < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++k)
            {
                phi_T_cell[q][k] = (cell_info.fe(t_rev.fetype_index)).shape_value(k,q);
                grad_phi_T_cell[q][k] = (cell_info.fe(t_rev.fetype_index)).shape_grad(k,q);
            }
    }
}

// ---                         ---
// --- explicit instantiations ---
// ---                         ---

template class NAME::LambdaTransportEquation<deal_II_dimension>;
