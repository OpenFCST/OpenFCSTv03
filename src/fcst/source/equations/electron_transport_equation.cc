// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: electron_transport_equation.cc
// - Description: This class describes steady-state electron transport using Ohm's law
// - Developers: Marc Secanell, Madhur Bhaiya, Valentin N. Zingan
//
// ----------------------------------------------------------------------------

#include "equations/electron_transport_equation.h"

namespace NAME = FuelCellShop::Equation;

// ---             ---
// --- Constructor ---
// ---             ---
template<int dim>
NAME::ElectronTransportEquation<dim>::ElectronTransportEquation(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data)
:
EquationBase<dim>(system_management,data)
{
    this->name_base_variable = "electronic_electrical_potential";
    this->equation_name = "Electron Transport Equation";
    FcstUtilities::log << "->FuelCellShop::Equation::ElectronTransportEquation" << std::endl;

    //----Initializing VariableInfo Structs------------------------------------------
    //----Setting indices_exist to false --------------------------------------------
    phi_s.indices_exist = false;

    //-- Counter to store if generic_data and constant_data need to be initialized.
    this->counter.resize(3, true);
}
// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::ElectronTransportEquation<dim>::~ElectronTransportEquation()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::ElectronTransportEquation<dim>::declare_parameters(ParameterHandler& param) const
{
    NAME::EquationBase<dim>::declare_parameters(param);
    
    param.enter_subsection("Equations");
    {
        param.enter_subsection(this->equation_name);
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
                param.declare_entry("Constant Electron Current Flux Boundary Conditions",
                                "",
                                Patterns::Map( Patterns::Integer(0), Patterns::Double() ),
                                "A comma-separated list of boundaries with Galvanostatic b.c., with prescribed values of electronic current fluxes [A/cm^2]."
                                "\n"
                                "Correct format is of a map, given as ''id1: value1, id2: value2, id3: value3, ...'' where "
                                "each id must be an unsigned integer, and each value can be a double (positive for electronic current flux leaving out and negative for current coming in).");
            }
            param.leave_subsection();

            ///////////////////
            // TO BE REMOVED //
            ///////////////////
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
NAME::ElectronTransportEquation<dim>::initialize(ParameterHandler& param)
{
    
    NAME::EquationBase<dim>::initialize(param);
    param.enter_subsection("Equations");
    {
        param.enter_subsection(this->equation_name);
        {
            param.enter_subsection("Boundary data");
            {
                if( !param.get("current_flux").empty() )
                {
                    // electron_current_flux_map = FcstUtilities::string_to_map<unsigned int, double>( param.get("current_flux") ); // will be uncommented later on
                }
            }
            param.leave_subsection();
            
            param.enter_subsection("Boundary conditions");
            {
                electron_current_flux_map = FcstUtilities::string_to_map<unsigned int, double>( Utilities::split_string_list( param.get("Constant Electron Current Flux Boundary Conditions")) );
            }
            param.leave_subsection();

            ///////////////////
            // TO BE REMOVED //
            ///////////////////
        }
        param.leave_subsection();
    }
    param.leave_subsection();

    // AssertThrow ( electron_current_flux_map.size() <= 1, ExcMessage("Currently, electronic current flux boundary condition can be applied atmost on one boundary only.") ); that's not true
    this->make_internal_cell_couplings();
}

// ---                           ---
// --- set_parameters ---
// ---                            ---
template<int dim>
void
NAME::ElectronTransportEquation<dim>::set_parameters(const std::vector<std::string>& name_dvar,
                                                     const std::vector<double>& value_dvar,
                                                     ParameterHandler& param)
{
    bool current_bdry_used(false);

    param.enter_subsection("Equations");
    {
        param.enter_subsection("Electron Transport Equation");
        {
            param.enter_subsection("Boundary conditions");
            {
                for (unsigned int i=0; i<name_dvar.size(); ++i)
                {
                    std::string temp = name_dvar.at(i);

                    if ( temp.find("current_bdryid_") != std::string::npos )
                    {
                        AssertThrow (!current_bdry_used, ExcMessage("current_bdryid_ can only be specified once in ElectronTransportEquation::set_parameters method."));
                        temp.erase(0,15);
                        unsigned int current_bdry_id = FcstUtilities::string_to_number<unsigned int>(temp);

                        std::stringstream s;
                        s << current_bdry_id << ":" << value_dvar.at(i);

                        param.set("Constant Electron Current Flux Boundary Conditions", s.str());

                        current_bdry_used = true;
                    }
                }
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

// ---                      ---
// --- assemble_cell_matrix ---
// ---                      ---

template<int dim>
void
NAME::ElectronTransportEquation<dim>::assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
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

    this->make_assemblers_cell_variable_data(cell_info, layer);


    //-------- Looping over Quadrature points ----------------------------
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
        //---------------LOOP over i -----------------------------------------------------------------
        for (unsigned int i=0; i < (cell_info.fe(phi_s.fetype_index)).dofs_per_cell; ++i)
        {
            //--------------LOOP(s) over j-------------------------------------------------------------

            //-----------Assembling Matrix for terms corresponding to "phi_s" BLOCK------------------------
            for (unsigned int j=0; j < (cell_info.fe(phi_s.fetype_index)).dofs_per_cell; ++j)
            {
                cell_matrices[phi_s.block_index].matrix(i,j) += grad_phi_phiS_cell[q][i] * sigmaSeff_cell * grad_phi_phiS_cell[q][j] * this->JxW_cell[q];
            }
        }
    }
}

// ---                        ---
// --- assemble_cell_residual ---
// ---                        ---
template<int dim>
void
NAME::ElectronTransportEquation<dim>::assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
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

    this->make_assemblers_cell_variable_data(cell_info, layer);


    for (unsigned int q=0; q < this->n_q_points_cell; ++q)
    {
        for (unsigned int i=0; i < (cell_info.fe(phi_s.fetype_index)).dofs_per_cell; ++i)
        {
            cell_residual.block(phi_s.solution_index)(i) += grad_phi_phiS_cell[q][i] * sigmaSeff_cell * cell_info.gradients[last_iter_cell][phi_s.solution_index][q] * this->JxW_cell[q];
        }
    }

}


// ---                      ---
// --- assemble_bdry_matrix ---
// ---                      ---
template<int dim>
void
NAME::ElectronTransportEquation<dim>::assemble_bdry_matrix(FuelCell::ApplicationCore::MatrixVector&                                 bdry_matrices,
                                                           const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                           FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
       // The boundary integral " q * [sigma_s * grad_delta_phi_s * n] " is always zero because of the boundary conditions we use:
       // - if phi_s = known then q = 0
       // - if [sigma_s * grad_phi_s * n] = known then [sigma_s * grad_delta_phi_s * n] = 0
}

// ---                        ---
// --- assemble_bdry_residual ---
// ---                        ---

template<int dim>
void
NAME::ElectronTransportEquation<dim>::assemble_bdry_residual(FuelCell::ApplicationCore::FEVector&                                     bdry_residual,
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

    //------ Constant Electron Current Flux boundaries --------------------------------------------------------------
    std::map<unsigned int, double>::const_iterator iter = electron_current_flux_map.find( bdry_info.dof_face->boundary_indicator() );
    if (iter != electron_current_flux_map.end() )
    {
        const double& electron_current_value = iter->second;

        this->make_assemblers_bdry_variable_data(bdry_info, layer);

        //-------- Looping over Quadrature points ----------------------------
        for (unsigned int q = 0; q < this->n_q_points_bdry; ++q)
        {
            for (unsigned int i = 0; i < (bdry_info.fe(phi_s.fetype_index)).dofs_per_cell; ++i)
            {
                bdry_residual.block(phi_s.solution_index)(i) += - phi_phiS_bdry[q][i] * electron_current_value * this->JxW_bdry[q];
            }
        }
    }
}

// ---                              ---
// --- make_internal_cell_couplings ---
// ---                              ---
template<int dim>
void
NAME::ElectronTransportEquation<dim>::make_internal_cell_couplings()
{
    AssertThrow(this->system_management->solution_in_userlist(this->name_base_variable), VariableShouldExistForEquation(this->name_base_variable, this->equation_name) );
    AssertThrow(this->system_management->equation_name_to_index(this->equation_name) == this->system_management->solution_name_to_index(this->name_base_variable),
              IndexDoNotMatch(this->name_base_variable, this->equation_name) );

    std::map< std::string, DoFTools::Coupling > tmp;

    std::vector< std::string> sol_names = this->system_management->get_solution_names();

    for (unsigned int i=0; i < sol_names.size(); ++i)
    {
        if (sol_names[i] == this->name_base_variable)
            tmp[this->name_base_variable] = DoFTools::always;

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
NAME::ElectronTransportEquation<dim>::print_equation_info() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "-------------------------------------------------------------------------------" << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "INTERNAL CELL COUPLINGS FOR \" Electron Transport Equation \":" << std::endl;
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
NAME::ElectronTransportEquation<dim>::make_assemblers_generic_constant_data()
{
    //-----------Filling VariableInfo structures------------------------------------------
    phi_s.solution_index = this->system_management->solution_name_to_index("electronic_electrical_potential");
    phi_s.block_index = this->system_management->matrix_block_index(this->equation_name, "electronic_electrical_potential");
    phi_s.fetype_index = this->system_management->block_info->base_element[phi_s.solution_index];
    phi_s.indices_exist = true;
}

// ---                                    ---
// --- make_assemblers_cell_constant_data ---
// ---                                    ---

template<int dim>
void
NAME::ElectronTransportEquation<dim>::make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
    Assert( phi_s.indices_exist, ExcMessage("make_assemblers_generic_constant_data function not called before.") );

    this->n_q_points_cell = (cell_info.fe(phi_s.fetype_index)).n_quadrature_points;
    last_iter_cell = cell_info.global_data->find_vector(this->solution_vector_name);

    //-------------Allocation------------------------------------------
    grad_phi_phiS_cell.resize( this->n_q_points_cell, std::vector< Tensor<1,dim> >( (cell_info.fe(phi_s.fetype_index)).dofs_per_cell ) );

    //-----------------------------------------------------------------
    this->JxW_cell.resize(this->n_q_points_cell);
}

// ---                                    ---
// --- make_assemblers_bdry_constant_data ---
// ---                                    ---

template<int dim>
void
NAME::ElectronTransportEquation<dim>::make_assemblers_bdry_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info)
{
    Assert( phi_s.indices_exist, ExcMessage("make_assemblers_generic_constant_data function not called before.") );

    this->n_q_points_bdry = (bdry_info.fe(phi_s.fetype_index)).n_quadrature_points;
    last_iter_bdry = bdry_info.global_data->find_vector(this->solution_vector_name);

    //-------------Allocation------------------------------------------
    phi_phiS_bdry.resize( this->n_q_points_bdry, std::vector<double>( (bdry_info.fe(phi_s.fetype_index)).dofs_per_cell ) );

    //-----------------------------------------------------------------
    this->JxW_bdry.resize(this->n_q_points_bdry);
}

// ---                                    ---
// --- make_assemblers_cell_variable_data ---
// ---                                    ---

template<int dim>
void
NAME::ElectronTransportEquation<dim>::make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                         FuelCellShop::Layer::BaseLayer<dim>* const layer)
{
    Assert( this->n_q_points_cell != 0, ExcMessage("make_assemblers_cell_constant_data function not called before.") );

    //---------------Effective Transport Properties---------------------------------------------------------------
    // ----- type infos -------------
    const std::type_info& Channel           = typeid(FuelCellShop::Layer::Channel<dim>);
    const std::type_info& GasDiffusionLayer = typeid(FuelCellShop::Layer::GasDiffusionLayer<dim>);
    const std::type_info& MicroPorousLayer  = typeid(FuelCellShop::Layer::MicroPorousLayer<dim>);
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
    const std::type_info& SolidLayer = typeid(FuelCellShop::Layer::SolidLayer<dim>);

    const std::type_info& base_layer = layer->get_base_type();

    // ----- Clearing and initializing the containers -----------------------------------------
    // ----- All containers intialized to zero by default -------------------------------------

    sigmaSeff_cell.clear(); // This automatically initializes it to zero Tensor.

    // ----- dynamic cast and filling the containers -----------------
    try
    {
        if( base_layer == Channel )
        {
            FuelCellShop::Layer::Channel<dim>* ptr = dynamic_cast< FuelCellShop::Layer::Channel<dim>* >(layer);
            sigmaSeff_cell  = ptr->get_effective_electronic_conductivity();
        }

        else if (base_layer == GasDiffusionLayer)
        {
            FuelCellShop::Layer::GasDiffusionLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::GasDiffusionLayer<dim>* >(layer);
            ptr->effective_electron_conductivity(sigmaSeff_cell);
        }

        else if (base_layer == MicroPorousLayer)
        {
            FuelCellShop::Layer::MicroPorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::MicroPorousLayer<dim>* >(layer);
            ptr->effective_electron_conductivity(sigmaSeff_cell);
            //std::cout<<"Effective electron conductivity in the MPL is : "<<sigmaSeff_cell<<std::endl;
        }

        else if (base_layer == CatalystLayer)
        {
            FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);
            ptr->effective_electron_conductivity(sigmaSeff_cell);
            //std::cout<<"Effective electron conductivity in the CL is : "<<sigmaSeff_cell<<std::endl;
        }

        else if (base_layer == SolidLayer)
        {
            FuelCellShop::Layer::SolidLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::SolidLayer<dim>* >(layer);
            ptr->effective_electron_conductivity(sigmaSeff_cell);
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
        this->JxW_cell[q] = (cell_info.fe(phi_s.fetype_index)).JxW(q);

        //------ Filling shape functions etc ----------------------------------------------------------------------
        //------ This avoids recalculating shape functions etc for efficiency -------------------------------------
        for (unsigned int k=0; k < (cell_info.fe(phi_s.fetype_index)).dofs_per_cell; ++k)
        {
            grad_phi_phiS_cell[q][k] = (cell_info.fe(phi_s.fetype_index)).shape_grad(k,q);
        }
    }
}

// ---                                    ---
// --- make_assemblers_bdry_variable_data ---
// ---                                    ---

template<int dim>
void
NAME::ElectronTransportEquation<dim>::make_assemblers_bdry_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                            FuelCellShop::Layer::BaseLayer<dim>* const layer)
{
    Assert( this->n_q_points_bdry != 0, ExcMessage("make_assemblers_bdry_constant_data function not called before.") );
    //---------------------------------------------------------------------------------------------------------------
    //------------Looping over quadrature points in the cell --------------------------------------------------------
    for (unsigned int q=0; q < this->n_q_points_bdry; ++q)
    {
        //-------JxW & normal_vectors----------------------------------------------------------------------------
        this->JxW_bdry[q] = (bdry_info.fe(phi_s.fetype_index)).JxW(q);

        //------ Filling shape functions etc --------------------------------------------------------------------
        //------ This avoids recalculating shape functions etc for efficiency -----------------------------------

        for (unsigned int k=0; k < (bdry_info.fe(phi_s.fetype_index)).dofs_per_cell; ++k)
        {
            phi_phiS_bdry[q][k] = (bdry_info.fe(phi_s.fetype_index)).shape_value(k,q);
        }
    }


}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
// CLASS TEST
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////


template<int dim>
void
NAME::ElectronTransportEquation<dim>::class_test()
{

    FuelCell::ApplicationCore::BlockInfo block_info;
    Table< 2, DoFTools::Coupling > cell_couplings;
    Table< 2, DoFTools::Coupling > flux_couplings;
    FuelCell::SystemManagement sys(block_info, cell_couplings, flux_couplings);
    FcstUtilities::log<<"Create object ElectronTransportEquation:"<<std::endl;
    FuelCellShop::Equation::ElectronTransportEquation<dim> test(sys);
    test.print_equation_info();

}

// ---                         ---
// --- explicit instantiations ---
// ---                         ---

template class NAME::ElectronTransportEquation<deal_II_dimension>;