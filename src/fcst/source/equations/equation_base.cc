// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: equation_base.cc
// - Description: This is a base class for all available openFCST equations
// - Developers: Valentin N. Zingan,    University of Alberta
//               Marc Secanell Gallart, University of Alberta
//               Aslan Kosakian,        University of Alberta
//
// ----------------------------------------------------------------------------

#include "equations/equation_base.h"

namespace NAME = FuelCellShop::Equation;

//////////////////////////////////////////////////
// CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
//////////////////////////////////////////////////
// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::EquationBase<dim>::EquationBase(FuelCell::SystemManagement& sys_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data)
:
system_management(&sys_management),
data(data),
solution_vector_name(this->data->get_solution_vector_name(this->data->get_nonlinear_solver())),
residual_vector_name(this->data->get_residual_vector_name(this->data->get_nonlinear_solver()))
{ 
    assemble_flags.assemble_generic_data = true;
    assemble_flags.assemble_cell_constant_data = true;
    assemble_flags.assemble_cell_variable_data_matrix = true;
}

// ---            ---
// --- Destructor ---
// ---            ---
template<int dim>
NAME::EquationBase<dim>::~EquationBase()
{ }

////////////////////////////
// GENERIC INITIALIZATION //
////////////////////////////
template<int dim>
void
NAME::EquationBase<dim>::declare_parameters(ParameterHandler& param) const
{
param.enter_subsection("Equations");
    {
        param.enter_subsection(this->equation_name);
        {
            //--
            param.enter_subsection("Initial data");
            {
                   param.declare_entry("Variable initial data",
                                       "false",
                                        Patterns::Bool(),
                                       "True, if initial data is NOT constant at least in one sub-domain. False, if initial data is constant or piece-wise constant.");

                   param.declare_entry( name_base_variable,
                                       """",
                                        Patterns::Map(   Patterns::Integer(0) , Patterns::Double()   ),
                                       "This information is used to setup an initial solution. Note this can be overwritten by Operating Conditions. \n"
                                       "Enter the mole fraction that you would like to use as initial value for each material in your mesh. \n"
                                       "The format should be as follows: material_id1:value1, material_id2:value2. For example, using "
                                       "2:0.1, 4:0.5 will set the initial mole fraction in material 2 to 0.1 and in material 4 to 0.5");
            }
            param.leave_subsection();
            //--
            param.enter_subsection("Boundary data");
            {
                   param.declare_entry("Variable boundary data",
                                       "false",
                                        Patterns::Bool(),
                                       "True, if Dirichlet boundary data is NOT constant at least at one sub-boundary. False, if Dirichlet boundary data is constant or piece-wise constant.");

                   param.declare_entry( name_base_variable,
                                       """",
                                        Patterns::Map(   Patterns::Integer(0) , Patterns::Double()   ),
                                       "This information is used to setup an initial solution. Note this can be overwritten by Operating Conditions.\n"
                                       "Enter the value that you would like to use as the Dirichlet boundary condition for each"
                                       "boundary_id. The format should be as follows: boundary_id1:value1, boundary_id2:value2. \n"
                                       "For example, using 3:0.1, 41:0.25 will set the value in boundary 3 to 0.1 and in boundary 41 to 0.25");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//--------------------------------------------------------------------
template<int dim>
void
NAME::EquationBase<dim>::initialize(ParameterHandler& param)
{
    param.enter_subsection("Equations");
    {
        param.enter_subsection(this->equation_name);
        {
            param.enter_subsection("Initial data");
            {
                   variable_initial_data = param.get_bool("Variable initial data");

                   if( !param.get(name_base_variable).empty() )
                   {
                          const std::map<types::material_id, double> tmp = FcstUtilities::string_to_map<types::material_id, double>( param.get(name_base_variable) );
                          this->component_materialID_value[name_base_variable] = tmp;
                   }
            }
            param.leave_subsection();

            param.enter_subsection("Boundary data");
            {
                   variable_boundary_data = param.get_bool("Variable boundary data");

                   if( !param.get(name_base_variable).empty() )
                   {
                          const std::map<types::boundary_id, double> tmp = FcstUtilities::string_to_map<types::boundary_id, double>( param.get(name_base_variable) );
                          this->component_boundaryID_value[name_base_variable] = tmp;
                   }
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();

}

//--------------------------------------------------------------------
template<int dim>
void
NAME::EquationBase<dim>::assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                              const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                              FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
    if ((this->data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NEWTONBASIC) || 
        (this->data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NEWTON3PP) || 
        (this->data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NEWTONLINESEARCH))
    {
        this->assemble_cell_Jacobian_matrix(cell_matrices, cell_info, layer);
    }
    else if ((this->data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::PICARD) || 
        (this->data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NONE))
    {
        this->assemble_cell_linear_matrix(cell_matrices, cell_info, layer);
    }
    else
    {
        FcstUtilities::log<<"Option not implemented"<<std::endl;
        exit(-1);
    }
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
template<int dim>
void
NAME::EquationBase<dim>::assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
    if ((this->data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NEWTONBASIC) || 
        (this->data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NEWTON3PP) || 
        (this->data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NEWTONLINESEARCH))
    {
        this->assemble_cell_residual_rhs(cell_residual, cell_info, layer);
    }
    else if ((this->data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::PICARD) || 
        (this->data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NONE))
    {
        this->assemble_cell_linear_rhs(cell_residual, cell_info, layer);
    }
    else
    {
        FcstUtilities::log<<"Option not implemented"<<std::endl;
        exit(-1);
    }
}    

////////////////
// CONVERTERS //
////////////////
// ---                        ---
// --- standard_to_block_wise ---
// ---                        ---
template<int dim>
void
NAME::EquationBase<dim>::standard_to_block_wise(FullMatrix<double>& target) const
{
    const std::vector<unsigned int> local_renumbering( this->system_management->block_info->local_renumbering );

    FullMatrix<double> tmp(target.m(), target.n());

    for(unsigned int i = 0; i < target.m(); ++i)
        for(unsigned int j = 0; j < target.n(); ++j)
            tmp(local_renumbering[i], local_renumbering[j]) = target(i,j);

    target = tmp;
}

// ---                        ---
// --- standard_to_block_wise ---
// ---                        ---
template<int dim>
void
NAME::EquationBase<dim>::standard_to_block_wise(Vector<double>& target) const
{
    const std::vector<unsigned int> local_renumbering( this->system_management->block_info->local_renumbering );

    Vector<double> tmp(target.size());

    for(unsigned int i = 0; i < target.size(); ++i)
        tmp(local_renumbering[i]) = target(i);

    target = tmp;
}

// ---                    ---
// --- dealII_to_appframe ---
// ---                    ---
template<int dim>
void
NAME::EquationBase<dim>::dealII_to_appframe(FuelCell::ApplicationCore::MatrixVector& dst,
                                            const FullMatrix<double>&                src,
                                            const std::vector<unsigned int>&         matrix_block_indices) const
{
    FullMatrix<double> block_wise(src);
    this->standard_to_block_wise(block_wise);
    
    BlockIndices local = this->system_management->block_info->local;
        
    unsigned int number = 0, i_from = 0; //First row index to access for each iteration of for loop
    
    const unsigned int size = local.size();
    for(unsigned int line = 0; line < size; ++line)
    {
        unsigned int j_from = 0; //First column index to access for each iteration of for loop
        for(unsigned int column = 0; column < size; ++column)
        {
            if( (*(this->system_management->cell_couplings))(line,column) != 0 )
            {
                const unsigned int n_lines   = local.block_size(line);
                const unsigned int n_columns = local.block_size(column);
                
                for(unsigned int i = 0; i < n_lines; ++i)
                    for(unsigned int j = 0; j < n_columns; ++j)
                        dst[number].matrix(i,j) += block_wise(i+i_from, j+j_from);
                
                number++;
            }
            j_from += local.block_size(column);
        }
         i_from += local.block_size(line);
    }
}

// ---                    ---
// --- dealII_to_appframe ---
// ---                    ---

template<int dim>
void
NAME::EquationBase<dim>::dealII_to_appframe(FuelCell::ApplicationCore::FEVector& dst,
                                            const Vector<double>&                src,
                                            const std::vector<unsigned int>&     residual_indices) const
{
    AssertThrow( dst.size() == src.size() , ExcDimensionMismatch( dst.size() , src.size() ) );
    
    Vector<double> block_wise(src);
    this->standard_to_block_wise(block_wise);
    
    BlockIndices local = this->system_management->block_info->local;
    
    unsigned int line_index = 0;
    
    const unsigned int size = local.size();
    for(unsigned int line = 0; line < size; ++line)
    {
        const unsigned int n_lines = local.block_size(line);
        for(unsigned int i = 0; i < n_lines; ++i, ++line_index)
            dst.block(line)(i) += block_wise(line_index);
    }
}

/////////////////////
// MINOR FUNCTIONS //
/////////////////////
// ---                   ---
// --- print_caller_name ---
// ---                   ---
template<int dim>
void
NAME::EquationBase<dim>::print_caller_name(const std::string& caller_name) const
{
    const std::type_info& info = typeid(*this);
    FcstUtilities::log << "Pure function " << caller_name << " called in Class " << info.name() << std::endl;
}

/////////////////////////////
// EXPLICIT INSTANTIATIONS //
/////////////////////////////
template class NAME::EquationBase<deal_II_dimension>;