// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: app_step3.cc
// - Description: Test for application framework. Solves deal.II step-3 tutorial
// - Developers: Marc Secanell, University of Alberta
// - $Id: app_step3.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <applications/app_step3.h>

namespace NAME = FuelCell::Application;

// ----------------------------------------------------------------------------
template <int dim>
NAME::AppStep3<dim>::AppStep3()
{
    FcstUtilities::log << "->Step3-" << dim;
    this->boundary_fluxes = false;
    this->interior_fluxes = false;
}

// ----------------------------------------------------------------------------
template <int dim>
void
NAME::AppStep3<dim>::initialize(ParameterHandler& param)
{
    
    // Initialize Mesh generator:
    param.enter_subsection("Grid generation");
    {
        param.set("Type of mesh","HyperCube");
        param.set("Initial refinement","3");
        
        param.enter_subsection("Internal mesh generator parameters");
        {
            param.enter_subsection("Dimensions");
            {
                param.set("HyperCube edge", "1.0");
            }
            param.leave_subsection(); 
        }
        param.leave_subsection(); 
    }
    param.leave_subsection(); 

    // Initialize System management:    
    param.enter_subsection("System management");
    {
        param.set("Number of solution variables","1");
        param.enter_subsection("Solution variables");
        {
            param.set("Solution variable 1","test_var");
        }
        param.leave_subsection();   
        param.enter_subsection("Equations");
        {
        param.set("Equation 1","Test Equation");
        }
    param.leave_subsection();    
    }
    param.leave_subsection();
    
    BlockMatrixApplication<dim>::initialize(param);
    
    // Specify the coupling between equations and variables
    this->cell_couplings.reinit(1,1);
    this->cell_couplings(0,0) = DoFTools::always;
    this->interior_fluxes = false;
    this->boundary_fluxes = false;
    
    // Initialize matrices and spartisity pattern for the whole system
    this->remesh_matrices();
}

// ----------------------------------------------------------------------------
template <int dim>
void NAME::AppStep3<dim>::cell_matrix(MatrixVector& cell_matrices,
                                const typename DoFApplication<dim>::CellInfo& info)
{
    Assert(cell_matrices.size() == 1,
           ExcDimensionMismatch(cell_matrices.size(), 1));
    Assert(cell_matrices[0].row == 0,
           ExcDimensionMismatch(cell_matrices[0].row, 0));
    Assert(cell_matrices[0].column == 0,
           ExcDimensionMismatch(cell_matrices[0].column, 0));
    
    const FEValuesBase<dim>& fe =  info.fe(0);  
    const unsigned int n_dofs = fe.dofs_per_cell;
    
    double factor(1.0);
    
    for (unsigned k=0;k<fe.n_quadrature_points;++k)
        for (unsigned i=0;i<n_dofs;++i)
            for (unsigned j=0;j<n_dofs;++j)
            {
                cell_matrices[0].matrix(i,j) += factor * (fe.shape_grad(j,k) * fe.shape_grad(i,k)) * fe.JxW(k);
            }
}

// ----------------------------------------------------------------------------
template <int dim>
void NAME::AppStep3<dim>::cell_residual(FEVector& cell_vector,
                                  const typename DoFApplication<dim>::CellInfo& info)
{
    const FEValuesBase<dim>& fe = info.fe(0);
    
    double factor(1.0);    
    
    const unsigned int n_dofs = fe.dofs_per_cell;
    
    for (unsigned k=0;k<fe.n_quadrature_points;++k)
        for (unsigned i=0;i<n_dofs;++i)
            cell_vector.block(0)(i) +=  factor * fe.shape_value(i,k) * fe.JxW(k);
}

// ----------------------------------------------------------------------------
template <int dim>
void NAME::AppStep3<dim>::dirichlet_bc(std::map<unsigned int, double>& boundary_values) const
{
    VectorTools::interpolate_boundary_values (*this->dof,
                                              0,
                                              ConstantFunction<dim>(0.0),
                                              boundary_values);
}

// ----------------------------------------------------------------------------
template <int dim>
double NAME::AppStep3<dim>::evaluate(const FEVectors&)
{
    return 0;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
template class NAME::AppStep3<deal_II_dimension>;
