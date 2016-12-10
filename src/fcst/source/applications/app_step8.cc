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
// - $Id: app_step8.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <applications/app_step8.h>

namespace NAME = FuelCell::Application;

// ----------------------------------------------------------------------------
template <int dim>
NAME::AppStep8<dim>::AppStep8()
{
    FcstUtilities::log << "->Step8-" << dim;
    this->boundary_fluxes = false;
    this->interior_fluxes = false;
}

// ----------------------------------------------------------------------------
template <int dim>
void
NAME::AppStep8<dim>::initialize(ParameterHandler& param)
{
    // Initialize Mesh generator:
    param.enter_subsection("Grid generation");
    {
        param.set("Type of mesh","HyperCube");
        param.set("Initial refinement","2");
        
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
        param.set("Number of solution variables","2");
        param.enter_subsection("Solution variables");
        {
            param.set("Solution variable 1","test_var");
            param.set("Solution variable 2","test_var");
        }
        param.leave_subsection();   
        param.enter_subsection("Equations");
        {
        param.set("Equation 1","Test Equation");
        param.set("Equation 2","Test Equation");
        }
    param.leave_subsection();    
    }
    param.leave_subsection();
    
    // Initialize blockmatrix
    BlockMatrixApplication<dim>::initialize(param);
  
    // Specify the coupling between equations and variables
    this->cell_couplings.reinit(2,2);
    this->cell_couplings(0,0) = DoFTools::always;
    this->cell_couplings(1,1) = DoFTools::always;
    this->cell_couplings(0,1) = DoFTools::always;
    this->cell_couplings(1,0) = DoFTools::always;
    this->interior_fluxes = false;
    this->boundary_fluxes = false;    
    
    // Initialize matrices and spartisity pattern for the whole system
    this->remesh_matrices();
}

// ----------------------------------------------------------------------------
template <int dim>
void NAME::AppStep8<dim>::cell_matrix(MatrixVector& cell_matrices,
                                const typename DoFApplication<dim>::CellInfo& info)
{
    //Compute cell_matrix for the first equation:
    const FEValuesBase<dim>& fe = info.fe(0);
    const unsigned int n_dofs = fe.dofs_per_cell;
    
    std::vector<double> lambda_values(fe.n_quadrature_points, 1.0); //creates vector with n_quadrature_point components
    std::vector<double> mu_values(fe.n_quadrature_points, 1.0); //creates vector with n_quadrature_point components
    
    //FcstUtilities::log<<"DOF PER CELL: "<<n_dofs<<std::endl;
    
    for (unsigned int m=0; m<cell_matrices.size(); ++m)
        for (unsigned i=0;i<n_dofs;++i)
            for (unsigned j=0;j<n_dofs;++j)
                for (unsigned k=0;k<fe.n_quadrature_points;++k)
                {
                    cell_matrices[m].matrix(i,j) += (
                        (fe.shape_grad(i,k)[cell_matrices[m].row]*   //in this case where n_comp == dim, cell_matrices[m].row == component_i
                        fe.shape_grad(j,k)[cell_matrices[m].column]* //in this case where n_comp == dim, cell_matrices[m].column == component_j
                        lambda_values[k]) 
                        +
                        (fe.shape_grad(i,k)[cell_matrices[m].column]*
                        fe.shape_grad(j,k)[cell_matrices[m].row]*
                        mu_values[k]) 
                        +
                        ((cell_matrices[m].row == cell_matrices[m].column) ?
                        (fe.shape_grad(i,k)*
                        fe.shape_grad(j,k)*
                        mu_values[k]) :
                        0)
                    )*fe.JxW(k);
                }       
}

// ----------------------------------------------------------------------------
template <int dim>
void NAME::AppStep8<dim>::cell_residual(FEVector& cell_vector,
                                  const typename DoFApplication<dim>::CellInfo& info)
{
    
    const FEValuesBase<dim>& fe = info.fe(0);
    const unsigned int n_quad = fe.n_quadrature_points;
    const unsigned int n_dofs = fe.dofs_per_cell;
    
    std::vector<Point<dim> > point(n_quad); //creates vector with n_quadrature_point components
    point = fe.get_quadrature_points();
    std::vector<double> values_0(n_quad);
    std::vector<double> values_1(n_quad);
    
    
    for (unsigned int p=0; p<n_quad; ++p)
    {
        Point<dim> point_1, point_2;
        point_1(0) = 0.5;
        point_2(0) = -0.5;
        
        if (((point[p]-point_1).norm_square() < 0.2*0.2) ||
            ((point[p]-point_2).norm_square() < 0.2*0.2))
            values_0[p] = 1;
        else
            values_0[p] = 0;
        
        if (point[p].square() < 0.2*0.2)
            values_1[p] = 1;
        else
            values_1[p] = 0;
    }
    
    for (unsigned k=0;k<fe.n_quadrature_points;++k)
        for (unsigned i=0;i<n_dofs;++i)
        {
            cell_vector.block(0)(i) +=  values_0[k] * fe.shape_value(i,k) * fe.JxW(k);
            cell_vector.block(1)(i) +=  values_1[k] * fe.shape_value(i,k) * fe.JxW(k);
        }
}

// ----------------------------------------------------------------------------
template <int dim>
void NAME::AppStep8<dim>::dirichlet_bc(std::map<unsigned int, double>& boundary_values) const
{
    VectorTools::interpolate_boundary_values (*this->dof,
                                              0,
                                              ZeroFunction<dim>(this->element->n_blocks()),
                                              boundary_values);
}

// ----------------------------------------------------------------------------
template <int dim>
double NAME::AppStep8<dim>::estimate(const FEVector& sol)
{
    
    this->cell_errors.reinit(this->tr->n_active_cells());
    typename FunctionMap<dim>::type newmann_boundary;
    KellyErrorEstimator<dim>::estimate(*this->dof,
                                       QGauss<dim-1>(3),
                                       newmann_boundary,
                                       sol,
                                       this->cell_errors);
    
    return 0;
}

// ----------------------------------------------------------------------------
template <int dim>
double NAME::AppStep8<dim>::evaluate(const FEVectors&)
{
    return 0;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
template class NAME::AppStep8<deal_II_dimension>;
