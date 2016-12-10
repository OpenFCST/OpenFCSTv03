// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: app_block_matrix_step8.cc
// - Description: Test for application framework. Solves deal.II step-8 tutorial
// - Developers: Marc Secanell, University of Alberta
// - $Id: app_step8.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

//-- OpenFCST
#include <application_core/block_matrix_application.h>

using namespace FuelCell::ApplicationCore;

namespace FuelCell
{
    namespace Application
    {
        
        /**
         * This class is used to test the application framework by solving Step-3 of the deal.II finite element libraries.
         * The problem solved is a simple Poisson equation with a forcing term equal to one.
         * 
         * This class, together with the code in FuelCell::UnitTest::ApplicationStep3Test also show a simple
         * application and how it can be solved by developing your own main file.
         *  
         * @note Please do not modify this application as it is used to test the code. Modifying the application will
         * result in unit test FuelCell::UnitTest::ApplicationStep3Test failure.
         *
         * @author Marc Secanell, 2014
         */
        template <int dim>
        class AppStep8 :
        public BlockMatrixApplication<dim>
        {
        public:
            /**
             * Constructor
             */
            AppStep8 ();
              
            /**
             * Initialize variables specific to the application
             */
            void initialize(ParameterHandler& param);
            
            /**
             * Implement the element-wise system matrix. In this case, the routine implements
             * the weak form of a Laplace operator
             */
            virtual void cell_matrix(MatrixVector& cell_matrices,
                                     const typename DoFApplication<dim>::CellInfo& cell);  
            
            /**
             * Implements the element-wise right-hand side. In this case, the routine
             * implements a unit vector
             */
            virtual void cell_residual(FEVector& cell_vector,
                                       const typename DoFApplication<dim>::CellInfo& cell);
            
            /**
             * Apply any Dirichlet boundary conditions. In this case all boundaries are set to zero.
             */
            virtual void dirichlet_bc(std::map<unsigned int, double>& boundary_values) const;
            
            /**
             * Estimate error in each cell
             */
            virtual double estimate(const FEVector& sol);
            
            /**
             * Evaluate any data necessary for testing.
             */
            virtual double evaluate (const FEVectors& src);
        };
    }
}
