//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: app_thermal_test.h 21-05-2013
//    - Description: Class designed to test thermal transport equation.
//    - Developers: Madhur Bhaiya
//    - Id: $Id: app_thermal_test.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELL__APP_THERMAL_TEST__H
#define _FUELCELL__APP_THERMAL_TEST__H

//-- OpenFCST
#include <application_core/optimization_block_matrix_application.h>
#include <layers/gas_diffusion_layer.h>
#include <equations/thermal_transport_equation.h>

namespace FuelCell
{
    namespace InitialSolution
    {
        /**
        * This class is used when solving the problem using Newton's method to provide an initial solution.
        * This function is called in VectorTools::interpolate(..,..,InitialSolution<dim> marc,...)
        * It provides a solution that satisfies Dirichlet boundaries and has a gradient.
        */
        template <int dim>
        class AppThermalTestIC
        :
        public Function<dim>
        {
        public:
            /**
            * Constructor
            */
            AppThermalTestIC (const std::string& type);
            /**
            * Destructor
            */
            virtual ~AppThermalTestIC (){};

            /**
            * This is the member function that computes the value of the initial
            * solution for a given point.
            */
            virtual double value (const Point<dim> &p,
                                  const unsigned int component = 0) const;
                                  
        private:
            
            std::string case_type;
        };
    } //end namespace InitialSolution
  

    namespace Application
    {
        //---------------------------------------------------------------------------
        //---------------------------------------------------------------------------
        //---------------------------------------------------------------------------
        /**
        * 
        * This class is used to test thermal transport equation (fourier equation) 
        * against analytical solutions.
        * 
        * @author Madhur Bhaiya; 2013
        */
        template <int dim>
        class AppThermalTest
        :
        public FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>
        {
        public:
            
            ///@name Constructors, destructor, and initalization
            //@{
    
            /**
            * Constructor.
            * @note the pointer data is initialized to boost::shared_ptr<> (), this means that
            * the pointer is empty and when we do data.get() it will return 0. This is good because at ApplicationBase
            * constructor an ApplicationData will be constructed.
            */
            AppThermalTest (boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data = boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> ());
        
            /**
            * Destructor
            */
            virtual ~AppThermalTest (){};
            /**
            * Declare all parameters that are needed for: 
            *   - the computation of the equation coefficients
            *   - the control of the linear system solution
            *   - ...
            */
            virtual void declare_parameters(ParameterHandler& param);
        
            /**
            * Function called by optimization loop in order to set the values in the 
            * ParameterHandler to the new design parameters.
            * Since ParameterHandler depends on the problem we are solving, set_parameters() is set
            * at the most inner loop of the application.
            */
            virtual void set_parameters(const std::vector<std::string>& name_dvar,
                                        const std::vector<double>& value_dvar,
                                        ParameterHandler& param){};    
        
            /**
            * Initilize and create the grid
            * Determines whether reading from a grid file or creating new
            */
            void initialize_triangulation(ParameterHandler& param);
        
            /**
            * Call the other initialize routines from the inherited classes
            */
            virtual void initialize(ParameterHandler& param);
      
            /**
            * Initialize nonlinear solution
            */  
            virtual void initialize_solution (FuelCell::ApplicationCore::FEVector& initial_guess,
                                              std::shared_ptr<Function<dim> > initial_function = std::shared_ptr<Function<dim> >());
            //@}
            
            ///@name Local CG FEM based assemblers
            //@{
            /**
            * Integration of local bilinear form. Here we loop over the quadrature 
            * points and over degrees of freedom in order to compute the matrix for the cell
            * This routine depends on the problem at hand and is called by assemble() in DoF_Handler
            * class.
            */
            virtual void cell_matrix(FuelCell::ApplicationCore::MatrixVector& cell_matrices,
                                     const typename DoFApplication<dim>::CellInfo& cell);  
            /**
            * Integration of the rhs of the equations. Here we loop over the quadrature 
            * points and over degrees of freedom in order to compute the right
            * hand side for each cell
            * This routine depends on the problem at hand and is called by residual() in DoF_Handler
            * class
            * @note This function is called residual because in the case of nonlinear systems
            * the rhs is equivalent to the residual
            */
            virtual void cell_residual(FuelCell::ApplicationCore::FEVector& cell_vector,
                                       const typename DoFApplication<dim>::CellInfo& cell);
            
            /**
            * Assemble local boundary matrix.
            */
            virtual void bdry_matrix(FuelCell::ApplicationCore::MatrixVector& bdry_matrices,
                                     const typename DoFApplication<dim>::FaceInfo& bdry_info);
            
            /**
            * Assemble local boundary residual.
            */
            virtual void bdry_residual(FuelCell::ApplicationCore::FEVector& bdry_vector,
                                       const typename DoFApplication<dim>::FaceInfo& bdry_info);
            //@}
        
            /**
            * Member function used to set dirichlet boundary conditions.
            * This function is application specific and it only computes the boundary_value
            * values that are used to constraint the linear system of equations that is being
            * solved
            */
            virtual void dirichlet_bc(std::map<unsigned int, double>& boundary_values) const;
      
            /**
            * Estimate error per cell
            */
            virtual double estimate(const FuelCell::ApplicationCore::FEVectors& sol);

            
            /**
            * Reimplementation of the routine in the base class BaseApplication in namespace AppFrame so
            * that the right labels are outputed and so that I can compute and output the source term.
            */
            virtual void data_out(const std::string& basename,
                                  const FuelCell::ApplicationCore::FEVectors& vectors);
            
            /**
            * Compute the value of all objective function and constraints 
            */
            virtual void cell_responses (std::vector<double>& resp,
                                        const typename DoFApplication<dim>::CellInfo& info,
                                        const FuelCell::ApplicationCore::FEVector& sol){};
            /**
            * This class is used to evaluate all responses that do not require looping over cells.
            * An example of one of this types of constraints is the solid volume fraction.
            */
            virtual void global_responses (std::vector<double>& resp,
                                            const FuelCell::ApplicationCore::FEVector& sol){};
                                      
                                      
        protected:
            

            ///@name MEA layers
            //@{
            boost::shared_ptr< FuelCellShop::Layer::GasDiffusionLayer<dim> > test_layer;
            //@}
        
            ///@name Physics equations
            //@{
            /**
            * ThermalTransportEquation object
            */
            FuelCellShop::Equation::ThermalTransportEquation<dim> thermal_transport;
            //@}
            
            /** String to determine the case type */
            std::string case_type;
            

        };
    }
}

#endif //_FUELCELL__AppThermalTest_H
