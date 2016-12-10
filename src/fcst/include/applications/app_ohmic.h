//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: app_ohmic.cc
//    - Description:
//    - Developers: M. Sabharwal, M. Secanell
//
//---------------------------------------------------------------------------

#ifndef _FUELCELL__APP_OHMIC__H
#define _FUELCELL__APP_OHMIC__H

//-- OpenFCST
#include <application_core/optimization_block_matrix_application.h>
#include <utils/operating_conditions.h>
#include <materials/PureGas.h>
#include <materials/GasMixture.h>
#include <equations/electron_transport_equation.h>
#include <postprocessing/data_out.h>

// Use namespace of deal.II
using namespace dealii;

namespace FuelCell
{

    namespace Application
    {
        //---------------------------------------------------------------------------
        //---------------------------------------------------------------------------
        //---------------------------------------------------------------------------
        /**
         * The application can be used to simulate a electron transport through a porous media. 
         * This class is used to solve a electronic transport problem using Ohm's Law. 
         *
         * The equation solved is written as follows:
         * \f$ \qquad \mathbf{\nabla} \cdot \left( \hat{\sigma}_{s,eff} \cdot \mathbf{\nabla} \phi_s \right) = 0 \quad \in \quad \Omega \f$
         *
         * - \f$ \hat{\sigma}_{s,eff} \f$ is effective electron conductivity tensor [\p S/cm].
         *
         * The governing equation in the weak form is actually linear and can be solved directly using a linear solver like UMFPACK or GMRES.
         * The main file for this application should look like this:
         * @code
         * subsection Simulator
         *  set simulator name = ohmic
         *  set simulator parameter file name = data.prm
         *  set solver name = Linear 
         *  set Dakota direct = false 
         * end
         * @endcode
         * NOTE: The solver name is set to Linear instead of the Newton solvers
         * The lines below indicate how to setup your boundary and initial conditions.
         * [Although initial conditions are not required for this problem it is advisable to use them just to see that the solver is working.]
         *
         * @code
         * subsection Equations
         * 
         * subsection Electron Transport Equation
         * 
         *  subsection Initial data
         *   set electronic_electrical_potential = 0: 0.1 # where 0 indicates the material_id setup in the grid and 0.1 is the electronic potential in V
         *  end
         * 
         *  subsection Boundary data
         *   set electronic_electrical_potential = 1: 0.4, 2:0.01 #where 1 & 2 denote the boundary ids and 0.4 and 0.01 are the electronic potentials in V at the respective boundary
         *  end
         * 
         * end
         * @endcode
         *
         * @author Mayank Sabharwal & Marc Secanell, 2014
         *
         */
        template<int dim>
        class AppOhmic : public OptimizationBlockMatrixApplication<dim>
        {
        public:
            
            ///@name Constructors, destructor, and initialization
            //@{
            
            /**
             * Constructor.
             */
            AppOhmic( boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = 
            boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >() );
            
            /**
             * Destructor.
             */
            ~AppOhmic();
            
            /**
             * Declare parameters.
             */
            virtual void declare_parameters(ParameterHandler& param);
            
            /**
             * Initialize parameters.
             */
            virtual void initialize(ParameterHandler& param);
            
            /**
             * The initial guess along with the appropriate BCs is formed here.
             */
              virtual void initialize_solution (FEVector& initial_guess,
                                                std::shared_ptr<Function<dim> > initial_function = std::shared_ptr<Function<dim> >());

            //@}
            
            ///@name Local CG FEM based assemblers
            //@{
            
            /**
             * Assemble local cell matrix.
             */
            virtual void cell_matrix(MatrixVector&                                 cell_matrices,
                                     const typename DoFApplication<dim>::CellInfo& cell_info);
            
            /**
             * Assemble local cell residual.
             */
            virtual void cell_residual(FuelCell::ApplicationCore::FEVector&          cell_res,
                                       const typename DoFApplication<dim>::CellInfo& cell_info);
            
            //@}
            
            ///@name Other functions
            //@{
            
            /**
             * BCs.
             */
            virtual void dirichlet_bc(std::map<unsigned int, double>& boundary_values) const;

            /**
             * Evaluate the response functions: Not implemented in this class; Returns 0.0
             */
            virtual double evaluate (const FuelCell::ApplicationCore::FEVectors& src);
            
            /**
             * Output results.
             */
            virtual void data_out(const std::string&         filename,
                                  const FEVectors& src);
            
            //@}
            
            ///@name Post-processing
            //@{
            
            /**
             * Compute some functionals.
             */
            virtual void bdry_responses(std::vector<double>&                                                     dst,
                                        const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                        const FuelCell::ApplicationCore::FEVector& src);
            //@}
            
        protected:            
            ///@name System matrix and boundary condition objects
            //@{

            /**
             *
             * component_boundaryID_value_maps info:
             *
             * subsection Initial data
             * set electronic_electrical_potential = 0: 0.1
             * end
             * subsection Boundary data
             * set electronic_electrical_potential = 1: 0.4, 2:0.01
             * end
             * end
             */
            
            //@}
            
            ///@name Operating conditions data
            //@{       
            /**
             * Operating conditions.
             */
            FuelCell::OperatingConditions OC;
	    
            FuelCellShop::Material::Oxygen oxygen;
	    
	    FuelCellShop::Material::Nitrogen nitrogen;
                                    
            //@}
            
            ///@name Layer objects
            //@{
            /**
             * Cathode GDL.
             */
            boost::shared_ptr< FuelCellShop::Layer::GasDiffusionLayer<dim> > CGDL;
            
            //@}
            
            ///@name Equation objects
            //@{
            /**
             * This object describes
             * the equations that we are going to
             * solve here.
             */
            FuelCellShop::Equation::ElectronTransportEquation<dim> electron_transport_equation;
            
            /**
             * This object describes
             * the equations that we are going to
             * solve here.
             */
            

        private:
            /**
             * Compute some functionals that are not needed for most applications (this section is not necessary in most cases.)
             */
            
            
        };
        
    } // Application
    
} // FuelCell

#endif
