//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: app_diffusion.h
//    - Description:
//    - Developers: M. Sabharwal, M. Secanell, A. Kosakian
//
//---------------------------------------------------------------------------

#ifndef _FUELCELL__APP_DIFFUSION__H
#define _FUELCELL__APP_DIFFUSION__H

//-- OpenFCST
#include <application_core/optimization_block_matrix_application.h>
#include <utils/operating_conditions.h>
#include <materials/PureGas.h>
#include <materials/GasMixture.h>
/**#include <materials/platinum.h>
#include <materials/nafion.h>
#include <materials/carbon.h>*/
#include <layers/gas_diffusion_layer.h>
#include <equations/ficks_transport_equation.h>
#include <equations/equation_auxiliaries.h>
#include <postprocessing/data_out.h>


// Use namespace of deal.II
using namespace dealii;
using namespace FuelCellShop::Equation::DebugTools;

namespace FuelCell
{
   
    namespace Application
    {
        //---------------------------------------------------------------------------
        //---------------------------------------------------------------------------
        //---------------------------------------------------------------------------
        /**
         * The application can be used to simulate a gas flow through a porous media. 
         * This class is used to solve a diffusion problem using Fick's law. 
         *
         * The equation solved is given by
         * \f[
         * R_1(x_{O_2}) = \nabla \cdot (c_{total}D^{eff}_{O_2} \nabla x_{O_2} ) = 0, \qquad x_{O_2} \in \Omega.
         * \f]
         * 
         * The governing equation in the weak form is actually linear and can be solved directly using an iterative solver, such as CG, BiStabCG or GMRes. 
         * You can also use direct solvers, such as UMFPACK, and parallel iterative or direct solvers.
         * 
         * The main file for this application should look like this:
         * @code
         * subsection Simulator
         *  set simulator name = diffusion
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
         * subsection Ficks Transport Equation - oxygen
         * 
         *  subsection Initial data
         *   set oxygen_molar_fraction = 0: 100 # where 0 indicates the material_id setup in the grid and 100 is the concentration of solute in mol/cm^3
         *  end
         * 
         *  subsection Boundary data
         *   set oxygen_molar_fraction = 1: 0.4, 2:0.01 #where 1 & 2 denote the boundary ids and 0.4 and 0.01 are the concentrations of solute in mol/cm^3 at the respective boundary
         *  end
         * 
         * end
         * @endcode
         *
         * @author Marc Secanell, 2013
         * @author Mayank Sabharwal, 2014
         *
         */
        template<int dim>
        class AppDiffusion : public OptimizationBlockMatrixApplication<dim>
        {
        public:
            
            ///@name Constructors, destructor, and initialization
            //@{
            
            /**
             * Constructor.
             */
            AppDiffusion( boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = 
            boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >() );
            
            /**
             * Destructor.
             */
            ~AppDiffusion();
            
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
            /**
             * Assemble local bdry matrix.
             */
            virtual void bdry_matrix(MatrixVector&                                 bdry_matrices,
                                     const typename DoFApplication<dim>::FaceInfo&  bdry_info);
            
            /**
             * Assemble local bdry residual.
             */
            virtual void bdry_residual(FuelCell::ApplicationCore::FEVector&          bdry_res,
                                       const typename DoFApplication<dim>::FaceInfo& bdry_info);
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
            
            ///@name Debug output
            //@{
            
            /**
             * Returns solution index.
             */
            unsigned int get_solution_index();
            //@}               
            
        protected:            
            ///@name System matrix and boundary condition objects
            //@{

            /**
             *
             * component_boundaryID_value_maps info:
             *
             * subsection Initial data
             * set oxygen_molar_fraction = 0: 100
             * end
             * subsection Boundary data
             * set oxygen_molar_fraction = 1: 0.4, 2:0.01
             * end
             * end
             */
            
            //@}
            
            ///@name Gases and operating conditions data
            //@{       
            /**
             * Operating conditions.
             */
            FuelCell::OperatingConditions OC;
            
            /**
             * Solute.
             */
            boost::shared_ptr<FuelCellShop::Material::PureGas> solute;
            
            
            
            /**
             * Solvent.
             */
            boost::shared_ptr<FuelCellShop::Material::PureGas> solvent;
            
                        
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
            FuelCellShop::Equation::FicksTransportEquation<dim> ficks_transport_equation;
            //@}

        private:
            
            ///@name Debug output
            //@{
            /**
             * Returns solution variable info object.
             */
            FuelCellShop::Equation::VariableInfo get_xi();
            
            /**
             * Object for debug output purposes.
             */
            DebugOutput<dim> equation_debug_output;              
            //@}               
            
            /**
             * Compute some functionals that are not needed for most applications (this section is not necessary in most cases.)
             */
        };
        
    } // Application
    
} // FuelCell

#endif
