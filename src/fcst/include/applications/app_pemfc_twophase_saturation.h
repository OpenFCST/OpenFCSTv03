//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: app_pemfc_nonisothermal.h 21-05-2013
//    - Description: Class designed to solve a twophase saturation PEMFC model.
//    - Developers: Jie Zhou
//    - $$
//
//---------------------------------------------------------------------------

#ifndef _FUELCELL__APP_PEMFC_TWOPHASE_SATURATION__H
#define _FUELCELL__APP_PEMFC_TWOPHASE_SATURATION__H

//-- OpenFCST 
#include <application_core/optimization_block_matrix_application.h>
#include <utils/operating_conditions.h>

#include <layers/gas_diffusion_layer.h>
#include <layers/micro_porous_layer.h>
#include <layers/catalyst_layer.h>
#include <layers/membrane_layer.h>
#include <materials/PureGas.h>

#include <equations/thermal_transport_equation.h>
#include <equations/reaction_source_terms.h>
#include <equations/proton_transport_equation.h>
#include <equations/lambda_transport_equation.h>
#include <equations/electron_transport_equation.h>
#include <equations/sorption_source_terms.h>
#include <equations/ficks_transport_equation.h>
#include <equations/saturation_transport_equation.h>

#include <postprocessing/data_out.h>
#include <postprocessing/response_current_density.h>
#include <postprocessing/response_ohmic_heat.h>
#include <postprocessing/response_sorption_heat.h>
#include <postprocessing/response_reaction_heat.h>
#include <postprocessing/response_water_sorption.h>

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
        class AppPemfcTPSaturationIC
        :
        public Function<dim>
        {
        public:
            /**
            * Constructor
            */
            AppPemfcTPSaturationIC (FuelCell::OperatingConditions* OC,
                                 boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid,
                                 FuelCell::SystemManagement* system_mgmt);
            /**
            * Destructor
            */
            ~AppPemfcTPSaturationIC ();

            /**
            * This is the member function that computes the value of the initial
            * solution for a given point.
            */
            void vector_value (const Point<dim> &p,
                               Vector<double> &v) const;
      
        private:
            /** Operating conditions class object*/
            FuelCell::OperatingConditions* OC;
            
            /** Geometry class object */
            boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid;
            
            FuelCell::SystemManagement* system;
        };
    } //end namespace InitialSolution
  

    namespace Application
    {
        //---------------------------------------------------------------------------
        //---------------------------------------------------------------------------
        //---------------------------------------------------------------------------
        /**
        * 
        * This class is used to solve the physical pheonoma on a complete
        * membrane electrode assembly, including non-isothermal and two phase effects.
        * This application model is detailed in the journal article: 
        *
        * \remarks
        * - This application can be used to model homogeneous, agglomerate or graded catalyst layers, utilizing different kinetic models.
        * - This application can also be used to simulate a MEA without MPL.
        *
        * For more information on the discretization and solution methodology, please read the following reference:
        * - M. Secanell et al. "Numerical Optimization of Proton Exchange Membrane Fuel Cell Cathode Electrodes", Electrochimica Acta, 52(7):2668-2682, February 2007.
        *
        * If using this application, please cite the article above as references.
        * 
        * @author Jie Zhou (bhaiya@ualberta.ca); 2013-14
        */
        template <int dim>
        class AppPemfcTPSaturation
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
            AppPemfcTPSaturation (boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data = 
            boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> ());
        
            /**
            * Destructor
            */
            ~AppPemfcTPSaturation ();
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
                                        ParameterHandler& param) {};    
            
            /**
            * Set up how many equations are needed and 
            * read in parameters for the parameter handler in order to initialize data
            */
            void _initialize(ParameterHandler& param);
        
            /**
            * Call the other initialize routines from the inherited classes
            */
            virtual void initialize(ParameterHandler& param);
      
            /**
            * Initialize nonlinear solution
            */
            virtual void initialize_solution (FEVector& initial_guess,
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
            * This class is called by responses to make sure that all responses requested are
            * implemented in either cell_responses, global_responses or face_responses.
            * @note Every time we add a new response that we can calculate we need to update this
            * file.
            */
            virtual void check_responses();
            
            /**
            * Compute the value of all objective function and constraints 
            */
            virtual void cell_responses (std::vector<double>& resp,
                                        const typename DoFApplication<dim>::CellInfo& info,
                                        const FuelCell::ApplicationCore::FEVector& sol);
            /**
            * This class is used to evaluate all responses that do not require looping over cells.
            * An example of one of this types of constraints is the solid volume fraction.
            */
            virtual void global_responses (std::vector<double>& resp,
                                            const FuelCell::ApplicationCore::FEVector& sol);
            
            /**
            * This class is used to evaluate the derivative of all the functionals that require looping over cells
            * with respect to the design variables.
            * This class is called by responses to evaluate the response at each cell.
            */
            virtual void cell_dresponses_dl(std::vector<std::vector<double> >& /*cell_df_dl*/,
                                            const typename DoFApplication<dim>::CellInfo& /*info*/,
                                            const FuelCell::ApplicationCore::FEVector& /*sol*/) {};
                                        
            /**
            * This class is used to evaluate the sensitivities of all responses that do not require looping over cells
            * with respect of the design variables.
            * An example of one of this types of constraints is the solid volume fraction.
            */
            virtual void global_dresponses_dl(std::vector<std::vector<double> >& df_dl,
                                            const FuelCell::ApplicationCore::FEVector& sol);
            /**
            * This class is used to evaluate the derivative of all the functionals that require looping over cells
            * with respect of the unknowns of the system of governing equations.
            * This class is called by responses to evaluate the response at each cell.
            */
            virtual void cell_dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector >& /*cell_df_du*/,
            const typename DoFApplication<dim>::CellInfo& /*info*/,
            std::vector<std::vector<double> >& /*src*/) {};
                                        
            /**
            * This class is used to evaluate the sensitivities of all responses that do not require looping over cells
            * with respecto of the unknowns of the system of governing equations.
            * An example of one of this types of constraints is the solid volume fraction.
            */
            virtual void global_dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector >& df_du,
                                            const FuelCell::ApplicationCore::FEVector& src);
            
            
            /**
            * Post-processing. Evaluate a functional such as the objective function of an
            * optimization problem
            */
            virtual double evaluate (const FuelCell::ApplicationCore::FEVectors& src);
            
            /**
            * Reimplementation of the routine in the base class BaseApplication in namespace AppFrame so
            * that the right labels are outputed and so that I can compute and output the source term.
            */
            virtual void data_out(const std::string& filename, 
                                  const FuelCell::ApplicationCore::FEVectors& src);                                      
                                      
        protected:
            
           ///@name Pre-processor and operating condition classes
            //@{
            
            /** Operating conditions */
            FuelCell::OperatingConditions OC;
            //@} 
      
            ///@name Material classes
            //@{
            
            /** The cathode/anode contains water vapour, so we need to create an object water
            * in order to compute viscosity, density, etc. for water
            */
            FuelCellShop::Material::WaterVapor water;
            
            /** The cathode contains oxygen, so we need to create an object oxygen
            * in order to compute viscosity, density, etc. for oxygen
            */
            FuelCellShop::Material::Oxygen oxygen;
            
            /** The cathode contains nitrogen as solvent, so we need to create an object nitrogen
            * in order to compute viscosity, density, etc. for nitrogen
            */
            FuelCellShop::Material::Nitrogen nitrogen;
            
            /** The anode contains hydrogen, so we need to create an object hydrogen
            * in order to compute viscosity, density, etc. for hydrogen
            */
            FuelCellShop::Material::Hydrogen hydrogen;
            //@}
            
            ///@name MEA layers
            //@{
            boost::shared_ptr< FuelCellShop::Layer::GasDiffusionLayer<dim> > AGDL; ///< Anode GDL Layer.
            boost::shared_ptr< FuelCellShop::Layer::GasDiffusionLayer<dim> > CGDL; ///< Cathhode GDL Layer.
            boost::shared_ptr< FuelCellShop::Layer::MicroPorousLayer<dim> > AMPL; ///< Anode MPL Layer.
            boost::shared_ptr< FuelCellShop::Layer::MicroPorousLayer<dim> > CMPL; ///< Cathode MPL Layer.
            boost::shared_ptr< FuelCellShop::Layer::CatalystLayer<dim> > ACL; ///< Anode Catalyst Layer.
            boost::shared_ptr< FuelCellShop::Layer::CatalystLayer<dim> > CCL; ///< Cathode Catalyst Layer.
            boost::shared_ptr< FuelCellShop::Layer::MembraneLayer<dim> > ML; ///< Membrane Layer.
            //@}           
        
            ///@name Physics equations
            //@{
            /**
            * ThermalTransportEquation object
            */
            FuelCellShop::Equation::ThermalTransportEquation<dim> thermal_transport;
	    
	    /**
	     * SaturationTransportEquation object
	     */
	    FuelCellShop::Equation::SaturationTransportEquation<dim> saturation_transport;
	    
	    /**
	     * ProtonTransportEquation object
	     */
            FuelCellShop::Equation::ProtonTransportEquation<dim> proton_transport;
            
            /**
            * LambdaTransportEquation object
            */
            FuelCellShop::Equation::LambdaTransportEquation<dim> lambda_transport;
            
            /**
            * ElectronTransportEquation object
            */
            FuelCellShop::Equation::ElectronTransportEquation<dim> electron_transport;
            
            /**
            * ReactionSourceTerms object
            */
            FuelCellShop::Equation::ReactionSourceTerms<dim> reaction_source_terms;
            
            /**
            * SorptionSourceTerms object
            */
            FuelCellShop::Equation::SorptionSourceTerms<dim> sorption_source_terms;
            
            FuelCellShop::Equation::FicksTransportEquation<dim> ficks_oxygen_nitrogen;
            
            FuelCellShop::Equation::FicksTransportEquation<dim> ficks_water_nitrogen;
            
            FuelCellShop::Equation::FicksTransportEquation<dim> ficks_water_hydrogen;
	    
	    
            //@}
            
            ///@name Post-processing objects (Functional evaluation)
            //@{
            FuelCellShop::PostProcessing::ORRCurrentDensityResponse<dim> ORRCurrent;
            
            FuelCellShop::PostProcessing::HORCurrentDensityResponse<dim> HORCurrent;
            
            FuelCellShop::PostProcessing::ProtonOhmicHeatResponse<dim> protonOhmicHeat;
            
            FuelCellShop::PostProcessing::ElectronOhmicHeatResponse<dim> electronOhmicHeat;
            
            FuelCellShop::PostProcessing::SorptionHeatResponse<dim> sorptionHeat;
            
            FuelCellShop::PostProcessing::ORRReactionHeatResponse<dim> catReactionHeat;
            
            FuelCellShop::PostProcessing::HORReactionHeatResponse<dim> anReactionHeat;
            
            FuelCellShop::PostProcessing::WaterSorptionResponse<dim> waterSorption;
            //@}
        
            /** Stores the design variable names so that the name can be appended to the .vtk file name. */
            std::vector<std::string> design_var;
            
            /** Stores the values of the design variables so that the number can be appended to the .vtk file name. */
            std::vector<double> design_var_value;
      
        private:
      
            double l_channel; ///< Width of the channel.
            
            double l_land; ///< Width of the landing.
            
            /**
            * Time constant for sorption isotherm [\p 1/s]
            */
            double time_k;
        };
    }
}

#endif //_FUELCELL__AppPemfcNIThermal_H
