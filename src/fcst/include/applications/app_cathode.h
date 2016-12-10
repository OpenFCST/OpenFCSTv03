// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: app_cathode.h
// - Description: This class describes diffusion in fuel cell cathodes
//                Ficks, one gas
// - Developers: Marc Secanell,      University of Alberta
//               Valentin N. Zingan, University of Alberta
//               
// ----------------------------------------------------------------------------

#ifndef _FCST_APPLICATION_APP_CATHODE_H_
#define _FCST_APPLICATION_APP_CATHODE_H_

//-- OpenFCST
#include <application_core/optimization_block_matrix_application.h>
#include <application_core/system_management.h>
#include <utils/operating_conditions.h>

#include <materials/PureGas.h>
#include <materials/GasMixture.h>
#include <materials/platinum.h>
#include <materials/nafion.h>
#include <materials/carbon.h>

#include <layers/gas_diffusion_layer.h>
#include <layers/micro_porous_layer.h>
#include <layers/catalyst_layer.h>

#include <equations/ficks_transport_equation.h>
#include <equations/electron_transport_equation.h>
#include <equations/proton_transport_equation.h>
#include <equations/reaction_source_terms.h>
#include <postprocessing/data_out.h>
#include <postprocessing/response_current_density.h>

using namespace dealii;
using namespace FuelCell::ApplicationCore;

namespace FuelCell
{
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    namespace Application
    {
        /**
         * 
         * This class is used to solve a system of equations similar to the one presented in the journal article
         * M. Secanell et al. "Numerical Optimization of Proton Exchange Membrane Fuel Cell Cathode Electrodes",
         * Electrochimica Acta, 52(7):2668-2682, February 2007 and also presented at the ECCM 2006 conference under
         * the abstract titled M. Secanell et al, "A PEM fuel cell electrode model for gradient-based optimization", June 2006.
         *
         *
         * The cathode catalyst layer and the gas diffusion layer are modelled using 
         * a set of equations found in the journal article: 
         * 
         * - M. Secanell et al. "Numerical Optimization of Proton Exchange Membrane Fuel Cell Cathode Electrodes", 
         * Electrochimica Acta, 52(7):2668-2682, February 2007. 
         * 
         * The simulation is isothermal, isotropic, models
         * the electrochemical reactions using tafel kinetics, oxygen diffusion using Ficks Law and electronic
         * and protonic transport using Ohms Law. 
         * 
         * The equations solved are written as follows:
         * \f[
         * R_1(\vec{u}) = \nabla \cdot (c_{total}D^{eff}_{O_2} \nabla x_{O_2} ) = \frac{1}{nF}\nabla\cdot i
         * \f]
         * \f[
         * R_2(\vec{u}) = \nabla \cdot (\sigma^{eff}_{m}\nabla\phi_m) = \nabla \cdot i
         * \f]
         * \f[
         * R_3(\vec{u}) = \nabla \cdot (\sigma^{eff}_{s}\nabla\phi_s) = -\nabla \cdot i
         * \f]
         * where, for the case of a macro-homogeneous catalyst layer,
         * \f[
         * \nabla \cdot i = A_v i^{ref}_0 \left( \frac{c_{O_2}} {c^{ref}_{O_2}} \right) \mbox{exp} 
         * \left( \frac{\alpha F}{RT}(\phi_m - \phi_s) \right)
         * \f]
         * The solution variables, \f$ \vec{u} \f$ are the protonic potential, \f$\phi_m\f$, the electronic potential, 
         * \f$\phi_s\f$ and, instead of the oxygen concentration, we use the oxygen molar fraction, 
         * \f$x_{O_2}\f$, that also accounts for the oxygen dissolving in the ionomer by using Henrys Law. 
         * 
         * The govering equations above are nonlinear and therefore they cannot be solved directly.
         * 
         * In OpenFCST, we have decided to solve the system of equations using a nonlinear Newton solver. Therefore,
         * instead of implemeting the equations above, #cell_matrix and #cell_residual
         * implement a linearization of the governing equations above, i.e.
         * \f[
         * \frac{\partial R_i(\vec{u^{n-1}})}{\partial u_j} u_j = R_i(\vec{u^{n-1}})
         * \f]
         * 
         * For more information on the governing equations, discretization and solution methodology
         * please read the following references:
         * - M. Secanell et al. "Numerical Optimization of Proton Exchange Membrane Fuel Cell Cathode Electrodes", 
         * Electrochimica Acta, 52(7):2668-2682, February 2007. 
         * 
         * Also, the folder /data/cathode/Secanell_EA07_Numerical_Optimization_PEMFC_Cathode_Electrodes contains
         * a sample datafile to run this application with a macro-homogeneous model. 
         * 
         * The application can be used to solve a cathode with and without an MPL. You can select the appropriate
         * geometry by changing the following in the input file:
         * 
         * @code
         * subsection Grid generation
         *   set Type of mesh = Cathode OR CathodeMPL
         *   subsection Internal mesh generator parameters
         *      subsection Dimensions
         *       set Cathode current collector width [cm] = 0.1 # Thickness of the rib of the bipolar plates (BPP)
         *       set Cathode channel width [cm] = 0.1           # Thickness of the channel on the BPP
         *       set Cathode CL thickness [cm] = 1E-3           # Thickness of the cathode catalyst layer [cm]
         *       set Cathode GDL thickness [cm]  = 2E-2         # Thickness of the cathode gas diffusion layer [cm]
         *       set Cathode MPL thickness [cm] = 2E-3          # Thickness of the cathode microporous layer [cm]
         *      end
         *   end
         * end
         * @endcode
         * 
         * The model can also be used to solve an agglomerate catalyst layer model. The governing equations are similar to the
         * ones outlined above, however, the volumetric current density source, i.e. \f$ \nabla \cdot i \f$ is obtained as
         * specified in the following article:
         * -  M. Secanell, K. Karan, A. Suleman and N. Djilali, "Multi-Variable Optimization of PEMFC Cathodes using 
         * an Agglomerate Model ", Electrochimica Acta, 52(22):6318-6337, 2007.
         * 
         * A sample dataset used to solve an agglomerate cathode catalyst layer model is given in 
         * data/cathode/Secanell_EA07_MultiVariable_Optimization_PEMFC_Cathodes_Agglomerate_Model
         * 
         * If using this function please cite the articles above as references.
         * 
         * @author Marc Secanell and Valentin N. Zingan, 2009-2014
         */
        
        template<int dim>
        class AppCathode : public OptimizationBlockMatrixApplication<dim>
        {
        public:
            
            ///@name Constructors, destructor, and initialization
            //@{
            
            /**
             * Constructor.
             */
            AppCathode( boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = 
            boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >() );
            
            /**
             * Destructor.
             */
            ~AppCathode();
            
            /**
             * Declare parameters.
             */
            virtual void declare_parameters(ParameterHandler& param);
            
            /**
             * Initialize parameters.
             */
            virtual void initialize(ParameterHandler& param);
            
            /**
             * The nonlinear solution initial guess
             * along with the appropriate BCs
             * is formed here.
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
             * Returns the current density
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
            virtual void cell_responses(std::vector<double>&                          dst,
                                        const typename DoFApplication<dim>::CellInfo& cell_info,
                                        const FEVector&);
            /**
             * 
             */
            void global_responses(std::vector<double>& resp,
                                  const FuelCell::ApplicationCore::FEVector& /*src*/);

            //@}
            
        protected:
            ///@name Pre-processor object
            //@{
            /**
             * Grid.
             */
            boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid;
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
            FuelCellShop::Material::Oxygen solute;
            FuelCellShop::Material::Hydrogen solute_anode;
            /**
             * Solvent.
             */
            FuelCellShop::Material::Nitrogen solvent;
            //@}
            
            ///@name Layer objects
            //@{
            /**
             * Cathode GDL.
             */
            boost::shared_ptr< FuelCellShop::Layer::GasDiffusionLayer<dim> > CGDL;
            
            /**
             * Cathode MPL.
             */
            boost::shared_ptr< FuelCellShop::Layer::MicroPorousLayer<dim> > CMPL;
            
            /**
             * Cathode CL.
             */
            boost::shared_ptr< FuelCellShop::Layer::CatalystLayer<dim> > CCL;
            
            /**
             * Run an anode instead of a cathode:
             */
            bool anode;
            //@}
            
            ///@name Equation objects
            //@{
            /**
             * This object describes
             * the equations that we are going to
             * solve here.
             */
            FuelCellShop::Equation::FicksTransportEquation<dim>* ficks_transport_equation;
            FuelCellShop::Equation::FicksTransportEquation<dim> ficks_transport_equation_cathode;
            FuelCellShop::Equation::FicksTransportEquation<dim> ficks_transport_equation_anode;
            
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
            FuelCellShop::Equation::ProtonTransportEquation<dim> proton_transport_equation;
            
            /**
             * The reaction source terms for
             * all underlying equations.
             */
            FuelCellShop::Equation::ReactionSourceTerms<dim> reaction_source_terms;
            //@}
            
            ///@name Post-processing objects (Functional evaluation)
            //@{
            FuelCellShop::PostProcessing::ORRCurrentDensityResponse<dim> ORRCurrent;
            FuelCellShop::PostProcessing::HORCurrentDensityResponse<dim> HORCurrent;
            //@}
            
            

        private:
            /**
             * Compute some functionals that are not needed for most applications (this section is not necessary in most cases.)
             */
            virtual void cell_responses_aux(std::vector<double>&                          dst,
                                            const typename DoFApplication<dim>::CellInfo& cell_info,
                                            const FEVector&);
            
            /**
             * Function to modify the default values of the data file in order to make sure that the equations
             * match those needed in the application.
             * 
             */
            void set_default_parameters_for_application(ParameterHandler &param)
            {
                param.enter_subsection("System management");
                {
                    param.set("Number of solution variables","3");
                    param.enter_subsection("Solution variables");
                    {
                        param.set("Solution variable 1","oxygen_molar_fraction");
                        param.set("Solution variable 2","protonic_electrical_potential");
                        param.set("Solution variable 3","electronic_electrical_potential");
                    }
                    param.leave_subsection();
                    
                    param.enter_subsection("Equations");
                    {
                        param.set("Equation 1","Ficks Transport Equation - oxygen");
                        param.set("Equation 2","Proton Transport Equation");
                        param.set("Equation 3","Electron Transport Equation");
                    }
                    param.leave_subsection();
                }
                param.leave_subsection();
                param.enter_subsection("Discretization");
                {
                    param.set("Element","FESystem[FE_Q(1)^3]");
                }
                param.leave_subsection();
            }
            
        };
        
    } // Application
    
} // FuelCell

#endif
