// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: electron_transport_equation.h
// - Description: This class describes steady-state electron transport using Ohm's law
// - Developers: Marc Secanell, Madhur Bhaiya, Valentin N. Zingan
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_EQUATION_ELECTRON_TRANSPORT_EQUATION_H_
#define _FCST_FUELCELLSHOP_EQUATION_ELECTRON_TRANSPORT_EQUATION_H_

// FCST includes
#include <equations/equation_base.h>
#include <layers/channel.h>
#include <layers/gas_diffusion_layer.h>
#include <layers/micro_porous_layer.h>
#include <layers/catalyst_layer.h>
#include <layers/solid_layer.h>

// STD
#include <sstream>
#include <string>

namespace FuelCellShop
{
    namespace Equation
    {
        /**
         * This class deals with <b>Electron Transport Equation</b>. This equation solves a steady-state <b>Ohm's law</b> for electron transport inside the layers.
         *
         * It is solved with respect to:
         * - \f$ \mathbf{\phi_{s}} \f$ \p(electronic_electrical_potential \p)
         *
         * where, \f$ \mathbf{\phi_{s}} \f$ is in Voltages.
         *
         * This equation can be written as:
         *
         * \f$ \qquad \mathbf{\nabla} \cdot \left( \hat{\sigma}_{s,eff} \cdot \mathbf{\nabla} \phi_s \right) = 0 \quad \in \quad \Omega \f$
         *
         * - \f$ \hat{\sigma}_{s,eff} \f$ is effective electron conductivity tensor [\p S/cm].
         *
         * To be well-posed, these equations are equipped with the appropriate boundary conditions. Dirichlet boundary conditions can be mentioned in the
         * subsection <b>"Boundary data"</b>. However, some other boundary types require additional information, which can also be provided by the parameter file.
         * We consider several types of boundary conditions:
         *
         * - <b>Galvanostatic (Constant Electron Current Flux):</b> A constant electronic current flux, \f$ C ~ \f$ [\p A/cm^2] is provided at the boundary. This value is provided using
         * parameter file. If the value is positive, it implies that electronic current flux is leaving out of the boundary, otherwise negative implies entering into the boundary.
         *
         * \f$ \qquad \mathbf{n} \cdot \left( \hat{\sigma}_{s,eff} \cdot \mathbf{\nabla} \phi_s  \right) = C \f$
         *
         * These are specified in the parameter file under subsection <b>"Boundary conditions"</b>, as a list of comma-separated \p map-like values.
         *
         * \em e.g. @code set Constant Electron Current Flux Boundary Conditions = 1: 2 , 4: -0.5 @endcode
         * where, \p boundary_id \b "1" has constant electron current flux value of \b 2 [\p A/cm^2] (leaving out of the boundary) and \p boundary_id \b "4"
         * has constant electron current flux value of \b -0.5 [\p A/cm^2] (entering into the boundary).
         *
         * - <b>No current flux (Insulated) / Symmetric:</b> A particular case of \p Neumann boundary condition.
         *
         * \f$ \qquad \mathbf{n} \cdot \left( \hat{\sigma}_{s,eff} \cdot \mathbf{\nabla} \phi_s \right) = 0 \f$
         *
         * \remarks
         * - There is no provision to specify boundary indicators for \p No \p current \p flux or \p Symmetric boundary conditions, as FEM
         * formulation automatically implies a particular boundary is one of these cases, by default.
         * - This class works with the following layer classes only:
         *      - FuelCellShop::Layer::Channel<dim>
         *      - FuelCellShop::Layer::GasDiffusionLayer<dim>
         *      - FuelCellShop::Layer::MicroPorousLayer<dim>
         *      - FuelCellShop::Layer::CatalystLayer<dim>
         *      - FuelCellShop::Layer::SolidLayer<dim>
         *
         * We solve the whole problem by linearizing the governing equation at each Newton iteration with subsequent
         * CG FEM discretization in space. The class contains the necessary class members to add the necessary contributions to cell_matrix
         * and cell_residual to the governing equations used to analyze electron transport by Ohm's law,
         *
         * <h3>Usage Details:</h3>
         *
         * @code
         * // Creating Equation object (in Application Header file)
         * FuelCellShop::Equation::ElectronTransportEquation<dim> electron_transport;
         *
         * // Declare parameters in application
         * electron_transport.declare_parameters(param);
         *
         * // Initialize in application
         * electron_transport.initialize(param);
         *
         * // Create a temporary vector in the application for storing couplings_map from all the equation used in the application.
         * std::vector<couplings_map> tmp;
         * ... // other equations
         * tmp.push_back( electron_transport.get_internal_cell_couplings() );
         *
         * // Look at ReactionSourceTerms class here, if source terms due to electrochemical reaction are to be considered.
         * // Making cell couplings using SystemManagement object created in the application
         * system_management.make_cell_couplings(tmp);
         *
         * // cell_matrix in application
         * // Do a check against layer and it should match with the layers currently working for this equation class.
         * // for eg: CGDL is FuelCellShop::Layer::DesignFibrousGDL<dim> object.
         * electron_transport.assemble_cell_matrix(cell_matrices, cell_info, &CGDL);
         *
         * // cell_residual in application
         * electron_transport.assemble_cell_residual(cell_vector, cell_info, &CGDL);
         * @endcode
         *
         *
         * \note This class doesn't assemble for reaction source terms; that is taken care off by \p ReactionSourceTerms class.
         * Please read the documentation of ReactionSourceTerms class, for additional methods to be implemented in the application.
         *
         * \warning If current source terms due to reaction are being considered, it's very important to use \p adjust_internal_cell_couplings
         * member function of ReactionSourceTerms class, before using  \p make_cell_couplings of \b SystemManagement at the application level.
         *
         * TODO Old Boundary conditions including Dirichlet Boundary Conditions and Constant Electron Current Flux Boundary Conditions are supposed to be replaced with the new subsections,
         *      see TO BE REMOVED comments in .cc file.
         *
         * \author Marc Secanell,      2013
         * \author Madhur Bhaiya,      2013
         * \author Valentin N. Zingan, 2012-2014 - afterward improvements, optimization, checkings, CG FEM bug fixings
         */

        template<int dim>
        class ElectronTransportEquation : public EquationBase<dim>
        {
        public:

            ///@name Constructors, destructor, and initalization
            //@{

            /**
             * Constructor.
             */
            ElectronTransportEquation(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = 
            boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >());

            /**
             * Destructor.
             */
            virtual ~ElectronTransportEquation();

            /**
             * Declare parameters.
             */
            virtual void declare_parameters(ParameterHandler& param) const;

            /**
            * Set parameters using the parameter file, in order to run parametric/optimization studies.
            *
            * Current list of available parameters:
            *   - \b current_bdryid_i : This is used to specify electronic current flux [\p A/cm^2] boundary condition at the specified boundary id <b>"i"</b>. <b>Please note</b> that \p boundary_id <b>"i"</b>
            *    should be an unsigned int number. For eg: In order to run parameteric study on boundary id \b 4, the argument in the parameter file should be \b current_bdryid_4.
            */
            virtual void set_parameters(const std::vector<std::string>& name_dvar,
                                        const std::vector<double>& value_dvar,
                                        ParameterHandler& param);

            /**
             * Initialize parameters.
             */
            virtual void initialize(ParameterHandler& param);

            //@}
            ///@name Local CG FEM based assemblers

            //@{
            /**
             * Assemble local cell matrix.
             * Computes the contribution of Ohm's law to the cell_matrix:
             * \f[
             *     -\frac{\partial R}{\partial u} \delta u
             * \f]
             */
            virtual void assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                              const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                              FuelCellShop::Layer::BaseLayer<dim>* const              layer);

            /**
             * Assemble local cell residual.
             */
            virtual void assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_rhs,
                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                FuelCellShop::Layer::BaseLayer<dim>* const              layer);

            /**
             * Assemble local boundary matrix.
             */
            virtual void assemble_bdry_matrix(FuelCell::ApplicationCore::MatrixVector&                                 bdry_matrices,
                                              const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                              FuelCellShop::Layer::BaseLayer<dim>* const              layer);

            /**
             * Assemble local boundary residual.
             */
            virtual void assemble_bdry_residual(FuelCell::ApplicationCore::FEVector&                                     bdry_rhs,
                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                FuelCellShop::Layer::BaseLayer<dim>* const              layer);

            //@}

            ///@name Accessors and info
            //@{

            /**
             * The function printing out
             * the equations info.
             */
            virtual void print_equation_info() const;

            /**
             * This member function creates an object of its own type and runs test to diagnose if there are any problems
             * with the routines that you have created.
             */
            void class_test();
            //@}

        protected:
            ///@name Boundary conditions
             //@{

             /**
              * \p std::map< \p unsigned \p int, \p double \p> container for details regarding <b>Galvanostatic</b> boundary conditions.
              * Here, \p Key (unsigned int) represents the \p boundary_id and \p Value (double) represents the constant electronic current flux values [\p A/\p cm^2].
              */
             std::map<unsigned int, double> electron_current_flux_map;


            ///@name Local CG FEM based assemblers - make_ functions
            //@{

            /**
             * This function computes Local CG FEM based
             * assemblers - constant data (generic).
             */
            virtual void make_assemblers_generic_constant_data();

            /**
             * This function computes
             * <b> Local CG FEM based assemblers - constant data (cell) </b>
             * and allocates the memory for
             * \p shape \p functions, \p shape \p function \p gradients, and
             * #JxW_cell in
             * <b> Local CG FEM based assemblers - variable data (cell) </b>.
             */
            virtual void make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info);

            /**
             * This function computes
             * <b> Local CG FEM based assemblers - constant data (boundary) </b>
             * and allocates the memory for
             * \p shape \p functions, #normal_vectors, and
             * #JxW_bdry in
             * <b> Local CG FEM based assemblers - variable data (boundary) </b>.
             */
            virtual void make_assemblers_bdry_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info);

            /**
             * This function computes
             * <b> Local CG FEM based assemblers - variable data (cell) </b>.
             */
            virtual void make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                            FuelCellShop::Layer::BaseLayer<dim>* const layer);

            /**
             * This function computes
             * <b> Local CG FEM based assemblers - variable data (boundary) </b>.
             */
            virtual void make_assemblers_bdry_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                            FuelCellShop::Layer::BaseLayer<dim>* const layer);

            //@}

            ///@name Other make_ functions
            //@{

            /**
             * This function fills out
             * \p internal_cell_couplings.
             */
            virtual void make_internal_cell_couplings();

            //@}

            ///@name Generic Constant Data
            //@{

            /**
             * VariableInfo structure corresponding to \p "electronic_electrical_potential".
             */
            VariableInfo phi_s;

            //@}

            ///@name Local CG FEM based assemblers - variable data (cell)
            //@{

            /**
             * Effective electronic conductivity,
             * [\p S/cm], of the cell.
             */
            Tensor<2,dim> sigmaSeff_cell;

            /**
             * \f$ \mathbf{\phi_s} \f$ shape function gradients.
             *
             * \p grad_phi_phiS_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{\phi_s} \f$ shape function gradient
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< Tensor<1,dim> > > grad_phi_phiS_cell;

            //@}

            ///@name Local CG FEM based assemblers - variable data (boundary)
            //@{

            /**
             * \f$ \mathbf{\phi_s} \f$ shape functions.
             *
             * \p phi_phiS_bdry \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{\phi_s} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the boundary.
             */
            std::vector< std::vector<double> > phi_phiS_bdry;

            //@}

            /**
             * Variable used to store the index in cell_info->global_data of the previous Newton solution
             * The solution at the previous iteration is used to compute cell_matrix and cell_residual
             */
            unsigned int last_iter_cell;

            /**
             * Variable used to store the index in bdry_info->global_data of the previous Newton solution
             * The solution at the previous iteration is used to compute bdry_matrix and bdry_residual
             */
            unsigned int last_iter_bdry;

        };

    } // Equation

} // FuelCellShop

#endif
