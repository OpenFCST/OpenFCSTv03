//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: proton_transport_equation.h
//    - Description: Equation class for Proton Transport (using Ohm's law)
//    - Developers: Madhur Bhaiya, M. Secanell, Valentin N. Zingan
//
//---------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_EQUATION_PROTON_TRANSPORT_EQUATION_H_
#define _FCST_FUELCELLSHOP_EQUATION_PROTON_TRANSPORT_EQUATION_H_

#include <equations/equation_base.h>

#include <layers/catalyst_layer.h>
#include <layers/membrane_layer.h>

// STD
#include <sstream>
#include <string>

namespace FuelCellShop
{
    namespace Equation
    {
        /**
         * This class deals with <b>Proton Transport Equation</b>.
         *
         * This equation solves a steady-state <b>Ohm's law</b> for proton transport inside the polymer electrolyte. Note that <b>Ohm's law</b> can be used
         * in the polymer electrolyte because electro-neutrality is assumed due to the constant number of negative ionic charges that are attached to the electrolyte.
         * Due to electro-neutrality, convective and diffusive current terms in the \p Nernst-Planck \p Equation become negligible, hence it simplifies to <b>Ohm's Law</b>.
         *
         * It is solved with respect to:
         * - \f$ \mathbf{\phi_{m}} \f$ \p(protonic_electrical_potential \p)
         *
         * where, \f$ \mathbf{\phi_{m}} \f$ is in Voltages.
         *
         * This equation can be written as:
         *
         * \f$ \qquad \mathbf{\nabla} \cdot \left( \hat{\sigma}_{m,eff} \cdot \mathbf{\nabla} \phi_m \right) = 0 \quad \in \quad \Omega \f$
         *
         * - \f$ \hat{\sigma}_{m,eff} \f$ is effective proton conductivity tensor [\p S/cm], which can be function of other variables,
         *  \em viz., \f$ \lambda \f$ \p(membrane_water_content \p) and \f$ T \f$ \p(temperature_of_REV \p).
         *
         * To be well-posed, these equations are equipped with the appropriate boundary conditions. All the boundary conditions can be described by
         * \p boundary_id \p(s \p) and \p boundary_type. Besides, this some boundary types require additional information, which can also be provided by the parameter file.
         * We consider several types of boundary conditions:
         *
         * - \b Dirichlet: At this boundary, we specify the \f$ \mathbf{\phi_{m}} \f$ values using parameter file.  It is to be noted that these
         * dirichlet values at the  boundary should be supplied to Initial Solution methods used in a particular application.
         * These are specified in the parameter file under subsection <b>"Dirichlet Boundary Conditions"</b>, as a list of comma-separated \p map-like values.
         *
         *  \em e.g. @code set Dirichlet Boundary Conditions = 2: 0.5, 4: -0.7 @endcode
         *
         * - <b>Galvanostatic (Constant Proton Current Flux):</b> A constant protonic current flux, \f$ C ~ \f$ [\p A/cm^2] is provided at the boundary. This value is provided using
         * parameter file. If the value is positive, it implies that protonic current flux is leaving out of the boundary, otherwise negative implies entering into the boundary.
         *
         * \f$ \qquad - \mathbf{n} \cdot \left( \hat{\sigma}_{m,eff} \cdot \mathbf{\nabla} \phi_m  \right) = C \f$
         *
         * These are specified in the parameter file under subsection <b>"Constant Proton Current Flux Boundary Conditions"</b>, as a list of comma-separated \p map-like values.
         *
         * \em e.g. @code set Constant Proton Current Flux Boundary Conditions = 1: 2 , 4: -0.5 @endcode
         * where, \p boundary_id \b "1" has constant proton current flux value of \b 2 [\p A/cm^2] (leaving out of the boundary) and \p boundary_id \b "4"
         * has constant proton current flux value of \b -0.5 [\p A/cm^2] (entering into the boundary).
         *
         * - <b>No current flux (Insulated) / Symmetric:</b> A particular case of \p Neumann boundary condition.
         *
         * \f$ \qquad - \mathbf{n} \cdot \left( \hat{\sigma}_{m,eff} \cdot \mathbf{\nabla} \phi_m \right) = 0 \f$
         *
         * \remarks
         * - There is no provision to specify boundary indicators for \p No \p current \p flux or \p Symmetric boundary conditions, as FEM
         * formulation automatically implies a particular boundary is one of these cases, by default.
         * - This class works with the following layer classes only:
         *      - FuelCellShop::Layer::CatalystLayer<dim>
         *      - FuelCellShop::Layer::MembraneLayer<dim>
         *
         * We solve the whole problem by linearizing the governing equation at each Newton iteration with subsequent
         * CG FEM discretization in space. The class contains the necessary class members to add the necessary contributions to cell_matrix
         * and cell_residual to the governing equations used to analyze proton transport by Ohm's law,
         *
         * <h3>Usage Details:</h3>
         *
         * @code
         * // Creating Equation object (in Application Header file)
         * FuelCellShop::Equation::ProtonTransportEquation<dim> proton_transport;
         *
         * // Declare parameters in application
         * proton_transport.declare_parameters(param);
         *
         * // Initialize in application
         * proton_transport.initialize(param);
         *
         * // Create a temporary vector in the application for storing couplings_map from all the equation used in the application.
         * std::vector<couplings_map> tmp;
         * ... // other equations
         * tmp.push_back( proton_transport.get_internal_cell_couplings() );
         *
         * // Look at ReactionSourceTerms class here, if source terms due to electrochemical reaction are to be considered.
         * // Making cell couplings using SystemManagement object created in the application
         * system_management.make_cell_couplings(tmp);
         *
         * // cell_matrix in application
         * // Do a check against layer and it should match with the layers currently working for this equation class.
         * // for eg: CCL is FuelCellShop::Layer::HomogeneousCL<dim> object.
         * proton_transport.assemble_cell_matrix(cell_matrices, cell_info, &CCL);
         *
         * // cell_residual in application
         * proton_transport.assemble_cell_residual(cell_vector, cell_info, &CCL);
         * @endcode
         *
         *
         * \note This class doesn't assemble for reaction source terms; that is taken care off by \p ReactionSourceTerms class.
         * Please read the documentation of ReactionSourceTerms class, for additional methods to be implemented in the application.
         *
         * \warning If current source terms due to reaction are being considered, it's very important to use \p adjust_internal_cell_couplings
         * member function of ReactionSourceTerms class, before using  \p make_cell_couplings of \b SystemManagement at the application level.
         *
         * TODO Old Boundary conditions including Dirichlet Boundary Conditions and Constant Proton Current Flux Boundary Conditions are supposed to be replaced with the new subsections,
         *      see TO BE REMOVED comments in .cc file.
         *
         * \author Madhur Bhaiya,      2013
         * \author Valentin N. Zingan, 2012-2014 - afterward improvements, optimization, checkings, CG FEM bug fixings
         */

        template<int dim>
        class ProtonTransportEquation : public EquationBase<dim>
        {
        public:

            ///@name Constructors, destructor, and initalization
            //@{

            /**
             * Constructor.
             */
            ProtonTransportEquation(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = 
            boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >());

            /**
             * Destructor.
             */
            virtual ~ProtonTransportEquation();

            /**
             * Declare parameters.
             */
            virtual void declare_parameters(ParameterHandler& param) const;

            /**
             * Initialize parameters. This class will call #make_internal_cell_couplings and #make_boundary_types.
             */
            virtual void initialize(ParameterHandler& param);

            //@}
            ///@name Local CG FEM based assemblers

            //@{
            /**
             * Assemble local cell matrix.
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
             * The function prints out the equation's info.
             */
            virtual void print_equation_info() const;

            /**
             * Method to provide access to Dirichlet boundary conditions map filled using the parameter file.
             * This function is useful in accessing the dircihlet values to be set in initial conditions methods of the application.
             */
            inline std::map<unsigned int, double> get_dirichlet_bdry_map() const
            {
                return dirichlet_bdry_map;
            }

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
              * \p std::map< \p unsigned \p int, \p double \p> container for Dirichlet boundary conditions.
              * Here, \p Key \p (unsigned \p int \p) represents the \p boundary_id and \p Value (double) represents the potential \p[V] .
              */
             std::map<unsigned int, double> dirichlet_bdry_map;

             /**
              * \p std::map< \p unsigned \p int, \p double \p> container for details regarding <b>Galvanostatic</b> boundary conditions.
              * Here, \p Key (unsigned int) represents the \p boundary_id and \p Value (double) represents the constant protonic current flux values [\p A/\p cm^2].
              */
             std::map<unsigned int, double> proton_current_flux_map;

            //@}

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
             * \p JxW_cell in
             * <b> Local CG FEM based assemblers - variable data (cell) </b>.
             */
            virtual void make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info);

            /**
             * This function computes
             * <b> Local CG FEM based assemblers - constant data (boundary) </b>
             * and allocates the memory for
             * \p shape \p functions, #normal_vectors, and
             * \p JxW_bdry in
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

            /**
             * This function fills out
             * \p boundary_types.
             */
            virtual void make_boundary_types();

            /**
             * This function fills out
             * \p output_types.
             */
            virtual void make_output_types()
            {};

            //@}

            ///@name Generic Constant Data
            //@{
            /**
             * VariableInfo structure corresponding to \p "protonic_electrical_potential".
             */
            VariableInfo phi_m;

            /**
             * VariableInfo structure corresponding to \p "membrane_water_content".
             */
            VariableInfo lambda;

            /**
             * VariableInfo structure corresponding to \p "temperature_of_REV".
             */
            VariableInfo t_rev;

            //@}

            ///@name Local CG FEM based assemblers - variable data (cell)
            //@{

            /**
             * Effective proton conductivity,
             * [\p S/cm], at all quadrature points of the cell.
             */
            std::vector<double> sigmaMeff_cell;

            /**
             * Derivative of effective protonic conductivity
             * w.r.t. \p "membrane_water_content",
             * at all quadrature points in the cell.
             */
            std::vector<double> dsigmaMeff_dlambda_cell;

            /**
             * Derivative of effective protonic conductivity
             * w.r.t \p "temperature_of_REV",
             * at all quadrature points in the cell.
             */
            std::vector<double> dsigmaMeff_dT_cell;

            /**
             * \f$ \mathbf{\phi_m} \f$ shape function gradients.
             *
             * \p grad_phi_phiM_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{\phi_m} \f$ shape function gradient
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< Tensor<1,dim> > > grad_phi_phiM_cell;

            /**
             * \f$ \mathbf{\lambda} \f$ shape functions.
             *
             * \p phi_lambda_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{\lambda} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_lambda_cell;

            /**
             * \f$ \mathbf{T} \f$ shape functions.
             *
             * \p phi_T_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{T} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_T_cell;

            //@}

            ///@name Local CG FEM based assemblers - variable data (boundary)
            //@{

            /**
             * \f$ \mathbf{\phi_{M}} \f$ shape functions.
             *
             * \p phi_phiM_bdry \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{\phi_{M}} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the boundary.
             */
            std::vector< std::vector<double> > phi_phiM_bdry;

            //@}

            /**
             * Counter set to \em TRUE when \p cell_residual is being assembled.
             * This ensures that only effective protonic conductivity is being calculated, not its derivatives. (improves speed)
             */
            bool cell_residual_counter;

            /**
             * Counter set to \em TRUE when \p bdry_residual is being assembled.
             * This ensures that only data required for bdry_residual is computed. (improves speed)
             */
            bool bdry_residual_counter;

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