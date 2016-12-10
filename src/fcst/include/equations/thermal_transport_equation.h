//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: thermal_transport_equation.h 18-04-2013
//    - Description: Equation class for Thermal transport
//    - Developers: Madhur Bhaiya
//
//---------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_EQUATION_THERMAL_TRANSPORT_EQUATION_H_
#define _FCST_FUELCELLSHOP_EQUATION_THERMAL_TRANSPORT_EQUATION_H_

#include <utils/fcst_units.h>
#include <equations/equation_base.h>
#include <layers/gas_diffusion_layer.h>
#include <layers/micro_porous_layer.h>
#include <layers/catalyst_layer.h>
#include <layers/membrane_layer.h>

//STD
#include <string>
#include <sstream>


namespace FuelCellShop
{
    namespace Equation
    {
        ///@name Exceptions
        //@{
        /**
         * Exception thrown when ohmic heating (\em Protonic or \em Electronic) is set to \p TRUE,
         * but respective solution variables, \f$ \phi_m \f$ or \f$ \phi_s \f$, are not
         * declared as user solution variables.
         */
        DeclException2(VariableNotFoundForOhmicHeat,
                       std::string,
                       std::string,
                       << "For "  << arg1 << " ohmic heating source term set as True, \"" << arg2 << "\" is not found as one of the solution variables.");
        //@}

        /**
         * This class deals with <b>Thermal Transport Equation</b>.
         *
         * This equation solves a steady-state, heat transfer equation (<b>Fourier's law of conduction</b>), alongwith enthalpy transport due to species diffusion in multi-component mixtures. Heat transfer
         * via convection is generally very little (roughly \b 6 orders of magnitude smaller) as compared to conductive heat transport in a PEMFC, hence neglected here. Also currently, <b>Fick's law of diffusion</b>
         * is being considered to account for enthalpy transport due to inter-diffusion amongst species.
         *
         * - It is solved with respect to
         *   \f$ \mathbf{T} \f$ \p(temperature_of_REV\p)
         *
         * where, \f$ \mathbf{T} \f$ is in Kelvins.
         *
         * This equation can be written as:
         *
         * \f$ \qquad -\mathbf{\nabla} \cdot \left( \hat{k}_{eff} \cdot \mathbf{\nabla} T \right) + \mathbf{\nabla} \cdot \left( \bar{H}_i \vec{J}_i \right) = \mathbf{S_{thermal}} \quad \in \quad \Omega \f$
         *
         * where
         *
         * \f$
         *    \qquad
         *    \mathbf{S_{thermal}} =
         *    \begin{cases}
         *    \hat{\sigma}_{s,eff} \cdot \left( \mathbf{\nabla} \phi_s \otimes  \mathbf{\nabla} \phi_s \right) \quad - \quad \text{Electronic ohmic heating in all layers except Membrane layer} \\ ~ \\
         *    \hat{\sigma}_{m,eff} \cdot \left( \mathbf{\nabla} \phi_m \otimes  \mathbf{\nabla} \phi_m \right) \quad - \quad \text{Protonic ohmic heating in Catalyst layers and Membrane layer}
         *    \end{cases}
         * \f$
         *
         * - \f$ \hat{k}_{eff} \f$ is Thermal conductivity tensor [\p W/\p(cm-K\p)], which can be a function of other variables, \em e.g. \f$ T \f$ (Temperature).
         * - \f$ \hat{\sigma}_{s,eff} \f$ is Electronic conductivity tensor [\p S/cm].
         * - \f$ \hat{\sigma}_{m,eff} \f$ is Protonic conductivity tensor [\p S/cm].
         * - \f$ \mathbf{S_{thermal}} \f$ is in [\p W/\p(cm^3\p)].
         * - \f$ \bar{H}_i \f$ is molar enthalpy [\p J/\p(mol\p)] of the diffusing specie \p "i".
         * - \f$ \vec{J}_i \f$ is molar diffusion flux [\p mol/\p(cm^2-s\p)] of the diffusing specie \p "i", for gaseous species it is given as:  \f$ \vec{J}_i = -\hat{D}_{i,eff} c_T \mathbf{\nabla} x_i \f$, where:
         *   - \f$ \hat{D}_{i,eff} \f$ is effective diffusivity tensor [\p cm^2/s] of the diffusing gas \p "i".
         *   - \f$ c_T \f$ is total concentration [\p mol/cm^3] of all the gaseous species involved in the cell.
         *   - \f$ x_i \f$ is molar fraction of the diffusing gas \p "i".
         * - Diffusing specie \p "i" can also be sorbed water, \f$ \lambda \f$.
         *
         * To be well-posed, these equations are equipped with the appropriate boundary conditions. Dirichlet boundary conditions can be mentioned in the
         * subsection <b>"Boundary data"</b>. However, some other boundary types require additional information, which can also be provided by the parameter file.
         * We consider several types of boundary conditions:
         *
         * - <b>Constant Heat Flux:</b> A constant heat flux, \f$ C ~ \f$ [\p W/\p cm^2] is provided at a boundary. This value is provided using
         * parameter file. If the value is positive, it implies that heat flux is leaving out of the boundary, otherwise negative implies entering into the boundary.
         *
         * \f$ \qquad -\mathbf{n} \cdot \left( \hat{k}_{eff} \cdot \mathbf{\nabla} T \right) = C \f$
         *
         * These are specified in the parameter file under subsection <b>"Boundary conditions"</b>, as a list of comma-separated \p map-like values.
         *
         * \em e.g. @code set Constant Heat Flux Boundary Conditions = 1: 34.5 , 4: -23.5 @endcode
         * where, \p boundary_id \b "1" has constant heat flux value of \b 34.5 [\p W/\p cm^2] (coming out of the boundary) and \p boundary_id \b "4"
         * has constant heat flux value of \b -23.5 [\p W/\p cm^2] (going into the boundary).
         *
         * - <b>Convective Heat Flux:</b> Convective heat flux based boundary condition is applied at the boundary. It requires \p heat \p transfer \p coefficient,
         * \f$ h \f$ [\p W/\p(cm^2-K \p)] and \p ambient \p temperature, \f$ T_{\infty} \f$ [\p K] at the boundary.
         *
         * \f$ \qquad -\mathbf{n} \cdot \left( \hat{k}_{eff} \cdot \mathbf{\nabla} T \right) = h(T - T_{\infty}) \f$
         *
         * These are specified in the parameter file under subsection <b>"Boundary conditions"</b>, as a list of comma-separated \p map-like values.
         *
         * \em e.g. @code set Convective Heat Flux Boundary Conditions = 1: 1.5;340, 3: 0.4;230 @endcode
         * where, \p boundary_id \b "1" has convective heat transfer coefficient of \b 1.5 [\p W/\p(cm^2-K \p)] and ambient temperature of \b 340 [\p K]; and
         * \p boundary_id \b "3" has convective heat transfer coefficient of \b 0.4 [\p K] and ambient temperature of \b 230 [\p K]. <b>Please note</b> that
         * <b>";"</b> is necessary to provide both of the required parameters, \em viz., \f$ h \f$ and \f$ T_{\infty} \f$ at a particular boundary.
         *
         * - <b>No Heat Flux (Insulated) / Symmetric:</b> A particular case of \p Neumann boundary condition.
         *
         * \f$ \qquad -\mathbf{n} \cdot \left( \hat{k}_{eff} \cdot \mathbf{\nabla} T \right) = 0 \f$
         *
         * \remarks
         * - \b Ohmic heating source terms (\p electronic and \p protonic) can be turned \em ON or \em OFF, using flags defined under subsection <b>"Boolean flags</b> in the parameter file.
         * - In order to have electronic ohmic heating source term, \f$ \phi_s \f$ [\p Volts ], should be one of the solution variables,
         * otherwise this class will throw an exception. Similarly, for protonic ohmic heating, \f$ \phi_m \f$ [ \p Volts ], should be one of the solution variables.
         * - Enthalpy tranport due to Fickian diffusion of species is defaulted to \em OFF. It can be turned \em ON by setting following paramter entry to \p TRUE, under subsection <b>"Boolean flags</b> :
         * @code set Enthalpy transport due to fickian diffusion of gases = true @endcode
         * - <b>Please note</b> that while considering enthalpy transport due to fickian diffusion, it is very necessary to set the \b solvent gas for a particular layer as the last entry in the input vector to
         * FuelCellShop::Layer::PorousLayer::set_gases method, in the initialization section of the application.
         * - Enthalpy transport due to transport of sorbed water is defaulted to \em OFF. It can be turned \em ON by setting following parameter entry to \p TRUE, under subsection <b>"Boolean flags</b> :
         * @code set Enthalpy transport associated with lambda transport = true @endcode
         * - <b>Please note</b> that while considering enthalpy transport due to lambda transport, it is very necessary to include \p membrane_water_content, \f$ \lambda \f$ as one of
         * the solution variables in the application.
         * - Lambda transport accounts for \b diffusion, \b electro-osmosis if \p protonic_electrical_potential, \f$ \phi_m \f$ is one of the solution variables and \b thermo-osmosis if
         * it is set to \p TRUE in parameter entry for <b>Membrane Water Content Transport Equation</b>.
         * - There is no provision to specify boundary indicators for \p No \p Heat \p Flux or \p Symmetric boundary conditions, as FEM
         * formulation automatically implies a particular boundary is one of these cases, by default.
         * - This class currently works with the following layer classes:
         *      - FuelCellShop::Layer::GasDiffusionLayer<dim>
         *      - FuelCellShop::Layer::MicroPorousLayer<dim>
         *      - FuelCellShop::Layer::CatalystLayer<dim>
         *      - FuelCellShop::Layer::MembraneLayer<dim>
         *
         * We solve the whole problem by linearizing the governing equation at each Newton iteration with subsequent
         * CG FEM discretization in space.
         *
         * <h3>Usage Details:</h3>
         *
         * @code
         * // Creating Equation object (in Application Header file)
         * FuelCellShop::Equation::ThermalTransportEquation<dim> thermal_transport;
         *
         * // Declare parameters in application
         * thermal_transport.declare_parameters(param);
         *
         * // Initialize in application
         * thermal_transport.initialize(param);
         *
         * // Create a temporary vector in the application for storing couplings_map from all the equation used in the application.
         * std::vector<couplings_map> tmp;
         * ... // other equations
         * tmp.push_back( thermal_transport.get_internal_cell_couplings() );
         *
         * // Look at ReactionSourceTerms class here, if source terms due to electrochemical reaction are to be considered.
         * // Making cell couplings using SystemManagement object created in the application
         * system_management.make_cell_couplings(tmp);
         *
         * // cell_matrix in application
         * // Do a check against layer and it should match with the layers currently working for this equation class.
         * // for eg: CGDL is FuelCellShop::Layer::DesignFibrousGDL<dim> object.
         * thermal_transport.assemble_cell_matrix(cell_matrices, cell_info, &CGDL);
         *
         * // cell_residual in application
         * thermal_transport.assemble_cell_residual(cell_vector, cell_info, &CGDL);
         *
         * // For bdry_matrices and bdry_vector; check for the layer present at that boundary. For eg: CGDL
         * thermal_transport.assemble_bdry_matrix(bdry_matrices, bdry_info, &CGDL);
         * thermal_transport.assemble_bdry_residual(bdry_vector, bdry_info, &CGDL);
         * @endcode
         *
         *
         * \note This class doesn't assemble for reaction heat source terms; that is taken care off by ReactionSourceTerms class.
         * Please read the documentation of ReactionSourceTerms class, for additional methods to be implemented in the application.
         *
         * \warning If heat source terms due to reaction are being considered, it's very important to use ReactionSourceTerms::adjust_internal_cell_couplings
         * method, before using FuelCell::SystemManagement::make_cell_couplings at the application level.
         *
         * \author Madhur Bhaiya, 2013
         */

        template<int dim>
        class ThermalTransportEquation : public EquationBase<dim>
        {
        public:

            ///@name Constructors, destructor, and initalization
            //@{

            /**
             * Constructor.
             */
            ThermalTransportEquation(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = 
            boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >());

            /**
             * Destructor.
             */
            virtual ~ThermalTransportEquation();

            /**
             * Declare parameters.
             */
            virtual void declare_parameters(ParameterHandler& param) const;

            /**
             * Initialize parameters.
             *
             * This class will call  #make_internal_cell_couplings.
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
            virtual void assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
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
            virtual void assemble_bdry_residual(FuelCell::ApplicationCore::FEVector&                                     bdry_residual,
                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                FuelCellShop::Layer::BaseLayer<dim>* const              layer);

            //@}

            ///@name Accessors and info
            //@{

            /**
             * This function prints out
             * the equations info.
             */
            virtual void print_equation_info() const;

            /**
             * Inline method to get whether the electronic ohmic heating in Gas diffusion layer is \p ON or \p OFF (#electron_ohmic_heat_gdl).
             * The typical use of this method is for the \p \b post-processing routines in the application.
             */
            inline bool get_electron_ohmic_heat_gdl() const
            {	return electron_ohmic_heat_gdl;}

            /**
             * Inline method to get whether the electronic ohmic heating in Microporous layer is \p ON or \p OFF (#electron_ohmic_heat_mpl).
             * The typical use of this method is for the \p \b post-processing routines in the application.
             */
            inline bool get_electron_ohmic_heat_mpl() const
            {	return electron_ohmic_heat_mpl;}

            /**
             * Inline method to get whether the electronic ohmic heating in Catalyst layer is \p ON or \p OFF (#electron_ohmic_heat_cl).
             * The typical use of this method is for the \p \b post-processing routines in the application.
             */
            inline bool get_electron_ohmic_heat_cl() const
            {	return electron_ohmic_heat_cl;}

            /**
             * Inline method to get whether the protonic ohmic heating in Catalyst layer is \p ON or \p OFF (#proton_ohmic_heat_cl).
             * The typical use of this method is for the \p \b post-processing routines in the application.
             */
            inline bool get_proton_ohmic_heat_cl() const
            {	return proton_ohmic_heat_cl;}

            /**
             * Inline method to get whether the protonic ohmic heating in Membrane layer is \p ON or \p OFF (#proton_ohmic_heat_ml).
             * The typical use of this method is for the \p \b post-processing routines in the application.
             */
            inline bool get_proton_ohmic_heat_ml() const
            {	return proton_ohmic_heat_ml;}

            //@}

        protected:

            ///@name Boolean flags
            //@{

            /**
             * This boolean data member indicates that
             * the electronic ohmic heating in Gas diffusion layer is \p ON or \p OFF.
             */
            bool electron_ohmic_heat_gdl;

            /**
             * This boolean data member indicates that
             * the electronic ohmic heating in Microporous layer is \p ON or \p OFF.
             */
            bool electron_ohmic_heat_mpl;

            /**
             * This boolean data member indicates that
             * the electronic ohmic heating in Catalyst layer is \p ON or \p OFF.
             */
            bool electron_ohmic_heat_cl;

            /**
             * This boolean data member indicates that
             * the protonic ohmic heating in Catalyst layer is \p ON or \p OFF.
             */
            bool proton_ohmic_heat_cl;

            /**
             * This boolean data member indicates that
             * the protonic ohmic heating in Membrane layer is \p ON or \p OFF.
             */
            bool proton_ohmic_heat_ml;

            /**
             * This boolean data member indicates that the enthalpy transport via Fickian diffusion is \p ON or \p OFF.
             */
            bool enthalpy_fickian_transport;

            /**
             * This boolean data member indicates that the enthalpy transport associated with lambda (sorbed water) transport is \p ON or \p OFF.
             */
            bool enthalpy_lambda_transport;

            /**
             * This flag indicates that lambda transport due to thermo-osmotic diffusion is \p ON or \p OFF.
             */
            bool flag_thermoosmosis;

            //@}

            ///@name Boundary conditions
            //@{

            /**
             * \p std::map< \p unsigned \p int, \p double \p> container for details regarding <b>Constant Heat Flux</b> boundaries.
             * Here, \p Key (unsigned int) represents the \p boundary_id and \p Value (double) represents the constant heat flux values [\p W/\p cm^2].
             * \b Positive heat flux represents heat leaving out of the boundary and \b Negative coming into the boundary.
             */
            std::map<unsigned int, double> const_heat_flux_map;

            /**
             * Container for details regarding <b>Convective Heat Flux</b> boundaries. Here, \p Key represents the \p boundary_id, and
             * \p Value vector stores the heat transfer coefficient [\p W/\p(cm^2-K \p)] and ambient temperature [\p K] respectively.
             */
            std::map<unsigned int, std::vector<double>> conv_heat_flux_map;

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


            /**
             * This function computes <b>Enthalpy transport assemblers - variable data (cell)</b> corresponding to enthalpy transport via Fickian diffusion.
             */
            template <typename porelayer>
            void gas_enthalpy_transport_assemblers_cell_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                             porelayer* const layer);

            /**
             * This function computes <b>Enthalpy transport assemblers - variable data (cell)</b> corresponding to enthalpy transport associated with lambda transport.
             */
            template <typename polymer>
            void lambda_enthalpy_transport_assemblers_cell_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                const polymer* const layer);

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
             * VariableInfo structure corresponding to base variable
             * of this equation class, \p "temperature_of_REV".
             */
            VariableInfo t_rev;

            /**
             * VariableInfo structure corresponding to \p "electronic_electrical_potential".
             */
            VariableInfo phi_s;

            /**
             * VariableInfo structure corresponding to \p "protonic_electrical_potential".
             */
            VariableInfo phi_m;

            /**
             * VariableInfo structure corresponding to \p "membrane_water_content".
             */
            VariableInfo lambda;

            /**
             * VariableInfo structure corresponding to \p "liquid_water_saturation".
             */
            VariableInfo s_liquid_water;
            
            VariableInfo p_liquid_water;

            /**
             * Map of VariableInfo structures of various gaseous species, whose diffusion is being considered in the application. For instance, in case of oxygen gas: String \p Key is \b "oxygen_molar_fraction",
             * while \p Value is VariableInfo structure corresponding to \p "oxygen_molar_fraction".
             */
            std::map< std::string, VariableInfo > xi_map;

            //@}

            ///@name Local CG FEM based assemblers - variable data (cell)
            //@{

            /**
             * Effective thermal conductivity, [\p W/\p(cm-k \p)],
             * at all quadrature points of the cell.
             */
            std::vector< Tensor<2,dim> > keff_cell;

            /**
             * Derivative of effective thermal conductivity
             * w.r.t \p "temperature_of_REV",
             * at all quadrature points of the cell.
             */
            std::vector< Tensor<2,dim> > dkeff_dT_cell;

            /**
             * Effective electronic conductivity,
             * [\p S/cm], of the cell.
             */
            Tensor<2,dim> sigmaSeff_cell;

            /**
             * Effective protonic conductivity,
             * [\p S/cm],
             * at all quadrature points of the cell.
             */
            std::vector<double> sigmaMeff_cell;

            /**
             * Derivative of effective protonic conductivity
             * w.r.t \p "temperature_of_REV",
             * at all quadrature points in the cell.
             */
            std::vector<double> dsigmaMeff_dT_cell;

            /**
             * Derivative of effective protonic conductivity
             * w.r.t \p "membrane_water_content",
             * at all quadrature points in the cell.
             */
            std::vector<double> dsigmaMeff_dlambda_cell;


            /**
             * \f$ \mathbf{T} \f$ shape functions.
             *
             * \p phi_T_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{T} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_T_cell;

            /**
             * \f$ \mathbf{T} \f$ shape function gradients.
             *
             * \p grad_phi_T_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{T} \f$ shape function gradient
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< Tensor<1,dim> > > grad_phi_T_cell;

            /**
             * \f$ \mathbf{\phi_s} \f$ shape function gradients.
             *
             * \p grad_phi_phiS_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{\phi_s} \f$ shape function gradient
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< Tensor<1,dim> > > grad_phi_phiS_cell;

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
             * \f$ \mathbf{\lambda} \f$ shape function gradients.
             *
             * \p grad_phi_lambda_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{\lambda} \f$ shape function gradient
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< Tensor<1,dim> > > grad_phi_lambda_cell;

            /**
             * \f$ \mathbf{s} \f$ shape functions.
             *
             * \p phi_s_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{s} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_s_cell;
            
            std::vector< std::vector<double> > phi_p_cell;


            /**
             * Gradient of \p "protonic_electrical_potential" at a previous Newton iteration,
             * at all quadrature points in the cell.
             */
            std::vector< Tensor<1,dim> > grad_phiM_cell_old;

            /**
             * Gradient of \p "electronic_electrical_potential" at a previous Newton iteration,
             * at all quadrature points in the cell.
             */
            std::vector< Tensor<1,dim> > grad_phiS_cell_old;

            /**
             * Gradient of \p "membrane_water_content" at a previous Newton iteration,
             * at all quadrature points in the cell.
             */
            std::vector< Tensor<1,dim> > grad_lambda_cell_old;

            //@}

            ///@name Enthalpy transport assemblers - variable data (cell)
            //@{

            /**
             * \p Value in the map represents \f$ \left( c_T \cdot \hat{D}_{i,eff} \cdot \frac{\partial \bar{H}_i}{\partial T} \right) ~\f$ [\p W/\p(cm-K\p)], at previous iteration step at all
             * quadrature points in the cell; corresponding to string \p Key representing for \em e.g \b "oxygen_molar_fraction".
             */
            std::map< std::string, std::vector< Tensor<2,dim> > > conc_Deff_dHdT_map;

            /**
             * \p Value in the map represents \f$ \frac{\partial \left( c_T \cdot \hat{D}_{i,eff} \cdot \frac{\partial \bar{H}_i}{\partial T} \right) }{\partial T} ~\f$ [\p W/\p(cm-K^2\p)], at previous
             * iteration step at all quadrature points in the cell; corresponding to string \p Key representing for \em e.g \b "oxygen_molar_fraction".
             */
            std::map< std::string, std::vector< Tensor<2,dim> > > dT_concDeffdHdT_map;

            /**
             * \p Value in the map represents \f$ \frac{\partial \left( c_T \cdot \hat{D}_{i,eff} \cdot \frac{\partial \bar{H}_i}{\partial T} \right) }{\partial s} ~\f$ [\p W/\p(cm-K\p)], at previous
             * iteration step at all quadrature points in the cell; corresponding to string \p Key representing for \em e.g \b "oxygen_molar_fraction".
             */
            std::map< std::string, std::vector< Tensor<2,dim> > > ds_concDeffdHdT_map;

            std::map< std::string, std::vector< Tensor<2,dim> > > dp_concDeffdHdT_map;
            /**
             * Map of \f$ \mathbf{x_i} \f$ shape function gradients.
             *
             * \p Value in the \p grad_phi_xi_map \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{x_i} \f$ shape function gradient
             * computed in \f$ q \f$-th quadrature point of the cell, corresponding to string \p Key representing for \em e.g \b "oxygen_molar_fraction".
             */
            std::map< std::string, std::vector< std::vector< Tensor<1,dim> > > > grad_phi_xi_map;

            /**
             * \f$ \left(\frac{n_d \sigma_{m,eff}}{F} \right) \frac{\partial \bar{H}_{\lambda}}{\partial T} \f$, at all quadrature points in the cell.
             */
            std::vector<double> electroosmotic_dhdT;

            /**
             * \f$ \frac{ \partial \left[ \left(\frac{n_d \sigma_{m,eff}}{F} \right) \frac{\partial \bar{H}_{\lambda}}{\partial T} \right]}{\partial \lambda} \f$, at all quadrature points in the cell.
             */
            std::vector<double> delectroosmotic_dhdT_dlambda;

            /**
             * \f$ \frac{ \partial \left[ \left(\frac{n_d \sigma_{m,eff}}{F} \right) \frac{\partial \bar{H}_{\lambda}}{\partial T} \right]}{\partial T} \f$, at all quadrature points in the cell.
             */
            std::vector<double> delectroosmotic_dhdT_dT;

            /**
             * \f$ \left(\frac{\rho_{dry}}{EW} \right) D_{\lambda ,eff} \frac{\partial \bar{H}_{\lambda}}{\partial T} \f$, at all quadrature points in the cell.
             */
            std::vector<double> backdiff_dhdT;

            /**
             * \f$ \frac{ \partial \left[ \left(\frac{\rho_{dry}}{EW} \right) D_{\lambda ,eff} \frac{\partial \bar{H}_{\lambda}}{\partial T} \right]}{\partial \lambda} \f$, at all quadrature points in the cell.
             */
            std::vector<double> dbackdiff_dhdT_dlambda;

            /**
             * \f$ \frac{ \partial \left[ \left(\frac{\rho_{dry}}{EW} \right) D_{\lambda ,eff} \frac{\partial \bar{H}_{\lambda}}{\partial T} \right]}{\partial T} \f$, at all quadrature points in the cell.
             */
            std::vector<double> dbackdiff_dhdT_dT;

            /**
             * \f$ \left(\frac{D_{T,eff}}{M_{water}} \right) \frac{\partial \bar{H}_{\lambda}}{\partial T} \f$, at all quadrature points in the cell.
             */
            std::vector<double> thermoosmotic_dhdT;

            /**
             * \f$ \frac{ \partial \left[ \left(\frac{D_{T,eff}}{M_{water}} \right) \frac{\partial \bar{H}_{\lambda}}{\partial T} \right]}{\partial T} \f$, at all quadrature points in the cell.
             */
            std::vector<double> dthermoosmotic_dhdT_dT;

            //@}

            ///@name Local CG FEM based assemblers - variable data (boundary)
            //@{

            /**
             * Derivative of effective thermal conductivity
             * w.r.t \p "temperature_of_REV",
             * at all quadrature points of the boundary.
             */
            std::vector< Tensor<2,dim> > dkeff_dT_bdry;

            /**
             * \f$ \mathbf{T} \f$ shape functions.
             *
             * \p phi_T_bdry \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{T} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the boundary.
             */
            std::vector< std::vector<double> > phi_T_bdry;

            //@}

            /**
             * Counter set to \em TRUE when \p cell_residual is being assembled.
             * This ensures that only effective transport properties are calculated, not their derivatives. (improves speed)
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
