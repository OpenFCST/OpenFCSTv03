//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: ficks_transport_equation.h
//    - Description: Header file of Ficks diffusion equation class.
//    - Developers: M. Secanell, M. Bhaiya, V. N. Zingan, A. Kosakian, M. Sabharwal
//
//---------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_EQUATION_NEW_FICKS_EQUATION_H_
#define _FCST_FUELCELLSHOP_EQUATION_NEW_FICKS_EQUATION_H_

// FCST includes
#include <utils/fcst_units.h>
#include <utils/fcst_constants.h>
#include <equations/equation_base.h>
#include <layers/gas_diffusion_layer.h>
#include <layers/micro_porous_layer.h>
#include <layers/catalyst_layer.h>
#include <equations/equation_auxiliaries.h>

// STD
#include <sstream>
#include <string>


namespace FuelCellShop
{
    namespace Equation
    {
        /**
         * This class deals with <b>Fick's Transport Equation</b>.
         *
         * This equation class solves for Fick's law of diffusion inside the porous layers. The solute gas and solvent gas
         * are normally passed inside the constructor of this equation class.
         *
         * It is solved with respect to:
         * - \f$ \mathbf{x_i} \f$ \p(solutegas_molar_fraction \p)
         *
         * This equation can be written as:
         *
         *  \f$ \qquad \mathbf{\nabla} \cdot \left( C_T \hat{D}_{i,eff} \mathbf{\nabla} x_i \right) = 0, \quad  x_i \in \Omega \f$
         *
         * - where, \f$ \hat{D}_{i,eff} \f$ is tensor of effective molecular diffusivity of solute gas in solvent gas [\p cm^2/s],
         * which can be function of \f$ T \f$ \p(temperature_of_REV \p) and \f$ s \f$ \p(liquid_water_saturation \p).
         * - \f$ C_T \f$ is total concentration [\p mol/cm^3].
         *
         * To be well-posed, these equations are equipped with the appropriate boundary conditions. All the boundary conditions can be described by
         * \p boundary_id \p(s \p) and \p boundary_type. Besides, this some boundary types require additional information, which can also be provided by the parameter file.
         * We consider following types of boundary conditions:
         *
         * - \b Dirichlet: At this boundary, we specify the \f$ x_i \f$ values. It is to be noted that this equation class is only used to specify whether
         * a particular boundary is Dirichlet or not. Molar fraction values at that boundary should be taken care by Initial Solution methods for a particular application.
         *
         * The \p boundary_ids are specified in the parameter file under subsection <b>"Dirichlet Boundary Indicators"</b>, as a list of comma-separated values.
         *
         *  \em e.g. @code set Dirichlet Boundary Indicators = 3, 4, 8 @endcode
         *
         * - <b>No gas flux / Symmetric:</b> A particular case of \p Neumann boundary condition.
         *
         * \f$ \qquad -\mathbf{n} \cdot \left( C_T \hat{D}_{i,eff} \cdot \mathbf{\nabla} x_i \right) = 0 \f$
         *
         * \remarks
         * - There is no provision to specify boundary indicators for \p No \p gas \p flux or \p Symmetric boundary conditions, as FEM
         * formulation automatically implies a particular boundary is one of these cases, by default.
         * - This class currently works with the following layer classes:
         *      - FuelCellShop::Layer::GasDiffusionLayer<dim>
         *      - FuelCellShop::Layer::MicroPorousLayer<dim>
         *      - FuelCellShop::Layer::CatalystLayer<dim>
         * - In the case of \b isothermal applications, it is necessary to use FuelCellShop::Layer::PorousLayer<dim>::set_gases_and_compute method
         * in the initialization of the application. This method sets the gases to be solved inside the layer and also computes
         * the isobaric isothermal bulk molecular diffusion coefficients.
         * - In the case of \b nonisothermal applications, it is necessary to use FuelCellShop::Layer::PorousLayer<dim>::set_gases method
         * in the initialization of the application. This method sets the gases to be solved inside the layer. It is recommended to set the solvent gas
         * as the last gas in the input vector for \p set_gases method.
         *
         * We solve the whole problem by linearizing the governing equation at each Newton iteration with subsequent
         * CG FEM discretization in space. The class contains the necessary class members to add the necessary contributions to cell_matrix
         * and cell_residual to the governing equations used to analyze gas transport via ficks diffusion model,
         *
         * <h3>Usage Details</h3>
         *
         * @code
         * // Creating Equation object (in Application Header file)
         * FuelCellShop::Equation::FicksTransportEquation<dim> oxygen_transport;
         *
         * // Declare parameters in application
         * oxygen_transport.declare_parameters(param);
         *
         * // Initialize in application
         * oxygen_transport.initialize(param);
         *
         * // Create a temporary vector in the application for storing couplings_map from all the equation used in the application.
         * std::vector<couplings_map> tmp;
         * ... // other equations
         * tmp.push_back( oxygen_transport.get_internal_cell_couplings() );
         *
         * // Look at ReactionSourceTerms class here, if source terms due to current production/consumption are to be considered.
         * // Making cell couplings using SystemManagement object created in the application
         * system_management.make_cell_couplings(tmp);
         *
         * // cell_matrix in application
         * // Do a check against layer and it should match with the layers currently working for this equation class.
         * // for eg: CCL is FuelCellShop::Layer::HomogeneousCL<dim> object.
         * oxygen_transport.assemble_cell_matrix(cell_matrices, cell_info, &CCL);
         *
         * // cell_residual in application
         * oxygen_transport.assemble_cell_residual(cell_vector, cell_info, &CCL);
         * @endcode
         *
         * \note This class doesn't assemble for current production/consumption source terms; that is taken care off by \p ReactionSourceTerms class.
         * Please read the documentation of ReactionSourceTerms class, for additional methods to be implemented in the application.
         *
         * \warning If current production/consumption source terms are being considered, it's very important to use \p adjust_internal_cell_couplings
         * member function of ReactionSourceTerms class, before using  \p make_cell_couplings of \b SystemManagement at the application level.
         *
         * TODO Old Boundary conditions including Dirichlet Boundary Indicators are supposed to be replaced with the new subsections,
         *      see TO BE REMOVED comments in .cc file.
         *
         * \author Marc Secanell,      2006-13
         * \author Madhur Bhaiya,      2013
         * \author Valentin N. Zingan, 2012-2014 - afterward improvements, optimization, checkings, CG FEM bug fixings
         */

        template<int dim>
        class FicksTransportEquation
        :
        public EquationBase<dim>
        {
        public:
            ///@name Constructors, destructor, and initalization
            //@{
            /**
             * Constructor.
             * @note This is the recommended construtor. The solute and solvent are passed immediately,
             * so that the class can setup its equation name based on the solute being considered.
             * The name of the equation is necessary before initialize() is called.
             */
            FicksTransportEquation(FuelCell::SystemManagement& system_management,
                                      FuelCellShop::Material::PureGas* solute,
                                      FuelCellShop::Material::PureGas* solvent,
                                      boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >()
                                     );

            /**
             * Constructor.
             * @warning This constructor is not recommened. If using this constructor, use member function set_solute_and_solvent to setup
             * the gases and name of the equation before calling initialize.
             */
            FicksTransportEquation(FuelCell::SystemManagement& system_management, std::string& name_section, boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >());

            /**
             * Destructor.
             */
            virtual ~FicksTransportEquation();

            /**
             * Declare parameters.
             */
            virtual void declare_parameters(ParameterHandler& param);

            /**
             * Initialize parameters.
             */
            virtual void initialize(ParameterHandler& param);

            /**
             * Method to set solute and solve if other constructor (not passing solute and solvent in the constructor)
             * is being used.
             * It will also setup name_equation and name_solution
             *
             */
            inline void set_solute_and_solvent( FuelCellShop::Material::PureGas* solute, FuelCellShop::Material::PureGas* solvent, ParameterHandler& param )
            {
                this->gas = solute;
                this->solvent = solvent;

                std::stringstream s;
                s <<"Ficks Transport Equation - "<<gas->name_material();
                this->equation_name = s.str();

                std::stringstream ss;
                ss <<gas->name_material()<<"_molar_fraction";
                this->name_base_variable = ss.str();
            }

            //@}

            ///@name Local CG FEM based assemblers
            //@{

            
            /**
             * Assemble local cell residual for nonlinear problems.
             */
            virtual void assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_rhs,
                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                FuelCellShop::Layer::BaseLayer<dim>* const                               layer
                                               );            
            
            /**
             * Assemble local boundary matrix.
             */
            virtual void assemble_bdry_matrix(FuelCell::ApplicationCore::MatrixVector&                                 bdry_matrices,
                                              const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                              FuelCellShop::Layer::BaseLayer<dim>* const                               layer
                                             );

            /**
             * Assemble local boundary residual.
             */
            virtual void assemble_bdry_residual(FuelCell::ApplicationCore::FEVector&                                     bdry_rhs,
                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                FuelCellShop::Layer::BaseLayer<dim>* const                               layer
                                               );

            //@}

            ///@name Accessors and info
            //@{
            /**
             * The function printing out
             * the equations info.
             */
            virtual void print_equation_info() const;

            /**
             * Member function used to test the functionality of the class.
             * It create an object of this class and test functionality.
             *
             */
            void class_test();
            
            /**
             * Member function used to compute the effective Knudsen diffusivity using the Knudsen radius entered using the field data in VTK mesh.
             * Cell index is used as the key to reference the local cell radius which is used to compute the Knudsen diffusion coefficient.
             */
            double compute_Knudsen_diffusivity(int cell_index);
            
            /**
             * Member function used to compute effective diffusion coefficient using the Bosanquet approximation by combining the Knudsen and Bulk diffusion coefficients.
             */
            Tensor<2,dim,double> effective_diffusion_coefficient(Tensor<2,dim,double>& D_bulk,
                                                                        double D_Knud);
            //@}

        protected:
            ///@name Accessors and info
            //@{
            /**
             * Computes the current using oxygen dissolution in the ionomer film by solving a Newton loop. 
             * The ICCP model is used to solve the current for the microstructure. 
             * This function runs the Newton loop to solve for the oxygen concentration at the ionomer-solid interface.
             * The function returns the computed current density. 
             * The coupled equations solved are:
             * 
             * \f$ \mathbf{N_{O_2}}\cdot \mathbf{n} = -k \left(c_{g|f}-c^{eq}_{g|f}\right) = \frac{D^{film}}{\delta} \left(c_{g|f}-c_{react}\right) = \frac{j(c_{react})}{4F} \f$
             * 
             * 
             * where \f$ k \f$ is the ionomer dissolution rate, \f$ D^{film} \f$ is the diffusion rate coefficient for oxygen in nafion, \f$ \delta \f$ is the ionomer film thickness, \f$ c_{g|f}\f$ is the oxygen concentration at gas-ionomer interface, \f$ c^{eq}_{g|f}\f$ is the equilibrium oxygen concentration computed using Henry's Law, \f$ c_{react} \f$ is the reactant concentration
             * at the ionomer-solid interface and \f$ j \f$ is the current density. The residual for the Newton loops is implemented in the function \p get_ICCP_residual. 
             * The derivative of the residual wrt to \f$ c_{react} \f$ is approximated using the central difference method. 
             */
            double compute_current(double x_gas);
            
            /**
             * Computes the residual for the compute_current function. The residual (\f$ R \f$) is defined as:
             * \f$ R = \frac{j}{4F} - \frac{\frac{kD^{film}}{\delta}}{k + \frac{D^{film}}{\delta}}\left(c^{eq}_{g|f}-c_{react}\right)
             * where \f$ k \f$ is the ionomer dissolution rate, \f$ D^{film} \f$ is the diffusion rate coefficient for oxygen in nafion, \f$ \delta \f$ is the ionomer film thickness, \f$ c^{eq}_{g|f}\f$ is the equilibrium oxygen concentration computed using Henry's Law, \f$ c_{react} \f$ is the reactant concentration
             * at the ionomer-solid interface.
             */
            double get_ICCP_residual(double c_react, double c_eq);
            //@}
            
            ///@name Local CG FEM based assemblers
            //@{            
            /**
             * Assembles linear steady-state cell matrix.
             */
            void assemble_cell_linear_matrix(FuelCell::ApplicationCore::MatrixVector& cell_matrices,
                                             const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                             FuelCellShop::Layer::BaseLayer<dim>* const layer );
            
            /**
             * Assembles system matrix for nonlinear case. 
             */            
            void assemble_cell_Jacobian_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                               const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                               FuelCellShop::Layer::BaseLayer<dim>* const                               layer );
                    
            //@}
            
            ///@name Boundary conditions
             //@{

             /**
              * Container of \p boundary_id \p(s \p) for Dirichlet boundaries
              */
             std::vector<unsigned int> dirichlet_bdry_ids;

             /**
              * \p std::map< \p unsigned \p int, \p double \p> container for details regarding <b>Neumann</b> boundary conditions.
              * Here, \p Key (unsigned int) represents the \p boundary_id and \p Value (double) represents the constant species flux values [mol/(cm^2 sec)].
              */
             std::map<unsigned int, double> species_flux;

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

            //@}

            ///@name Generic Constant Data
            //@{

            /**
             * Gas for which the equation is setup.
             */
            FuelCellShop::Material::PureGas* gas;

            /**
             * Solvent gas for which the equation is setup. The diffusion coefficient that we would use is given by
             * D_gas,solvent.
             */
            FuelCellShop::Material::PureGas* solvent;

            /**
              * VariableInfo structure corresponding to base variable
              * of this equation class, \p "solute.name()_molar_fraction".
              */
            VariableInfo xi;

            /**
             * VariableInfo structure corresponding to \p "temperature_of_REV".
             */
            VariableInfo t_rev;

            /**
             * VariableInfo structure corresponding to \p "liquid_water_saturation".
             */
            VariableInfo s_liquid_water;
            
            /**
             * VariableInfo structure corresponding to \p "capillary_pressure".
             */
            VariableInfo p_liquid_water;

            //@}

            ///@name Local CG FEM based assemblers - variable data (cell)
            //@{
            
            /**
             * Tensor with concentration times effective molecular diffusivity [\p mol/\p(cm-s \p)],
             * at all quadrature points in the cell.
             */
            std::vector< Tensor< 2, dim > > conc_Deff_cell;

            /**
             * Derivative of concentration times effective molecular diffusivity
             * w.r.t \p "temperature_of_REV" [\p mol/\p(cm-s-K \p)],
             * at all quadrature points in the cell.
             */
            std::vector< Tensor<2,dim> > dconc_Deff_dT_cell;

            /**
             * Derivative of concentration times effective molecular diffusivity
             * w.r.t \p "liquid_water_saturation" [\p mol/\p(cm-s \p)],
             * at all quadrature points in the cell.
             */
            std::vector< Tensor<2,dim> > dconc_Deff_ds_cell;
            
            
                        /**
             * Derivative of concentration times effective molecular diffusivity
             * w.r.t \p "liquid_water_saturation" [\p mol/\p(cm-s \p)],
             * at all quadrature points in the cell.
             */
            std::vector< Tensor<2,dim> > dconc_Deff_dp_cell;

            /**
             * \f$ \mathbf{x_{i}} \f$ shape function
             *
             * \p phi_xi_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{x_{i}} \f$ shape function
             * computed at \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< double > >  phi_xi_cell;
            
            /**
             * \f$ \mathbf{x_{i}} \f$ shape function gradients.
             *
             * \p grad_phi_xi_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{x_{i}} \f$ shape function gradient
             * computed at \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< Tensor<1,dim> > > grad_phi_xi_cell;

            /**
             * \f$ \mathbf{T} \f$ shape functions.
             *
             * \p phi_T_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{T} \f$ shape function
             * computed at \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< double > > phi_T_cell;

            /**
             * \f$ \mathbf{s} \f$ shape functions.
             *
             * \p phi_s_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{s} \f$ shape function
             * computed at \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< double > > phi_s_cell;
            
            /**
             * \f$ \mathbf{s} \f$ shape functions.
             *
             * \p phi_s_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{s} \f$ shape function
             * computed at \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector< double > > phi_p_cell;

            //@}

            ///@name Local CG FEM based assemblers - variable data (boundary)
            //@{

            /**
             * \f$ \mathbf{x_i} \f$ shape functions.
             *
             * \p phi_xi_bdry \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{x_i} \f$ shape function
             * computed at \f$ q \f$-th quadrature point of the boundary.
             */
            std::vector< std::vector<double> > phi_xi_bdry;

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
            

            
        private:
            
            /**
             * Shared pointer for the kinetics class used when the reaction at the boundaries is switched on
             */
            boost::shared_ptr<FuelCellShop::Kinetics::BaseKinetics> kinetics;
            
            /**
             * Shared pointer for the electrolyte class used when the reaction on the boundaries is switched on
             */
            boost::shared_ptr<FuelCellShop::Material::PolymerElectrolyteBase> electrolyte;
            
            /**
             * Shared pointer for the catalyst class used when the reaction on the boundaries is switched on
             */
            boost::shared_ptr<FuelCellShop::Material::CatalystBase> catalyst;
            
            /**
             * Concentration of the gas [in mol/cm3] computed using ideal gas law
             */
            double concentration;
            
            /**
             * Electronic potential [in V] to be used for the reaction term 
             */
            double electronic_pot;
            
            /**
             * Protonic potential [in V] to be used for the reaction term 
             */
            double protonic_pot;
            
            /**
             * Temperature [in K] of the layer 
             */
            double temperature;
            
            /**
             * Oxygen dissolution rate constant in ionomer film
             */
            double k_ionomer;
            
            /**
             * Ionomer film thickness around the catalyst
             */
            double delta;
            
            /**
             * Ratio of the active platinum area to the solid-gas interface in the reconstruction
             */
            double prefactor;
            
    };

} // Equation

} // FuelCellShop

#endif