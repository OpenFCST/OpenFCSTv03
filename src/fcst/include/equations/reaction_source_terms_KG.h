// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: reaction_source_terms_KG.h
// - Description: This class is used to assemble both cell matrix and cell residual
//                for reaction source terms in the catalyst layers for various equation classes
// - Developers: Marc Secanell, 2015
//
// ----------------------------------------------------------------------------
#ifndef _FCST_FUELCELLSHOP_EQUATION_REACTION_SOURCE_TERMS_KG_H_
#define _FCST_FUELCELLSHOP_EQUATION_REACTION_SOURCE_TERMS_KG_H_

#include <equations/equation_base.h>
#include <equations/reaction_heat.h>
#include <equations/reaction_source_terms_base.h>

namespace FuelCellShop
{
    namespace Equation
    {
        /**
         * This class assembles the reaction source terms for all other transport equations, if there's any.
         *
         * \author Marc Secanell,      2015
         * \author Valentin N. Zingan, 2013 - all couplings with fluid transport equations
         */

        template<int dim>
        class ReactionSourceTermsKG : public ReactionSourceTermsBase<dim>
        {
        public:

            ///@name Constructors, destructor, and initalization
            //@{

            /**
             * Constructor.
             */
            ReactionSourceTermsKG(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = 
            boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >());

            /**
             * Destructor.
             */
            virtual ~ReactionSourceTermsKG();

            /**
             * Declare parameters.
             */
            virtual void declare_parameters(ParameterHandler& param) const;

            /**
             * Initialize parameters.
             */
            virtual void initialize(ParameterHandler& param);

            /**
             * Set the pointer to cathode kinetics in the object.
             * \note If an application uses cathode kinetics model, this method should be called neccessarily before initializing this object (using #initialize method).
             */
            inline void set_cathode_kinetics(FuelCellShop::Kinetics::BaseKinetics* kinetics)
            {
                cathode_kinetics = kinetics;
            }

            /**
             * Set the pointer to anode kinetics in the object.
             * \note If an application uses anode kinetics model, this method should be called neccessarily before initializing this object (using #initialize method).
             */
            inline void set_anode_kinetics(FuelCellShop::Kinetics::BaseKinetics* kinetics)
            {
                anode_kinetics = kinetics;
            }

            //@}

            ///@name Local CG FEM based assemblers
            //@{

            /**
             * Assemble local cell matrix.
             */
            virtual void assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                              const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                              FuelCellShop::Layer::BaseLayer<dim>* const                               layer);

            /**
             * Assemble local cell residual.
             */
            virtual void assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                FuelCellShop::Layer::BaseLayer<dim>* const                               layer);

            //@}

            ///@name Accessors & Info
            //@{

            /**
             * This function is used to adjust \p std::vector \p < \p internal_cell_couplings \p >, which is generated after getting \p internal_cell_couplings
             * from all the equations being used in the application.
             *
             * \note It's very important to use this function, if we are considering source terms due to reaction in the catalyst layers. It's also noteworthy that
             * this function should be used after whole vector of \p internal_cell_couplings is created, and \b BEFORE \p make_cell_couplings of \b SystemManagement
             * is called.
             */
            virtual void adjust_internal_cell_couplings( std::vector< couplings_map >& dst,
                                                         const std::vector< FuelCellShop::Material::PureGas* >& gases = std::vector< FuelCellShop::Material::PureGas* >() );

            /**
             * This function prints out
             * the \p info for this class.
             */
            virtual void print_equation_info() const;

            /**
             * Accessor for cathode catalyst layer FuelCellShop::Kinetics::BaseKinetics pointer.
             */
            inline FuelCellShop::Kinetics::BaseKinetics* get_cathode_kinetics() const
            {
                return cathode_kinetics;
            }

            /**
             * Accessor for anode catalyst layer FuelCellShop::Kinetics::BaseKinetics pointer.
             */
            inline FuelCellShop::Kinetics::BaseKinetics* get_anode_kinetics() const
            {
                return anode_kinetics;
            }

            //@}

        protected:

            ///@name Local CG FEM based assemblers - make_ functions
            //@{

            /**
             * This function computes Local CG FEM based
             * assemblers - constant data (generic).
             *
             * \warning <b>For the DEVELOPERS:</b> \p block_index for VariableInfo are not filled here, but \p indices_exist flag is set to \b TRUE (as it
             * is needed at other places before \p cell_matrix or \p cell_residual assembly).
             * \p Developers need to be wary of this fact that the \p block_index are still not filled yet, and they need to <b>fill them before</b> doing \p cell_matrix assembly.
             */
            virtual void make_assemblers_generic_constant_data();

            /**
             * This function computes <p> Local CG FEM based assemblers - constant data (cell) </p>
             * and allocates the memory for shape functions, shape function gradients, and
             * \p JxW_cell in <p> Local CG FEM based assemblers - variable data (cell) </p>.
             */
            virtual void make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info);

            /**
             * This function computes <p> Local CG FEM based assemblers - variable data (cell) </p>.
             */
            virtual void make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                            FuelCellShop::Layer::BaseLayer<dim>* const layer);

            //@}

            ///@name Kinetics Objects
            //@{
            /**
             * Pointer to Cathode Kinetics object, initialized in the constructor.
             */
            FuelCellShop::Kinetics::BaseKinetics* cathode_kinetics;

            /**
             * Pointer to Anode Kinetics object, initialized in the constructor.
             */
            FuelCellShop::Kinetics::BaseKinetics* anode_kinetics;

            //@}

           ///@name Counters
            //@{

            /**
             * Variable used to store the index in cell_info->global_data of the previous Newton solution
             * The solution at the previous iteration is used to compute cell_matrix and cell_residual
             */
            unsigned int last_iter_cell;
            
            //@}
            
            ///@name Methods and Data related to fluid transport equations
            //@{
            
            /**
             * This function fills out the
             * max number of \p matrix_block_indices which need to be updated due to the source terms.
             * The real number might be less depending on in which CL we are.
             */
            virtual void make_matrix_block_indices();
            
            /**
             * This function fills out the
             * max number of \p residual_indices which need to be updated due to the source terms.
             * The real number might be less depending on in which CL we are.
             */
            virtual void make_residual_indices();
            
            /**
             * Number of species, \f$ N \f$.
             */
            unsigned int n_species;
            
            /**
             * Keeps oxygen set of fluid transport equations.
             */
            unsigned int indexO2;
            
            /**
             * Keeps hydrogen set of fluid transport equations.
             */
            unsigned int indexH2;
            
            /**
             * Keeps water vapor set of fluid transport equations.
             */
            unsigned int indexH2O;
            
            /**
             * Oxygen source multiplier,      \f$ \quad \left[ \frac{\text{g}}{A sec} \right] \f$.
             */
            double multiplierO2;
            
            /**
             * Hydrogen source multiplier,    \f$ \quad \left[ \frac{\text{g}}{A sec} \right] \f$.
             */
            double multiplierH2;
            
            /**
             * Water vapor source multiplier, \f$ \quad \left[ \frac{\text{g}}{A sec} \right] \f$.
             */
            double multiplierH2O;
            
            /**
             * Molar mass of pure gas, \f$ M_i \quad \left[ \frac{\text{g}}{\text{mol}} \right] \f$.
             */
            std::vector<double> molar_mass;
            
            /**
             * Density extractors.
             */
            std::vector< FEValuesExtractors::Scalar > density_extractors;
            
            /**
             * Electronic electrical potential extractor.
             */
            FEValuesExtractors::Scalar electronic_electrical_potential_extractor;
            
            /**
             * Protonic electrical potential extractor.
             */
            FEValuesExtractors::Scalar protonic_electrical_potential_extractor;
            
            /**
             * Constant temperature of species mixture
             * in the quadrature points of a cell,
             * \f$ T_{\text{mixture}}^{\text{const}} \quad \left[ \text{K} \right] \f$.
             */
            std::vector<double> T_mixture;
            
            /**
             * ORR current density
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector<double> ORR_current_density;
            
            /**
             * ORR current density derivative
             * with respect to oxygen concentration (gas, NOT gas-liquid)
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector<double> DORR_current_density_Doxygen_concentration;
            
            /**
             * ORR current density derivative
             * with respect to electronic electrical potential
             * in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector<double> DORR_current_density_Delectronic_electrical_potential;
            
            /**
             * ORR current density derivative with respect to protonic electrical potential
             * in the quadrature points of a cell at a previous Newton iteration.
             */
            std::vector<double> DORR_current_density_Dprotonic_electrical_potential;
            
            /**
             * HOR current density in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector<double> HOR_current_density;
            
            /**
             * HOR current density derivative with respect to hydrogen concentration (gas, NOT gas-liquid)
             * in the quadrature points of a cell at a previous Newton iteration.
             */
            std::vector<double> DHOR_current_density_Dhydrogen_concentration;
            
            /**
             * HOR current density derivative with respect to electronic electrical potential
             * in the quadrature points of a cell at a previous Newton iteration.
             */
            std::vector<double> DHOR_current_density_Delectronic_electrical_potential;
            
            /**
             * HOR current density derivative with respect to protonic electrical potential
             * in the quadrature points of a cell at a previous Newton iteration.
             */
            std::vector<double> DHOR_current_density_Dprotonic_electrical_potential;
            
            /**
             * Density of each species in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector< std::vector<double> > density_old;
            
            /**
             * Electronic electrical potential in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector<double> electronic_electrical_potential_old;
            
            /**
             * Protonic electrical potential in the quadrature points of a cell
             * at a previous Newton iteration.
             */
            std::vector<double> protonic_electrical_potential_old;
            
            /**
             * Density shape functions.
             *
             * \p phi_density \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th density shape function
             * computed in \f$ q \f$-th quadrature point of a cell
             * for species \f$ s \f$.
             */
            std::vector< std::vector< std::vector<double> > > phi_density;
            
            /**
             * Electronic electrical potential shape functions.
             *
             * \p phi_electronic_electrical_potential \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th electronic electrical potential shape function
             * computed in \f$ q \f$-th quadrature point of a cell.
             */
            std::vector< std::vector<double> > phi_electronic_electrical_potential;
            
            /**
             * Protonic electrical potential shape functions.
             *
             * \p phi_protonic_electrical_potential \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th protonic electrical potential shape function
             * computed in \f$ q \f$-th quadrature point of a cell.
             */
            std::vector< std::vector<double> > phi_protonic_electrical_potential;
            
            /**
             * For internal use only.
             */
            std::string eq_generic_prefix;
            
            /**
             * For internal use only.
             */
            std::vector<std::string> eq_postfixes;
            
            /**
             * For internal use only.
             */
            std::vector<std::string> var_postfixes;
            
            /**
             * For internal use only.
             */
            std::string eq_name;
            
            /**
             * For internal use only.
             */
            std::string var_name;
            
            //@}
            /**
             * Function used to set the appropriate couplings between equations. In this case, given the reacting species equation number
             * and the species equation number, the function sets up the coupling with electrical and protonic potential in the case
             * of electrochemical reactions. It also checks to make sure the coupling could be set.
             * 
             * Used in adjusnt_internal_cell_couplings
             */
            inline void set_species_couplings (unsigned int reacting_species_equation_number, unsigned int species_number, std::vector< couplings_map >& equation_map)
            {
                unsigned int index = reacting_species_equation_number;
                unsigned int g = species_number;
                
                eq_name = "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - mass conservation - " + eq_postfixes[g-1];
                          
                couplings_map::iterator iter = equation_map[index].find(eq_name);
                
                if( iter != equation_map[index].end() )
                {
                    iter->second["electronic_electrical_potential"] = DoFTools::always;
                    iter->second["protonic_electrical_potential"]   = DoFTools::always;
                }
                else
                    AssertThrow( false, ExcInternalError() );
            }
            
            inline void set_density_couplings (unsigned int species_number, couplings_map::iterator iter)
            {

                if( species_number != -1 )
                {
                    var_name = "density_" + var_postfixes[species_number-1];
                    iter->second[var_name] = DoFTools::always;
                }
                else
                    AssertThrow( false, ExcNotImplemented() );
            }
        };
             
    } // Equation
             
} // FuelCellShop
        
#endif