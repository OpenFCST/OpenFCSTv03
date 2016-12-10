// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: reaction_source_terms_base.h
// - Description: This class is used to assemble both cell matrix and cell residual
//                for reaction source terms in the catalyst layers for various equation classes
// - Developers: Marc Secanell, 2015
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_EQUATION_REACTION_SOURCE_TERMS_BASE_H_
#define _FCST_FUELCELLSHOP_EQUATION_REACTION_SOURCE_TERMS_BASE_H_

#include <utils/fcst_constants.h>

#include <equations/equation_base.h>
#include <equations/reaction_heat.h>

#include <layers/catalyst_layer.h>

#include <reactions/tafel_kinetics.h>
#include <reactions/double_trap_kinetics.h>
#include <reactions/dual_path_kinetics.h>

namespace FuelCellShop
{
    namespace Equation
    {

        ///@name Exceptions
        //@{
        /**
         * Exception thrown when a particular "Required" solution variable
         * is not found for \p cathode or \p anode kinetics,
         * in \p user-defined solution variables.
         */
        DeclException2(VariableNotFoundForKinetics,
                       std::string,
                       std::string,
                       << "For " << arg1 << " kinetics source terms, \"" << arg2 << "\" is not found as one of the solution variables.");
        //@}

        /**
         * This class assembles the reaction source terms for all other transport equations, if there's any.
         *
         * \author Madhur Bhaiya,      2013
         * \author Valentin N. Zingan, 2013 - all couplings with fluid transport equations
         */

        template<int dim>
        class ReactionSourceTermsBase : public EquationBase<dim>
        {
        public:

            ///@name Constructors, destructor, and initalization
            //@{

            /**
             * Constructor.
             */
            ReactionSourceTermsBase(FuelCell::SystemManagement& system_management,boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = 
            boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >());

            /**
             * Destructor.
             */
            virtual ~ReactionSourceTermsBase();

            /**
             * Declare parameters.
             */
            virtual void declare_parameters(ParameterHandler& param) const;

            /**
             * Initialize parameters.
             */
            virtual void initialize(ParameterHandler& param){};

            //@}

            ///@name Local CG FEM based assemblers
            //@{

            /**
             * Assemble local cell matrix.
             */
            virtual void assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                              const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                              FuelCellShop::Layer::BaseLayer<dim>* const                               layer){};

            /**
             * Assemble local cell residual.
             */
            virtual void assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                FuelCellShop::Layer::BaseLayer<dim>* const                               layer){};

            //@}

            ///@name Accessors & Info
            //@{

            /**
             * This function prints out
             * the \p info for this class.
             */
            virtual void print_equation_info() const {};

            //@}

        protected:

            
             };
             
        } // Equation
             
} // FuelCellShop
        
#endif