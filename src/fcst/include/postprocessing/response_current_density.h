// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: responses.h
// - Description: This is a file containing the declaration of several FCST response evaluators
// - Developers: Marc Secanell Gallart,    University of Alberta
// - $Id: response_current_density.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__RESPONSE_CURRENT_DENSITY_H
#define _FUELCELLSHOP__RESPONSE_CURRENT_DENSITY_H

//Include STL
#include <exception>      // std::exception

// Include OpenFCST routines:
#include <layers/catalyst_layer.h>

#include "postprocessing/base_response.h"

using namespace dealii;

namespace FuelCellShop
{

    namespace PostProcessing
    {

        /**
         * Class used to calculate the ORR current density and coverages (if provided in the kinetic model)
         * by the catalyst layer.
         * 
         * Since the surface area and volume of the catalyst layer is not known, the layer return the
         * total current density produced, i.e.
         * \f[
         * \int_{V_{el}} i dV_{el}
         * \f]
         * and the sum of all coverages at the element, i.e.
         * \f[
         * \int_{V_{el}} \theta_i dV_{el}
         * \f]
         * In order for the class to provide current density and average coverage, the returned quantity
         * has to be divided by surface area and catalyst layer volume respectively. If the surface
         * area and volume are known at the input file, they can be specified and the quantities will
         * automatically be corrected.
         * 
         * Note that this class will only work if an object of type CatalystLayer is passed as an argument.
         * 
         * \author M. Secanell, 2014
         */
        template <int dim>
        class ORRCurrentDensityResponse : public BaseResponse<dim>
        {
        public:            
            ///@name Constructor, declaration and initialization
            //@{
            ORRCurrentDensityResponse(const FuelCell::SystemManagement& sm)
            :
            BaseResponse<dim>(sm)
            {
                oxygen_density_name = "";
            }
            
            ~ORRCurrentDensityResponse() {}
            
            /**
             * 
             */
            void declare_parameters(ParameterHandler& param) const;
            /**
             * Initialize class parameters.
             */
            void initialize(ParameterHandler& param);
            //@}
            ///@name Compute functional
            //@{
            /**
             * This member function computes the volumetric current density produced inside the electrode
             * due to the oxygen reduction reaction. In order to compute the current density, the following functional is
             * evaluated
             * \f[
             * I = int_{\Omega} i dV
             * \f]
             * 
             * In addition, if the catalyst layer uses a kinetic model that provides coverages those values are
             * also provided in the response parameter, i.e., resp.
             * 
             * In order to access the current produced at the electrode use:
             * @code
             * ORR_responses[FuelCellShop::PostProcessing::ResponsesNames::current]
             * @endcode
             * 
             * The coverages are accessed using FuelCellShop::PostProcessing::ResponsesNames::O_coverage 
             * and FuelCellShop::PostProcessing::ResponsesNames::OH_coverage
             */
            void compute_responses(const typename DoFApplication<dim>::CellInfo& info,
                                   FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                   std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& resp) const;
            /**
             * Routine used in order to compute the response with a modified solution (not the one stored
             * in info)
             * 
             */                       
            void compute_responses(std::vector< FuelCellShop::SolutionVariable > solution_variables,
                                   const typename DoFApplication<dim>::CellInfo& info,
                                   FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                   std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& resp) const;
            
            /**
             * For the case of multi-component solvers, we need to specify which one of the variables contains the oxygen concentration,
             * e.g., density_species_1, density_species_2 and so on. Provide this value here.
             */
            void set_oxygen_density_name(std::string name)
            {
                oxygen_density_name = name;
            }
        private:
            /**
             * Surface area
             */
            double S_CL;
            /**
             * Volume
             */
            double V_CL;
            /**
             * VariableInfo structure corresponding to the
             * reactant_molar_fraction.
             */
            FuelCellShop::Equation::VariableInfo xi;
            /**
             * VariableInfo structure corresponding to the
             * electrolyte potential.
             */
            FuelCellShop::Equation::VariableInfo phiM;
            /**
             * VariableInfo structure corresponding to the
             * electronic phase potential.
             */            
            FuelCellShop::Equation::VariableInfo phiS;
            /**
             * Set oxygen concentration name
             */
            std::string oxygen_density_name;
            
        };
        
        
        //===============================================================================================================
        //===============================================================================================================
        //===============================================================================================================
        //===============================================================================================================
        
        /**
         * Class used to calculate the current density at the anode catalyst layer.
         * 
         * Since the surface area and volume of the catalyst layer is not known, the layer return the
         * total current density produced, i.e.
         * \f[
         * \int_{V_{el}} i dV_{el}
         * \f]
         * 
         * In order for the class to provide current density and average coverage, the returned quantity
         * has to be divided by surface area and catalyst layer volume respectively. If the surface
         * area and volume are known at the input file, they can be specified and the quantities will
         * automatically be corrected.
         * 
         * Note that this class will only work if an object of type CatalystLayer is passed as an argument.
         * 
         * \author M. Secanell, 2014
         */
        template <int dim>
        class HORCurrentDensityResponse : public BaseResponse<dim>
        {
        public:            
            ///@name Constructor, declaration and initialization
            //@{
            HORCurrentDensityResponse(const FuelCell::SystemManagement& sm)
            :
            BaseResponse<dim>(sm)
            {}
            
            ~HORCurrentDensityResponse() {}
            
            /**
             * 
             */
            void declare_parameters(ParameterHandler& param) const;
            /**
             * Initialize class parameters.
             */
            void initialize(ParameterHandler& param);
            //@}
            ///@name Compute functional
            //@{
            /**
             * This member function computes the volumetric current density produced inside the electrode
             * due to the oxygen reduction reaction. In order to compute the current density, the following functional is
             * evaluated
             * \f[
             * I = int_{\Omega} i dV
             * \f]
             * 
             * In addition, if the catalyst layer uses a kinetic model that provides coverages those values are
             * also provided in the response parameter, i.e., resp.
             * 
             * In order to access the current produced at the electrode use:
             * @code
             * ORR_responses[FuelCellShop::PostProcessing::ResponsesNames::current]
             * @endcode
             * 
             * The coverages are accessed using FuelCellShop::PostProcessing::ResponsesNames::O_coverage 
             * and FuelCellShop::PostProcessing::ResponsesNames::OH_coverage
             */
            void compute_responses(const typename DoFApplication<dim>::CellInfo& info,
                                   FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                   std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& resp) const;
            /**
             * Routine used in order to compute the response with a modified solution (not the one stored
             * in info)
             * 
             */                       
            void compute_responses(std::vector< FuelCellShop::SolutionVariable > solution_variables,
                                   const typename DoFApplication<dim>::CellInfo& info,
                                   FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                   std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& resp) const;
            //
            //
        private:
            /**
             * Surface area
             */
            double S_CL;
            /**
             * Volume
             */
            double V_CL;
            /**
             * VariableInfo structure corresponding to the
             * reactant_molar_fraction.
             */
            //FuelCellShop::Equation::VariableInfo x_H2;
            /**
             * VariableInfo structure corresponding to the
             * electrolyte potential.
             */
            FuelCellShop::Equation::VariableInfo phiM;
            /**
             * VariableInfo structure corresponding to the
             * electronic phase potential.
             */            
            FuelCellShop::Equation::VariableInfo phiS;
        };
        
    }   
}

#endif