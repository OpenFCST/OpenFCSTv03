// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: response_sorption_heat.h
// - Description: This is header file for response class computing heat generated due to sorption of water.
// - Developers: Madhur Bhaiya
// - $Id: 
//
// ----------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__RESPONSE_SORPTION_HEAT_H
#define _FUELCELLSHOP__RESPONSE_SORPTION_HEAT_H

//Include STL
#include <exception>      // std::exception

// Include OpenFCST routines:
#include <equations/sorption_source_terms.h>
#include <layers/catalyst_layer.h>
#include <postprocessing/base_response.h>

using namespace dealii;

namespace FuelCellShop
{

    namespace PostProcessing
    {

        /**
         * Class used to calculate the heat generated due to sorption of water inside the catalyst layer.
         * 
         * This class returns the total heat generated due to sorption inside the layer, \em i.e.
         * \f[
         * \int_{\Omega} k h_{sorp} \frac{\rho_{dry}}{EW} (\lambda_{eq} - \lambda) d\Omega
         * \f]
         *
         * Applications can divide the obtained total quantity by surface area, volume etc. of the layer, in order to determine 
         * surface density, volumetric density etc.
         *
         * \remarks
         * - This class will not work if anything other than the CL layer object is passed as argument.
         * - FuelCellShop::Equation::SorptionSourceTerms object reference is required in the constructor, besides the FuelCell::SystemManagement object reference.
         * - \f$ \lambda \f$, \f$ x_{H_2O} \f$ and \f$ T_{rev} \f$ should be solved for in the application.
         * - Flag parameters from the SorptionSourceTerms, eg, 'Heat source/sink due to sorption/desorption' is Enabled or Not, are used by this 
         * class as well. So for eg, if 'Heat source/sink due to sorption/desorption' is set to false, computing heat due to sorption in CL will return zero.
         *
         * <h3> Usage </h3>
         * In order to use this classes, first add them to your application as an object.
         * @code
         * #include "postprocessing/response_sorption_heat.h"
         * 
         * (...)
         * 
         *  FuelCellShop::PostProcessing::SorptionHeatResponse<dim> sorptionHeat;
         * @endcode
         * 
         * In the class constructor, construct the class passing SystemManagement and SorptionSourceTerms object. 
         * These objects are used in order to find the solution variables as appropriate, and flag parameters.
         * @code
         * //---------------------------------------------------------------------------
         * template <int dim>
         * NAME::AppPemfcNIThermal<dim>::AppPemfcNIThermal(boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data)
         * :
         * OptimizationBlockMatrixApplication<dim>(data),
         * system_management(this->block_info, this->cell_couplings, this->flux_couplings),
         * sorption_source_terms(system_management)
         * sorptionHeat(system_management, &sorption_source_terms)
         * {}
         * @endcode
         * 
         * Note that unlike some response classes, this class is not required to declare parameters,
         * as it is already deriving flag parameters from SorptionSourceTerms object.
         
         * Next, the object has to be initialized once SystemManagement and SorptionSourceTerms objects have
         * already been initialized:
         * @code
         * // Initialize parameters in system management:
         * system_management.initialize(param);
         * sorption_source_terms.initialize(param);
         * //Initialize post-processing routines:
         * sorptionHeat.initialize(param);
         * @endcode
         * 
         * Finally, the object is ready for use in cell_responses in your application:
         * @code
         * // Compute Total heat due to sorption generated in the CCL
         *    std::map<FuelCellShop::PostProcessing::ResponsesNames, double> respMap;
         * if (CCL->belongs_to_material(material_id)) //the material is the cathode catalyst layer
         *        sorptionHeat.compute_responses(cellInfo, CCL.get(), respMap);            
         * @endcode
         * 
         * \author Madhur Bhaiya, 2014
         */
        template <int dim>
        class SorptionHeatResponse : public BaseResponse<dim>
        {
        public:            
            ///@name Constructor, declaration and initialization
            //@{
            SorptionHeatResponse(const FuelCell::SystemManagement& sm,
                                 const FuelCellShop::Equation::SorptionSourceTerms<dim>* sst)
            :
            BaseResponse<dim>(sm),
            sorption_source(sst)
            {}
            
            ~SorptionHeatResponse() {}
            
            /**
             * Initialize class.
             * This function basically checks for various pre-requisites such as, certain solution 
             * variables are solved for in the application.
             */
            void initialize(ParameterHandler& param);
            
            //@}
            ///@name Compute functional
            //@{
            /**
             * This member function computes the heat generated due to sorption inside the layer.
             *
             * In order to access the sorption heat generated at the layer use:
             * @code
             * respMap[FuelCellShop::PostProcessing::ResponsesNames::sorption_heat]
             * @endcode
             *
             */
            void compute_responses(const typename DoFApplication<dim>::CellInfo& info,
                                   FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                   std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& respMap) const;
                                   
            /**
             * Routine used in order to compute the response with a modified solution (not the one stored
             * in CellInfo object.)
             *
             * \note Currently NOT IMPLEMENTED.
             */                       
            void compute_responses(std::vector< FuelCellShop::SolutionVariable > solution_variables,
                                   const typename DoFApplication<dim>::CellInfo& info,
                                   FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                   std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& respMap) const
            {
                throw std::runtime_error("SorptionHeatResponse::compute_responses(solution_variables, info, layer, respMap) not implemented");
            }
            //@}
            //
        private:
            /**
             * Pointer to SorptionSourceTerms object.
             */
            const FuelCellShop::Equation::SorptionSourceTerms<dim>* sorption_source;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "water_molar_fraction".
             */            
            FuelCellShop::Equation::VariableInfo xWater;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "membrane_water_content".
             */            
            FuelCellShop::Equation::VariableInfo lambda;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "temperature_of_REV".
             */            
            FuelCellShop::Equation::VariableInfo tRev;
            
            /**
             * Time constant, \f$ k \f$ [\p 1/s]
             */
            double time_constant;
        };
    }   
}

#endif