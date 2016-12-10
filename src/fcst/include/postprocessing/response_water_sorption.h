// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: response_water_sorption.h
// - Description: This file defines response class computing amount of water sorbed in the catalyst layers.
// - Developers: Marc Secanell Gallart and Madhur Bhaiya
// - $Id:
//
// ----------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__RESPONSE_WATER_SORPTION_H
#define _FUELCELLSHOP__RESPONSE_WATER_SORPTION_H

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
         * Class used to calculate the amount of water sorbed inside the catalyst layer. The amount is reported in
         * mol H2O/(s). Most applications then divide by the cell active area based on the given geometry, so if not 
         * stated, the value should be assumed to be in mol H2O/(s cm2).
         * 
         * This class returns the total water amount sorbed inside the layer, \em i.e.
         * \f[
         * \int_{\Omega} k \frac{\rho_{dry}}{EW} (\lambda_{eq} - \lambda) d\Omega
         * \f]
         *
         * If the response value is negative, it implies that water is desorbed from the catalyst layer. Applications
         * can divide the obtained total quantity by surface area, volume etc. of the layer, in order to determine 
         * surface density, volumetric density etc.
         *
         * \remarks
         * - This class will not work if anything other than the CL layer object is passed as argument.
         * - FuelCellShop::Equation::SorptionSourceTerms object reference is required in the constructor, besides
         * the FuelCell::SystemManagement object reference.
         * - \f$ \lambda \f$ and \f$ x_{H_2O} \f$, should be solved for in the application.
         *
         *
         * <h3> Usage </h3>
         * In order to use this classes, first add them to your application as an object.
         * @code
         * #include "postprocessing/response_water_sorption.h"
         * 
         * (...)
         * 
         *  FuelCellShop::PostProcessing::WaterSorptionResponse<dim> waterSorption;
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
         * waterSorption(system_management, &sorption_source_terms)
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
         * waterSorption.initialize(param);
         * @endcode
         * 
         * Finally, the object is ready for use in cell_responses in your application:
         * @code
         * // Compute Total amount of water sorbed in the CCL
         *    std::map<FuelCellShop::PostProcessing::ResponsesNames, double> respMap;
         * if (CCL->belongs_to_material(material_id)) //the material is the cathode catalyst layer
         *        waterSorption.compute_responses(cellInfo, CCL.get(), respMap);            
         * @endcode
         * 
         * \author M. Secanell, 2014
         * \author Madhur Bhaiya, 2014
         */
        template <int dim>
        class WaterSorptionResponse : public BaseResponse<dim>
        {
        public:            
            ///@name Constructor, declaration and initialization
            //@{
            WaterSorptionResponse(const FuelCell::SystemManagement& sm,
                                  const FuelCellShop::Equation::SorptionSourceTerms<dim>* sst)
            :
            BaseResponse<dim>(sm),
            sorption_source(sst)
            {}
            
            ~WaterSorptionResponse() {}
            
            /**
             * Initialize class parameters.
             */
            void initialize(ParameterHandler& param) 
            {                
                if ( this->system_management->solution_in_userlist("membrane_water_content") )
                {
                    lambda.solution_index = this->system_management->solution_name_to_index("membrane_water_content"); 
                    lambda.fetype_index = this->system_management->block_info->base_element[lambda.solution_index];
                    lambda.indices_exist = true;
                }
                else
                    throw std::runtime_error("membrane_water_content variable required for WaterSorptionResponse");
                
                if ( this->system_management->solution_in_userlist("water_molar_fraction") )
                {
                    x_w.solution_index = this->system_management->solution_name_to_index("water_molar_fraction"); 
                    x_w.fetype_index = this->system_management->block_info->base_element[x_w.solution_index];
                    x_w.indices_exist = true;
                }
                else
                    throw std::runtime_error("water_molar_fraction variable required for WaterSorptionResponse");
                    
                if ( this->system_management->solution_in_userlist("temperature_of_REV") )
                {
                    tRev.solution_index = this->system_management->solution_name_to_index("temperature_of_REV"); 
                    tRev.fetype_index = this->system_management->block_info->base_element[tRev.solution_index];
                    tRev.indices_exist = true;
                }
                
                time_k = sorption_source->get_time_constant();
            }
            //@}
            
            ///@name Compute functional
            //@{
            /**
             * This member function computes the water adsorbed/desorbed from the electrolyte in the
             * catalyst layer. In order to compute the water sorption, the following functional is
             * evaluated
             * \f[
             * I = int_{\Omega} (\frac{(k_t*\rho_{Dry})}{EW})(\lambda_{eq}(a_w) - \lambda) d \Omega
             * \f]
             *             
             * In order to access the sorbed water at the electrode use:
             * @code
             * ORR_responses[FuelCellShop::PostProcessing::ResponsesNames::sorbed_water]
             * @endcode
             * 
             */
            void compute_responses(const typename DoFApplication<dim>::CellInfo& info,
                                   FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                   std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& respMap) const
            {
                respMap.clear();
                
                // Make sure you are in a CL
                const std::type_info& base_layer = layer->get_base_type();
                const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
                AssertThrow(base_layer == CatalystLayer,
                            ExcMessage("WaterSorptionResponse can only be used with a CatalystLayer object"));
                                
                // Get number of quadrature points in the cell:
                unsigned int n_q_points_cell = (info.fe(lambda.fetype_index)).n_quadrature_points;
                
                // Find where the solution is stored in info:
                unsigned int solIndex = info.global_data->find_vector("Solution");
                
                // Create CL:
                FuelCellShop::Layer::CatalystLayer<dim>* catalyst_layer = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);
                
                double rhoDry = catalyst_layer->get_electrolyte()->get_density();
                double EW = catalyst_layer->get_electrolyte()->get_EW();
                catalyst_layer->get_electrolyte()->set_water_molar_fraction( FuelCellShop::SolutionVariable(&info.values[solIndex][x_w.solution_index], water_molar_fraction) );
                if (tRev.indices_exist)
                    catalyst_layer->get_electrolyte()->set_temperature( FuelCellShop::SolutionVariable(&info.values[solIndex][tRev.solution_index], temperature_of_REV) );
                
                std::vector<double> lambdaValue( info.values[solIndex][lambda.solution_index] );
                
                std::vector<double> lambdaEq(n_q_points_cell, 0.0);
                catalyst_layer->get_electrolyte()->sorption_isotherm(lambdaEq);
                
                for (unsigned int q = 0; q < n_q_points_cell; ++q){
                    double JxW = info.fe(lambda.fetype_index).JxW(q);
                    respMap[FuelCellShop::PostProcessing::ResponsesNames::sorbed_water] += ( ( ((time_k*rhoDry)/EW)*(lambdaEq[q] - lambdaValue[q]) ) * JxW );       
                }
            }
            /**
             * Routine used in order to compute the response with a modified solution (not the one stored
             * in info)
             * 
             */                       
            void compute_responses(std::vector< FuelCellShop::SolutionVariable > solution_variables,
                                   const typename DoFApplication<dim>::CellInfo& info,
                                   FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                   std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& resp) const
            {
                throw std::runtime_error("WaterSorptionResponse::compute_responses(solution_variables, info, layer, resp) not implemented");
            }
            //
        private:
            /**
             * Pointer to SorptionSourceTerms object.
             */
            const FuelCellShop::Equation::SorptionSourceTerms<dim>* sorption_source;
            
            /**
             * Rate of sorption/desorption. Specified by SorptionSourceTerms in parameter file section
             * Sorption Source Terms >> Water soption time constant [1/s]
             */
            double time_k;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "membrane_water_content".
             */
            FuelCellShop::Equation::VariableInfo lambda;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "water_molar_fraction".
             */
            FuelCellShop::Equation::VariableInfo x_w;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "temperature_of_REV".
             */            
            FuelCellShop::Equation::VariableInfo tRev;
        };
        
    }   
}

#endif