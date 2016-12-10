// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: response_ohmic_heat.h
// - Description: This is header file for electronic and protonic ohmic heat response evaluator classes.
// - Developers: Madhur Bhaiya
// - $Id: 
//
// ----------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__RESPONSE_OHMIC_HEAT_H
#define _FUELCELLSHOP__RESPONSE_OHMIC_HEAT_H

//Include STL
#include <exception>      // std::exception

// Include OpenFCST routines:
#include <equations/thermal_transport_equation.h>
#include <postprocessing/base_response.h>

using namespace dealii;

namespace FuelCellShop
{

    namespace PostProcessing
    {

        /**
         * Class used to calculate the electronic ohmic heat generated in the electron conducting layers, \em viz., \b GDL, \b MPL and \b CL.
         * 
         * This class returns the total electronic ohmic heat produced in the layer, \em i.e.
         * \f[
         * \int_{\Omega} \hat{\sigma}_{s,eff} \cdot \left( \mathbf{\nabla} \phi_s \otimes  \mathbf{\nabla} \phi_s \right) d\Omega
         * \f]
         *
         * Applications can divide the obtained total quantity by surface area, volume etc. of the layer, in order to determine 
         * surface density, volumetric density etc.
         *
         * \remarks
         * - This class will not work if anything other than the GDL, MPL or CL layer objects are passed as argument.
         * - This class will not work if thermal transport equation is not being solved for. FuelCellShop::Equation::ThermalTransportEquation 
         * object reference is required in the constructor, besides the FuelCell::SystemManagement object reference.
         * - Electron transport should also be solved for, \em i.e. , \f$ \phi_s \f$, should be a solution variable.
         * - Flag parameters from the ThermalTransportEquation, for instance, 'Electronic ohmic heat in GDL' is Enabled or Not, are used by this 
         * class as well. So for eg, if 'Electronic ohmic heat in GDL' is set to false, computing electronic ohmic heat in GDL will return zero.
         *
         * <h3> Usage </h3>
         * In order to use this classes, first add them to your application as an object.
         * @code
         * #include "postprocessing/response_ohmic_heat.h"
         * 
         * (...)
         * 
         *  FuelCellShop::PostProcessing::ElectronOhmicHeatResponse<dim> electronOhmicHeat;
         * @endcode
         * 
         * In the class constructor, construct the class passing SystemManagement and ThermalTransportEquation object. 
         * These objects are used in order to find the solution variables as appropriate, and flag parameters.
         * @code
         * //---------------------------------------------------------------------------
         * template <int dim>
         * NAME::AppPemfcNIThermal<dim>::AppPemfcNIThermal(boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data)
         * :
         * OptimizationBlockMatrixApplication<dim>(data),
         * system_management(this->block_info, this->cell_couplings, this->flux_couplings),
         * thermal_equation(system_management)
         * electronOhmicHeat(system_management, &thermal_equation)
         * {}
         * @endcode
         * 
         * Note that unlike some response classes, this class is not required to declare parameters,
         * as it is already deriving flag parameters from ThermalTransportEquation object.
         
         * Next, the object has to be initialized once SystemManagement and ThermalTransportEquation objects have
         * already been initialized:
         * @code
         * // Initialize parameters in system management:
         * system_management.initialize(param);
         * thermal_equation.initialize(param);
         * //Initialize post-processing routines:
         * electronOhmicHeat.initialize(param);
         * @endcode
         * 
         * Finally, the object is ready for use in cell_responses in your application:
         * @code
         * // Compute Total Electronic Ohmic Heat generated in the CGDL
         *    std::map<FuelCellShop::PostProcessing::ResponsesNames, double> respMap;
         * if (CGDL->belongs_to_material(material_id)) //the material is the cathode gas diffusion layer 
         *        electronOhmicHeat.compute_responses(cellInfo, CGDL.get(), respMap);            
         * @endcode
         * 
         * \author Madhur Bhaiya, 2014
         */
        template <int dim>
        class ElectronOhmicHeatResponse : public BaseResponse<dim>
        {
        public:            
            ///@name Constructor, declaration and initialization
            //@{
            ElectronOhmicHeatResponse(const FuelCell::SystemManagement& sm,
                                      const FuelCellShop::Equation::ThermalTransportEquation<dim>* tte)
            :
            BaseResponse<dim>(sm),
            thermal_equation(tte),
            factor_GDL(0),
            factor_MPL(0),
            factor_CL(0)
            {}
            
            ~ElectronOhmicHeatResponse() {}
            
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
             * This member function computes the ohmic heat generated due to electron transport inside the layer.
             *
             * In order to access the electronic ohmic heat at the layer use:
             * @code
             * respMap[FuelCellShop::PostProcessing::ResponsesNames::electron_ohmic_heat]
             * @endcode
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
                throw std::runtime_error("ElectronOhmicHeatResponse::compute_responses(solution_variables, info, layer, respMap) not implemented");
            }
            //@}
            //
        private:
            /**
             * Pointer to ThermalTransportEquation object.
             */
            const FuelCellShop::Equation::ThermalTransportEquation<dim>* thermal_equation;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "electronic_electrical_potential".
             */            
            FuelCellShop::Equation::VariableInfo phiS;
            
            /**
             * Factor is 1 if electronic ohmic heating in GDL is enabled, else 0.
             */
            unsigned int factor_GDL;
            
            /**
             * Factor is 1 if electronic ohmic heating in MPL is enabled, else 0.
             */
            unsigned int factor_MPL;
            
            /**
             * Factor is 1 if electronic ohmic heating in CL is enabled, else 0.
             */
            unsigned int factor_CL;
        };
        
        
        //===============================================================================================================
        //===============================================================================================================
        //===============================================================================================================
        //===============================================================================================================
        
        
        /**
         * Class used to calculate the protonic ohmic heat generated in the proton conducting layers, \em viz., \b Membrane and \b CL.
         * 
         * This class returns the total protonic ohmic heat produced in the layer, \em i.e.
         * \f[
         * \int_{\Omega} \hat{\sigma}_{m,eff} \cdot \left( \mathbf{\nabla} \phi_m \otimes  \mathbf{\nabla} \phi_m \right) d\Omega
         * \f]
         *
         * Applications can divide the obtained total quantity by surface area, volume etc. of the layer, in order to determine 
         * surface density, volumetric density etc.
         *
         * \remarks
         * - This class will not work if anything other than the Membrane or CL layer objects are passed as argument.
         * - This class will not work if thermal transport equation is not being solved for. FuelCellShop::Equation::ThermalTransportEquation 
         * object reference is required in the constructor, besides the FuelCell::SystemManagement object reference.
         * - Proton transport should also be solved for, \em i.e. , \f$ \phi_m \f$, should be a solution variable.
         * - Flag parameters from the ThermalTransportEquation, for instance, 'Protonic ohmic heat in CL' is Enabled or Not, are used by this 
         * class as well. So for eg, if 'Protonic ohmic heat in CL' is set to false, computing protonic ohmic heat in CL will return zero.
         *
         * <h3> Usage </h3>
         * In order to use this classes, first add them to your application as an object.
         * @code
         * #include "postprocessing/response_ohmic_heat.h"
         * 
         * (...)
         * 
         *  FuelCellShop::PostProcessing::ProtonOhmicHeatResponse<dim> protonOhmicHeat;
         * @endcode
         * 
         * In the class constructor, construct the class passing SystemManagement and ThermalTransportEquation object. 
         * These objects are used in order to find the solution variables as appropriate, and flag parameters.
         * @code
         * //---------------------------------------------------------------------------
         * template <int dim>
         * NAME::AppPemfcNIThermal<dim>::AppPemfcNIThermal(boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data)
         * :
         * OptimizationBlockMatrixApplication<dim>(data),
         * system_management(this->block_info, this->cell_couplings, this->flux_couplings),
         * thermal_equation(system_management)
         * protonOhmicHeat(system_management, &thermal_equation)
         * {}
         * @endcode
         * 
         * Note that unlike some response classes, this class is not required to declare parameters,
         * as it is already deriving flag parameters from ThermalTransportEquation object.
         
         * Next, the object has to be initialized once SystemManagement and ThermalTransportEquation objects have 
         * already been initialized:
         * @code
         * // Initialize parameters in system management:
         * system_management.initialize(param);
         * thermal_equation.initialize(param);
         * //Initialize post-processing routines:
         * protonOhmicHeat.initialize(param);
         * @endcode
         * 
         * Finally, the object is ready for use in cell_responses in your application:
         * @code
         * // Compute Total Protonic Ohmic Heat generated in the CCL
         *    std::map<FuelCellShop::PostProcessing::ResponsesNames, double> respMap;
         * if (CCL->belongs_to_material(material_id)) //the material is the cathode catalyst layer 
         *        protonOhmicHeat.compute_responses(cellInfo, CCL.get(), respMap);            
         * @endcode
         * 
         * \author Madhur Bhaiya, 2014
         */
        template <int dim>
        class ProtonOhmicHeatResponse : public BaseResponse<dim>
        {
        public:            
            ///@name Constructor, declaration and initialization
            //@{
            ProtonOhmicHeatResponse(const FuelCell::SystemManagement& sm,
                                    const FuelCellShop::Equation::ThermalTransportEquation<dim>* tte)
            :
            BaseResponse<dim>(sm),
            thermal_equation(tte),
            factor_CL(0),
            factor_ML(0)
            {}
            
            ~ProtonOhmicHeatResponse() {}
            
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
             * This member function computes the ohmic heat generated due to proton transport inside the layer.
             *
             * In order to access the protonic ohmic heat at the layer use:
             * @code
             * respMap[FuelCellShop::PostProcessing::ResponsesNames::proton_ohmic_heat]
             * @endcode
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
                throw std::runtime_error("ProtonOhmicHeatResponse::compute_responses(solution_variables, info, layer, respMap) not implemented");
            }
            //@}
            //
        private:
            /**
             * Pointer to ThermalTransportEquation object.
             */
            const FuelCellShop::Equation::ThermalTransportEquation<dim>* thermal_equation;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "protonic_electrical_potential".
             */            
            FuelCellShop::Equation::VariableInfo phiM;
            
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
             * Factor is 1 if protonic ohmic heating in CL is enabled, else 0.
             */
            unsigned int factor_CL;
            
            /**
             * Factor is 1 if protonic ohmic heating in Membrane is enabled, else 0.
             */
            unsigned int factor_ML;
        };
    }   
}

#endif