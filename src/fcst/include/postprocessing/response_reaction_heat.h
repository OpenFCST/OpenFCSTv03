// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: response_reaction_heat.h
// - Description: This is header file for response class computing heat generated due to electrochemical reactions.
// - Developers: Madhur Bhaiya
// - $Id: 
//
// ----------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__RESPONSE_REACTION_HEAT_H
#define _FUELCELLSHOP__RESPONSE_REACTION_HEAT_H

//Include STL
#include <exception>      // std::exception

// Include OpenFCST routines:
#include <equations/reaction_source_terms.h>
#include <layers/catalyst_layer.h>
#include <postprocessing/base_response.h>

using namespace dealii;

namespace FuelCellShop
{

    namespace PostProcessing
    {

        /**
         * Class used to calculate the heat generated due to ORR inside the cathode catalyst layer.
         * 
         * For mathematical formulation details of various heat terms, \em viz., reversible, irreversible and 
         * water vaporization sink, please look at documentation of FuelCellShop::Equation::ReactionHeat.
         * Also refer: M. Bhaiya, A. Putz and M. Secanell, "Analysis of non-isothermal effects on polymer electrolyte fuel cell electrode assemblies", Electrochimica Acta, 147C:294-309, 2014.
         * DOI: http://dx.doi.org/10.1016/j.electacta.2014.09.051
         * 
         * Applications can divide the obtained total quantity by surface area, volume etc. of the layer, in order to determine 
         * surface density, volumetric density etc.
         *
         * \remarks
         * - This class will not work if anything other than the CL layer object is passed as argument.
         * - FuelCellShop::Equation::ReactionSourceTerms object reference is required in the constructor, besides the FuelCell::SystemManagement object reference.
         * - \f$ \phi_s \f$, \f$ \phi_m \f$, \f$ x_{O_2} \f$ and \f$ T_{rev} \f$ should be solved for in the application.
         * - Flag parameters from the ReactionSourceTerms, eg, 'Irreversible heat source due to ORR' is Enabled or Not, are used by this 
         * class as well. So for eg, if 'Irreversible heat source due to ORR' is set to false, computing irreversible heat due to ORR will return zero.
         *
         * <h3> Usage </h3>
         * In order to use this classes, first add them to your application as an object.
         * @code
         * #include "postprocessing/response_reaction_heat.h"
         * 
         * (...)
         * 
         *  FuelCellShop::PostProcessing::ORRReactionHeatResponse<dim> catReactionHeat;
         * @endcode
         * 
         * In the class constructor, construct the class passing SystemManagement and ReactionSourceTerms object. 
         * These objects are used in order to find the solution variables as appropriate, and flag parameters.
         * @code
         * //---------------------------------------------------------------------------
         * template <int dim>
         * NAME::AppPemfcNIThermal<dim>::AppPemfcNIThermal(boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data)
         * :
         * OptimizationBlockMatrixApplication<dim>(data),
         * system_management(this->block_info, this->cell_couplings, this->flux_couplings),
         * reaction_source_terms(system_management)
         * catReactionHeat(system_management, &reaction_source_terms)
         * {}
         * @endcode
         * 
         * Note that unlike some response classes, this class is not required to declare parameters,
         * as it is already deriving flag parameters from ReactionSourceTerms object.
         
         * Next, the object has to be initialized once SystemManagement and ReactionSourceTerms objects have
         * already been initialized:
         * @code
         * // Initialize parameters in system management:
         * system_management.initialize(param);
         * reaction_source_terms.set_cathode_kinetics(CCL->get_kinetics());
         * reaction_source_terms.initialize(param);
         * //Initialize post-processing routines:
         * catReactionHeat.initialize(param);
         * @endcode
         * 
         * Finally, the object is ready for use in cell_responses in your application:
         * @code
         * // Compute reaction heat terms due to ORR generated in the CCL
         *    std::map<FuelCellShop::PostProcessing::ResponsesNames, double> respMap;
         * if (CCL->belongs_to_material(material_id)) //the material is the cathode catalyst layer
         *        catReactionHeat.compute_responses(cellInfo, CCL.get(), respMap);            
         * @endcode
         * 
         * \author Madhur Bhaiya, 2014
         */
        template <int dim>
        class ORRReactionHeatResponse : public BaseResponse<dim>
        {
        public:            
            ///@name Constructor, declaration and initialization
            //@{
            ORRReactionHeatResponse(const FuelCell::SystemManagement& sm,
                                    const FuelCellShop::Equation::ReactionSourceTerms<dim>* rst)
            :
            BaseResponse<dim>(sm),
            reaction_source(rst),
            reaction_heat(NULL)
            {}
            
            ~ORRReactionHeatResponse()
            {
                delete reaction_heat;
            }
            
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
             * This member function computes the heat generated due to ORR inside the cathode catalyst layer.
             *
             * In order to access various ORR reactions heat terms:
             * @code
             * // To access total reaction heat generated due to ORR
             * respMap[FuelCellShop::PostProcessing::ResponsesNames::ORR_reaction_heat]
             *
             * // To access irreversible heat generated due to ORR
             * respMap[FuelCellShop::PostProcessing::ResponsesNames::ORR_irrev_heat]
             *
             * // To access reversible heat generated due to ORR
             * respMap[FuelCellShop::PostProcessing::ResponsesNames::ORR_rev_heat]
             *
             * // To access water vaporization heat sink due to ORR
             * respMap[FuelCellShop::PostProcessing::ResponsesNames::ORR_watervap_heat]
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
                throw std::runtime_error("ORRReactionHeatResponse::compute_responses(solution_variables, info, layer, respMap) not implemented");
            }
            //@}
            //
        private:
            /**
             * Pointer to ReactionSourceTerms object.
             */
            const FuelCellShop::Equation::ReactionSourceTerms<dim>* reaction_source;
            
            /**
             * Pointer to BaseKinetics object.
             */
            FuelCellShop::Kinetics::BaseKinetics* kinetics;
            
            /**
             * Pointer to ReactionHeat object. This is initialized automatically
             * inside this class.
             */
            FuelCellShop::Equation::ReactionHeat* reaction_heat;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "oxygen_molar_fraction".
             */            
            FuelCellShop::Equation::VariableInfo xOxygen;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "electronic_electrical_potential".
             */            
            FuelCellShop::Equation::VariableInfo phiS;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "protonic_electrical_potential".
             */            
            FuelCellShop::Equation::VariableInfo phiM;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "temperature_of_REV".
             */            
            FuelCellShop::Equation::VariableInfo tRev;
        };
        
        
        //===============================================================================================================
        //===============================================================================================================
        //===============================================================================================================
        //===============================================================================================================
        
        
        /**
         * Class used to calculate the heat generated due to HOR inside the anode catalyst layer.
         * 
         * For mathematical formulation details of various heat terms, \em viz., reversible and irreversible, please look at documentation of FuelCellShop::Equation::ReactionHeat.
         * Also refer: M. Bhaiya, A. Putz and M. Secanell, "Analysis of non-isothermal effects on polymer electrolyte fuel cell electrode assemblies", Electrochimica Acta, 147C:294-309, 2014.
         * DOI: http://dx.doi.org/10.1016/j.electacta.2014.09.051
         * 
         * Applications can divide the obtained total quantity by surface area, volume etc. of the layer, in order to determine 
         * surface density, volumetric density etc.
         *
         * \remarks
         * - This class will not work if anything other than the CL layer object is passed as argument.
         * - FuelCellShop::Equation::ReactionSourceTerms object reference is required in the constructor, besides the FuelCell::SystemManagement object reference.
         * - \f$ \phi_s \f$, \f$ \phi_m \f$, \f$ T_{rev} \f$, and (\f$ x_{H_2} \f$ OR \f$ x_{H_2O} \f$) should be solved for in the application.
         * - Flag parameters from the ReactionSourceTerms, eg, 'Irreversible heat source due to HOR' is Enabled or Not, are used by this 
         * class as well. So for eg, if 'Irreversible heat source due to HOR' is set to false, computing irreversible heat due to HOR will return zero.
         *
         * <h3> Usage </h3>
         * In order to use this classes, first add them to your application as an object.
         * @code
         * #include "postprocessing/response_reaction_heat.h"
         * 
         * (...)
         * 
         *  FuelCellShop::PostProcessing::HORReactionHeatResponse<dim> anReactionHeat;
         * @endcode
         * 
         * In the class constructor, construct the class passing SystemManagement and ReactionSourceTerms object. 
         * These objects are used in order to find the solution variables as appropriate, and flag parameters.
         * @code
         * //---------------------------------------------------------------------------
         * template <int dim>
         * NAME::AppPemfcNIThermal<dim>::AppPemfcNIThermal(boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data)
         * :
         * OptimizationBlockMatrixApplication<dim>(data),
         * system_management(this->block_info, this->cell_couplings, this->flux_couplings),
         * reaction_source_terms(system_management)
         * anReactionHeat(system_management, &reaction_source_terms)
         * {}
         * @endcode
         * 
         * Note that unlike some response classes, this class is not required to declare parameters,
         * as it is already deriving flag parameters from ReactionSourceTerms object.
         
         * Next, the object has to be initialized once SystemManagement and ReactionSourceTerms objects have
         * already been initialized:
         * @code
         * // Initialize parameters in system management:
         * system_management.initialize(param);
         * reaction_source_terms.set_anode_kinetics(ACL->get_kinetics());
         * reaction_source_terms.initialize(param);
         * //Initialize post-processing routines:
         * anReactionHeat.initialize(param);
         * @endcode
         * 
         * Finally, the object is ready for use in cell_responses in your application:
         * @code
         * // Compute reaction heat terms due to HOR generated in the CCL
         *    std::map<FuelCellShop::PostProcessing::ResponsesNames, double> respMap;
         * if (ACL->belongs_to_material(material_id)) //the material is the anode catalyst layer
         *        anReactionHeat.compute_responses(cellInfo, ACL.get(), respMap);            
         * @endcode
         * 
         * \author Madhur Bhaiya, 2014
         */
        template <int dim>
        class HORReactionHeatResponse : public BaseResponse<dim>
        {
        public:            
            ///@name Constructor, declaration and initialization
            //@{
            HORReactionHeatResponse(const FuelCell::SystemManagement& sm,
                                    const FuelCellShop::Equation::ReactionSourceTerms<dim>* rst)
            :
            BaseResponse<dim>(sm),
            reaction_source(rst),
            reaction_heat(NULL)
            {}
            
            ~HORReactionHeatResponse()
            {
                delete reaction_heat;
            }
            
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
             * This member function computes the heat generated due to HOR inside the anode catalyst layer.
             *
             * In order to access various HOR reactions heat terms:
             * @code
             * // To access total reaction heat generated due to HOR
             * respMap[FuelCellShop::PostProcessing::ResponsesNames::HOR_reaction_heat]
             *
             * // To access irreversible heat generated due to HOR
             * respMap[FuelCellShop::PostProcessing::ResponsesNames::HOR_irrev_heat]
             *
             * // To access reversible heat generated due to HOR
             * respMap[FuelCellShop::PostProcessing::ResponsesNames::HOR_rev_heat]
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
                throw std::runtime_error("HORReactionHeatResponse::compute_responses(solution_variables, info, layer, respMap) not implemented");
            }
            //@}
            //
        private:
            /**
             * Pointer to ReactionSourceTerms object.
             */
            const FuelCellShop::Equation::ReactionSourceTerms<dim>* reaction_source;
            
            /**
             * Pointer to BaseKinetics object.
             */
            FuelCellShop::Kinetics::BaseKinetics* kinetics;
            
            /**
             * Pointer to ReactionHeat object. This is initialized automatically
             * inside this class.
             */
            FuelCellShop::Equation::ReactionHeat* reaction_heat;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "hydrogen_molar_fraction".
             */            
            FuelCellShop::Equation::VariableInfo xHydrogen;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "water_molar_fraction".
             */            
            FuelCellShop::Equation::VariableInfo xWater;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "electronic_electrical_potential".
             */            
            FuelCellShop::Equation::VariableInfo phiS;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "protonic_electrical_potential".
             */            
            FuelCellShop::Equation::VariableInfo phiM;
            
            /**
             * VariableInfo structure corresponding to the
             * \p "temperature_of_REV".
             */            
            FuelCellShop::Equation::VariableInfo tRev;
        };
    }   
}

#endif