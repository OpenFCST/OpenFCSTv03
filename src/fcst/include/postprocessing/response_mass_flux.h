// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: response_mass_flux.h
// - Description: This is header file for electronic and protonic ohmic heat response evaluator classes.
// - Developers: Chad Balen
// - $Id: 
//
// ----------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__RESPONSE_MASS_FLUX_H
#define _FUELCELLSHOP__RESPONSE_MASS_FLUX_H

#include <postprocessing/base_response.h>

using namespace dealii;

namespace FuelCellShop
{
    namespace PostProcessing
    {

        /**
         * Class used to calculate the total mass flux at a boundary.
         * 
         * This class returns the total mass flux at specified boundary IDs, \em i.e.
         * \f[
         * \oint_{\Gamma} \rho_i \boldsymbol{\hat{u}}_i \cdot \boldsymbol{n} \, \mathrm{d}\Gamma
         * \f]
         *
         * User can then divide the obtained total quantity by the length (2D simulation) or surface area (3D simulation) in order to determine 
         * average mass flux. 
         * 
         * \b Note: be careful if combining multiple boundary ids into one. For example if you have two outlets that have the same BCs and you call them by
         * the same boundary ID then calculating the average mass flux gets a bit more tricky, as you need to normalize the integral over each outlet. If the
         * outlets are the same dimensions then just divide by the length of one of the outlets (NOT the combined length of the outlets).
         *
         * \remarks
         * - This class will not work if density and velocity are not being solved for.
         * - At this moment only implemented for bdry_responses() in application class
         * 
         * <h3> Usage </h3>
         * In the input data file, the following parameters can be specified (see declare_parameters ):
         * @code
         * subsection Output Variables
         *     set Compute boundary responses = true
         *     set num_output_vars = 3
         *     set Output boundary id = 1, 2
         *     set Output_var_0 = total_mass_flux_species_1
         *     set Output_var_2 = total_mass_flux_species_2
         *     set Output_var_3 = total_mass_flux_species_3
         * end
         * @endcode
         * \note: you can specify up to 5 different total_mass_flux_species. If you require more add more to application_core/optimization_block_matrix_application.cc at the end of the file
         * 
         * <h3> Implementation in Applications </h3>
         * In order to use this classes, first add it to your application as an object.
         * @code
           #include "postprocessing/response_mass_flux.h"
           
           (...)
           
           FuelCellShop::PostProcessing::MassFluxResponse<dim> mass_flux_response;
           @endcode
         * 
         * In the class constructor, construct the class passing SystemManagement. 
         * This object is used in order to find the solution variables as appropriate, and flag parameters.
         * @code
           //---------------------------------------------------------------------------
           template<int dim>
           NAME::AppCompressibleFlows<dim>::AppCompressibleFlows( boost::shared_ptr<ApplicationData> data )
           :
           OptimizationBlockMatrixApplication<dim>(data),
           fluid_transport_equations(this->system_management, data),
           mass_flux_response(this->system_management),
           fluid("fluid"),
           CChannel("channel", fluid)
           { }
           @endcode
         * 
         * Note that unlike some response classes, this class does not require declare parameters,
         * as information is declared by the OptimizationBlockMatrixApplication class. initialize() is used
         * so the parameters do not need to be passed to MassFluxResponse.
         * 
         * Next, the object has to be initialized once SystemManagement and ThermalTransportEquation objects have
         * already been initialized:
         * @code
           template<int dim>
           void
           NAME::AppCompressibleFlows<dim>::initialize(ParameterHandler& param)
           {
                // Initialize post-processing routines:
                mass_flux_response.initialize(param);
           }
           @endcode
         * 
         * Finally, the object is ready for use in bdry_responses in your application:
         * @code
           template<int dim>
           void
           NAME::AppCompressibleFlows<dim>::bdry_responses(std::vector<double>&                                                     dst,
                                                           const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                           const FuelCell::ApplicationCore::FEVector&                               src)
           {
                mass_flux_response.bdry_responses(dst, bdry_info);
           }
           @endcode
         * 
         * \author Chad Balen, 2016
         */
        template <int dim>
        class MassFluxResponse : public BaseResponse<dim>
        {
        public:            
            ///@name Constructor, declaration and initialization
            //@{
            MassFluxResponse(const FuelCell::SystemManagement& sm)
            :
            BaseResponse<dim>(sm)
            {}
            
            ~MassFluxResponse() {}
            
            //@}
            ///@name Compute functional
            //@{
            /**
             * This member function computes the mass flux in a cell.
             *
             * \note Currently NOT IMPLEMENTED.
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
                                   std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& respMap) const;
            
            /**
             * This member function computes the total mass flux along specified boundaries.
             *
             * Call this function in the application class in the bdry_responses() function:
             * @code
               template<int dim>
               void
               NAME::AppCompressibleFlows<dim>::bdry_responses(std::vector<double>&                                                     dst,
                                                               const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                               const FuelCell::ApplicationCore::FEVector&                               src)
               {
                    mass_flux_response.bdry_responses(dst, bdry_info);
               }
               @endcode
             */
            void bdry_responses(std::vector<double>&                                                     dst,
                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                std::vector<std::string>&                                                responseNames,
                                std::vector<unsigned int>&                                               bdryIDs) const;
            
            //@}
            //
        private:
            ///@name Compute parameters
            //@{
            
            /**
             * private variable for storing the number of species to calculate the mass flux for. Set in initialize().
             */
            unsigned int numOutputVars;
            
            /**
             * private variable for storing a vector of the boundary ids to calculate mass flux along. Set in initalize();
             */
            std::vector<unsigned int> bdryIDs;
            
            /**
             * private variable for storing the names of species to calculate the mass flux for. This is important in case the user
             * has three species but only wants the mass flux for two of them. As well, this way the user can specify what order they
             * want the output in. Format for strings should be: "total_mass_flux_species_" + integer.
             */
            std::vector<std::string> nameOutputVars;
            
            //@}
        };
    }
}
#endif