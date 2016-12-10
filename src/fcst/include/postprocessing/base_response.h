// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: base_response.h
// - Description: This is a base class for all available FCST response evaluators
// - Developers: Marc Secanell Gallart,    University of Alberta
// - $Id: base_response.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _FUELCELLSHOP__BASE__RESPONSE_H
#define _FUELCELLSHOP__BASE__RESPONSE_H

// Include deal.II classes
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_values.h>


//Include STL
#include <cmath>
#include <iostream>

// Include OpenFCST routines:
#include <application_core/fcst_variables.h>
#include <application_core/system_management.h>
#include <equations/equation_base.h>

using namespace dealii;

namespace FuelCellShop
{
    /**
     * Namespace used for all classes use for post-processing either for evaluating
     * a new quantity at a quadrature point for outputing with the solution,
     * for the evaluation of functions or for the evaluation of functionals
     * such as the current density.
     * 
     */
    namespace PostProcessing
    {
        /**
         * Enumeration with names for different responses:
         */
        enum ResponsesNames
        {
            nothing = 0,
            
            ORR_current,
            HOR_current,
            O_coverage,
            OH_coverage,
            sorbed_water,
            water_cathode,
            water_anode,
            
            ORR_reaction_heat,
            HOR_reaction_heat,
            ORR_irrev_heat,
            HOR_irrev_heat,
            ORR_rev_heat,
            HOR_rev_heat,
            ORR_watervap_heat,
            sorption_heat,
            electron_ohmic_heat,
            proton_ohmic_heat,
            PEM_proton_ohmic_heat
        };
        
        /**
         * Virtual class used to develop a common interface to a set of functions used
         * to evaluate functionals that are obtained at postprocessing.
         * 
         * Examples of functionals would be the current density, the water sorbed in the 
         * catalyst layer by the ionomer, the friction factor at the channel, etc.
         * 
         * <h3> Usage </h3>
         * In order to use this classes, first add them to your application as an object.
         * @code
         * #include "postprocessing/responses.h"
         * 
         * (...)
         * 
         *  FuelCellShop::PostProcessing::ORRCurrentDensityResponse<dim> ORRCurrent;
         * @endcode
         * 
         * In the class constructor, construct the class passing an object SystemManagement. This object is used in 
         * order to find the solution variables as appropriate.
         * @code
         * //---------------------------------------------------------------------------
         * template <int dim>
         * NAME::AppCathode<dim>::AppCathode(boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data)
         * :
         * OptimizationBlockMatrixApplication<dim>(data),
         * system_management(this->block_info, this->cell_couplings, this->flux_couplings),
         * ORRCurrent(system_management)
         * {}
         * @endcode
         * 
         * Next, the object has to be initialized once SystemManagement has already been initialized
         * correctly:
         * @code
         * // Initialize parameters in system management:
         * system_management.initialize(param);
         * //Initialize post-processing routines:
         * ORRCurrent.initialize(param);
         * system_management.make_cell_couplings(tmp);
         * this->remesh_dofs();
         * this->remesh_matrices();
         * @endcode
         * 
         * Finally, the object is ready for use in cell_responses in your application:
         * @code
         * / Compute ORR responses in the CL
         *    std::map<FuelCellShop::PostProcessing::ResponsesNames, double> ORR_responses;
         * if (CCL->belongs_to_material(material_id)) //the material is the catalyst layer 
         *        ORRCurrent.compute_responses(info, CCL.get(), ORR_responses);            
         * @endcode
         * 
         * \author M. Secanell, 2014
         */
        template <int dim>
        class BaseResponse : public Subscriptor
        {
        
        protected:
            /**
             * Constructor 
             * 
             * Usually, all constructors will contain a pointer to SystemManagement and a pointer to a child
             * of BaseLayer
             */
            BaseResponse(const FuelCell::SystemManagement& sm)
            :
            system_management(&sm)
            {}          
            
            /**
             *  Destructor
             */
            virtual ~BaseResponse(){}          
            
            ///@name Initalization
            //@{
            /** 
             * Declare any necessary parameters to compute the functional.
             */
            virtual void declare_parameters(ParameterHandler& param) const {}
            
            /**
             * Initialize class parameters.
             * 
             * @warning This routine has to be called once SystemManagement has been initialized and
             * all blocks in the matrix have been declared, i.e. after calling app->remesh_dofs();
             */
            virtual void initialize(ParameterHandler& param) {}
            //@}
            
            ///@name Compute functional
            //@{ 
            /**
             * Compute several functionals (response) at a given finite element using the information in CellInfo
             * in layer specified in parameter @param layer.
             * 
             * Responses are stored in resp which return a map with the name of the response and its value
             * at the given either cell or boundary.
             */
            virtual void compute_responses(const typename DoFApplication<dim>::CellInfo& info,
                                           FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                           std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& resp) const = 0;
           /**
            * In rare instances, the user would like to specify a given solution vector that is not the final
            * solution to the problem. This member function is built for such instances. Instead of using CellInfo to
            * create automatically the SolutionVariable vector, in this case the user is allowed to provide this 
            * function.
            */ 
           virtual void compute_responses(std::vector< FuelCellShop::SolutionVariable > solution_variables,
                                          const typename DoFApplication<dim>::CellInfo& info,
                                           FuelCellShop::Layer::BaseLayer<dim>* const layer, 
                                           std::map<FuelCellShop::PostProcessing::ResponsesNames, double>& resp) const = 0;
            //@}
                                           
            /** Pointer to system management */
            const FuelCell::SystemManagement* system_management;   

        };
        
    }   
}

#endif