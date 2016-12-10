// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: simulation_selector.h
// - Description: This class selects an openFCST application which will run
// - Developers: P. Dobson,
//               M. Secanell,
//               A. Koupaei,
//               V. Zingan,
//               M. Bhaiya,
//               M. Sabharwal
//
// ----------------------------------------------------------------------------

#ifndef _SIMULATION_SELECTOR_H_
#define _SIMULATION_SELECTOR_H_

#include <string>
#include <iostream>

#include <boost/shared_ptr.hpp>
#include <application_core/application_data.h>
#include <application_core/optimization_block_matrix_application.h>

#include <solvers/adaptive_refinement.h>
#include <solvers/newton_basic.h>
#include <solvers/newton_w_line_search.h>
#include <solvers/newton_w_3pp.h>
#include <solvers/picard.h>

/////////////////////////////////////////////////////////////////
// FUEL CELL APPLICATIONS WITH DIFFUSION-FICKS-BASED TRANSPORT //
/////////////////////////////////////////////////////////////////

#include <applications/app_cathode.h>
#include <applications/app_pemfc.h>
#include <applications/app_pemfc_nonisothermal.h>
#include <applications/app_pemfc_twophase_saturation.h>
#include <applications/app_thermal_test.h>
#include <applications/app_test.h>
#include <applications/app_diffusion.h>
#include <applications/app_ohmic.h>

////////////////////////
// OTHER APPLICATIONS //
////////////////////////

#include <applications/app_read_mesh.h>

using namespace boost;

/**
 * This class selects an openFCST application which will run.
 *
 * - Add the *.h file of your app,
 * - Add the name of your app to the list in \p get_simulator_names() function,
 * - Add the specification of your app to the list in \p get_simulator_specifications() function (optional),
 * - Add the \p shared_ptr to your app to the list of \p shared_ptrs in \p select_application() function using the same pattern.
 *
 * @author P. Dobson, M. Secanell, A. Koupaei, V. Zingan, M. Bhaiya, ESDLab 2006-2015
 */

template<int dim>
class SimulationSelector
{
public:

       /**
        * Constructor.
        */
       SimulationSelector(boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = 
            boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >());

       /**
        * Destructor.
        */
      ~SimulationSelector();

       /**
        * Declare parameters.
        */
       void declare_parameters(ParameterHandler& param) const;

       /**
        * Initialize parameters.
        */
       void initialize(ParameterHandler& param);

       /**
        * Assign the application you would like to solve. All applications should be inherited from
        * \p FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>.
        *
        * If you develop a new application, you would assign the application an application name and add
        * the following to this section:
        * \code
        *         else if(name_application.compare("app_cathode") == 0)
        *         {
        *                FcstUtilities::log << "YOU ARE CURRENTLY SOLVING A CATHODE MODEL" << std::endl;
        *                return shared_ptr< FuelCell::Application::AppCathode<dim> > (new FuelCell::Application::AppCathode<dim>);
        *         }
        * \endcode
        */
       boost::shared_ptr< FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim> > select_application();

       /**
        * Select solver.
        */
       boost::shared_ptr< FuelCell::ApplicationCore::ApplicationWrapper > select_solver(FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>* app_lin);

       /**
        * Select the solution method.
        * Only adaptive refinement is available.
        */
       boost::shared_ptr< FuelCell::ApplicationCore::AdaptiveRefinement<dim> > select_solver_method(FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>* app_lin,
                                                                                                    FuelCell::ApplicationCore::ApplicationWrapper*                      newton_solver,
                                                                                                    const FuelCell::ApplicationCore::FEVector&                          solution = FuelCell::ApplicationCore::FEVector());

protected:

       /**
        * This function forms the string of names.
        */
       const std::string get_simulator_names() const
       {
              std::stringstream result;

              result << "cathode"
                     << " | "
                     << "anode"
                     << " | "
                     << "cathode_MPL"
                     << " | "
                     << "cathodeNIT"
                     << " | "
                     << "MEA"
                     << " | "
                     << "meaNIT"
                     << " | "
                     << "meaTwoPhaseSaturationNIT"
                     << " | "
                     << "test"
                     << " | "
                     << "thermalTest"
                     << " | "
                     << "test_mesh"
                     << " | "
                     << "diffusion"
                     << " | "
                     << "ohmic"
                     << " | "
                     << "meaTwoPhaseNITcapillary"
                     << " | "
                     << "Capillary_Testing";

              return result.str();
       }

       /**
        * This function forms the string of names.
        */
       const std::string get_simulator_specifications() const
       {
              std::stringstream result;

              result << "None"
                     << " | "
                     << "sphere"
                     << " | "
                     << "with_channel"
                     << " | "
                     << "without_channel"
                     << " | "
                     << "reaction"
                     << "|"
                     << "knudsen"
                     << "|"
                     << "reaction_and_knudsen";

              return result.str();
       }

       /**
        * This function forms the string of names.
        */
       const std::string get_nonlinear_solver_names() const
       {
              std::stringstream result;

              result << "None"
                     << " | "
                     << "NewtonBasic"
                     << " | "
                     << "NewtonLineSearch"
                     << " | "
                     << "Newton3pp"
		     << " | "
                     << "Picard";

              return result.str();
       }

       /**
        * This function forms the string of names.
        */
       const std::string get_refinement_methods() const
       {
              std::stringstream result;

              result << "AdaptiveRefinement";

              return result.str();
       }
       

       //////////
       // DATA //
       //////////
       
       /** 
        * Data structure storing information to be shared between applications 
        */
       boost::shared_ptr <FuelCell::ApplicationCore::ApplicationData> data;

       /**
        * Name of application.
        */
       std::string name_application;

       /**
        * Variable storing the name of the concrete application to be solved
        * from the broader class of applications defined by \p name_application.
        * Example: \p name_application = \p app_incompressible_flows and \p app_specification = \p Poiseuille.
        */
       std::string app_specification;

       /**
        * The name of a nonlinear solver.
        */
       std::string name_nonlinear_solver;

       /**
        * Name of refinement method.
        */
       std::string name_refinement_method;
};

#endif