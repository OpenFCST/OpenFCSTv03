// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: simulation_selector.cc
// - Description: This class selects an openFCST application which will run
// - Developers: P. Dobson,
//               M. Secanell,
//               A. Koupaei,
//               V. Zingan,
//               M. Bhaiya,
//               M. Sabharwal,
//               J. Zhou,
//               A. Kosakian
//
// ----------------------------------------------------------------------------

#include "utils/simulation_selector.h"

// ----------------------------------------------------------------------------
template<int dim>
SimulationSelector<dim>::SimulationSelector(boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data)
:
data(data)
{ }

// ----------------------------------------------------------------------------
template<int dim>
SimulationSelector<dim>::~SimulationSelector()
{ }

// ----------------------------------------------------------------------------
template<int dim>
void
SimulationSelector<dim>::declare_parameters(ParameterHandler& param) const
{
     param.enter_subsection("Simulator");
     {
          param.declare_entry("simulator name",
                              "cathode",
                               Patterns::Selection( get_simulator_names() ),
                              "Name of the application that you would like to run ");

          param.declare_entry("simulator specification",
                              "None",
                               Patterns::Selection( get_simulator_specifications() ),
                              "Select a sub-application ");

          param.declare_entry("nonlinear solver name",
                              "NewtonLineSearch",
                               Patterns::Selection( get_nonlinear_solver_names() ),
                              "Select the type of non-linear solver ");
          param.declare_entry("linear solver name",
          #ifdef OPENFCST_WITH_PETSC
                              "MUMPS",
                              Patterns::Selection("MUMPS|CG|Bicgstab|ILU-GMRES|UMFPACK"),
          #else
                              "UMFPACK",
                              Patterns::Selection("UMFPACK|CG|Bicgstab|ILU-GMRES|MUMPS"),
          #endif
                              "Select the linear solver you would like to use to solve the problem.");

          param.declare_entry("refinement method",
                              "AdaptiveRefinement",
                               Patterns::Selection( get_refinement_methods() ),
                              " ");
     }
     param.leave_subsection();
}

// ----------------------------------------------------------------------------
template<int dim>
void
SimulationSelector<dim>::initialize(ParameterHandler& param)
{
     param.enter_subsection("Simulator");
     {
          name_application  = param.get("simulator name");
          app_specification = param.get("simulator specification");
          data->set_nonlinear_solver(param.get("nonlinear solver name"));
          data->set_linear_solver(param.get("linear solver name"));
          data->set_refinement_solver(param.get("refinement method"));
     }
     param.leave_subsection();
}

// ----------------------------------------------------------------------------
template<int dim>
shared_ptr< FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim> >
SimulationSelector<dim>::select_application()
{
    data->enter_flag(app_specification, true);

    //-- Cathode without convection, pseudo-homogeneous
    if ( (name_application.compare("cathode") == 0) || (name_application.compare("anode") == 0))
    {
        FcstUtilities::log << "YOU ARE CURRENTLY SOLVING A CATHODE MODEL" << std::endl;
        return shared_ptr<FuelCell::Application::AppCathode<dim> > (new FuelCell::Application::AppCathode<dim>(data));
    }
    //-- Complete MEA without convection
    else if (name_application.compare("MEA") == 0)
    {
        FcstUtilities::log << "YOU ARE CURRENTLY SOLVING AN MEA MODEL" << std::endl;
        return shared_ptr<FuelCell::Application::AppPemfc<dim> > (new FuelCell::Application::AppPemfc<dim>(data));

    }
    //Run application defined in app_thermal_test.h
    else if (name_application.compare("thermalTest") == 0)
    {
        FcstUtilities::log << "YOU ARE CURRENTLY SOLVING A Thermal Transport Equation TEST CASE" << std::endl;
        return shared_ptr<FuelCell::Application::AppThermalTest<dim>> (new FuelCell::Application::AppThermalTest<dim>(data));
    }
    //Run application defined in app_pemfc_nonisothermal.h
    else if (name_application.compare("meaNIT") == 0)
    {
        FcstUtilities::log << "YOU ARE CURRENTLY SOLVING A NON-ISOTHERMAL MEA MODEL" << std::endl;
        return shared_ptr<FuelCell::Application::AppPemfcNIThermal<dim>> (new FuelCell::Application::AppPemfcNIThermal<dim>(data));
    }
    
    else if (name_application.compare("meaTwoPhaseSaturationNIT") == 0)
    {
        FcstUtilities::log << "YOU ARE CURRENTLY SOLVING A TWO-PHASE, NON-ISOTHERMAL MEA MODEL" << std::endl;
        return shared_ptr<FuelCell::Application::AppPemfcTPSaturation<dim>> (new FuelCell::Application::AppPemfcTPSaturation<dim>(data));
    }   
    //Run application defined in app_pemfc_twophase_nonisothermal.h
    //run mesh test application
    else if (name_application.compare("test_mesh") == 0)
    {
        FcstUtilities::log << "YOU ARE CURRENTLY RUNNING THE TEST MESH APPLICATION" << std::endl;
        return shared_ptr<FuelCell::Application::AppReadMesh<dim> > (new FuelCell::Application::AppReadMesh<dim>(data));
    }
    else if (name_application.compare("diffusion") == 0)
    {
        if (app_specification.compare("reaction") == 0)
        {
            FcstUtilities::log << "YOU ARE CURRENTLY RUNNING THE DIFFUSION APPLICATION WITH REACTION" << std::endl;
            data->enter_flag("reaction",true);
        }
        else if (app_specification.compare("knudsen") == 0)
        {
            FcstUtilities::log << "YOU ARE CURRENTLY RUNNING THE DIFFUSION APPLICATION WITH LOCAL KNUDSEN EFFECTS" << std::endl;
            data->enter_flag("knudsen",true);
        }
        else if (app_specification.compare("reaction_and_knudsen")==0)
        {
            FcstUtilities::log << "YOU ARE CURRENTLY RUNNING THE DIFFUSION APPLICATION WITH REACTION AND LOCAL KNUDSEN EFFECTS" << std::endl;
            data->enter_flag("reaction",true);
            data->enter_flag("knudsen",true);
        }
        else
            FcstUtilities::log << "YOU ARE CURRENTLY RUNNING DIFFUSION APPLICATION" << std::endl;
        
        return shared_ptr<FuelCell::Application::AppDiffusion<dim> > (new FuelCell::Application::AppDiffusion<dim>(data));
    }
    else if (name_application.compare("ohmic") == 0)
    {
        FcstUtilities::log << "YOU ARE CURRENTLY RUNNING THE OHMIC ELECTRON TRANSPORT APPLICATION" << std::endl;
        return shared_ptr<FuelCell::Application::AppOhmic<dim> > (new FuelCell::Application::AppOhmic<dim>(data));
    }
    else
    {
        FcstUtilities::log << "Application not found. See " << __FILE__ << " and " << __FUNCTION__ << std::endl;
        Assert(false, ExcNotImplemented());
        FcstUtilities::log << "Application not implemented. See " << __FILE__ << std::endl;
    }
}

// ----------------------------------------------------------------------------
template<int dim>
shared_ptr< FuelCell::ApplicationCore::ApplicationWrapper >
SimulationSelector<dim>::select_solver(FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>* app_lin)
{
    if  (data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NONE)
    {
        FcstUtilities::log << "YOUR PROBLEM IS STEADY-STATE AND LINEAR" << std::endl;
        return shared_ptr<FuelCell::ApplicationCore::ApplicationWrapper> (new FuelCell::ApplicationCore::ApplicationWrapper(*app_lin));
    }
    else if (data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::PICARD )
    {
        FcstUtilities::log << "YOU ARE USING PICARD SOLVER FOR STEADY-STATE PROBLEM" << std::endl;
        return shared_ptr<FuelCell::ApplicationCore::Picard> (new FuelCell::ApplicationCore::Picard(*app_lin));
    }
    else if( data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NEWTONBASIC )
    {
        FcstUtilities::log << "YOU ARE USING NewtonBasic NEWTON SOLVER FOR STEADY-STATE PROBLEM" << std::endl;
        return shared_ptr<FuelCell::ApplicationCore::NewtonBasic> ( new FuelCell::ApplicationCore::NewtonBasic(*app_lin) );
    }
    else if( data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NEWTON3PP )
    {
        FcstUtilities::log << "YOU ARE USING Newton3pp NEWTON SOLVER FOR STEADY-STATE PROBLEM" << std::endl;
        return shared_ptr<FuelCell::ApplicationCore::Newton3pp> (new FuelCell::ApplicationCore::Newton3pp(*app_lin));
    }
    else if( data->get_nonlinear_solver() == FuelCell::ApplicationCore::NonLinearSolver::NEWTONLINESEARCH )
    {
        FcstUtilities::log << "YOU ARE USING NewtonLineSearch NEWTON SOLVER FOR STEADY-STATE PROBLEM" << std::endl;
        return shared_ptr<FuelCell::ApplicationCore::NewtonLineSearch> (new FuelCell::ApplicationCore::NewtonLineSearch(*app_lin));
    }
    
    else
    {
        Assert(false, ExcNotImplemented());
        FcstUtilities::log << "This steady-state nonlinear solver is not implemented yet! See " << __FILE__ << std::endl;
    }
}

// ----------------------------------------------------------------------------
template<int dim>
shared_ptr< FuelCell::ApplicationCore::AdaptiveRefinement<dim> >
SimulationSelector<dim>::select_solver_method(FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>* app_lin,
                                              FuelCell::ApplicationCore::ApplicationWrapper*                      newton_solver,
                                              const FuelCell::ApplicationCore::FEVector&                          solution)
{
    if (data->get_refinement_solver() == FuelCell::ApplicationCore::RefinementSolver::ADAPTIVE)
    {
        return shared_ptr<FuelCell::ApplicationCore::AdaptiveRefinement<dim> > (new FuelCell::ApplicationCore::AdaptiveRefinement<dim> (*app_lin, *newton_solver, solution));
    }
    else
    {
        Assert(false, ExcNotImplemented());
        FcstUtilities::log << "This refinement method is not implemented yet! See " << __FILE__ << std::endl;
    }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
template class SimulationSelector<deal_II_dimension>;

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
