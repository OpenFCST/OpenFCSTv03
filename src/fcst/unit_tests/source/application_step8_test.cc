//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2014 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: application_step3_test.cc
//    - Description: Test for application. This application solves Step-3 in the deal.II tutorial
//    - Developers: Marc Secanell,    University of Alberta 
//    - $Id: application_step8_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include <application_step8_test.h>

namespace NAME = FuelCell::UnitTest;

//---------------------------------------------
void 
NAME::ApplicationStep8Test::runApplication()
{
    ParameterHandler param;
    FuelCell::Application::AppStep8<deal_II_dimension> app;
        
    app.declare_parameters(param);
    param.enter_subsection("Discretization");{
        param.set("Element","FESystem[FE_Q(1)^2]");
    }
    param.leave_subsection();

    param.enter_subsection("Adaptive refinement");{
        param.set("Refinement","adaptive");  //#{global | steps | adaptive}
    }
    param.leave_subsection();

    param.enter_subsection("Output");{
        param.enter_subsection("Data");{
            param.set("Output format","vtk");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    app.initialize(param);
    app.get_data()->set_nonlinear_solver("None");
    app.get_data()->set_linear_solver("CG");
    FEVector u;
    FEVector rhs;
    
    
    
    for (unsigned int cycle=0; cycle<2; ++cycle)
    {

        FcstUtilities::log<< "Cycle "<<cycle<<" : "<<std::endl;

        std::string filename = "step-8-grid-"+std::to_string(cycle);
        std::string sol_filename = "sol-step-8-grid-"+std::to_string(cycle);
    
        if (cycle == 0)
        {
            app.remesh_dofs();
            app.remesh_matrices();
        }
        else
        {

            app.estimate(u);
            app.remesh();
        }
        
        // initialize vectors with new grid
        app.init_vector(u);
        app.init_vector(rhs);
        FEVectors data;

        data.add_vector(rhs, "residual");
        data.add_vector(rhs, "Solution");

        // build rhs
        app.residual(rhs, data);
        
        // notify to assemble before solving
        app.notify(Event::assign("LinearAssembly"));

        app.solve(u, data);
        app.grid_out(filename);
        FcstUtilities::log << "u:" << u.n_blocks() << ':';
        for (unsigned int b=0;b<u.n_blocks();++b)
            FcstUtilities::log << ' ' << u.block(b).size();
        FcstUtilities::log << std::endl;

        //app.data_out(sol_filename,data);

    }
    
    // ========================================
    // USE FOR DEBUGGING ONLY:
    // ========================================
    /*
    AppFrame::FEVectors vectors;
    vectors.add_vector(u,"Solution");
    app.data_out("sol-step-8",vectors);
    */
    /*
    for (unsigned int i=0; i<u.block(0).size(); ++i)
    {
        FcstUtilities::log<<"Displacement x-dir, node "<<i<<" : "<<u.block(0)(i)<<std::endl;       
        FcstUtilities::log<<"Displacement y-dir, node "<<i<<" : "<<u.block(1)(i)<<std::endl;        
    }
    */
    // ========================================
    
    // Test the solution at a given node (Assuming 3 global refinements:)
    // First test BC
    TEST_ASSERT_DELTA_MSG(0.0, u(0), 1e-06, "ApplicationStep8Test::runApplication failed");
    TEST_ASSERT_DELTA_MSG(0.0, u(0), 1e-06, "ApplicationStep8Test::runApplication failed");
    // Then test somewhere else (here assumed cycle = 2)
    double expectedAnswer0(0.00522218);
    double expectedAnswer1(0.00351314);
    TEST_ASSERT_DELTA_MSG(expectedAnswer0, u.block(0)(20), 1e-06, "ApplicationStep8Test::runApplication failed");
    TEST_ASSERT_DELTA_MSG(expectedAnswer1, u.block(1)(20), 1e-06, "ApplicationStep8Test::runApplication failed");
    

}
