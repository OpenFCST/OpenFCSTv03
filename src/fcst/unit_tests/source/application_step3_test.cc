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
//    - $Id: application_step3_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include <application_step3_test.h>

namespace NAME = FuelCell::UnitTest;

//---------------------------------------------
void 
NAME::ApplicationStep3Test::runApplication()
{
    ParameterHandler param;
    FuelCell::Application::AppStep3<deal_II_dimension> app;
        
    app.declare_parameters(param);
    param.enter_subsection("Output");{
        param.enter_subsection("Data");{
            param.set("Output format","vtk");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    app.initialize(param);
    app.get_data()->set_nonlinear_solver("None");
    #ifdef OPENFCST_WITH_PETSC
    app.get_data()->set_linear_solver("MUMPS");
    #else
    app.get_data()->set_linear_solver("UMFPACK");
    #endif

    
    
    FEVector u;
    FEVector rhs;
    app.init_vector(u);
    app.init_vector(rhs);
    
    FEVectors data;
    data.add_vector(rhs, "Solution");
    data.add_vector(rhs, "residual");

    // build rhs
    app.residual(rhs, data);
    
    
    app.solve(u, data);


    // ========================================
    // USE FOR DEBUGGING ONLY:
    // ========================================
    /*
    AppFrame::FEVectors vectors;
    vectors.add_vector(u,"Solution");
    app.data_out("test.vtk",vectors);
    */
    // ========================================
    
    // Test the solution at a given node (Assuming 3 global refinements:)
    double expectedAnswer(0.298393);
    TEST_ASSERT_DELTA_MSG(expectedAnswer, u(24), 1e-06, "ApplicationStep3Test::runApplication failed");
}
