//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: picard_base.cc
//    - Description: Base class for all variants of the Picard solver
//    - Developers: Mayank Sabharwal
//
//---------------------------------------------------------------------------

#include <solvers/picard_base.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/lac/block_vector.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace FuelCell::ApplicationCore;

//---------------------------------------------------------------------------

PicardBase::PicardBase(ApplicationBase& app)
: ApplicationWrapper(app),
assemble_now(false),
abs_tolerance(0.),
rel_tolerance(0.),
maxsteps(100),
debug_solution(false),
debug_residual(false),
debug(0)
{

}

//---------------------------------------------------------------------------
void
PicardBase::declare_parameters(ParameterHandler& param)
{
    param.enter_subsection("Picard");
    {
        param.declare_entry("Max steps",
			    "100",
			    Patterns::Integer(),
			    "Maximum number of steps for the Picard iterations");
	param.declare_entry("Absolute tolerance",
			    "1e-20",
			    Patterns::Double(),
			    "Absolute Tolerance for convergence of solution");
        param.declare_entry("Relative tolerance",
                            "1e-5",
                            Patterns::Double(),
                            "Relative Tolerance for convergence of solution");
	param.declare_entry("Debug level", 
                            "0", 
                            Patterns::Integer(),
                            "Write debug output to #FcstUtilities::log; the higher the number, the more output.");
        param.declare_entry("Debug solution", 
                            "false", 
                            Patterns::Bool(),
                            "Output the solution at every Picard iteration.");
        param.declare_entry("Debug residual", 
                            "false", 
                            Patterns::Bool(),
                            "Output the residual at every Picard iteration.");
    }
    param.leave_subsection();
    
    ApplicationWrapper::declare_parameters(param);
}

//---------------------------------------------------------------------------
void
PicardBase::initialize(ParameterHandler& param)
{
    ApplicationWrapper::initialize(param);
    param.enter_subsection("Picard");
    {
           maxsteps = param.get_integer("Max steps");
           abs_tolerance = param.get_double("Absolute tolerance");
           rel_tolerance = param.get_double("Relative tolerance");
           debug = param.get_integer("Debug level");
           debug_solution = param.get_bool("Debug solution");
           debug_residual = param.get_bool("Debug residual");
    }
    param.leave_subsection();
}



//---------------------------------------------------------------------------
void
PicardBase::assemble()
{
    assemble_now = true;
}

//---------------------------------------------------------------------------
double
PicardBase::residual (FuelCell::ApplicationCore::FEVector& dst, const FuelCell::ApplicationCore::FEVectors& src)
{
    //FcstUtilities::log <<"Compute Residual!!\n";
    const double result = app->residual(dst, src);
    return result;
}

//---------------------------------------------------------------------------
void
PicardBase::debug_output(const FEVector& sol,
                         const FEVector& update,
                         const FEVector& residual) const
{
    if (debug_solution)
    {
        std::ostringstream streamOut;
        streamOut << step;
        std::string name = "Picard_" + streamOut.str();
        FEVectors out;
        out.add_vector(sol, "Solution");
        app->data_out(name, out);
    }
    else if (debug_residual)
    {
        std::ostringstream streamOut;
        streamOut << step;
        std::string name = "Picard_" + streamOut.str();
        FEVectors out;
        out.add_vector(residual, "residual");
        app->data_out(name, out);
    }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------