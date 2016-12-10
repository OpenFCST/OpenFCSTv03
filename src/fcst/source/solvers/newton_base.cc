//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: newton_base.cc
//    - Description: Base class for all variants of the Newton-Raphson solver
//    - Developers: Jason Boisvert, M. Secanell, V. Zingan
//
//---------------------------------------------------------------------------

#include <solvers/newton_base.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/lac/block_vector.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace FuelCell::ApplicationCore;

const Event newtonBase::bad_derivative = Event::assign("Newton");

//---------------------------------------------------------------------------

newtonBase::newtonBase(ApplicationBase& app)
: ApplicationWrapper(app),
assemble_now(false),
assemble_threshold(0.),
debug_solution(false),
debug_update(false),
debug_residual(false),
debug(0)
{
    //FcstUtilities::log << "->Newton";
    //FcstUtilities::log <<"In Newton Base!\n";

}

//---------------------------------------------------------------------------
void
newtonBase::declare_parameters(ParameterHandler& param)
{
    param.enter_subsection("Newton");
    {
        // Specialization of SolverControl which returns success if either the specified tolerance is achieved or 
        // if the initial residual (or whatever criterion was chosen by the solver class) is reduced by a given factor. 
        // This is useful in cases where you don't want to solve exactly, but rather want to gain two digits or if the 
        // maximal number of iterations is achieved. For example: The maximal number of iterations is 20, the reduction factor is
        // 1% und the tolerance is 0.1%. The initial residual is 2.5. The process will break if 20 iteration are comleted or 
        // the new residual is less then 2.5*1% or if it is less then 0.1%.
        ReductionControl::declare_parameters (param);
        
        //
        param.declare_entry("Assemble threshold", 
                            "0.", 
                            Patterns::Double(),
                            "");
        param.declare_entry("Debug level", 
                            "0", 
                            Patterns::Integer(),
                            "Write debug output to #FcstUtilities::log; the higher the number, the more output.");
        param.declare_entry("Debug solution", 
                            "false", 
                            Patterns::Bool(),
                            "Output the solution at every Newton iteration.");
        param.declare_entry("Debug update", 
                            "false", 
                            Patterns::Bool(),
                            "Output the solution update at every Newton iteration.");
        param.declare_entry("Debug residual", 
                            "false", 
                            Patterns::Bool(),
                            "Output the residual at every Newton iteration.");
        param.declare_entry("Total number of special blocks",
                            "0",
                            Patterns::Integer(),
                            "For internal use only.");
        
        for(unsigned int index = 1; index <= 10; ++index)
        {
            std::ostringstream streamOut;
            streamOut << index;
            const std::string name = "special_block_" + streamOut.str();
            
            param.declare_entry(name.c_str(),
                                "0",
                                Patterns::Integer(),
                                " ");
        }
    }
    param.leave_subsection();
    
    ApplicationWrapper::declare_parameters(param);
}

//---------------------------------------------------------------------------
void
newtonBase::_initialize(ParameterHandler& param)
{
    param.enter_subsection("Newton");
    {
           control.parse_parameters (param);
           
           // Control uses deallog which is not parallel to output information. 
           // We need to deactivate printing to screen so that we do not get repeated messages:      
           control.log_history(false);
           control.log_result(false);
           
           assemble_threshold = param.get_double("Assemble threshold");
           debug = param.get_integer("Debug level");
           debug_solution = param.get_bool("Debug solution");
           debug_update   = param.get_bool("Debug update");
           debug_residual = param.get_bool("Debug residual");
           n_blocks       = param.get_integer("Total number of special blocks");

           for(unsigned int index = 1; index <= n_blocks; ++index)
           {
                  std::ostringstream streamOut;
                  streamOut << index;
                  const std::string name = "special_block_" + streamOut.str();

                  blocks.push_back( param.get_integer(name.c_str()) );
           }
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
void
newtonBase::initialize(ParameterHandler& param)
{
    ApplicationWrapper::initialize(param);
    _initialize(param);
}

//---------------------------------------------------------------------------
double
newtonBase::threshold (double value)
{
    const double t = assemble_threshold;
    assemble_threshold = value;
    return t;
}

//---------------------------------------------------------------------------
void
newtonBase::assemble()
{
    assemble_now = true;
}

//---------------------------------------------------------------------------
double
newtonBase::residual (FuelCell::ApplicationCore::FEVector& dst, const FuelCell::ApplicationCore::FEVectors& src)
{
    //FcstUtilities::log <<"Compute Residual!!\n";
    const double result = app->residual(dst, src);
    return result;
}

//---------------------------------------------------------------------------
void
newtonBase::debug_output(const FEVector& sol,
                         const FEVector& update,
                         const FEVector& residual) const
{
    if (debug_solution)
    {
        std::ostringstream streamOut;
        streamOut << step;
        std::string name = "Newton_" + streamOut.str();
        FEVectors out;
        out.add_vector(sol, "Solution");
        app->data_out(name, out);
    }
    else if (debug_update)
    {
        std::ostringstream streamOut;
        streamOut << step;
        std::string name = "Newton_" + streamOut.str();
        FEVectors out;
        out.add_vector(update, "update");
        app->data_out(name, out);
    }
    else if (debug_residual)
    {
        std::ostringstream streamOut;
        streamOut << step;
        std::string name = "Newton_" + streamOut.str();
        FEVectors out;
        out.add_vector(residual, "residual");
        app->data_out(name, out);
    }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------