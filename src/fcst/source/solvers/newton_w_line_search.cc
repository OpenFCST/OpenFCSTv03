//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: newton_w_line_search.cc
//    - Description: Variant of the Newton-Raphson solver
//    - Developers: Marc Secanell
//
//---------------------------------------------------------------------------

#include <solvers/newton_w_line_search.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/lac/block_vector.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

using namespace FuelCell::ApplicationCore;

//---------------------------------------------------------------------------
const Event NewtonLineSearch::bad_derivative = Event::assign("Newton");

//---------------------------------------------------------------------------
NewtonLineSearch::NewtonLineSearch(ApplicationBase& app)
    : newtonBase(app)
{
  FcstUtilities::log << "->NewtonLineSearch";
}

//---------------------------------------------------------------------------
void
NewtonLineSearch::declare_parameters(ParameterHandler& param)
{
    param.enter_subsection("Newton");
    {
        param.declare_entry("Line search",
                            "false", 
                            Patterns::Bool(),
                            "A line search is performed at each step.");
        param.declare_entry("Initial Overrelaxation", 
                            "1.", 
                            Patterns::Double(),
                            "This value multiplies the step size in the first iteration, e.g., if the Newton update is 0.25, "
                            "the applied update would be the value in Initial Overrelaxation times 0.25");
        param.declare_entry("Number of iterations with overrelaxation",
                            "1",
                            Patterns::Integer(),
                            "Step during which overrelaxation is applied. At each step the Netwon step is multiplied by"
                            "the step times Initial overrelaxation, e.g., step one is 0.1, then step to is 0.2 and so on until 1 (no overrelaxation).");
        param.declare_entry("Solution variable not allowed to be negative",
                            "0", 
                            Patterns::Integer(),
                            "This enforces a given solution variable to be positive. If it goes negative, step size is reduced.");        
    }
    param.leave_subsection();

    newtonBase::declare_parameters(param);
}

//---------------------------------------------------------------------------
void
NewtonLineSearch::_initialize (ParameterHandler& param)
{
    param.enter_subsection("Newton");
    line_search = param.get_bool("Line search");
    overrelax = param.get_double("Initial Overrelaxation");
    overrelax_steps = param.get_integer("Number of iterations with overrelaxation");
    block_to_fix = param.get_integer("Solution variable not allowed to be negative");    
    param.leave_subsection ();
}

//---------------------------------------------------------------------------
void
NewtonLineSearch::initialize (ParameterHandler& param)
{
    newtonBase::initialize(param);
    _initialize(param);
}

//---------------------------------------------------------------------------
void
NewtonLineSearch::solve (FuelCell::ApplicationCore::FEVector& u, const FuelCell::ApplicationCore::FEVectors& in_vectors)
{
    this->step = 0;

    if (debug>2)
        FcstUtilities::log << "u: " << u.l2_norm() << std::endl;

    FEVector Du;
    FEVector res;

    res.reinit(u);
    Du.reinit(u);
    FEVectors src1;
    FEVectors src2;
    src1.add_vector(u, "Newton iterate");
    src1.merge(in_vectors);
    src2.add_vector(res, "Newton residual");
    src2.merge(src1);

    // fill res with (f(u), v)
    double residual = app->residual(res, src1);
    FcstUtilities::log << "Overall residual at iteration "<<this->step<<" = " << residual << std::endl;
    double old_residual = residual;

    // Output the solution at the Newton iteration if residual debug is on
    this->debug_output(u, Du, res);

    //
    while (control.check(this->step++, residual) == SolverControl::iterate)
    {
        // assemble (Df(u), v)
        if (residual/old_residual >= assemble_threshold)
            app->notify (bad_derivative);

        Du.reinit(u);

        try
        {
            app->solve (Du, src2);
        }
        catch (SolverControl::NoConvergence& e)
        {
            FcstUtilities::log << "Inner iteration failed after "
            << e.last_step << " steps with residual "
            << e.last_residual << std::endl;
        }

        ////////////////////////////////////////////////////////////
        //-- Step size control
        ////////////////////////////////////////////////////////////
        if (line_search)
        {
            FcstUtilities::log << "Line Search" << std::endl;
            // -- Use a line search to increase robustness
            double min_residual = residual;
            double opt_alpha = -0.25;
            unsigned int num_points = 1;
            const double* ref = get_data()->scalar("Refinement");
            if (this->step < 10 &&  *ref < 1)
                num_points = 20;
            else if ((this->step >=10 && this->step < 20) &&  *ref < 1)
                num_points = 10;
            else
                num_points = 5;

            // ------------ GEOMETRIC LINE SEARCH

            // Scale: I will look at 2^0, 2^-1, 2^-2, ..., 2^{-num_points}
            double scale = 2;
            for (unsigned int alpha = 1; alpha < num_points + 1; ++alpha)
            {
                // scale u
                u.add(-pow(scale,-double(alpha-1)), Du);
                residual = app->residual(res, src1);

                if ((residual < min_residual) && !(std::isnan(residual)))
                {
                    opt_alpha = -pow(scale,-double(alpha-1));
                    min_residual = residual;
                }
                //  FcstUtilities::log<<"Residual for "<<pow(scale,-double(alpha-1))<<" is: "<<residual<<std::endl;
                //
                //  scale back u
                u.add(pow(scale,-double(alpha-1)), Du);
            }
            /*
            // ------ LINEAR LINE SEARCH (Not tested) -------------
            for (unsigned int alpha = 1; alpha < num_points; ++alpha)
            {
                u.add(-1.0/num_points, *Du);
                residual = app->residual(*res, u);
                if (residual < min_residual && !std::isnan(residual))
                {
                    opt_alpha = alpha*(1.0/num_points);
                    min_residual = residual;
                }
            }
            // unscale u
            u.add(1.0, *Du);
            */
            //FcstUtilities::log<<"Norm of change in solution: "<<Du->l2_norm()<<std::endl;
            //FcstUtilities::log<<"Line search with "<<num_points<<" number of points. Optimal alpha: "<<opt_alpha<<std::endl;
            FcstUtilities::log << "Step size = " << -opt_alpha << std::endl;
            u.add(opt_alpha, Du);
            residual = app->residual(res, src1);

        }
        else // No line search
        {
            ////////////////////////////////////////////////////////////
            // a) First, use overrelaxation during the initial
            // iterations because Newton's algorithm tends to overshoot:
            double alpha = 1.0;
            const double* ref = get_data()->scalar("Refinement");

            if (this->step < overrelax_steps &&  *ref < 1)
            {
                alpha = overrelax*this->step;
                alpha = (alpha <= 1) ? alpha : 1;
                FcstUtilities::log << "Step size = " << alpha << std::endl;
            }


            // Update solution
            u.add(-alpha,Du);

            // Update residual
            old_residual = residual;
            residual = app->residual(res, src1);

            // b) Next, make sure that the chosen step reduces the residual
            // overwise reduce the residual further:
            unsigned int step_size = 0;
            
            //-- Make sure none of the values are negative:
            bool negative_values = find_negative_values(u);

            //--
            while (residual >= old_residual || std::isnan(residual) || negative_values)
            {                
                ++step_size;

                u.add(alpha*pow(2.0,-double(step_size)), Du);
                residual = app->residual(res, src1);
                
                //-- Make sure none of the values are negative:
                negative_values = find_negative_values(u);
                
                // -- Some output:                    
                if (control.log_history())
                    FcstUtilities::log << "Step size: 2^{-" << step_size << "}" << std::endl;
                if (step_size > 15)
                {
                    if (control.log_history())
                        FcstUtilities::log << "Step size too small!" << std::endl;
                    break;
                }
            }
        }

        // Output the global residual and the equation specific residual:
        for (unsigned int i = 0; i<res.n_blocks(); i++)
            FcstUtilities::log << "Residual for equation "<<i<<" is: "<<res.block(i).l2_norm() << std::endl;
        FcstUtilities::log << "Overall residual at iteration "<<this->step<<" = " << residual << std::endl;

        // Debug output options:
        this->debug_output(u, Du, res);

    }
    
    // in case of failure: throw exception
    if (control.last_check() != SolverControl::success)
        AssertThrow (false, ExcMessage ("No convergence in Newton solver"));

}

//---------------------------------------------------------------------------
bool
NewtonLineSearch::find_negative_values(const FEVector& u)
{
    bool negative_values = false;
    /*
    if (block_to_fix >= -1.0e-5)
    {
        for(unsigned int i = 0; i < u.block(block_to_fix).size(); ++i)
            if( u.block(block_to_fix)(i) < -1.0e-6)
            {
                negative_values = true;
                FcstUtilities::log << "Negative value... u("<<i<<")="<<u.block(block_to_fix)(i)<<"resize..."<<std::endl;
                break;                   
            }
    }
    */
    return negative_values;
    
}