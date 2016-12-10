//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: newton_w_3pp.cc
//    - Description: Base class for all variants of the Newton-Raphson solver
//    - Developers: Jason Boisvert and M. Secanell
//
//---------------------------------------------------------------------------


#include <solvers/newton_w_3pp.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/lac/block_vector.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace FuelCell::ApplicationCore;

const Event Newton3pp::bad_derivative = Event::assign("Newton");

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Newton3pp::Newton3pp(ApplicationBase& app)
: newtonBase(app)

{
    FcstUtilities::log << "->Newton3pp";
}

//---------------------------------------------------------------------------
void
Newton3pp::declare_parameters(ParameterHandler& param)
{
  newtonBase::declare_parameters(param);

  param.enter_subsection("Newton");
  {
         param.declare_entry("Include step size prediction",
                             "true",
                              Patterns::Bool(),
                             " ");
          param.declare_entry("Use predefined parameters",
                              "true",
                              Patterns::Bool(),
                             " ");
         param.declare_entry("Max. reduction factor",
                             "0.9",
                              Patterns::Double(),
                             " ");
         param.declare_entry("Min. reduction factor",
                             "0.5",
                              Patterns::Double(),
                             " ");
         param.declare_entry("Line search convergence, alpha [0-1]",
                             "0.001",
                              Patterns::Double(),
                             " ");
  }
  param.leave_subsection();
}

//---------------------------------------------------------------------------
void
Newton3pp::initialize(ParameterHandler& param)
{
  newtonBase::initialize(param);

  param.enter_subsection("Newton");
  {
         with_prediction = param.get_bool("Include step size prediction");
         use_predefined = param.get_bool("Use predefined parameters");
         max_reduction_factor = param.get_double("Max. reduction factor");
         min_reduction_factor = param.get_double("Min. reduction factor");
         alpha= param.get_double("Line search convergence, alpha [0-1]");
  }
  param.leave_subsection();
}

//---------------------------------------------------------------------------
void
Newton3pp::solve(FuelCell::ApplicationCore::FEVector& u,
                 const FuelCell::ApplicationCore::FEVectors& in_vectors)
{

    // Initialize:
    this->step = 0;

    int line_search_iterations = 0;
    int max_line_search_reductions =  25;

    // lambda == h_j in the theory
    double lambda_trial = 1.0;
    double lambda_previous  = 1.0;        //h_{j-2}
    double lambda_current =lambda_trial;  //h_{j-1}
    double lambda_old = 1.0;
    double lambda_2 = 0;

    //   double residual_trial_norm = 0.0;
    //   double residual_norm = 0.0;
    //   double residual_hold_norm = 0.0;

    double mf_original = 0.0;
    double mf_current = 0.0;
    double mf_previous = 0.0;

    if (use_predefined)
        set_predefined_parameters();

    FcstUtilities::log.push ("Newton");
    if (debug > 2)
        FcstUtilities::log << "u: " << u.l2_norm() << std::endl;

    FEVector Du;
    FEVector res;

    res.reinit(u);
    FEVectors src1;
    FEVectors src2;
    src1.add_vector(u, "Newton iterate");
    src1.merge(in_vectors);
    src2.add_vector(res, "Newton residual");
    src2.merge(src1);

    get_data()->enter("Newton", u);

    // fill res with (f(u), v)
    double residual = app->residual(res, src1);
    double old_residual = residual;

    // Output the solution at the Newton iteration if residual debug is on
    this->debug_output(u, Du, res);

    //Begin Newton Iterations
    while (control.check(this->step++, residual) == SolverControl::iterate)
    {
        if( (this->step > 2) && with_prediction )
        {
            if(lambda_old < lambda_2* (1-this->alpha))
                lambda_current = lambda_old;
            else
                lambda_current = std::min(1.0, 2.0*lambda_old);
        }
        else lambda_current = 1.0;

        //FcstUtilities::log << "Initial lambda: " << lambda_current << "\n";
        //FcstUtilities::log << "residual: " << residual <<"\n";

        // assemble (Df(u), v)
        if(residual/old_residual >= assemble_threshold)
            app->notify (bad_derivative);
        Du.reinit(u);
        //Solver Linear System
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
        //-- Three Point Parabolic Line Search Method
        ////////////////////////////////////////////////////////////

        // get coefficients
        mf_original = sum_of_squares(old_residual);
        //start search

        //Try with full newton step;
        u.add(-lambda_current,Du);
        old_residual =residual;
        residual = app->residual(res, src1);

        //Get coefficients with full Newton step
        mf_current = sum_of_squares(residual);
        mf_previous = mf_current;

        while (std::isnan(residual) || !armijo(lambda_current,old_residual,residual) )
        {
            residual = app->residual(res, src1);
            //FcstUtilities::log << "Starting line search at: " <<  line_search_iterations << " for a max of: "<<max_line_search_reductions << "\n";
            //FcstUtilities::log<<"I am in search loop.\n";

            u.add(lambda_current,Du);
            //compute new step
            lambda_trial = three_point_step (line_search_iterations,lambda_current,lambda_previous,mf_original,mf_current,mf_previous);
	    
            //FcstUtilities::log<< "Computed the step: " << lambda_trial << "\n";

            if (std::isnan(lambda_trial))
                lambda_trial = lambda_current * (max_reduction_factor+min_reduction_factor)*0.5;
            //Now the safe Gurads.
            lambda_trial = std::max(lambda_trial,min_reduction_factor * lambda_current);
            lambda_trial = std::min(lambda_trial,max_reduction_factor * lambda_current);

            FcstUtilities::log << "Step after safe guards: "<< lambda_trial << std::endl;


            //update lambda
            lambda_previous = lambda_current;
            lambda_current = lambda_trial;

            //Try Stepsize
            u.add(-lambda_current,Du);

            //calcuate residual_trial
            residual = app->residual(res, src1);

            //new coefficents
            mf_previous = mf_current;
            mf_current = sum_of_squares(residual);

            //incrment line searc iterations and check for mac converg reach
            line_search_iterations++;

            //FcstUtilities::log<<"IN LOOP: "<< line_search_iterations <<  " \n";
            //FcstUtilities::log<<"Residual: "<<residual <<"\n";
            if(line_search_iterations > max_line_search_reductions)
            {
                FcstUtilities::log << "Reached maximum line search iterations, Exiting and restarting.\n";
                break;
            }
            //FcstUtilities::log << "residual:" << residual << "\n";
            //FcstUtilities::log << "old residual:" << old_residual << "\n";


            //Update lambda_old
        }
        if(this->step > 2)
            lambda_2= lambda_old;
        lambda_old = lambda_current;
        line_search_iterations = 0;


        // Output the global residual and the equation specific residual:
        FcstUtilities::log << "Overall residual "<<"at iteration "<<this->step<<" = " << residual << std::endl;
        for (unsigned int i = 0; i<res.n_blocks(); i++)
            FcstUtilities::log << "Residual for equation "<<i<<" is: "<<res.block(i).l2_norm() << std::endl;

        // Debug output options:
        this->debug_output(u, Du, res);
    }


    get_data()->erase_vector("Newton");
    //get_data()->block_vector_pool.free(u_trial);
    FcstUtilities::log.pop();

    // in case of failure: throw
    // exception
    if (control.last_check() != SolverControl::success)
        throw SolverControl::NoConvergence (control.last_step(),
                                            control.last_value());
    // otherwise exit as normal

    //step=6;
    numIter= this->step - 1;
    get_data()->enter("niters",numIter);
}

//---------------------------------------------------------------------------
double
Newton3pp::sum_of_squares (double residual_norm)
{
    return (0.5 * residual_norm * residual_norm);
}

//---------------------------------------------------------------------------
double
Newton3pp::three_point_step (int line_search_iterations,double lambda_current, double lambda_previous, double mf_original, double mf_current, double mf_previous)
{

    if (line_search_iterations == 0 ) return max_reduction_factor * lambda_current;

    // Compute coefficients of interpolation polynomial
    double c_2 = 0.0;
    c_2 = lambda_previous * (mf_current - mf_original) - lambda_current * (mf_previous - mf_original);
    if (c_2 >= 0) return (max_reduction_factor * lambda_current);
    // Equivalent to:
    //double ddp = c_2/((lambda_current-lambda_previous)*lambda_previous*lambda_current);

    //FcstUtilities::log<<"ddp: "<<ddp<<std::endl;
    //FcstUtilities::log<<"ce2: "<<c_2<<std::endl;

    //if (ddp < 0) return (max_reduction_factor * lambda_current);
    // Since (lambda_current-lambda_previous) < 0
    double c_1 = 0.0;
    c_1 = lambda_current * lambda_current * (mf_previous - mf_original) - lambda_previous * lambda_previous * (mf_current - mf_original);
    return -c_1 * 0.5 / c_2;

}

//---------------------------------------------------------------------------
bool
Newton3pp::armijo(double step_length,double residual_norm,double new_residual_norm)
{

    //set acceptance factor
    double fact = (1.0 - alpha*step_length);

    //FcstUtilities::log <<"alpha: " << alpha <<"acc: "<< fact << "step: " << step_length << "\n";

    return ( new_residual_norm <= ( fact * residual_norm ));

}
