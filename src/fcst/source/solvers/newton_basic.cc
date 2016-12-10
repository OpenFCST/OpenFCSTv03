// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: newton_basic.cc
// - Description: This class performs basic Newton iterations with a constant weight
// - Developers: Valentin N. Zingan, University of Alberta
//
// ----------------------------------------------------------------------------

#include <solvers/newton_basic.h>

const Event NewtonBasic::bad_derivative = Event::assign("Newton");

       //////////////////////////////////////////////////
       //////////////////////////////////////////////////
       // CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
       //////////////////////////////////////////////////
       //////////////////////////////////////////////////

// ---             ---
// --- Constructor ---
// ---             ---

NewtonBasic::NewtonBasic(ApplicationBase& app)
:
newtonBase(app)
{
    FcstUtilities::log << "->NewtonBasic";
}

// ---            ---
// --- Destructor ---
// ---            ---

NewtonBasic::~NewtonBasic()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

void
NewtonBasic::declare_parameters(ParameterHandler& param)
{
    newtonBase::declare_parameters(param);

    param.enter_subsection("Newton");
    {
        param.declare_entry("Basic Newton constant weight",
                            "1.0",
                            Patterns::Double(),
                            " ");
    }
    param.leave_subsection();
}

// ---            ---
// --- initialize ---
// ---            ---

void
NewtonBasic::initialize(ParameterHandler& param)
{
    newtonBase::initialize(param);

    param.enter_subsection("Newton");
    {
        weight = param.get_double("Basic Newton constant weight");
    }
    param.leave_subsection();
}

       ////////////////////
       ////////////////////
       // SOLVE FUNCTION //
       ////////////////////
       ////////////////////

// ---       ---
// --- solve ---
// ---       ---

void
NewtonBasic::solve(FEVector&        u,
                   const FEVectors& in_vectors)
{
    // initialize "step"
    this->step = 0;

    // output the L2 norm of the initial guess if needed
    if( debug > 2 )
        FcstUtilities::log << "L2 norm of the initial guess: " << u.l2_norm() << std::endl;

    // "Du" is the solution of a linear system
    FEVector Du;
    Du.reinit(u);

    // "res" is the residual of a linear system
    FEVector res;
    res.reinit(u);

    // "src1" contains "Solution" and "Newton iterate"
    FEVectors src1;
    src1.add_vector(u, "Newton iterate");
    src1.merge(in_vectors);

    // "src2" contains "Solution", "Newton iterate", and "Newton residual"
    FEVectors src2;
    src2.add_vector(res, "Newton residual");
    src2.merge(src1);

    // "ApplicationBase::data" contains "Newton"
    this->get_data()->enter("Newton", u);

    // fill "res" using the initial guess info
    double residual     = app->residual(res, src1);
    res *= -1.0;
    double old_residual = residual;

    FcstUtilities::log << "iter  = " << step     << std::endl;
    FcstUtilities::log << "error = " << residual << std::endl; // L2 norm of the residual

    // the basic Newton's loop with a constant weight
    while(control.check(this->step++, residual) == SolverControl::iterate)
    {
        // checking
        if(residual/old_residual >= assemble_threshold)
        app->notify(bad_derivative);

        // reset "Du"
        Du.reinit(u);

        // solve a linear system
        // we pass u^n, Du = 0, res^n
        try
        {
            app->solve(Du, src2);
        }
        catch(SolverControl::NoConvergence& e)
        {
            FcstUtilities::log << "Inner iteration failed after "
                               << e.last_step
                               << " steps with residual "
                               << e.last_residual
                               << std::endl;
        }

        // update "u"
        u.add(weight, Du);

        bool flag = true;
        double k  = 1.0;
        while(flag)
        {
            bool flag2 = false;
            for(unsigned int index = 0; index < this->blocks.size(); ++index)
            {
                const unsigned int no_block = this->blocks[index];

                for(unsigned int i = 0; i < u.block(no_block).size(); ++i)
                    if( u.block(no_block)(i) < 0.0 )
                    {
                        flag2 = true;
                        break;
                    }

                if(flag2)
                    break;
            }

            flag = flag2;

            if(flag)
            {
                FcstUtilities::log << "Reduction!!!" << std::endl;
                u.add(-1.0*k*weight, Du); // reset
                u.add( 0.5*k*weight, Du); // update
                k /= 2.0;                  // reduction
            }
        }

        // update "old_residual"
        old_residual = residual;

        // reset "res"
        res.reinit(u);

        // fill "res" using "u"
        residual = app->residual(res, src1);
        res *= -1.0;

        FcstUtilities::log << "iter  = " << step     << std::endl;
        FcstUtilities::log << "error = " << residual << std::endl; // L2 norm of the residual

    } // end while

    // deallocation
    get_data()->erase_vector("Newton");
    
    if( control.last_check() != SolverControl::success )
        throw SolverControl::NoConvergence(control.last_step(),
                                        control.last_value());
    numIter = this->step - 1;
    get_data()->enter("niters", numIter);
}