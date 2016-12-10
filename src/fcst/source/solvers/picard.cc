//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: picard.cc
//    - Description: Picard solver
//    - Developers: Mayank Sabharwal
//
//---------------------------------------------------------------------------

#include <solvers/picard.h> 
#include <deal.II/base/data_out_base.h>
#include <deal.II/lac/block_vector.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

using namespace FuelCell::ApplicationCore;


//---------------------------------------------------------------------------
Picard::Picard(ApplicationBase& app)
    : PicardBase(app),
    underrelaxation(false)
{
  FcstUtilities::log << "->Picard";
}

//---------------------------------------------------------------------------
void
Picard::declare_parameters(ParameterHandler& param)
{
    PicardBase::declare_parameters(param);
    param.enter_subsection("Picard");
    {
        param.declare_entry("Under-relaxation",
                            "false",
                            Patterns::Bool(),
                            "Use adaptive under-relaxation");
        param.declare_entry("Alpha", 
                            "4", 
                            Patterns::Double(),
                            "Alpha value for underrelaxation;Range of 2-6 works best");
        param.declare_entry("Gamma min", 
                            "0.6", 
                            Patterns::Double(),
                            "Gamma min value for underrelaxation;Range of 0.1-0.6 works best");
    }
    param.leave_subsection();      
}

//---------------------------------------------------------------------------
void
Picard::initialize (ParameterHandler& param)
{
    PicardBase::initialize(param);
    param.enter_subsection("Picard");
    {
        underrelaxation = param.get_bool("Under-relaxation");
        alpha = param.get_double("Alpha");
        gamma_min = param.get_double("Gamma min");
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
void
Picard::solve (FuelCell::ApplicationCore::FEVector& u, const FuelCell::ApplicationCore::FEVectors& in_vectors)
{
    this->step = 0;

    if (debug>2)
        FcstUtilities::log << "u: " << u.l2_norm() << std::endl;

    FEVector u_n;
    FEVector res;

    res.reinit(u);
    u_n.reinit(u);
    FEVectors src1;
    FEVectors src2;
    src1.add_vector(u, "Solution");
    src1.merge(in_vectors);
    src2.add_vector(res, "residual");
    src2.merge(src1);

    double residual = app->residual(res, src1);
    double old_residual = residual;
    
    double abs_error = 1.e6;
    double rel_error = 1.e6;
    
    this->debug_output(u, u_n, res);
    app->notify(Event::assign("Picard"));
    double old_error = 1e5;
    
    //
    while ((abs_error > this->abs_tolerance) && (rel_error > this->rel_tolerance) && (this->step < this->maxsteps))
    {
        u_n.reinit(u);
        app->solve (u_n, src2);


        abs_error = 0;
        rel_error = 0;
        FEVector error;
        error.reinit(u);
        
        double gamma, delta=0, gamma_min=this->gamma_min, alpha= this->alpha, epsilon=this->abs_tolerance;

        this->compute_errors(u,u_n,error,abs_error,rel_error,delta);
            
        // Adaptive under-relaxation scheme
        if (this->underrelaxation)
        {
            
            FcstUtilities::log<<"Delta: "<<delta<<std::endl;

            double rho=0.95;    //damping factor if solution starts to diverge
            
            // If solution diverges then gamma and alpha are further damped 
            while(old_error<abs_error && rho>0 && this->step>1)
            {
                FcstUtilities::log<<"New solution diverging! Damping gamma and alpha"<<std::endl;
                gamma_min*=rho;
                
                alpha = - log((gamma*rho - gamma_min)/(1-gamma_min))/(delta-epsilon);
                gamma = gamma_min+(1-gamma_min)*exp (-alpha*(delta-epsilon));
                
                for(unsigned int i=0; i<u_n.size();i++)
                    u(i)=u(i)+gamma*(u_n(i)-u(i));
                residual = app->residual(res, src1);
                app->notify(Event::assign("Picard"));
                u_n.reinit(u);
                app->solve(u_n,src2);
                abs_error = 0;
                rel_error = 0;
                
                FcstUtilities::log<<"Gamma: "<<gamma<<std::endl;
                double new_delta = 0;
                this->compute_errors(u,u_n,error,abs_error,rel_error,new_delta);

                FcstUtilities::log<<"Delta: "<<new_delta<<std::endl;
                if (old_error>abs_error)
                    delta = new_delta;
                rho-=0.1;
            }
            if (delta>epsilon)
                gamma = gamma_min+(1-gamma_min)*exp (-alpha*(delta-epsilon));
            else
                gamma = 1;
            FcstUtilities::log<<"Gamma: "<<gamma<<std::endl;
            for(unsigned int i=0; i<u_n.size();i++)
                u(i)=u(i)+gamma*(u_n(i)-u(i));
        }
        //Pure Picard scheme
        else
            u=u_n;
        
        // Changing type of solver in order to assemble the cell residual using the assemble_cell_residual function of equation class
        FEVector temp_res;
        FEVectors temp_src;
        temp_src.merge(src1);
        this->data->set_nonlinear_solver("NewtonLineSearch");
        double new_residual = app->residual(temp_res,temp_src,false);
        
        // Reset the non-linear solver to Picard
        this->data->set_nonlinear_solver("Picard");
        old_residual = residual;
        residual = app->residual(res, src1);

        // Output the global residual and the equation specific residual:
        FcstUtilities::log << "Overall absolute error at iteration "<<this->step<<" = " << abs_error << std::endl;
        FcstUtilities::log << "Overall relative error at iteration "<<this->step<<" = " << rel_error << std::endl;


        for (unsigned int i = 0; i<res.n_blocks(); i++)
            FcstUtilities::log << "Residual for equation "<<i<<" is: "<<temp_res.block(i).l2_norm() << std::endl;

        // Debug output options:
        this->debug_output(u, u_n, res);
        old_error=abs_error;
        this->step+=1;
    }

    AssertThrow((abs_error< this->abs_tolerance)||(rel_error< this->rel_tolerance),SolverControl::NoConvergence (this->step,abs_error));
   
}

void 
Picard::compute_errors ( FEVector &u, FEVector &u_n, FEVector &error, double &abs_error, double &rel_error, double &delta)
{
    int flag =0;
    double dofs=u.size();
    for(unsigned int i=0; i<u_n.size();i++)
    {
        if (u_n(i)<0)
        {
            flag=1;
            u_n(i)=0;
        }
        abs_error+=pow((u(i)-u_n(i)),2);
        error(i)=u(i)-u_n(i);
        if (fabs(u(i)-u_n(i))>delta)
            delta = fabs(u(i)-u_n(i));
    }
    abs_error = sqrt(abs_error)/dofs;
    rel_error = error.l2_norm()/u_n.l2_norm();
    if (flag)
        FcstUtilities::log<<"Negative values in solution were set to zero!!!"<<std::endl;
}