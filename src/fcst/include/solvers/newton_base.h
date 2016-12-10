//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: newton_base.h
//    - Description: Base class for all variants of the Newton-Raphson solver
//    - Developers: J. Boisvert, M. Secanell, V. Zingan
//
//---------------------------------------------------------------------------

#ifndef __deal2__appframe__newton_base_h
#define __deal2__appframe__newton_base_h

#include <deal.II/lac/solver_control.h>

#include <application_core/application_wrapper.h>

using namespace dealii;
using namespace FuelCell::ApplicationCore;

namespace FuelCell
{
namespace ApplicationCore
{
    /**
       Base class for all classes performing Newton's iteration. Children of this class mainly differ
       due to the global search method implemented as discussed below.
      
       \b Theory
      
       Newton's method approximates a value for \f$ \mathbf{x}^* \in R^n \f$ such that
       \f[
       \mathbf{F}(\mathbf{x}^*) = 0
       \f]
       where \f$ \mathbf{F} \f$ is a system of \f$ n \f$ nonlinear equations.
      
       Solution \f$ \mathbf{x}^* \f$  is iteratively determined by evaluating
       \f[
       \mathbf{x}_i = \mathbf{x}_{i-1} + \delta \mathbf{x}
       \f]
       where \f$ \mathbf{x}_i \f$ is ith Newton iterate starting from an initial guess \f$ \mathbf{x}_0 \f$ and
       \f[
       \delta \mathbf{x} = - \left( \frac{\partial \mathbf{F}(\mathbf{x}_{i-1})}{\partial \mathbf{x}} \right)^{-1} \mathbf{F}(\mathbf{x}_{i-1})
       \f]
      
       Matrix
       \f[
       \frac{\partial \mathbf{F}(\mathbf{x}_{i-1})}{\partial \mathbf{x}}
       \f]
       is Jacobian matrix evaluated at \f$ \mathbf{x}_{i-1} \f$. Most steady-state openFCST linear applications compute the Jacobian and solve the
       linear system above in order to provide the Newton solver with iteration update, \f$ \delta \mathbf{x} \f$, based on the system of PDEs
       that user wants to solve.
      
       Newton iterations are performed until the error between the iterate and the final solution is sufficiently small:
       \f[
       \| \mathbf{e}_i \| = \| \mathbf{x}^* - \mathbf{x}_i \|
       \f]
      
       The error \f$ \mathbf{e}_i \f$ is not known in practice. However, the norm of the
       residual \f$ \| \mathbf{F(\mathbf{x}_i)} \| \f$ proves to be an effective indicator of the rate of decay of the
       error. In the case of openFCST, the L2 norm of the residual is used as a termination criteria.
      
       openFCST checks if
       \f[
       \| \mathbf{F(\mathbf{x}_i)} \| \leq \tau_a
       \f]
       where \f$ \tau_a \f$ is user-supplied absolute tolerance. If the absolute tolerance is satisfied, the Newton
       iterations are terminated. This user-supplied parameter is an input parameter specified in "set Tolerance = 1e-5" as
       discussed below.
      
       Further, user can specify a maximum number of iteration in the Newton subsection of the input file using "set Max steps = 10".
       In this case, the solver would terminate when the max. number of iterations is reached regardless of the value
       of the residual, i.e. the nonlinear problem is not solved.
      
       A third option for stopping convergence, not used in openFCST, is to terminate when
       \f[
       \| \mathbf{F(\mathbf{x}_i)} \| \leq \tau_r \| \mathbf{F(\mathbf{x}_0)} \|
       \f]
       where  \f$ \tau_r \f$ is user-supplied relative tolerance.
      
       Newton’s method should converge to the solution as the number of Newton
       iterations becomes large. If the initial guess \f$ \mathbf{x}_0 \f$ is poor, however, convergence
       may fail to occur. Often the failure of convergence is a consequence of the
       Newton step, \f$ \delta \mathbf{x} \f$, overshooting the correct value of \f$ \mathbf{x}^* \f$.
      
       The problem of overshooting can usually be solved by using only a fraction
       of the Newton direction for each Newton iteration; see section 1.6 of [2] for
       additional details. In other words, a new Newton step, \f$ \delta \mathbf{x}^* \f$, is applied to each Newton
       iteration such that
       \f[
       \delta \mathbf{x}^* = h \delta \mathbf{x}
       \f]
       where \f$ h \f$ is known as the step size.
      
       The method used to determine \f$ h \f$ is referred to as the global-convergence
       strategy. The strategy can greatly affect both the efficiency and versatility of
       Newton’s method.
      
       Several Newton solvers have been implemented that use a different global-convergence
       strategy. In particular, user has three options:
       - NewtonBasic: No global-convergence strategy is implemented;
       - NewtonLineSearch: A line search is implemented as a global strategy;
       - Newton3pp: A three-point parabolic method is implemented as a global-strategy; two options are possible, with or without prediction
       of step size.
      
       <h3>Usage Details:</h3>
      
       @code
       FuelCell::Application::AppCathode<dim> app_linear;     // Create a linear application (computes the step size by solving linear problem)
       FuelCell::ApplicationCore::NewtonBasic(*app_linear);                    // Associate linear problem to nonlinear solver
       app_linear->init_solution(solution);                   // Initialize solution
       FuelCell::ApplicationCore::FEVector solution;                           // Create a vector where solution will be stored
       FuelCell::ApplicationCore::FEVectors vectors;                           // Pass any other constant data to Newton solver
       vectors.add_vector(solution,"Solution");
       FcstUtilities::log << "Solving..." << std::endl;
       app->solve(solution, vectors);                         // Solve
       app->data_out(sol_filename,vectors);                   // Output solution
       @endcode
      
       For a more detailed example, please see AdaptiveRefinement.
      
      
       <h3> Parameters </h3>
      
       The parameters that control all Newton solvers are defined in the subsection Newton in the data input file.
       The parameters are the following:
       \code
       subsection Newton
        set Max steps          = 100           # Maximum number of iterations
        set Tolerance          = 1.e-8         # Absolute tolerance
        set Reduction          = 1.e-20        # Maximum allowed reduction during line search. If the residual is not reduced by at least this number the iteration is discarded.
        set Assemble threshold = 0.0           # Used as a threshold to state if the Jacobian needs to be re-evaluated
        set Debug level        = 0
        set Debug residual     = false         # Would you like the code to output the residual at every Newton iteration?
        set Debug solution     = false         # Would you like the code to output the solution at every Newton iteration?
        set Debug update       = false
       end
       \endcode
      
       <h3> References </h3>
      
       [1] C. T. Kelley. Iterative methods for linear and nonlinear equations, volume 16 of Frontiers in Applied Mathematics. Society for Industrial and
       Applied Mathematics (SIAM), Philadelphia, PA, 1995. With separately available software.
      
       [2] C. T. Kelley. Solving nonlinear equations with Newton’s method. Fundamentals of Algorithms. Society for Industrial and Applied Mathematics
       (SIAM), Philadelphia, PA, 2003.
      
       @author Marc Secanell and Jason Boisvert, 2009-13
     */

    class newtonBase:public ApplicationWrapper
    {
    public:
        /**
           The Event set if
           convergence is becoming bad
           and a new matrix should be
           assembled.
         */
        static const FuelCell::ApplicationCore::Event bad_derivative;
        /**
           Constructor, receiving the application computing the residual and
           solving the linear problem.
         */
        newtonBase(ApplicationBase& app);

        /**
           Declare the input parameters. The following are the parameters used and their default values:
           In section: "Newton"
           - "Tolerance",
           - "Max steps",
           - "Assemble threshold", Default: "0."
           - "Debug level", Default: "0",
           - "Debug solution", Default: "false", Patterns::Bool()
           - "Debug update", Default: "false", Patterns::Bool()
           - "Debug residual", "false", Patterns::Bool()
         */
        virtual void declare_parameters (ParameterHandler& param);

        /**
           Read the parameters.
         */
        virtual void initialize (ParameterHandler& param);

        /**
           Read the parameters local to newtonBase.
         */
        void _initialize (ParameterHandler& param);

        /**
           Set the maximal residual reduction allowed without triggering
           assembling in the next step. Return the previous value.
          */
        double threshold(double new_value);

        /**
           Control object for the Newton iteration.
         */
        void initialize_initial_guess(BlockVector<double>& dst){};

        /**
           Instead of assembling, this function only sets a flag, such that
           the inner application will be required to assemble a new derivative
           matrix next time solve() is called.
         */
        virtual void assemble();

        /**
           Returns the L2-norm of the residual and the residual vector in "dst" using
           the residual function in the ApplicationBase used to initialize the application.
           The FEVectors rhs should contain at least the following vectors:
           -
          
           \note Finish this section
         */
        virtual double residual(FuelCell::ApplicationCore::FEVector& dst,
                                const FuelCell::ApplicationCore::FEVectors& rhs);

        /**
           The actual Newton solver.
          
           In this section, the Newton solver performs the iterations necessary to solve the nonlinear system.
           During the process it will call assemble() and residual() as needed. Once it is done, it will return
           the solution as parameter \b{u}.
          
         */
        virtual void solve(FuelCell::ApplicationCore::FEVector& u,
                           const FuelCell::ApplicationCore::FEVectors& in_vectors) = 0;

        /**
           Control object for the Newton iteration.
         */
        ReductionControl control;

        /**
           Number of Iterations;
         */
        double numIter;

    protected:
        /**
           Function used to output any necessary information at each Newton iteration
           for debugging purposes
         */
        void debug_output(const FEVector& sol,
                          const FEVector& update,
                          const FEVector& residual) const;

        /**
           This flag is set by the function assemble(),
           indicating that the matrix must be assembled anew upon
           start.
         */
        bool assemble_now;

        /**
           Threshold for re-assembling matrix.
          
           If the quotient of two consecutive residuals is
           smaller than this threshold, the system matrix is not
           assembled in this step.
          
           @note This parameter should be adjusted to the residual gain of the inner solver.
         */
        double assemble_threshold;
        /**
           Print updated solution after each step into file
           <tt>Newton_uNNN</tt>?
         */
        bool debug_solution;
        /**
           Print Newton update after each step into file
           <tt>Newton_dNNN</tt>?
         */
        bool debug_update;
        /**
           Print Newton residual after each step into file
           <tt>Newton_rNNN</tt>?
         */
        bool debug_residual;
        /**
           Write debug output to #FcstUtilities::log; the higher the
           number, the more output.
         */
        unsigned int debug;
        /**
           The number of a basic Newton iteration.
         */
        unsigned int step;

        /**
           This vector specifies the blocks
           of the global solution
           which are supposed to be treated
           specially.
          
           Such blocks usually include
           the density or molar fraction of oxygen (burns as a result of chemical reactions)
           and the density or molar fraction of hydrogen (burns as well as a result of chemical reactions).
          
           These blocks are defined in the parameter
           files manually. Therefore, you need to know them beforehand.
          
           Knowing these blocks, the non-linear solver checks against the negative values
           in those blocks and reduces the solution update if needed to keep the values
           positive.
         */
        std::vector<unsigned int> blocks;

        /**
           The total number of \p blocks.
         */
        unsigned int n_blocks;
    };
}
}

#endif