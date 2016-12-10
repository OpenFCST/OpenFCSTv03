//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: newton_w_line_search.h
//    - Description: Variant of the Newton-Raphson solver
//    - Developers: Marc Secanell
//
//---------------------------------------------------------------------------

#ifndef __deal2__appframe__newton_w_lines_search_h
#define __deal2__appframe__newton_w_lines_search_h

#include <solvers/newton_base.h>

namespace FuelCell
{
namespace ApplicationCore
{
    /**
     * Application class performing a Newton's iteration as described in newtonBase.
     *
     * The implemented algorithm attempts to find a step size length, \f$ h \f$,
     * such that the scalar function
     * \f[
     * f (h) = \| \mathbf{F}(\mathbf{x}_{i−1} + h \delta \mathbf{x} ) \|_2^2.
     * \f]
     * is reduced at every step. The meaning of \f$ \delta \mathbf{x} \f$ and how it is obtained
     * is described in the documentation for the base class, i.e. newtonBase.
     *
     * The method used to determine \f$ h \f$ is referred to as the global-convergence
     * strategy. The strategy can greatly affect both the efficiency and versatility of
     * Newton’s method.
     *
     * This class contains two methods to determine \f$ h \f$ . Two options are possible.
     *
     * <h3> Overrelaxation and step size control</h3>
     *
     * If, in the Newton subsection in the parameter file, the option "Line search" is set to false, i.e. set Line search = false,
     * then a normal Netwon iteration is performed where \f$ h = 1 \f$. However, since Newton methods are notorious
     * for overshooting the solution at the initial stages, the algorithm allows the user to specify an
     * overrelaxation constant, \f$ h \f$, for a given number of steps, i.e. \p overrelax_steps, such that the step size during the
     * initial stages of the algorithm, i.e. for iteration number < \p overrelax_steps, \f$ h \f$ is given as follows:
     * \f[
     * h = overrelaxation*step
     * \f]
     * where \p step is the iteration number in the Newton solver. The \p overrelaxation constant is specified in the
     * Newton subsection in the input file with set Initial Overrelaxation = 0.2. The parameter \p overrelax_steps is set in
     * the same section with set Number of iterations with overrelaxation = 1.
     *
     * In addition to the overrelaxation option, the algorithm makes sure that after every iteration the residual is reduced.
     * If the residual is not reduced, then the step size, \f$ h = 1 \f$, is discarded and updated with
     * \f[
     * h = 2^{-i}
     * \f]
     * where \f$ i \f$ is the number of failed attempts.
     *
     * <h3>Line search</h3>
     *
     * If, in the Newton subsection in the parameter file, the option "Line search" is set to true, i.e. set Line search = true,
     * then a line search is performed at every iteration in order to estimate the value of \f$ h \f$.
     *
     * In this case, the line search simply involves evaluating the residual at several points, by default 5 points, with
     * \f$ h_i = 2^{-i} \f$ and \f$ i = 1,...,5 \f$.
     *
     * At each iteration a line search is performed by evaluating the L2 norm of the residual for each
     * \f$ h_i \f$ and then the best step size is selected.
     *
     *
     * <h3> Parameters </h3>
     *
     * The parameters that control all Newton solvers are defined in the subsection Newton in the data input file.
     * The parameters are the following:
     * \code
     * subsection Newton
     *  //NewtonLineSearch specific parameters
     *  set Line search = false                # Specify if the code should perform a line search
     *  set Initial Overrelaxation = 1         # Specify the over-relaxation length, \alpha, for the initial step, i.e. u_i = u_{i-1} + \alpha \delta u
     *  set Number of iterations with overrelaxation = 1 # Specify for how many steps should the code overrelax the update ()
     *  // All the parameters below are from NetwonBase:
     *  set Max steps          = 100           # Maximum number of iterations
     *  set Tolerance          = 1.e-8         # Absolute tolerance
     *  set Reduction          = 1.e-20        # Maximum allowed reduction during line search. If the residual is not reduced by at least this number the iteration is discarded.
     *  set Assemble threshold = 0.0           # Used as a threshold to state if the Jacobian needs to be re-evaluated
     *  set Debug level        = 0
     *  set Debug residual     = false         # Would you like the code to output the residual at every Newton iteration?
     *  set Debug solution     = false         # Would you like the code to output the solution at every Newton iteration?
     *  set Debug update       = false
     * end
     * \endcode
     *
     * @author Marc Secanell, 2007
     */
    class NewtonLineSearch : public newtonBase
    {
    public:
        /**
         * The Event set by NewtonLineSearch if
         * convergence is becoming bad
         * and a new matrix should be
         * assembled.
         */
        static const FuelCell::ApplicationCore::Event bad_derivative;
        /**
         * Constructor, receiving the application computing the residual and
         * solving the linear problem.
         */
        NewtonLineSearch (ApplicationBase& app);

        /** Declare the input parameters. */
        virtual void declare_parameters (ParameterHandler& param);

        /** Read the parameters local to NewtonLineSearch. */
        void _initialize (ParameterHandler& param);

        /** Read the parameters */
        virtual void initialize (ParameterHandler& param);

        /**
         * The actual Newton solver.
         */
        virtual void solve(FuelCell::ApplicationCore::FEVector& u,
                           const FuelCell::ApplicationCore::FEVectors& in_vectors);


    private:

        /**
         * Value used to multiply \f$\delta u\f$, the solution update, from the Newton algorithm
         * during the first overrelax_steps iterations. We will multiply the first
         * iteration by \f$\alpha \delta u\f$, where \f$ \alpha \f$ is the value of overrelax and then
         * the following iterations by \f$it \alpha \delta u\f$ where \f$ it \f$ is the iteration number
         */
        double overrelax;
        /**
         * Number of iterations to which we apply overrelax
         */
        unsigned int overrelax_steps;
        /**
         * Do you want to use a line search at each step?
         */
        bool line_search;
        
        /**
         * Block not allowed to be negative:
         */
        unsigned int block_to_fix;
        
        /**
         * Routine used to find if there are any negative values on a given solution block. The block searched is given
         * in 
         * - Solution variable not allowed to be negative
         * 
         * If you do not want to search any block, provide a negative value.
         * 
         */
        bool find_negative_values(const FEVector& u);
    };
}
}

#endif