//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: newton_w_3pp.h
//    - Description: Base class for all variants of the Newton-Raphson solver
//    - Developers: Jason Boisvert and M. Secanell
//
//---------------------------------------------------------------------------

#ifndef __deal2__appframe__newton_w_3pp_h
#define __deal2__appframe__newton_w_3pp_h

#include <solvers/newton_base.h>

namespace FuelCell
{
namespace ApplicationCore
{
    /**
     * Application class performing Newton's iteration. At each iteration a line search is performed using
     * a three point parabolic algorithm to enhance robustness of the algorithm and to improve convergence to a global solution.
     * Either a standard three point parabolic algorithm or a three point parabolic algorithm with prediction can be used depending
     * on the input parameter "set Include step size prediction" in the Newton section of the input file. If this option is
     * set to true, the algorithm with prediction is used. Othewise, the algorithm without prediction is used.
     *
     * The three-point parabolic model, described in section 8.3.1 of [1], attempts
     * to find a step size length, \f$ h \f$, such that it minimizes the scalar function
     * \f[
     * f (h) = \| \mathbf{F}(\mathbf{x}_{i−1} + h \delta \mathbf{x} ) \|_2^2.
     * \f]
     * where the meaning of \f$ \delta \mathbf{x} \f$ and how it is obtain is described in the documentation for the base class, i.e. newtonBase.
     *
     * We estimate the norm above using the following three-point interpolating polynomial
     * \f[
     * p(h) = f (h) = f(0) + \frac{h}{h_{j-1} − h_{j−2}} \left( \frac{(h − h_{j−2}) ( f(h_{j−1}) − f(0) )}{h_{j−1}} + \frac{(h_{j−1} − h)(f(h_{j−2})-f(0)}{h_{j−2}} \right )
     * \f]
     * where 0, \f$ h_{j−1} \f$, \f$ h_{j−2} \f$ are the initial and the two lastly rejected step sizes in the interpolation.
     *
     * Using the equation above, and assuming that \f$ p > 0 \f$, the step size \f$ h_j \f$ is given by
     * \f[
     * h_j = − \frac{p'(0)}{p''(0)}
     * \f]
     * If \f$ p'' ≤ 0 \f$, \f$ h_j \f$ is reduced by a fixed amount.
     *
     * The three-point parabolic method is continued until a sufficient decrease in the residual \f$ \| \mathbf{F}(x_i) \| \f$ occurs based on the condition
     * \f[
     * \| \mathbf{F}(\mathbf{x_{i−1}} + h_j \delta \mathbf{x}) \| < (1 − \alpha h_j ) \| \mathbf{F}(x_{i−1}) \|.
     * \f]
     * where \f$ \alpha ∈ [0, 1] \f$. In practice, \f$ \alpha \f$ is typically \f$ 10^{−4}\f$ [2].
     *
     * When implementing this method, there are a few considerations to take into account.
     *
     * To prevent a reduction in the step size that is either too large or too small
     * to be effective, we allow the user to specify a maximum and minimum reduction rate, i.e. \f$ \sigma_0 , \sigma_1 \f$.
     * Therefore, \f$ h_j \f$ must equal \f$ \sigma h_{j−1} \f$, where \f$ \sigma ∈ [ \sigma_0 , \sigma_1] \f$.
     * This is guaranteed by making the following check after \f$ h_j \f$  is determined
     * \f[
     * h_j = \left\{
     * \begin{array}{ll}
     *    \sigma_0 h_{j−1} & \mbox{ if } h_j < \sigma_0 h_{j−1} \\
     *    \sigma_1 h_{j−1} & \mbox{ if } h_j > \sigma_1 h_{j−1} \\
     *    h_j  & \mbox{otherwise}.
     * \end{array}
     * \right.
     * \f]
     * The three-point interpolating polynomial also requires some initial values for \f$ h_0 \f$ and \f$ h_1 \f$.
     * If a maximum reduction factor is set, then it is convenient
     * to choose
     * \f[
     * h_0 = 1,
     * h_1 = \sigma_1
     * \f]
     *
     * <h3> Three point parabolic with prediction </h3>
     *
     * The value of \f$ h_0 \f$ can be seen as an initial guess for the three-point parabolic
     * method. In the previous section, we simply set \f$ h_0 = 1 \f$. Of course, it is very
     * likely that 1 is a poor choice for the initial guess. This results in additional
     * line-search iterations.
     *
     * In an effort to reduce the number of line-search iterations, we can use the
     * previous step size to more accurately choose an initial step size \f$ h_{0,i} \f$ for the
     * current Newton iteration \f$ i \f$. However, there is also a possibility that the previous
     * step size is a poor choice for an initial guess, therefore we use the condition
     * \f[
     * h_{0,i} = \left\{
     * \begin{array}{ll}
     * h_{i−1} & \mbox{ if } h_{i-1} < h_{i−2}(1-\alpha) \\
     * min(1, 2 h_{i-1}) & \mbox{otherwise}.
     * \end{array}
     * \right.
     * \f]
     *
     * The step sizes \f$ h_{i−1}\f$  and \f$ h_{i−2} \f$ are from previous two Newton iterations. Initially,
     * \f$ h_{0,1} = h_{0,2} = 1 \f$.
     *
     * Using this strategy results in a predictor-corrector approach to the three-point parabolic method.
     * The condition above predicts the value of the step size. After, the three-point parabolic method
     * simply corrects for the differences between \f$ h_{i−1} \f$ and \f$ h_i \f$.
     *
     * @note This is the method by default (set Include step size prediction = true)
     *
     * <h3> Parameters </h3>
     *
     * The parameters that control all Newton solvers are defined in the subsection Newton in the data input file.
     * The parameters are the following:
     * \code
     * subsection Newton
     *  set Include step size prediction = true       # Use prediction with the three-point parabolic method
     *  set Use predefined parameters = false             # Use optimized reduction factor and line search parameters.
     *  set Max. reduction factor = 0.9
     *  set Min. reduction factor = 0.5
     *  set Line search convergence, alpha [0-1] = 0.01
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
     *
     * <h3> References </h3>
     *
     * [1] C. T. Kelley. Iterative methods for linear and nonlinear equations, volume 16 of Frontiers in Applied Mathematics. Society for Industrial and
     * Applied Mathematics (SIAM), Philadelphia, PA, 1995. With separately available software.
     *
     * [2] C. T. Kelley. Solving nonlinear equations with Newton’s method. Fundamentals of Algorithms. Society for Industrial and Applied Mathematics
     * (SIAM), Philadelphia, PA, 2003.
     *
     * @author Jason Boisvert, Raymond Spiteri and Marc Secanell, 2009-13
     */

    class Newton3pp : public newtonBase
    {
    public:

        static const FuelCell::ApplicationCore::Event bad_derivative;

        /**
         * Constructor
         */
        Newton3pp(ApplicationBase& app);

        /**
         * Declare any additional parameters necessary for the solver. In this case the only additional parameter is a bool
         * parameter estating if we want to use a Newton step prediction, i.e.
         *
         * In subsection Newton you have the following option:
         * - set Include step size prediction = true (false)
         *
         * If true then a three-point parabolic method with prediction is used. Otherwise, the standard three-point parabolic is used.
         */
        void declare_parameters(ParameterHandler& param);

        /**
         * Initialize NewtonBase and any additional parameters.
         */
        void initialize(ParameterHandler& param);

        /**
         * The actual Newton solver.
         *
         * In this section, the Newton solver performs the iterations necessary to solve the nonlinear system.
         * During the process it will call assemble() and residual() as needed. Once it is done, it will return
         * the solution as parameter \b{u}.
         */
        virtual void solve(FuelCell::ApplicationCore::FEVector& u,
                           const FuelCell::ApplicationCore::FEVectors& in_vectors);

    private:
        /**
         * Define default parameters for convergence. These depend on refinement level.
         */
        inline void set_predefined_parameters()
        {
            // Depends in which refinement
            const double* ref = get_data()->scalar("Refinement");
            //FcstUtilities::log << "Ref: "<< *ref <<"\n";
            if (*ref == 0)
            {
                max_reduction_factor = 0.90;
                min_reduction_factor = 0.01;
                alpha=0.01;
            }
            else
            {
                max_reduction_factor = 0.90;
                min_reduction_factor = 0.50;
                alpha=0.001;
            }
        }

        /**
         * Sum of squares
         */
        double sum_of_squares(double residual_norm);

        /**
         * Merit Function
         */
        double three_point_step (int line_search_iterations,double lambda_current, double lambda_previous, double mf_original, double mf_current, double mf_previous);

        /**
         * Determines step size by three point parabolic function
         */
        bool armijo(double step_length,double residual_norm,double new_residual_norm);

        /**
         * Convergance criterion
         *
         * To prevent a reduction in the step size that is either too large or too small
         * to be effective, we allow the user to specify a maximum and minimum reduction rate, i.e. \f$ \sigma_0 , \sigma_1 \f$.
         * Therefore, \f$ h_j \f$ must equal \f$ \sigma h_{j−1} \f$, where \f$ \sigma ∈ [ \sigma_0 , \sigma_1] \f$.
         * This is guaranteed by making the following check after \f$ h_j \f$  is determined
         * \f[
         * h_j = \left\{
         * \begin{array}{ll}
         * \sigma_0 h_{j−1} & \mbox{ if } h_j < \sigma_0 h_{j−1} \\
         * \sigma_1 h_{j−1} & \mbox{ if } h_j > \sigma_1 h_{j−1} \\
         * h_j  & \mbox{otherwise}.
         * \end{array}
         * \right.
         * \f]
         *
         * \p max_reduction_factor corresponds to \f$ \sigma_1 \f$
         */
        double max_reduction_factor;
        /**
         * \p min_reduction_factor corresponds to \f$ \sigma_0 \f$
         */
        double min_reduction_factor;
        /**
         * This parameter is used in order to assess when the line search should be stopped. This value should be
         * between zero and one. In practice, \f$ \alpha \f$ is typically \f$ 10^{−4}\f$.
         *
         * The larger the value of \f$ \alpha \f$, the strictest the convergence criteria, i.e.
         * \f[
         * \| \mathbf{F}(\mathbf{x_{i−1}} + h_j \delta \mathbf{x}) \| < (1 − \alpha h_j ) \| \mathbf{F}(x_{i−1}) \|.
         * \f]
         * will be.
         *
         * This value can be modified in the input file in subsection Newton using "set alpha = 1e-4".
         */
        double alpha;

        /**
         * Boolean variable that specifies if the three-point method is used with (true) or without (false)
         * prediction, i.e. the first iteration step is based on previous iteration steps.
         */
        bool with_prediction;

        /**
         * Flag so that predefined parameters are used.
         */
        bool use_predefined;
    };
}
}

#endif