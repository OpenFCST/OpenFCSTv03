//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: picard_base.h
//    - Description: Base class for all variants of the Picard solver
//    - Developers: Mayank Sabharwal
//
//---------------------------------------------------------------------------

#ifndef __deal2__appframe__picard_base_h
#define __deal2__appframe__Picard_base_h

#include <deal.II/lac/solver_control.h>

#include <application_core/application_wrapper.h>

using namespace dealii;
using namespace FuelCell::ApplicationCore;

namespace FuelCell
{
namespace ApplicationCore
{
    /**
     * Base class for all classes performing Picard iterations. Children of this class mainly differ
     * in the formulation of the problem. Currently only one algorithm which inverts the complete stiffness
     * matrix is implemented.
     *
     * <h3> Parameters </h3>
     *
     * The parameters that control Picard solvers are defined in the subsection Picard in the data input file.
     * The parameters are the following:
     * \code
     * subsection Picard
     *  set Max steps          = 100           # Maximum number of iterations
     *  set Absolute tolerance = 1.e-8         # Absolute tolerance
     *  set Relative tolerance = 1.e-3         # Relative tolerance
     *  set Debug level        = 0
     *  set Debug residual     = false         # Would you like the code to output the residual at every Picard iteration?
     *  set Debug solution     = false         # Would you like the code to output the solution at every Picard iteration?
     * end
     * \endcode
     *
     * @author Mayank Sabharwal, 2015
     */

    class PicardBase:public ApplicationWrapper
    {
    public:

        /**
         * Constructor, receiving the application computing the residual and
         * solving the linear problem.
         */
        PicardBase(ApplicationBase& app);

        /**
         * Declare the input parameters. The following are the parameters used and their default values:
         * In section: "Picard"
         * - "Absolute Tolerance",
         * - "Relative Tolerance",
         * - "Max steps",
         */
        virtual void declare_parameters (ParameterHandler& param);

        /**
         * Read the parameters.
         */
        virtual void initialize (ParameterHandler& param);

	/**
         * Control object for the Picard iteration.
         */
        void initialize_initial_guess(BlockVector<double>& dst){};

        /**
         * Instead of assembling, this function only sets a flag, such that
         * the inner application will be required to assemble a new derivative
         * matrix next time solve() is called.
         */
        virtual void assemble();

        /**
         * Returns the L2-norm of the residual and the residual vector in "dst" using
         * the residual function in the ApplicationBase used to initialize the application.
         * The FEVectors rhs should contain at least the following vectors:
         * -
         *
         * \note Finish this section
         */
        virtual double residual(FuelCell::ApplicationCore::FEVector& dst,
                                const FuelCell::ApplicationCore::FEVectors& rhs);

        /**
         * The actual Picard solver.
         *
         * In this section, the Picard solver performs the iterations necessary to solve the nonlinear system.
         * During the process it will call assemble() and residual() as needed. Once it is done, it will return
         * the solution as parameter \b{u}.
         *
         */
        virtual void solve(FuelCell::ApplicationCore::FEVector& u,
                           const FuelCell::ApplicationCore::FEVectors& in_vectors) = 0;

        


    protected:
        /**
         * Function used to output any necessary information at each Picard iteration
         * for debugging purposes
         */
        void debug_output(const FEVector& sol,
                          const FEVector& update,
                          const FEVector& residual) const;

        /**
         * This flag is set by the function assemble(),
         * indicating that the matrix must be assembled anew upon
         * start.
         */
        bool assemble_now;

        /**
         * Absolute tolerance for the Picard convergence loop
         */
        double abs_tolerance;
        
        /**
         * Relative tolerance for the Picard convergence loop
         */
        double rel_tolerance;
        
	/**
	 * Maximu no. of steps for the Picard iterations"
	 */
	int maxsteps;
        /**
         * Print updated solution after each step into file
         */
        bool debug_solution;

        /**
         * Print Picard residual after each step into file
         * <tt>Picard_rNNN</tt>?
         */
        bool debug_residual;
        /**
         * Write debug output to #FcstUtilities::log; the higher the
         * number, the more output.
         */
        unsigned int debug;
        /**
         * The number of a basic Picard iteration.
         */
        unsigned int step;
        
        /**
         * Number of Iterations;
         * */
        double numIter;

        /**
         * This vector specifies the blocks of the global solution which are supposed to be treated specially.
         * Such blocks usually include the density or molar fraction of oxygen (burns as a result of chemical reactions)
         * and the density or molar fraction of hydrogen (burns as well as a result of chemical reactions).
         * These blocks are defined in the parameter files manually. Therefore, you need to know them beforehand.
         * Knowing these blocks, the non-linear solver checks against the negative values in those blocks and reduces the solution update if needed to keep the values positive.
         */
        std::vector<unsigned int> blocks;

        /**
         * The total number of \p blocks.
         */
        unsigned int n_blocks;
    };
}
}

#endif