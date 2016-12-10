//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: picard.h
//    - Description: Picard Solver for Non-linear problems
//    - Developers: Mayank Sabharwal
//
//---------------------------------------------------------------------------

#ifndef __deal2__appframe__picard_h
#define __deal2__appframe__picard_h

#include <solvers/picard_base.h>

namespace FuelCell
{
namespace ApplicationCore
{
    /**
     * Application class performing a Picard iteration as described in the PicardBase.
     * This class inverts the entire stiffness matrix based on Algorithm 1 described in J.N. Reddy.
     * The iterations for next guess are calculated as follows:
     * \f[
     * u^N = \left( K(u^{N-1}) \right)^{-1} F(u^{N-1})
     * \f]
     * 
     * <h3> Adaptive Underrelaxation </h3>
     * An adaptive underrelaxation method is available with the Picard solver as proposed by Durbin and Delemos[1]. 
     * The scheme involves modifying the underrelaxation factor &gamma; with the error. As the iterations proceed the value of &gamma; approaches unity. The expression used to evaluate &gamma; is,
     * \f[
     * \gamma = \gamma _{min} + (1-\gamma _{min})exp(-\alpha (\delta -\epsilon)) \qquad for\  \delta > \epsilon
     * \f] 
     * and 
     * \f[
     * \gamma = 1 \qquad for \ \delta<\epsilon
     * \f]
     * where &delta; is the L\_infinity norm of the change in the solution,
     * \f[
     * \delta = max_i\mid u^{n+1}_i - u^n_i \mid
     * \f]
     * In addition if the solution starts to diverge then the parameter &gamma;<SUB>min</SUB> is scaled down to minimize divergence.
     * 
     * <EM> [1] Durbin, Timothy, and David Delemos. "Adaptive underrelaxation of Picard iterations in ground water models." Groundwater 45.5 (2007): 648-651.</EM>
     * 
     * 
     * @author Mayank Sabharwal, 2015
     */
    class Picard : public PicardBase
    {
    public:
        /**
         * Constructor, receiving the application computing the residual and
         * solving the linear problem.
         */
        Picard (ApplicationBase& app);

        /** Declare the input parameters. */
        virtual void declare_parameters (ParameterHandler& param);

        /** Read the parameters */
        virtual void initialize (ParameterHandler& param);

        /**
         * The actual Picard solver.
         */
        virtual void solve(FuelCell::ApplicationCore::FEVector& u,
                           const FuelCell::ApplicationCore::FEVectors& in_vectors);
    private:
        
        /**
         * Flag for using adaptive under-relaxation
         */
        bool underrelaxation;
        
        /**
         * Alpha value to be used in the under-relaxation scheme proposed in Durbin and Delemos(2007)
         */
        double alpha;
        
        /**
         * Gamma_min value to be used in the under-relaxation scheme proposed in Durbin and Delemos(2007)
         */
        double gamma_min;
        
        void compute_errors ( FEVector &u, FEVector &u_n, FEVector &error, double &abs_error, double &rel_error, double &delta);


    };
}
}

#endif