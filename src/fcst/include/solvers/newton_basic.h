// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: newton_basic.h
// - Description: This class performs basic Newton iterations with a constant weight
// - Developers: Valentin N. Zingan, University of Alberta
// - Id: $Id: newton_basic.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _DEALII_APPFRAME_NEWTON_BASIC_H_
#define _DEALII_APPFRAME_NEWTON_BASIC_H_

#include <solvers/newton_base.h>

namespace FuelCell
{
    namespace ApplicationCore
    {

        /**
        * This class performs basic Newton iterations with a constant weight.
        *
        * \author Valentin N. Zingan, 2012
        */

        class NewtonBasic : public newtonBase
        {
            public:

                ///@name Constructors, destructor, and initialization
                //@{

                /**
                * Constructor.
                */
                NewtonBasic(ApplicationBase& app);

                /**
                * Destructor.
                */
                ~NewtonBasic();

                /**
                * Declare parameters.
                */
                virtual void declare_parameters(ParameterHandler& param);

                /**
                * Initialize parameters.
                */
                virtual void initialize(ParameterHandler& param);

                //@}

                ///@name Solve function
                //@{

                /**
                * This function implements
                * basic Newton iterations
                * with a constant weight.
                */
                virtual void solve(FEVector&        u,
                                    const FEVectors& in_vectors);

                //@}

                //////////
                // DATA //
                //////////

                ///@name Events
                //@{

                /**
                * This event is set
                * if the convergence
                * becomes bad.
                */
                static const Event bad_derivative;

                //@}

            private:

                //////////
                // DATA //
                //////////

                ///@name Basic Newton iteration parameters
                //@{

                /**
                * A constant weight is the same for all
                * basic Newton iterations.
                */
                double weight;

                //@}

        };

    }
}

#endif