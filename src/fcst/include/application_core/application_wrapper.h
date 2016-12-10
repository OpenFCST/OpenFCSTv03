// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2009 by Guido Kanschat
// Copyright (C) 2006-2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: application_wrapper.h
// - Description: This class implements either iterative or time-stepping
//                wrapper of applications
// - Developers: Guido Kanschat,     Texas A&M University
//               Valentin N. Zingan, University of Alberta
//               Marc Secanell,      University of Alberta
//
// ----------------------------------------------------------------------------

#ifndef _FUEL_CELL_APPLICATION_CORE_APPLICATION_WRAPPER_H_
#define _FUEL_CELL_APPLICATION_CORE_APPLICATION_WRAPPER_H_

#include <application_core/application_base.h>

using namespace dealii;

namespace FuelCell
{
    namespace ApplicationCore
    {
        
        /**
         * This class implements either iterative or time-stepping
         * wrapper of applications.
         *
         * @author Guido Kanschat
         */
        
        class ApplicationWrapper : public ApplicationBase
        {
        public:
            
            /**
             * Constructor for a derived application. A SmartPointer
             * to <tt>app</tt> is stored in the wrapper. Therefore,
             * <tt>app</tt> must live longer than the
             * ApplicationWrapper object.
             */
            ApplicationWrapper(ApplicationBase& app);
            
            /**
             * Destructor.
             */
            ~ApplicationWrapper();
            
            /**
             * Declare parameters for a parameter file.
             * 
             * Current implementation simply calls the declare_parameters 
             * of the derived application.
             */
            virtual void declare_parameters(ParameterHandler& param);
            
            virtual void initialize(ParameterHandler& param);
            
            virtual void remesh();
            
            virtual void init_vector(FEVector& dst) const;
            
            virtual double residual(FEVector&        dst,
                                    const FEVectors& src,
                                    bool             apply_boundaries = true);
            
            virtual void solve(FEVector&        dst,
                               const FEVectors& src);
            
            virtual void Tsolve(FEVector&        dst,
                                const FEVectors& src);
            
            virtual double estimate(const FEVectors& src);
            
            virtual double evaluate(const FEVectors& src);
            
            virtual void grid_out(const std::string& filename) const;
            
            virtual void data_out(const std::string& filename,
                                  const FEVectors&   src);
            
            virtual std::string id() const;
            
            virtual void notify(const Event& reason);
            
        protected:
            /**
             * Gain access to the inner application
             */
            SmartPointer<ApplicationBase> get_wrapped_application()
            {
                return app;
            }
            
            /**
             * Pointer to the application
             * this one depends upon.
             */
            SmartPointer<ApplicationBase> app;
        };
        
    } // ApplicationCore
    
} // FuelCell

#endif