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
// - Class: application_base.h
// - Description: This class implements the interface of applications
// - Developers: Guido Kanschat,     Texas A&M University
//               Valentin N. Zingan, University of Alberta
//               Marc Secanell,      University of Alberta
//
// ----------------------------------------------------------------------------

#ifndef _FUEL_CELL_APPLICATION_CORE_APPLICATION_BASE_H_
#define _FUEL_CELL_APPLICATION_CORE_APPLICATION_BASE_H_

// openFCST objects:
#include <application_core/application_data.h>
#include <application_core/fe_vectors.h>
#include <application_core/event.h>
#include <utils/logging.h>

// deal.II objects:
#include <deal.II/base/parameter_handler.h>

// C++ objects:
#include <boost/shared_ptr.hpp>

using namespace dealii;

namespace FuelCell
{
    /**
     * @brief Namespace containing the basic application framework used to loop over cells and create objects passed to FuelCellShop::Equation objects
     * 
     * This namespace encapsulates all the routines that are used to generate an application. ApplicationBase and ApplicationWrapper are declared
     * here. All applications inherit from the objects developed here.
     * 
     * @note This namespace should only be modified by advanced users as all applications depend on these routines
     */
    namespace ApplicationCore
    {
        /**
         * Base class for applications.
         *
         * <h3>Use of applications</h3>
         *
         * <h3>Application chains</h3>
         *
         * The idea of these application classes is the possibility to build a
         * chain out of these in order to have several predefined nested
         * solvers. For this, we distinguish applications roughly in two
         * classes, both derived from ApplicationBase:
         *
         * <ol>
         * <li> Linear applications 
         * <li> Wrapper applications 
         * </ol>
         * 
         * Linear applications implement the real finite element
         * code like computing residuals on mesh cells, assembling matrices
         * and solving linear systems. Linear application classes are the heart of openFCST. 
         * These classes use Equation, Layer, Kinetics, Material and Post-processing class objects 
         * in order to build and solve the underlying partial differential equations that describe 
         * the physical processes of the problem at hand, and output the solution in a convenient 
         * format that can be analyzed in ParaView, such as VTU format. 
         * 
         * openFCST uses the continuous Galerkin finite element method (CG FEM) to solve the fuel cell governing equations. 
         * This method for solving partial differential equations can be summarized in the following steps:
         * <ol>
         * <li> Discretization of the n-dimensional physical domain, \f$ \Omega \f$, into the set \f$ \mathbb{T}_h \f$ of 
         * disjoint non-overlapping elements \f$ K \f$, known commonly as either mesh or grid, 
         * <li> introduction of the finite element space \f$ \mathcal{V}_h^p \f$ consisting of continuous scalar-valued polynomial functions \f$ v_h \f$ 
         * of degree \f$ p > 0 \f$ such that for each element \f$ K \in \mathbb{T}_h: v_h = \sum_{i=1}^{N_s} v_{h,i} \psi_i(\mathbf{x}) \f$
         * where \f$ v_h \f$ is either the test or solution approximation function, \f$ N_s = (p+1)^{\text{dim}} \f$ is the number of support points, 
         * \f$ v_{h,i} \f$ are the values of \f$ v_h \f$ at the support points, and \f$ \psi_i(\mathbf{x}) \f$ are basis functions which 
         * are compactly supported on \f$ K \f$ and elements adjacent to \f$ K \f$,
         * <li> linearization of governing equations if needed,
         * <li> derivation, using the continuous Galerkin method, of the global and then element-wise weak formulation of the linearized problem, 
         * <li> assemble of the global linear system of equations using the finite element-wise weak form of the PDEs, 
         * <li> solution of the global linear system of equations using the numerical methods of linear algebra, and 
         * <li> post-processing of the solution.
         * </ol>
         * 
         * Linear application classes are directly responsible for steps 1, 2, 5 and 6. Steps 3 and 4 are performed by 
         * the Equation, Layer, Kinetics, Material objects inside the application. Step 7 is performed by Post-processing objects inside the linear application.
         * 
         * The execution flow of openFCST is as follows. First, declare_parameters is called to make sure the input file 
         * contains all the necessary information and to provide default data for any missing variables, then initialize()
         * is used to read the input file provided by the user and setup the problem computational domain, 
         * effective transport properties, and solver parameters. Then, the system of equations is constructed 
         * and solved using assembly, residual and solve. 
         *
         * Wrapper applications derived from ApplicationCopy; these usually implement a new solve() function as an iterative solver
         * around another application. Wrapper applications are used to implement adaptive refinement 
         * and nonlinear problems. Solving the aforementioned cases, generally involves solving a linear problem iteratively. 
         * Wrapper applications therefore use either another wrapper application that contains a linear application 
         * or a linear application directly in order to solve the required linear problem. 
         * They implement all functions of the ApplicationBase interface either by forwarding them to the next inner
         * application by ApplicationCopy or by providing their own implementation.
         *
         * @author Guido Kanschat, 2006
         * @author Marc Secanell, 2006-14
         * @author Valentin Zingan, 2012-14
         */
        
        
        class ApplicationBase : public Subscriptor
        {
        public:
            
            /**
             * Constructor for an application. One task of the constructor is declaring
             * parameters for an application. Since this is done consecutively for derived
             * classes, all parameters may be declared in a modular way. Derived classes should
             * create unique subsections in the ParameterHandler.
             *
             * @note Actually, this constructor does nothing, since the base class declares no parameters.
             */
            ApplicationBase(boost::shared_ptr<ApplicationData> data = boost::shared_ptr<ApplicationData>());
            
            /**
             * Copy constructor. Generate a new application sharing the same
             * ApplicationData with another one.
             */
            ApplicationBase(const ApplicationBase& other);
            
            /**
             * Virtual destructor.
             */
            virtual ~ApplicationBase();
            
            /**
             * Declare parameters for a parameter file. 
             * This member function is used to identify a parameter in the parameter file.
             *
             * @note Pure virtual function MUST be redeclared in derived classes.
             */
            virtual void declare_parameters(ParameterHandler& )
            {
                print_caller_name(__FUNCTION__);
            }
            
                        /**
             * Print default parameters for the application to a file.
             * Note that this member function MUST be called after declare_parameters.
             * 
             * Example:
             * \code
             * ApplicationBase app;
             * ParamterHandler param;
             * app.declare_parameters(param);
             * app.print_parameters_to_file(param, "default_parameters.txt", ParameterHandler::XML);
             * \endcode
             */
            void print_parameters_to_file(ParameterHandler& param,
                                          const std::string& file_name,
                                          const ParameterHandler::OutputStyle& style);
            /**
             * 
             * Pure virtual initialization function. This member function is used to read the parameters
             * declared in the constructor from ParameterHandler and then, used to initialize all object 
             * including the object that generates the mesh and Equation, Layer, Kinetics, Material objects.
             *
             * @note This function must be overloaded in any derived class using its own values
             * from the parameter file.
             *
             * @note Initialization functions of derived classes must make sure to call all functions
             * initialize() of the base classes.
             */
            virtual void initialize(ParameterHandler& param) = 0; // No implementation here
            
            /**
             * Generate the next mesh depending on the mesh generation parameters.
             *
             * Pure virtual.
             */
            virtual void remesh()
            {
                print_caller_name(__FUNCTION__);
            }
            
            /**
             * Initialize vector to problem size.
             *
             * Pure virtual.
             */
            virtual void init_vector(FEVector& ) const
            {
                print_caller_name(__FUNCTION__);
            }
            
            /**
             * Initialize vector to problem size.
             */
            virtual void start_vector(FEVector& dst, std::string ) const
            {
                print_caller_name(__FUNCTION__);
                init_vector(dst);
                
            };
            
            /**
             * Compute residual of <tt>src</tt> and store it
             * into <tt>dst</tt>. The bool variable is used
             * to specify if boundary conditions should be applied
             * using FilteredMatrix.
             *
             * Pure virtual.
             */
            virtual double residual(FEVector&               /*dst*/,
                                    const FEVectors&        /*src*/,
                                    bool apply_boundaries = true)
            {
                print_caller_name(__FUNCTION__);
                
                return -1;
            }
            
            /**
             * Solve the system assembled with right hand side in FEVectors <tt>src</tt> and return the
             * result in FEVector <tt>dst</tt>.
             *
             * Pure virtual.
             */
            virtual void solve(FEVector&        dst,
                               const FEVectors& src) = 0;
            
            /**
             * Solve the dual system assembled with right hand side <tt>rhs</tt> and return
             * the result in <tt>start</tt>.
             *
             * Pure virtual.
             */
            virtual void Tsolve(FEVector&        /*dst*/,
                                const FEVectors& /*src*/)
            {
                print_caller_name(__FUNCTION__);
            }
            
            /**
             * Estimate cell-wise errors.
             *
             * Pure virtual.
             */
            virtual double estimate(const FEVectors& )
            {
                print_caller_name(__FUNCTION__);                
                return -1;
            }
            
            /**
             * Evaluate a functional.
             *
             * Pure virtual.
             */
            virtual double evaluate(const FEVectors& )
            {
                print_caller_name(__FUNCTION__);
                return -1;
            }
            
            /**
             * Write the mesh in the format specified by the ParameterHandler.
             *
             * Pure virtual.
             */
            virtual void grid_out(const std::string& )
            {
                print_caller_name(__FUNCTION__);
            }
            
            /**
             * Write data in the format specified by the ParameterHandler.
             *
             * Pure virtual.
             */
            virtual void data_out(const std::string& filename,
                                  const FEVectors&   src)
            {
                print_caller_name(__FUNCTION__);
            }
            
            /**
             * Get access to the protected variable #data.
             */
            boost::shared_ptr<ApplicationData> get_data();
            
            /**
             * Get read-only access to the protected variable #data.
             */
            const boost::shared_ptr<ApplicationData> get_data() const;
            
            /**
             * Return a unique identification string for this application. In particular, it returns
             * the name of the class in which the object you requested id() from belongs to. 
             * 
             * For example:
             * \code
             * cBase* a = new cBase
             * cBase* b = new cDerived
             * \endcode
             * then,
             * \code
             * a->id() will return cBase and
             * b->id() will return cDerived
             * \endcode
             */
            virtual std::string id() const;
            
            /**
             * Add a reason for
             * assembling.
             */
            virtual void notify(const Event& reason);
            
            /**
             * All @p true in @p notifications.
             */
            virtual void clear();
            
            /**
             * All @p false in @p notifications.
             */
            virtual void clear_events(); 
            
            /**
             * Returns solution index. Pure virtual.
             */
            virtual unsigned int get_solution_index();
            
        protected:        
            
            /**
             * Print caller name.
             */
            void print_caller_name(const std::string& caller_name) const;
            
            /**
             * Object for auxiliary data. This object is always valid, since it is
             * either handed to the constructor or created there.
             */
            boost::shared_ptr<ApplicationData> data;
            
            /**
             * Accumulate reasons for
             * assembling here. If any of
             * those is set, the function
             * solve() of a terminal
             * application must take care
             * of assembling the matrix.
             */
            Event notifications;
            /**
             * Reasons for assembling
             */
            //std::map<std::string, bool> notifications;

            #ifdef OPENFCST_WITH_PETSC
            /**
             * MPI communication objects and processor identifiers
             */
            const unsigned int n_mpi_processes;
            const unsigned int this_mpi_process;
            MPI_Comm mpi_communicator;

            #endif
        };
    
    }// ApplicationCore

} // FuelCell

#endif
