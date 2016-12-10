//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 20014 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: app_test.cc
//    - Description: Simple application used to test different equation classes. 
//                   Specially useful to test derivatives.
//    - Developers: M. Secanell
//    - $Id: app_test.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------
// 
// #ifndef _FUELCELL__APP_TEST__H
// #define _FUELCELL__APP_TEST__H
// 
// // Include deal.II classes
// #include "base/parameter_handler.h"
// #include "base/function_lib.h"
// #include "base/function.h"
// #include "base/quadrature_lib.h"
// 
// #include "lac/block_vector.h"
// #include "lac/solver_cg.h"
// #include "lac/solver_gmres.h"
// #include "lac/precondition.h"
// #include "lac/precondition_block.h"
// #include "lac/block_matrix_array.h"
// #include "lac/filtered_matrix.h"
// #include "lac/sparse_ilu.h"
// #include "lac/sparse_direct.h"
// 
// #include "grid/grid_generator.h"
// #include "grid/tria_accessor.h"
// #include "grid/tria_iterator.h"
// #include "grid/tria_boundary_lib.h"
// 
// #include "fe/fe_values.h"
// 
// #include "numerics/vector_tools.h"
// #include "numerics/matrix_tools.h"
// #include "numerics/error_estimator.h"
// 
// #include "boost/shared_ptr.hpp"
// 
// // Include appframe classes
// #include "appframe/base.h"
// #include "block_matrix_application.h"
// #include "appframe/matrix_shop_cell.h"
// #include "appframe/residual_shop_cell.h"
// 
// // Include FuelCell classes
// #include "system_management.h"
// #include "optimization_block_matrix_application.h"
// #include "solver_utils.h"
// #include "operating_conditions.h"
// #include "dummy_GDL.h"
// #include "PureGas.h"
// #include "GasMixture.h"
// #include "geometries.h"
// #include "homogeneous_CL.h"
// #include "platinum.h"
// #include "tafel_kinetics.h"
// #include "fcst_constants.h"
// #include "ideal_gas.h"
// #include "SGL_24_BA.h"
// #include "design_fibrous_GDL.h"
// //Include STL
// #include "fstream"
// #include "iostream"
// #include "sstream"
// 
// // Use namespace of deal.II
// using namespace dealii;
// 
// namespace FuelCell
// {
//   /**
//    * This namespace is used for all those auxiliary classes that are used by AppTest.
//    * These are mainly classes that inherit Function<dim> and that are necessary in order to call
//    * some subroutine from deal.II. For example InitialSolution is created in order to use the
//    * deal.II class VectorInterpolate which in turn is used to set up the initial solution to the problem. 
//    */
//   namespace InitialSolution
//   {
//     /**
//      * This class is used when solving the problem using Newton's method to provide an initial solution.
//      * This function is called in VectorTools::interpolate(..,..,InitialSolution<dim> marc,...)
//      * It provides a solution that satisfies Dirichlet boundaries and has a gradient.
//      */
//     template <int dim>
//       //class InitialSolution
//       class AppTestIC
//       :
//     public Function<dim>
//     {
//     public:
//       /**
//        * Constructor
//        */
//       AppTestIC (FuelCell::OperatingConditions* );
//       /**
//        * Destructor
//        */
//       ~AppTestIC ();
//       /**
//        * 
//        */
//        double value (const Point<dim> &p, unsigned int) const;
//       /**
//        * This is the member function that computes the value of the initial
//        * solution for a given point.
//        */
//       void vector_value (const Point<dim> &p, 
// 			 Vector<double> &v) const;
//       
//     private:
//       /**
//       Operating conditions class object
//       */
//       FuelCell::OperatingConditions* OC;
//      
//     };
//   } //end namespace InitialSolution
//   
// 
// namespace Application
// {
//   //---------------------------------------------------------------------------
//   //---------------------------------------------------------------------------
//   //---------------------------------------------------------------------------
//   /**
//    * 
//    * This class is used to develop and test new applications. All functions are blank so all will
//    * have to be implemented.
//    */
//   template <int dim>
//     class AppTest
//     :
//   public FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>
//   {
//   public:
//     
//     /**
//      * Constructor.
//      * @note the pointer data is initialized to boost::shared_ptr<> (), this means that
//      * the pointer is empty and when we do data.get() it will return 0. This is good because at ApplicationBase
//      * constructor an ApplicationData will be constructed.
//      */
//     AppTest (boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data = 
// 			    boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> ());
//     
//     /**
//      * Destructor
//      */
//     ~AppTest();
//     /**
//      * Declare all parameters that are needed for: 
//      *   - the computation of the equation coefficients
//      *   - the control of the linear system solution
//      *   - ...
//      */
//     virtual void declare_parameters(ParameterHandler& param);
// 
//     /**
//      * Function called by optimization loop in order to set the values in the 
//      * ParameterHandler to the new design parameters.
//      * Since ParameterHandler depends on the problem we are solving, set_parameters() is set
//      * at the most inner loop of the application.
//      */
//     virtual void set_parameters(const std::vector<std::string>& name_dvar,
// 				const std::vector<double>& value_dvar,
// 				ParameterHandler& param);    
//     
//     /**
//      * Set up how many equations are needed and 
//      * read in parameters for the parameter handler in order to initialize data
//      */
//     void _initialize(ParameterHandler& param);
//     
//     /**
//      * Call the other initialize routines from the inherited classes
//      */
//     virtual void initialize(ParameterHandler& param);
//     
//     /**
//      * Initialize nonlinear solution
//      * Note that this routine is extremely important since it is used to resize the solution vector. Therefore,
//      * it always has to exist. Otherwise you will receive an error such as: 
//      * An error occurred in line <306> of file </home/secanell/Programs/fcst/contrib/deal.II/deal.II/include/deal.II/lac/block_indices.h> in function
//      * std::pair<unsigned int, unsigned int> dealii::BlockIndices::global_to_local(unsigned int) const
//      * The violated condition was: 
//      *    i<total_size()
//      *  The name and call sequence of the exception was: ExcIndexRange(i, 0, total_size())
//      * Additional Information: Index 0 is not in [0,0[
//      * To initialize the solution vector at least you need the following:
//      * src.reinit(this->block_info.global);
//      */
//     void init_solution(FuelCell::ApplicationCore::FEVector& src) ;
//     /**
//      * Integration of local bilinear form. Here we loop over the quadrature 
//      * points and over degrees of freedom in order to compute the matrix for the cell
//      * This routine depends on the problem at hand and is called by assemble() in DoF_Handler
//      * class
//      * The matrix to be assembled is:
//      \f[
//        \begin{array}{l}
//        M(i,j).block(0) = \int_{\Omega} a \nabla \phi_i \nabla \phi_j d\Omega + \int_{\Omega} \phi_i \frac{\partial f}{\partial u_0}\Big|_n \phi_j d\Omega \\
//        \end{array}
//      \f]
//      */
//     virtual void cell_matrix(MatrixVector& cell_matrices,
// 			     const typename DoFApplication<dim>::CellInfo& cell);  
//     /**
//      * Integration of the rhs of the equations. Here we loop over the quadrature 
//      * points and over degrees of freedom in order to compute the right
//      * hand side for each cell
//      * This routine depends on the problem at hand and is called by residual() in DoF_Handler
//      * class
//      * @note This function is called residual because in the case of nonlinear systems
//      * the rhs is equivalent to the residual
//      */
//     virtual void cell_residual(FuelCell::ApplicationCore::FEVector& cell_vector,
// 			       const typename DoFApplication<dim>::CellInfo& cell);
// 
//     /**
//      * Integration of local bilinear form.
//      */
//     virtual void cell_dresidual_dlambda(std::vector<FuelCell::ApplicationCore::FEVector >& cell_vector,
// 					const typename DoFApplication<dim>::CellInfo& cell,
// 					std::vector<std::vector<double> >& src);
//     
//     /**
//      * Member function used to set dirichlet boundary conditions.
//      * This function is application specific and it only computes the boundary_value
//      * values that are used to constraint the linear system of equations that is being
//      * solved
//      */
//     virtual void dirichlet_bc(std::map<unsigned int, double>& boundary_values) const;
//     
//     /**
//      * Solve the linear system of equations
//      */
//     virtual void solve(FuelCell::ApplicationCore::FEVector& start,
// 		       const FuelCell::ApplicationCore::FEVectors& rhs);
//     
//     /**
//      * Estimate error per cell
//      */
//     virtual double estimate(const FuelCell::ApplicationCore::FEVectors& sol);
//     
//     /**
//      * This class is called by responses to make sure that all responses requested are
//      * implemented in either cell_responses, global_responses or face_responses.
//      * @note Every time we add a new response that we can calculate we need to update this
//      * file.
//      */
//     virtual void check_responses();
//     
//     /**
//      * Compute the value of all objective function and constraints 
//      */
//     virtual void cell_responses (std::vector<double>& resp,
// 				 const typename DoFApplication<dim>::CellInfo& info,
// 				 const FuelCell::ApplicationCore::FEVector& sol);
//     /**
//      * This class is used to evaluate all responses that do not require looping over cells.
//      * An example of one of this types of constraints is the solid volume fraction.
//      */
//     virtual void global_responses (std::vector<double>& resp,
// 				   const FuelCell::ApplicationCore::FEVector& sol);
//     
//     /**
//      * This class is used to evaluate the derivative of all the functionals that require looping over cells
//      * with respect to the design variables.
//      * This class is called by responses to evaluate the response at each cell.
//      */
//     virtual void cell_dresponses_dl(std::vector<std::vector<double> >& cell_df_dl,
// 				    const typename DoFApplication<dim>::CellInfo& info,
// 				    const FuelCell::ApplicationCore::FEVector& sol);
//     
//     /**
//      * This class is used to evaluate the sensitivities of all responses that do not require looping over cells
//      * with respect of the design variables.
//      * An example of one of this types of constraints is the solid volume fraction.
//      */
//     virtual void global_dresponses_dl(std::vector<std::vector<double> >& df_dl,
// 				      const FuelCell::ApplicationCore::FEVector& sol);
//     /**
//      * This class is used to evaluate the derivative of all the functionals that require looping over cells
//      * with respect of the unknowns of the system of governing equations.
//      * This class is called by responses to evaluate the response at each cell.
//      */
//     virtual void cell_dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector >& cell_df_du,
// 				    const typename DoFApplication<dim>::CellInfo& info,
// 				    std::vector<std::vector<double> >& src);
// 
//     /**
//      * This class is used to evaluate the sensitivities of all responses that do not require looping over cells
//      * with respecto of the unknowns of the system of governing equations.
//      * An example of one of this types of constraints is the solid volume fraction.
//      */
//     virtual void global_dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector >& df_du,
// 				      const FuelCell::ApplicationCore::FEVector& src);
// 
//     /**
//      * Post-processing. Evaluate a functional such as the objective function of an
//      * optimization problem
//      */
//     virtual double evaluate (const FuelCell::ApplicationCore::FEVectors& src);
// 
//     /**
//      * Reimplementation of the routine in the base class BaseApplication in namespace AppFrame so
//      * that the right labels are outputed and so that I can compute and output the source term.
//      */
//     /*virtual void data_out(const std::string &basename, 
// 			  const FuelCell::ApplicationCore::FEVectors &src);
//     */
//   protected:
//     /**
//      * Structure where we store the problem we want to solve.
//      * Each vector component contains a string with the 
//      * name of the equation we want to solve
//      * Then, the number of components is equation_names.size()
//      */
//     std::vector<std::string> equation_names;
//     /**
//      * Structure where we store the name of each component
//      * in our problem. The component names are stored in the same
//      * way as they are stored in the solution.
//      */
//     std::vector<std::string> component_names;
//        
//     /** Operating conditions */
//     FuelCell::OperatingConditions OC;
//            
//     /** Create an object to generate the mesh*/
//     FuelCellShop::Geometry::Cathode<dim> grid;
//     
//     /** The object GDL layer will contain all the information relevant to the
//     * the GDL. We can request any effective property from this class
//     */
//     boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > CGDL;
// 
//     /** The object CCL layer will contain all the information relevant to the
//     * the catalyst layer. We can request any effective property from this class
//     */
//     FuelCellShop::Layer::HomogeneousCL<dim> CCL;  
//     
//     /** 
//       * List of solution variables that are required to solve the equations in the 
//       * kinetics class (e.g. oxygen concentration, protonic/electronic potential, 
//       * temperature etc.) Used by set_solution to ensure that all required variables are set.
//       */
//     std::vector<std::string> required_solution_names;    
//     
// 	FuelCell::SystemManagement sys_mgnt;
//   };
// 
// }
// }
// 
// #endif //_FUELCELL__APP_TEST_H
