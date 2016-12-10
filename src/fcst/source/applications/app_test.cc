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
//    - $Id: app_test.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

// #include <app_test.h>
// 
// namespace NAME = FuelCell::Application;
// using namespace FuelCell::ApplicationCore;

//---------------------------------------------------------------------------
// template <int dim>
// NAME::AppTest<dim>::AppTest(boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data)
//   :
//   OptimizationBlockMatrixApplication<dim>(data),  
//   CCL("Cathode catalyst layer"),
//   sys_mgnt(this->block_info, this->cell_couplings, this->flux_couplings)
// {
// 	//Move up: PhysicsBase1 pb1(sys_mgnt)
//   FcstUtilities::log << "->FuelCell::Application::App_test-" << dim <<"d"<<std::endl;
// }

// //---------------------------------------------------------------------------
// template <int dim>
// NAME::AppTest<dim>::~AppTest()
// {}

// //---------------------------------------------------------------------------
// template <int dim>
// void 
// NAME::AppTest<dim>::declare_parameters(ParameterHandler& param)
// {
//   OptimizationBlockMatrixApplication<dim>::declare_parameters(param);
// 
//   grid.declare_parameters(param);
//   FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Cathode gas diffusion layer", param);
//   CCL.declare_parameters(param);
//   
//   sys_mgnt.declare_parameters(param);
// }
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void 
// NAME::AppTest<dim>::set_parameters(const std::vector<std::string>& /*name_dvar*/,
// 						const std::vector<double>& /*value_dvar*/,
// 						ParameterHandler& /*param*/)
// {
// 
// }
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void
// NAME::AppTest<dim>::_initialize(ParameterHandler& param)
// {
// //Initialize the problem data:
//   // This specifies the type of problem we have
//   equation_names.clear();
//   equation_names.push_back("Thermal Transport");
//   component_names.clear();
//   component_names.push_back("Temperature");
//   sys_mgnt.initialize(param);
//   sys_mgnt.print_system_info();
//   
//   int index;
//   //FcstUtilities::log<<sys_mgnt.temperature(index)<<" with index :"<<index<<std::endl;
//   //FcstUtilities::log<<sys_mgnt.gas_pressure(index)<<" with index :"<<index<<std::endl;
//   
//   // Generate grid and refine
//   grid.initialize(param);
//   grid.generate_grid(*this->tr);
//   
//   CGDL = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Cathode gas diffusion layer",param);
//   
//   CCL.initialize(param);
//   // Specify the coupling between equations and variables and between fluxes:
//   // This is problem dependant:
//   // In this case the system is fully coupled since the source term of all equations depends on
//   // all parameters.
//   this->cell_couplings.reinit(this->element->n_blocks(), this->element->n_blocks());
//   this->cell_couplings(0,0) = DoFTools::always;
// 	
// 	/*
// 	//PhysicsBase1 pb1(sys_mgnt); //should be in constructor
// 	for (unsigned int i = 0; i<sys_mgnt.num_solution_variables(); ++i)
// 	{
// 		if (i ==1)
// 		{
// 			pb1.cell_couplings(this->cell_couplings(i)); //< in PhysicsBase Table< 1, DoFTools::Coupling > &
// 			pb2.cell_couplings(this->cell_couplings(i)); //< in PhysicsBase Table< 1, DoFTools::Coupling > &
// 		}
// 		else if (i == 2)
// 		{
// 			//Bla, bla, bla
// 		}
// 	}
//   */
//   // Initialize dof per cell for the system of equations and information on
//   //block sizes and indices
//   this->remesh_dofs();
//   // Initialize matrices and spartisity pattern for the whole system
//   this->remesh_matrices();
//   
//   // Print default parameter file:
//   //this->print_parameters_to_file(param,"default_parameters.xml",ParameterHandler::XML); 
// }
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void
// NAME::AppTest<dim>::initialize(ParameterHandler& param)
// {
//   OptimizationBlockMatrixApplication<dim>::initialize(param);
//   _initialize(param);
// }
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void
// NAME::AppTest<dim>::init_solution(FuelCell::ApplicationCore::FEVector& dst)
// {
//   // resize vector:
//   dst.reinit(this->block_info.global);
//   FuelCell::InitialSolution::AppTestIC<dim> initial_solution(&OC);
// 
//   VectorTools::interpolate (*this->dof, 
// 			    initial_solution, 
// 			    dst);  
// }
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void 
// NAME::AppTest<dim>::cell_matrix(MatrixVector& cell_matrices,
// 				       const typename DoFApplication<dim>::CellInfo& info)
// {
//   // -- Loop over each equation and set the matrix to the appropriate values:
//   for (unsigned int eq_index = 0; eq_index<this->element->n_blocks(); ++eq_index)
//     {
//       // -- Choose the correct finite element for the equation
//       unsigned int fe_type = this->element->block_to_base_index(eq_index).first;
//       const FEValuesBase<dim>& fe = info.fe(fe_type);
//       //const unsigned int n_quad = fe.n_quadrature_points;
//       
//       //In this case we only have one equation:
//        // -- Compute the coefficients for the equations:
//       if (equation_names[eq_index] == "Thermal Transport")
// 	  {
// 		  if (CGDL->belongs_to_material(info.cell->material_id()))
// 		  {
// 			  double k_T;
// 			  CGDL->effective_thermal_conductivity(k_T);
// 			  AppShop::Matrix::Cell::laplacian(cell_matrices[0].matrix, fe, 100);
// 		   }
// 		   else if (CCL.belongs_to_material(info.cell->material_id()))
// 		   {
// 			   AppShop::Matrix::Cell::laplacian(cell_matrices[0].matrix, fe, 1);  
// 		}
// 	}
//    }
// }
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void 
// NAME::AppTest<dim>::cell_residual(FuelCell::ApplicationCore::FEVector& cell_vector,
// 				  const typename DoFApplication<dim>::CellInfo& info)
// { 
//   // -- Loop over each equation and set the matrix to the appropriate values:
//   for (unsigned int eq_index = 0; eq_index<this->element->n_blocks(); ++eq_index)
//     {
//       //In this case we only have one equation:
//        // -- Compute the coefficients for the equations:
//       if (equation_names[eq_index] == "Thermal Transport")
// 	{
// 	  if (CGDL->belongs_to_material(info.cell->material_id()))
// 	  {
// 	    cell_vector.block(eq_index) = 0;
// 	  }
// 	  else if (CCL.belongs_to_material(info.cell->material_id()))
// 	  {
// 	    cell_vector.block(eq_index) = 0;
// 	  }
// 	}
// 	
//     }
// }
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void 
// NAME::AppTest<dim>::dirichlet_bc(std::map<unsigned int, double>& boundary_values) const
// {
//   std::vector<bool> comp_mask (this->element->n_blocks(), false);
//   comp_mask[0]=true;
//   
//   ConstantFunction<dim> PEM(250);
//   ConstantFunction<dim> Channel(200);
//   ConstantFunction<dim> Land(500);
// 
//   // -- boundary component 1: membrane/catalyst layer interface
//   VectorTools::interpolate_boundary_values (*this->dof,
// 					    grid.get_boundary_id("c_CL/Membrane"),
// 					    PEM,
// 					    boundary_values,
// 					    comp_mask);
//   // -- boundary component 2: area underneath ribs
//   VectorTools::interpolate_boundary_values (*this->dof,
// 					    grid.get_boundary_id("c_BPP/GDL"),
// 					    Land,
// 					    boundary_values,
// 					    comp_mask);
//  
//   // -- boundary component 3: area underneath gas channel
//   // No Dirichlet boundary conditions for phases but we have Dirichlet BC for gas composition. Since for the linear 
//   //system we solve for delta_u, the dirichlet BC for delta_u are to have the value to be zero 
//   VectorTools::interpolate_boundary_values (*this->dof,
// 					    grid.get_boundary_id("c_Ch/GDL"),
// 					    Channel,
// 					    boundary_values,
// 					    comp_mask);
// 
// } 
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void 
// NAME::AppTest<dim>::solve(FuelCell::ApplicationCore::FEVector& start,
// 				 const FuelCell::ApplicationCore::FEVectors& rhs)
//   
// { 
//   // Assemble system matrix:
//   this->assemble(rhs); //Note the second component of
// /* 
//   std::ofstream file_out;
//   file_out.open("system_matrix.txt");
//   this->matrix.print_formatted(file_out,3,true,0,0);
//   file_out.close();
//  */
//   //*********** UMFPACK *****
//   // Copy RHS to solution vector
//   FuelCell::ApplicationCore::FEVector residual;
//   residual.reinit(start);
//   std::ofstream fp_out;
//   double res = this->residual(residual, rhs, false);
//   fp_out.open("residual_vector.txt");
//   residual.print(fp_out);
//   fp_out.close();
//   
//   // Apply BC to matrix, rhs and sol
//   std::map<unsigned int, double> boundary_values;
//   this->dirichlet_bc(boundary_values);
//   MatrixTools::apply_boundary_values(boundary_values,
// 				     this->matrix,
// 				     start, 
// 				     residual);
// 
//   // UMFPACK can only receive a Vector, so transform the assign the BlockVector to a normal
//   // vector
//   Vector<double> start_copy;
//   start_copy = residual;
//   // Facorize matrix and solve:
//   SparseDirectUMFPACK solver;
//   solver.initialize(this->matrix);
//   solver.solve(start_copy);
//   start = start_copy;
// 
//   // Distribute solution:
//   //Return the solution to all the nodes taking into accout the hanging nodes
//   this->hanging_node_constraints.distribute(start);
// 
// }
// 
// //---------------------------------------------------------------------------
// template <int dim>
// double 
// NAME::AppTest<dim>::estimate(const FuelCell::ApplicationCore::FEVectors& vectors)
// { 
// //   // Reinitialize error estimator vector and set Newmann boundaries:
//   this->cell_errors.reinit(this->tr->n_active_cells());
// 
//   // Extract soluton from FEVectors
//   unsigned int ind = vectors.find_vector("Solution");
//   const BlockVector<double>& sol = vectors.vector(ind);
// 
//   //component_mask denote which components of the finite element space shall be used to 
//   //compute the error estimator. If it is left as specified by the default value (i.e. an empty array), all components are 
//   //used.
//   std::vector<bool> comp_mask (this->element->n_blocks());
//   
//   // Set refinement flags for material_id = 1
//   // -- material_id 1: catalyst layer
//   comp_mask[0] = true; // - comp 1 = x02
//       
//   this->cell_errors.reinit(this->tr->n_active_cells());
//   //typename FunctionMap<dim>::type newmann_boundary;
//   KellyErrorEstimator<dim>::estimate(*this->dof,
// 				     QGauss<dim-1>(3),
// 				     typename FunctionMap<dim>::type(),
// 				     sol,
// 				     this->cell_errors,
// 				     comp_mask);
//   
//   return 0.0;
// }
// 
// //---------------------------------------------------------------------------
// //---------------------------------------------------------------------------
// //---------------------------------------------------------------------------
// //-----------    POSTPROCESSING ROUTINES     --------------------------------
// //---------------------------------------------------------------------------
// //---------------------------------------------------------------------------
// //---------------------------------------------------------------------------
// 
// //---------------------------------------------------------------------------
// template <int dim>
// double 
// NAME::AppTest<dim>::evaluate (const FuelCell::ApplicationCore::FEVectors& vectors)
// { 
// 	return 0.0;
// }
// 
// //---------------------------------------------------------------------------
// /*template <int dim>
// void 
// NAME::AppTest<dim>::data_out(const std::string &basename, 
// 				    const FuelCell::ApplicationCore::FEVectors &vectors)
// {  
// 
// }
// */
// 
// //---------------------------------------------------------------------------
// //---------------------------------------------------------------------------
// //---------------------------------------------------------------------------
// //-----------    OPTIMIZATION ROUTINES     ----------------------------------
// //---------------------------------------------------------------------------
// //---------------------------------------------------------------------------
// //---------------------------------------------------------------------------
// //---------------------------------------------------------------------------
// template <int dim>
// void 
// NAME::AppTest<dim>::cell_dresidual_dlambda(std::vector<FuelCell::ApplicationCore::FEVector >&,
// 						  const typename DoFApplication<dim>::CellInfo&,
// 						  std::vector<std::vector<double> >& /*src*/)
// {
// 
// }
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void 
// NAME::AppTest<dim>::check_responses()
// {
// 
// }
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void 
// NAME::AppTest<dim>::cell_responses (std::vector<double>& resp,
// 					   const typename DoFApplication<dim>::CellInfo& info,
// 					   const FuelCell::ApplicationCore::FEVector& /*src*/)
// { 
// 
// }
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void 
// NAME::AppTest<dim>::global_responses(std::vector<double>&,
// 						  const FuelCell::ApplicationCore::FEVector& /*src*/)
// {
// 
// } 
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void 
// NAME::AppTest<dim>::cell_dresponses_dl(std::vector<std::vector<double> >&,
// 					      const typename DoFApplication<dim>::CellInfo&,
// 					      const FuelCell::ApplicationCore::FEVector& /*src*/)
// {
// 
// }
// 
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void 
// NAME::AppTest<dim>::global_dresponses_dl(std::vector<std::vector<double> >&,
// 						      const FuelCell::ApplicationCore::FEVector& )
// {
// 
// }
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void 
// NAME::AppTest<dim>::cell_dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector >&,
// 						    const typename DoFApplication<dim>::CellInfo&,
// 						    std::vector<std::vector<double> >& /*src*/)
// {
// 
// }
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void 
// NAME::AppTest<dim>::global_dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector >&,
// 						      const FuelCell::ApplicationCore::FEVector& /*src*/)
// {
// 
// }
// 
// 
// //---------------------------------------------------------------------------
// //---------------------------------------------------------------------------
// //---------------------------------------------------------------------------
// //---------------------------------------------------------------------------
// //---------------------------------------------------------------------------
// template <int dim>
// FuelCell::InitialSolution::AppTestIC<dim>::AppTestIC(FuelCell::OperatingConditions* OC)
//   :
//   Function<dim> (1),
//   OC(OC)
// {}
// 
// //---------------------------------------------------------------------------
// template <int dim>
// FuelCell::InitialSolution::AppTestIC<dim>::~AppTestIC()
// {}
// 
// //---------------------------------------------------------------------------
// template <int dim>
// double
// FuelCell::InitialSolution::AppTestIC<dim>::value (const Point<dim> &/*p*/, 
// 						  unsigned int /*i*/) const
// {
//   return 0;
// }
// 
// //---------------------------------------------------------------------------
// template <int dim>
// void
// FuelCell::InitialSolution::AppTestIC<dim>::vector_value (const Point<dim> &/*p*/, 
// 								Vector<double> &v) const
// {
//   for(unsigned int i = 0; i<v.size(); ++i)
//   {
//     v(i)=0; 
//   }
// }
// 
// //---------------------------------------------------------------------------
// //---------------------------------------------------------------------------
// //---------------------------------------------------------------------------
// // Explicit instantiations. 
// template class NAME::AppTest<deal_II_dimension>;
// template class FuelCell::InitialSolution::AppTestIC<deal_II_dimension>;