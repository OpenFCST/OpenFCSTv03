//---------------------------------------------------------------------------
//
//    Copyright (C) 2006, 2007, 2008, 2009 by Marc Secanell
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <applications/app_read_mesh.h>


namespace NAME = FuelCell::Application;
using namespace FuelCell::ApplicationCore;

//---------------------------------------------------------------------------
template <int dim>
NAME::AppReadMesh<dim>::AppReadMesh(boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data)
:
OptimizationBlockMatrixApplication<dim>(data)
{
	FcstUtilities::log << "->FuelCell::Application::AppReadMesh_test-" << dim <<"d"<<std::endl;
}

//---------------------------------------------------------------------------
template <int dim>
NAME::AppReadMesh<dim>::~AppReadMesh()
{}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppReadMesh<dim>::declare_parameters(ParameterHandler& param)
{
	OptimizationBlockMatrixApplication<dim>::declare_parameters(param);
}



//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppReadMesh<dim>::_initialize(ParameterHandler& param)
{
	//Initialize the problem data:
	// This specifies the type of problem we have
	equation_names.clear();
	equation_names.push_back("Oxygen Diffusion");

	component_names.clear();
	component_names.push_back("Oxygen molar fraction");

	// Make sure that the number of input finite elements matches with the problem:
	Assert (this->element->n_blocks() == equation_names.size(),
			  ExcDimensionMismatch(this->element->n_blocks(), equation_names.size()));
	Assert (this->element->n_blocks() == component_names.size(),
			  ExcDimensionMismatch(this->element->n_blocks(), component_names.size()));
    
	// Specify the coupling between equations and variables and between fluxes:
	// This is problem dependant:
	// In this case the system is fully coupled since the source term of all equations depends on
	// all parameters.
	this->cell_couplings.reinit(this->element->n_blocks(), this->element->n_blocks());
	this->cell_couplings(0,0) = DoFTools::always;

	// Initialize matrices and spartisity pattern for the whole system
	this->remesh_matrices();

}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppReadMesh<dim>::initialize(ParameterHandler& param)
{
	OptimizationBlockMatrixApplication<dim>::initialize(param);
	_initialize(param);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppReadMesh<dim>::initialize_solution(FuelCell::ApplicationCore::FEVector& dst,
                                            std::shared_ptr<Function<dim> > initial_function)
{ 
	bool solution_read_success = false;
	bool use_transfer_solution = false;
	// resize vector:
	dst.reinit(this->block_info.global);

	if (!use_transfer_solution)
	{
		FuelCell::InitialSolution::AppReadMeshIC<dim> initial_solution(component_names);
		FcstUtilities::log << "Generating default initial solution" << std::endl;
		dst.reinit(this->block_info.global);
		VectorTools::interpolate (*this->dof,
									initial_solution,
									dst);
	}
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::AppReadMesh<dim>::data_out(const std::string &basename,
								const FuelCell::ApplicationCore::FEVectors &vectors)
{
    // Read the solution vector out of the many vectors:
    unsigned int i = vectors.find_vector("Solution");
    const BlockVector<double>& src = vectors.vector(i);

    // Output solution with its corresponent names:
    std::vector<std::string> solution_names;
    solution_names.push_back("dummy");

    // --- solution interpretations ---

    this->solution_interpretations.clear();
    this->solution_interpretations.resize(this->element->n_blocks(),
                                          DataComponentInterpretation::component_is_scalar);

    FuelCell::ApplicationCore::DoFApplication<dim>::data_out(basename,
                                            src,
                                            solution_names);

	// Read the solution vector out of the many vectors:

	const std::string suffix = this->d_out.default_suffix();
	if (suffix == "") return;

	//construct the full filename, where basename is given by the adaptive refinement loop and suffix is the file extension (i.e. .vtk)
	const std::string filename = basename + suffix;
	FcstUtilities::log << "Datafile:" << filename << std::endl;

	std::ofstream of (filename.c_str ());
	this->d_out.build_patches();
	this->d_out.write (of);
	this->d_out.clear();

	FcstUtilities::log << std::endl;
	FcstUtilities::log << "----- Mesh information -----" << std::endl;
	FcstUtilities::log << "- Number of cells: " << this->tr->n_active_cells() << std::endl;
	FcstUtilities::log << std::endl;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
FuelCell::InitialSolution::AppReadMeshIC<dim>::AppReadMeshIC(std::vector<std::string> names)
:
Function<dim> (1)
{
	set_solution_names(names);
}

//---------------------------------------------------------------------------
template <int dim>
FuelCell::InitialSolution::AppReadMeshIC<dim>::~AppReadMeshIC()
{}

//---------------------------------------------------------------------------
template <int dim>
double
FuelCell::InitialSolution::AppReadMeshIC<dim>::value (const Point<dim> &/*p*/, const unsigned int) const
{
	// size checking for vectors
	double x_o2 = 0.0;

	for(unsigned int i = 0; i<component_names.size(); ++i)
	{
		if (component_names[i] == "Oxygen molar fraction")
			return x_o2;
		else
		{
			FcstUtilities::log << "Unknown solution name: "<< component_names[i] << " used in initial solution" << std::endl;
			exit(1);
		}
	}

}

//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::InitialSolution::AppReadMeshIC<dim>::vector_value (const Point<dim> &/*p*/,
																				Vector<double> &v) const
{
	// size checking for vectors
	Assert(v.size() == this->n_components,
				ExcDimensionMismatch (v.size(), this->n_components));
	;
	double x_o2 = 0.0;

	for(unsigned int i = 0; i<component_names.size(); ++i)
	{
		if (component_names[i] == "Oxygen molar fraction")
			v(i) = x_o2;
		else
		{
			FcstUtilities::log << "Unknown solution name: "<< component_names[i] << " used in initial solution" << std::endl;
			exit(1);
		}
	}

}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::AppReadMesh<deal_II_dimension>;
template class FuelCell::InitialSolution::AppReadMeshIC<deal_II_dimension>;
