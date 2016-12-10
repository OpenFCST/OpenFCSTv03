// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2009 by Guido Kanschat
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: application_data.h
// - Description: This class implements general data of applications
// - Developers: Guido Kanschat,     Texas A&M University
//               Valentin N. Zingan, University of Alberta
//               Marc Secanell,      University of Alberta
//               Mayank Sabharwal,   University of Alberta
//               Aslan Kosakian,     University of Alberta
//
// ----------------------------------------------------------------------------

#ifndef _FUEL_CELL_APPLICATION_CORE_APPLICATION_DATA_H_
#define _FUEL_CELL_APPLICATION_CORE_APPLICATION_DATA_H_
//-- dealII
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/base/parameter_handler.h>

//-- OpenFCST
#include <utils/logging.h>

//-- C++ Standard libraries
#include <algorithm>
#include <fstream>


using namespace dealii;

namespace FuelCell
{
namespace ApplicationCore
{

/**
 * The vector class used by applications.
 */
typedef BlockVector<double> FEVector;

/**
 * Enumeration class for Linear solvers
 */
enum class LinearSolver {UMFPACK,CG,ILU_GMRES,MUMPS,BICGSTAB};

/**
 * Enumeration class for Non-linear solvers
 */
enum class NonLinearSolver {NONE,NEWTONBASIC,NEWTON3PP,NEWTONLINESEARCH,PICARD};

/**
 * Enumeration class for refinement types
 */
enum class RefinementSolver {ADAPTIVE};

/**
 * Here we handle general data of applications. It is the idea of this
 * class that all applications inheriting from ApplicationBase share
 * the <b>same</b> ApplicationData object, thus accessing the same
 * vector pools and using the object as a means for communicating
 * global information.
 *
 * First purpose is providing vector memory objects for Vector and
 * BlockVector objects for applications. These objects, #vector_pool
 * and #block_vector_pool can be accessed directly.
 *
 * Furthermore, it provides a map which allows applications to
 * exchange data in named registers. To this end, the functions
 * enter_flag(std::string,const bool&), enter(std::string,const double&),
 * enter(std::string,const FEVector&), and enter(std::string,const std::vector<double>&)
 * can be used to enter boolean flag, scalar and vector data, respectively.
 * These data are entered as references, so
 * the controlling application can change them at any time. All other
 * applications sharing the same ApplicationData can access them read
 * only through the functions flag(), scalar(), vector(), and std_vector().
 *
 * @author Guido Kanschat
 * @author Valentin N. Zingan
 */

class ApplicationData
{
    
public:
    /**
     * Constructor.
     */
    ApplicationData();
    
    /**
     * Destructor, releasing all data.
     */
    ~ApplicationData();
    
    /**
     * Empty member function at the moment.
     */
    void declare_parameters(ParameterHandler& param) const;
    
    /**
     * Empty member function at the moment.
     */
    void initialize(ParameterHandler& param);
    
    
    /**
     * Register a named boolean flag.
     */
    void enter_flag(std::string name,
                    const bool& s);
    
    /**
     * Register a named scalar.
     */
    void enter(std::string   name,
               const double& s);
    
    /**
     * Register a named vector.
     */
    void enter(std::string     name,
               const FEVector& v);
    
    /**
     * Register a named std vector.
     */
    void enter(std::string                name,
               const std::vector<double>& v);
    
    /**
     * Register a named std vector of std vectors.
     */
    void enter(std::string                name,
               const std::vector< std::vector<double> >& v);    
    
    /**
     * Delete a previously registered boolean flag.
     */
    void erase_flag(std::string name);
    
    /**
     * Delete a previously registered scalar.
     */
    void erase_scalar(std::string name);
    
    /**
     * Delete a previously registered vector.
     */
    void erase_vector(std::string name);
    
    /**
     * Delete a previously registered std vector.
     */
    void erase_std_vector(std::string name);
    
    /**
     * Delete a previously registered std vector of std vectors.
     */
    void erase_std_vector_of_std_vectors(std::string name);
    
    /**
     * This function returns @p true
     * if a boolean flag exists (doesn't matter if the flag itself is @p true or @p false).
     * Otherwise returns @p false.
     */
    bool flag_exists(const std::string& name) const;
    
    /**
     * Get read-only access to a registered boolean flag.
     */
    bool flag(std::string name) const;
    
    /**
     * Get read-only access to a registered scalar. It only offers read access and
     * returns a null pointer, if the name has not been registered.
     */
    const double* scalar(std::string name) const;
    
    /**
     * Get read-only access to a registered vector. It only offers read access and
     * returns a null pointer, if the name has not been registered.
     */
    const FEVector* vector(std::string name) const;
    
    /**
     * Get read-only access to a registered std vector. It only offers read access and
     * returns a null pointer, if the name has not been registered.
     */
    const std::vector<double>* std_vector(std::string name) const;

    /**
     * Get read-only access to a registered std vector of std vectors. It only offers read access and
     * returns a null pointer, if the name has not been registered.
     */
    const std::vector< std::vector<double> >* std_vector_std_vector(std::string name) const;
    
    /**
     * List all stored objects to FcstUtilities::log.
     */
    void log() const;
    
    /**
     * VectorMemory object for simple vectors. All applications should allocate
     * their vectors here, so they can be reused and operating system memory management can
     * be avoided.
     */
    GrowingVectorMemory< Vector<double> > vector_pool;
    
    /**
     * VectorMemory object for block vectors. All applications should allocate
     * their vectors here, so they can be reused and operating
     * system memory management can be avoided.
     */
    GrowingVectorMemory<FEVector> block_vector_pool;
    
    /**
     * Function to return the linear solver type
     */
    FuelCell::ApplicationCore::LinearSolver get_linear_solver()
    {
        return this->lin_solver;
    }
    
    /**
     * Function to return the non-linear solver type
     */  
    FuelCell::ApplicationCore::NonLinearSolver get_nonlinear_solver()
    {
        return this->nonlin_solver;
    }   
    
    /**
     * Function to return refinement type
     */
    FuelCell::ApplicationCore::RefinementSolver get_refinement_solver()
    {
        return this->refinement_solver;
    }
    
    /**
     * Function to set the type of linear solver
     */ 
    void set_linear_solver(const std::string name)
    {
        if (name.compare("UMFPACK") == 0)
            this->lin_solver = FuelCell::ApplicationCore::LinearSolver::UMFPACK;
        else if ( name.compare("CG") == 0 )
            this->lin_solver = FuelCell::ApplicationCore::LinearSolver::CG;
        else if ( name.compare("ILU-GMRES") == 0 )
            this->lin_solver = FuelCell::ApplicationCore::LinearSolver::ILU_GMRES;
        else if ( name.compare("MUMPS") == 0 )
            this->lin_solver = FuelCell::ApplicationCore::LinearSolver::MUMPS;
        else if ( name.compare("Bicgstab") == 0)
            this->lin_solver = FuelCell::ApplicationCore::LinearSolver::BICGSTAB;
        else
            throw(ExcNotFound("linear_solver", name));
    }
    
    /**
     * Function to set the type of non-linear solver
     */ 
    void set_nonlinear_solver(const std::string name)
    {
        if (name.compare("None") == 0)
            this->nonlin_solver = FuelCell::ApplicationCore::NonLinearSolver::NONE;
        else if ( name.compare("NewtonBasic") == 0 )
            this->nonlin_solver = FuelCell::ApplicationCore::NonLinearSolver::NEWTONBASIC;
        else if ( name.compare("Newton3pp") == 0 )
            this->nonlin_solver = FuelCell::ApplicationCore::NonLinearSolver::NEWTON3PP;
        else if( name.compare("NewtonLineSearch") == 0)
            this->nonlin_solver = FuelCell::ApplicationCore::NonLinearSolver::NEWTONLINESEARCH;
        else if ( name.compare("Picard") == 0 )
            this->nonlin_solver = FuelCell::ApplicationCore::NonLinearSolver::PICARD;
        else
            throw(ExcNotFound("nonlinear_solver", name));
    }
    
    /**
     * Function to set the type of refinement
     */ 
    void set_refinement_solver(const std::string name)
    {
        if (name.compare("AdaptiveRefinement") == 0)
            this->refinement_solver = FuelCell::ApplicationCore::RefinementSolver::ADAPTIVE;
        else
            throw(ExcNotFound("refinement_solver", name));
    }
    
    /**
     * Function to return solution vector name in the FEVectors object
     */  
    std::string get_solution_vector_name(FuelCell::ApplicationCore::NonLinearSolver in)
    {
        if ((in == FuelCell::ApplicationCore::NonLinearSolver::NEWTONBASIC) || (in == FuelCell::ApplicationCore::NonLinearSolver::NEWTON3PP) || (in == FuelCell::ApplicationCore::NonLinearSolver::NEWTONLINESEARCH))
            return "Newton iterate";
        else if ((in == FuelCell::ApplicationCore::NonLinearSolver::PICARD) || (in == FuelCell::ApplicationCore::NonLinearSolver::NONE))
            return "Solution";         
    }
    
    /**
     * Function to return the residual vector name in the FEVectors object
     */
    std::string get_residual_vector_name(FuelCell::ApplicationCore::NonLinearSolver in)
    {
        if ((in == FuelCell::ApplicationCore::NonLinearSolver::NEWTONBASIC) || (in == FuelCell::ApplicationCore::NonLinearSolver::NEWTON3PP) || (in == FuelCell::ApplicationCore::NonLinearSolver::NEWTONLINESEARCH))
            return "Newton residual";
        else if ((in == FuelCell::ApplicationCore::NonLinearSolver::PICARD) || (in == FuelCell::ApplicationCore::NonLinearSolver::NONE))
            return "residual";         
    }
    
    /**
     * The typedef for the map of boolean flags.
     */
    typedef std::map< std::string, bool > flag_map;
    
    /**
     * The typedef for the map of scalars.
     */
    typedef std::map< std::string, const double* > scalar_map;
    
    /**
     * The typedef for the map of
     * vectors.
     */
    typedef std::map< std::string, const FEVector* > vector_map;
    
    /**
     * The typedef for the map of std vectors.
     */
    typedef std::map< std::string, const std::vector<double>* > std_vector_map;
    
    /**
     * The typedef for the map of std vectors of std vectors.
     */
    typedef std::map< std::string, const std::vector< std::vector<double> >* > std_vector_std_vector_map;
    
    /**
     * Map of field data with names corresponding to the physical meaning of data, which can be read from VTK files. 
     * This can be used for passing Knudsen raidus of the microstructures by mapping Knudsen radius to cell index or any such data which can be mapped to cell index.
     */ 
    std::map<std::string,std::map<int,double>> field_data;
    
    /**
     * Exception thrown when a named scalar or vector was searched but not found.
     */
    DeclException2(ExcNotFound,
                   char*,
                   std::string,
                   << "A " << arg1 << " with name " << arg2 << " was not stored in this data object");
    
private:
    
    /**
     * A map linking names of data to actual boolean flags.
     */
    flag_map named_flags;
    
    /**
     * A map linking names of data to actual scalars.
     */
    scalar_map named_scalars;
    
    /**
     * A map linking names of data to actual vectors.
     */
    vector_map named_vectors;
    
    /**
     * A map linking names of data to actual std vectors.
     */
    std_vector_map named_std_vectors;
    
    /**
     * A map linking names of data to actual std vectors of std vectors.
     */
    std_vector_std_vector_map named_std_vectors_of_std_vectors;    
    
    
    
protected:
    
    FuelCell::ApplicationCore::LinearSolver lin_solver;
    
    FuelCell::ApplicationCore::NonLinearSolver nonlin_solver;
    
    FuelCell::ApplicationCore::RefinementSolver refinement_solver;
    
};

} // ApplicationCore

} // FuelCell

#endif