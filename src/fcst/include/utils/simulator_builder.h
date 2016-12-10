//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: simulation_builder.cc
//    - Description:
//    - Developers: M. Secanell, P. Dobson, Valentin N. Zingan
//
//---------------------------------------------------------------------------

#ifndef SIMULATOR_BUILDER_H
#define SIMULATOR_BUILDER_H

// STL
#include <fstream>
#include <vector>
#include <ctime>
#include <stdlib.h>

// Boost
#include <boost/smart_ptr.hpp>
#include <boost/program_options.hpp>

// deal.II
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

// FCST
#include <utils/polarization_curve.h>
#include <utils/parametric_study.h>
#include <grid/geometries.h>
#include <grid/geometry.h>
#include <utils/simulation_selector.h>
#include <microscale/agglomerate_ionomer_1D.h>
#include <microscale/agglomerate_water_1D.h>
#include <reactions/tafel_kinetics.h>
#include "FCST_TEST_SUITE.h"
#include <utils/fcst_utilities.h>

#include "contribs/dakota_interface.h"
// These files can only be used if DAKOTA is linked to our package:
#ifdef _WITH_DAKOTA
#include "contribs/dakota_direct_interface.h"
#include "contribs/dakota_application.h"
#endif

//---------------------------------------------------------
/**
 * 
 */
enum FileConversionOption
{
    NONE = 0, 
    XML2PRM, 
    PRM2XML
};

/**
 * This class is used to output data or to initialize and launch simulations.
 * It reads from the command line either an input file or any other flag and 
 * acts accordingly. 
 * Based on the flags it will print instructions to screen or setup a simulation.
 *
 * @note This is the only class called in main and it manages the problem to be solved.
 *
 * @author M. Secanell, P. Wardlaw and P. Dobson.
 *
 * and the concept of sequential applications and simplified parameter study without Dakota by
 * @author Valentin N. Zingan
 */
template <int dim>
class SimulatorBuilder
{
public:

    ///@name Constructors, destructor, and initalization
    //@{
    /**
     * Constructor
     */
    SimulatorBuilder();

    /**
     * Destructor
     */
    ~SimulatorBuilder();

    /**
     * Member function used to read in the command line and initialize the necessary data
     */
    void parse_inputs(int argc, char *argv[]);

    /**
     * Member function used to call the declaration of the parameters, initialize data
     * and set the application and solver objects to those specified in the input file
     */
    void scan();
    //@}

    ///@name Execute
    //@{
    /**
     * Run the simulation.
     * Determines whether optimization is being used and runs a single simulation
     * or optimization routines
     */
    virtual void run();
    //@}

protected:
    ///@name Initialization (called by scan())
    //@{
    /**
     * Declare all necessary parameters to read the input files
     * 
     * The following parameters are set here:
     * - "Dakota direct": Set to true if you would like to solve an optimization problem using Dakota. 
     *  Transfer of information between the application and Dakota is via the class DakotaDirectInterface
     * - "simulator parameter file name": Specify the file that contains all the application 
     *  (See the application that you would like to run, e.g. AppCathode)
     * - "optimization parameter file name" : Name of the file that specifies all the information related to 
     *  the optimization class (See DakotaApplication for parameters)
     * - "Run tests" : Set to true if you would like to run a single function. This routine can be re-implemented 
     *  to test different objects without running a full simulation
     * - "save transfer files": If set to false the mesh and solution from a parametric study will be deleted.
     * 
     */
    void declare_parameters(ParameterHandler& param) const;

    /**
     * Initialize the local variables declared by the parameter handler
     */
    void initialize(ParameterHandler& param);

    /**
     * Member function used to attach the deal.II logfile to the application and provide the name of the file.
     * If no name is given, then the machine hostname is used.
     */
    void open_logfile (const std::string&);
    //@}
    ///@name Execution routines (called by run)
    //@{
    /**
     * Set up the required optimization objects (i.e.
     * - Dakota::ParallelLibrary: Is where you specify what type of Message Passing Interface (MPI) you would like to run. E.g. Running in "serial mode" single processor, or in "parallel" on multiple processors.
     * - Dakota::ProblemDescDB: Stores optimization problem description and optimization strategy
     * - DakotaApplication: Manages inputs and sends data to application via DakotaDirectInterface (See DakotaApplication for more detials)
     * - DirectApplicInterface: Communicates Dakota parameters to application, sets up application (manages input from Dakota to Application),
     * runs application and results from the Application to Dakota.
     * 
     * It must be a child of Dakota::DirectApplicInterface
     *
     *  http://dakota.sandia.gov/docs/dakota/5.2/html-dev/DakLibrary.html
     *
     * \todo3 Print results to file from direct interface?
     */
    void run_optimization();

    /**
     * This routine is mainly used for testing single member functions in classes.
     * It will run in isolation from the main simulation framework, so that you are not
     * solving a finite element problem.
     * */
    void run_test();
    
    ///@name XML converters:
    //@{
    /**
     * Outputs the SimulatorBuilder parameter options to file "main.xml" in current working directory.
     */
    void output_default_main();

    /**
     * Outputs the AdaptiveRefinement and DakotaApplication parameters options to files "data.xml" and "opt.xml" in current working directory.
     *
     * @note "opt.xml" is only produced if compiling with Dakota
     */
    void output_default_other(std::string main_file);

    /**
     * Opens project specified by main file @param main_file
     * and outputs corresponding "main" "data" and "opt" in current working directory.
     * 
     * Possible options for conversion are xml2prm and prm2xml.
     *
     * @note "opt.xml" is only produced if compiling with Dakota and if opt file is found
     */
    void convert_file(std::string main_file);
    //@}
    
    ///@name Print info
    //@{
    /**
     * Member function that is used to print information about the program. The info
     * in print logo appears everytime the program is executed.
     */
    void print_logo() const;

    /**
     * Print \p timer info.
     */
    void print_timer_info() const;
    //@}

    ///@name Constant parameters
    //@{
    /** Local varible for the program name */
    const std::string program_name;

    /** Local variable for the program version */
    const std::string program_version;
    //@}

    ///@name Objects storing data read from all data files, linear applications and ApplicationWrapper applications
    //@{
    /**
     * Parameter handler object that will be used to store all input
     * data for the application
     */
    ParameterHandler param;
    
    /** 
     * Data structure storing information to be shared between applications 
     */
    boost::shared_ptr <FuelCell::ApplicationCore::ApplicationData> data;

    /** Pointer where linear application is stored */
    boost::shared_ptr <FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim> >app_lin;

    /** Pointer to the non-linear solver */
    boost::shared_ptr <FuelCell::ApplicationCore::ApplicationWrapper> newton;

    /**
     * Pointer to a solver application which applies adaptive refinement to the grid
     */
    boost::shared_ptr <FuelCell::ApplicationCore::AdaptiveRefinement<dim> > solver;

    /**
     * Object which stores the name of the application and solver
     * Allows the user to select the application and set objects to the empty pointers
     */
    boost::shared_ptr< SimulationSelector<dim> > sim_selector;
    //@}
    ///@name File names and general options:
    //@{
    /** 
     * Stores the name of the input file with the application selector data.
     * 
     * This value is specified as the argument when you launch openFCST.
     */
    std::string input_file;

    /** 
     * Stores a pointer to the log file object.
     * 
     * This value is specified in the main.xml file under section
     * @code
     * subsection Log
     * 
     * end
     * @endcode
     */
    boost::shared_ptr<std::ofstream> log_file;

    /** 
     * Stores the name of the parameter file containing the physical data
     * for the simulation.  I.e. Fuel cell data.
     */
    std::string simulator_parameter_file_name;
    
    /**
     * Decision varible for selecting the type of analysis. 
     *
     * This variable is set in the input file in:
     * @code
     * subsection Simulator
     *    set Analysis type = 
     * end
     * @endcode
     *
     * Possible options are:
     * - Unit Tests
     * - Analysis
     * - Parametric Study
     * - Polarization Curve
     * - Optimization
     * where the options specify the type of study you would like to perform. The options are:
     * - Unit Tests: Performs internal tests to make sure openFCST is running correctly.
     * - Analysis: Performs an analysis run, i.e. uses the model in simulator name and the provided here to perform an analysis run.
     * - Parametric Study: Performs a parameteric study using same values as analysis run and modifying parameter specified in Parametric study.
     * - Polarization Curve: Computes a polarization curve using same values as analysis run.
     * - Optimization: Solves an optimization problem using the values in the analysis files and the data in Optimization study and within.
     */
    std::string analysis_type;
    /** 
     * Set to true if you want to keep the mesh and solution from the simulation. 
     */
    bool save_transfer_files;
    //@}
    
    ///@name Auxiliary data:
    //@{
    /** 
     * Object used to calculate the CPU and Run time for the simulation 
     */
    Timer timer;
    //@}
    
    ///@name Dakota dependent parameters
    //@{
    /**
     * Variable set in parse_inputs which determines whether the program is
     * being called from DAKOTA
     */
    bool dakota_use;

    /**
     * Decision varible for using the direct dakota interface
     */
    bool dakota_direct;

    /** Stores the name of the parameter file containing the optimization data
     */
    std::string optimization_parameter_file_name;

    /**
     * Variable set in parse_inputs which stores the name of the results (responses) file
     * if the simulation is being called from DAKOTA
     */
    std::string dakota_results;

    /**
     * Variable set in parse_inputs which stores the name of the parameter file
     * if the simulation is being called from DAKOTA
     */
    std::string dakota_parameters;
    //@}
    ///@name Polarization curve parameters:
    //@{
    /** Polarization curve object. */
    FuelCell::PolarizationCurve<dim> curve;
    //@}
    ///@name Parametric study section:
    //@{
    /** Polarization curve object. */
    FuelCell::ParametricStudy<dim> param_study;    
    //@}
};

#endif
