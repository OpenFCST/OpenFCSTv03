//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: dakota_direct_interface.h
//    - Description: Class used to create an interface between FCST and DAKOTA
//    - Developers: Malte Krack, Peter Dobson, Marc Secanell <secanell@ualberta.ca>
//
//---------------------------------------------------------------------------

#ifndef dakota_direct_interface_h
#define dakota_direct_interface_h

//STL libraries:
#include<string>
#include<fstream>
#include<iostream>
#include<vector>
#include<ctime>

//deal.II:
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/timer.h>

// Fuel Cells:
#include <application_core/optimization_block_matrix_application.h>
#include <contribs/experimental_data.h>
#include <utils/simulation_selector.h>
#include <solvers/adaptive_refinement.h>
#include <contribs/dakota_application.h>
#include <utils/fcst_utilities.h>
#include <application_core/application_data.h>


#ifdef _WITH_DAKOTA

#include <DirectApplicInterface.hpp>
#include <ApplicationInterface.hpp>

/**
 * \brief Namespace under Dakota used for interfacing simulations
 *
 * Classes inside this namespace are used to interface the fuel cell analysis code
 * with DAKOTA (an optimization toolbox).
 */
namespace SIM
{

    /**
     * This class is used to read the input file from Dakota in its original format, use
     * this inforamtion to launch the fuel cell simulator and finally from the data from the
     * fuel cell simulator, write the output to Dakota.
     *
     * This application needs to be inherited from Dakota::DirectApplicInterface which belongs to the
     * Dakota library.
     *
     * These applications are selected in SimulatorBuilder<dim>::run_optimization() and they are assigned to
     * Dakota in  DakotaApplication::assign_interface(Dakota::DirectApplicInterface* optimization_interface)
     *
     * For more information please visit the Dakota reference guide at http://dakota.sandia.gov/licensing/votd/html-dev/
     *
     * @author M. Krack, P. Dobson and M. Secanell, 2010-12
     */
    template <int dim>
    class DakotaDirectInterface : public Dakota::DirectApplicInterface
    {
    public:

        /**
         * Default constructor for an object of this class.
         */
        //DakotaDirectInterface(boost::shared_ptr<const Dakota::ProblemDescDB > problem_db);

        /**
         * Constructor for an object of this class.
         * Links the Dakota database with the linear applications
         * and the solver.
         */
        DakotaDirectInterface(DakotaApplication &fcst_interface,
                              boost::shared_ptr<const Dakota::ProblemDescDB > problem_db,
                              ParameterHandler& param,
                              boost::shared_ptr <FuelCell::ApplicationCore::ApplicationData> data,
                              boost::shared_ptr<SimulationSelector<dim> > sim_selector,
                              std::string &param_file);

        /** Destructor */
        ~DakotaDirectInterface()
        {
        };


    protected:
        /**
         * Function called by DAKOTA's iterators in order to run the simulation code and
         * assign the corresponding result values by using synchronize_variables().
         */
        int derived_map_ac(const Dakota::String& ac_name);

        /**
         * Member function that is used to get the design variable names and values
         * and assigns the values to the parameter handler object
         */
        void dakota_adopt_parameters(boost::shared_ptr<FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim> > app_linear);

        /**
         * Member function that is used to assign the application responses and gradients
         * to the result variables of DAKOTA.
         */
        void dakota_assign_results(const std::vector<double>& responses,
                                   const std::vector<std::vector<double> >& dresponses_dl);
        /**
         * Member function that is used to assign the application responses to the
         * result variables of DAKOTA
         */
        void dakota_assign_results(const std::vector<double>& responses);

        /**
         * Pointer to the Parameter Handler object.
         */
        ParameterHandler* param;

        /**
         * Pointer to the SimulationSelector object which gives access to the file names
         * and selector functions for linear applications (CATHODE MODEL, NON-ISOTHERMAL CATHODE MODEL,
         * MEA Model,...), newton solves (NewtonBasic, Newton3ppC, ...), and the meshing method that you
         * would like to use (Adaptive Refinement).
         */
        boost::shared_ptr<SimulationSelector<deal_II_dimension> > sim_selector;

        /**
         * Pointer to the DakotaApplication object, so that the design variables can be accessed with
         * synchronize_variables().
         */
        DakotaApplication* fcst_interface;

        /**
         * Pointer to the problem description data base
         */
        boost::shared_ptr<const Dakota::ProblemDescDB> problem_db;

        /**
         *
         */
        std::string simulator_parameter_file_name;

        /**
         * Object used to calculate the CPU and Run time for the optimization evaluation
         */
        Timer eval_timer;

        /**
         * Number of design variables
         */
        unsigned int n_dvar;

        /**
         * Number of responses
         */
        unsigned int n_resp;

        /**
         * Member that stores the name of the design variables
         */
        std::vector<std::string> name_design_var;

        /**
         * Member that stores the name of the responses
         */
        std::vector<std::string> name_responses;
        
        /** 
         * Data structure storing information to be shared between applications 
         */
        boost::shared_ptr <FuelCell::ApplicationCore::ApplicationData> data;

    };


    // Least-Squares Interface
    /**
     * This class is used to solve a least-squares parameter estimation problem.
     * Data is read in from file and the solve function loops over the data points and
     * obtains a solution from the fuel cell simulator.  The results are returned to DAKOTA
     * for the optimization loop.
     */
    template <int dim>
    class DakotaLeastSquaresInterface : public DakotaDirectInterface<dim>
    {
    public:
        /**
         * Constructor
         */
        DakotaLeastSquaresInterface(DakotaApplication &fcst_interface,
                                    boost::shared_ptr<const Dakota::ProblemDescDB > problem_db,
                                    ParameterHandler& param,
                                    boost::shared_ptr <FuelCell::ApplicationCore::ApplicationData> data,
                                    boost::shared_ptr<SimulationSelector<dim> > sim_selector,
                                    std::string &param_file);

        /**
         * Destructor
         */
        ~DakotaLeastSquaresInterface()
        {
        };

        /**
         * Local initialize function for this class
         */
        void _initialize(ParameterHandler& param);

    protected:
        /**
         * Function called by DAKOTA's iterators in order to run the simulation code and
         * assign the corresponding result values by using synchronize_variables().
         */
        int derived_map_ac(const Dakota::String& ac_name);

    private:

        /**
         * Experimental data file for NLS parameter estimation - Operating Conditions & output
         */
        std::string NLS_data_file;
        /**
         * Output vector storing residuals from experimental data and output.
         * Used for NLS optimization iteration
         */
        std::vector<double> NLS_residual;
        /**
         * String storing the option of absolute, weighted, or no residuals (for sensitivity)
         */
        std::string NLS_residual_option;
        /**
         * String storing the choice of returning the residual vector (for NLS methods) or the L2 norm (for all other methods)
         */
        std::string NLS_residual_value;
        /**
         * Vector to store experimental current density output values
         */
        std::vector<double> experimental_current;
        /**
         * Vector to store Operating Condition names (set in NLS_data_file)
         */
        std::vector<std::string> OC_names;
        /**
         * Array to store Operating Condition values (set in NLS_data_file)
         */
        std::vector<std::vector<double> > OC_values;
        /**
         * Vector to store the parameter names to set
         */
        std::vector<std::string> param_names;
        /**
         * Vector to store the parameter values to set
         */
        std::vector<double> param_values;

        /** Object used to calculate the CPU and Run time for the NLS loop */
        Timer nls_timer;
        
        /** 
         * Data structure storing information to be shared between applications 
         */
        boost::shared_ptr <FuelCell::ApplicationCore::ApplicationData> data;
    };

}

#endif  //Dakota


#endif
