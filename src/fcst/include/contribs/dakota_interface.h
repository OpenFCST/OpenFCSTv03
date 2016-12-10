//---------------------------------------------------------------------------
//    $Id: dakota_interface.h 2605 2014-08-15 03:36:44Z secanell $
//
//    Copyright (C) 2006 by Marc Secanell
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef dakota_interface_h
#define dakota_interface_h

//deal.II:
#include <deal.II/base/parameter_handler.h>

// Fuel Cells:
#include <application_core/optimization_block_matrix_application.h>
#include <utils/fcst_utilities.h>

//STL libraries:
#include<string>
#include<fstream>
#include<iostream>
#include<vector>

namespace SIM
{
     /**
      * Classes used to interface the fuel cell analysis code with DAKOTA (an optimization toolbox).
      * 
      * This class reads the input file from Dakota in its original format, use
      * this inforamtion to launch the fuel cell simulator and finally from the data from the
      * fuel cell simulator, write the output to Dakota.
      */
     template <int dim>
     class DakotaInterface
     {
     public:
          /**
           * Constructor for an object of this class. In order to be able to read and
           * write the objects needs to know where the Dakota files are. Furthermore,
           * the design variables modify the analysis parameters, so the input file location
           * is also necessary
           */
          DakotaInterface(const std::string input_file,
          				ParameterHandler& param,
						const std::string dakota_parameters,
						const std::string dakota_results,
						FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>& app_lin,
						FuelCell::ApplicationCore::ApplicationWrapper& app);

          /** Destructor */
          ~DakotaInterface()
          {
          };

          /**
           *
           */
          void run();

     private:
          /**
           * Declare all parameters that are needed for:
           *   - the computation of the equation coefficients
           *   - the control of the linear system solution
           *   - ...
           */
          void declare_parameters(ParameterHandler& param);

          /**
           * Set up how many equations are needed and
           * read in parameters for the parameter handler in order to initialize data
           */
          void initialize(ParameterHandler& param);

          /**
           * Member function that is used to read the file from Dakota and perform the necessary changes on the
           * ParameterHandler file
           */
          void DakotaReadIn(const std::string dakota_parameters,
                            FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>& app,
                            ParameterHandler &param);

          /**
           * Member function that is used to write the Dakota output file
           */
          void DakotaWriteOut(const std::string dakota_results,
                              const std::vector<double>& responses,
                              const std::vector<std::vector<double> >& dresponses_dl);

          /**
           * Gradients necessary?
           */
          bool gradients;
          /**
           * Number of refinements
           */
          unsigned int n_ref;

          /**
           * Dakota version used
           */
          std::string dakota_version;

          /**
           * Name of the analysis file
           */
          const std::string input_file;

          /**
           * Name of the parameters file from DAKOTA
           */
          const std::string dakota_parameters;

          /**
           * Name of the results file to DAKOTA
           */
          const std::string dakota_results;

          /**
           *
           */
          std::vector<int> ASV;

          /** 
           * Pointer to application
           */
          FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim> *app_linear;

          /**
           * Pointer to nonlinear application
           */
          FuelCell::ApplicationCore::ApplicationWrapper *app;
		
		/** Pointer to parameter handler object */
		ParameterHandler *param;
     };

}

#endif
