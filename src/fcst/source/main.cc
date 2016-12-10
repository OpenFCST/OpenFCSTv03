//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: Main program. This program is used select the appropriate application that we want to run.
//    - Description: Implementation code for the Iterator class
//    - Developers: M. Secanell, Peter Dobson and Michael Moore
//    - Id: $Id: main.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include "simulator_builder.h"

int main (int argc, char *argv[])
{
    int exit_code = 0;
    
    try
    {
        // Initialize MPI
        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
        {
            SimulatorBuilder<deal_II_dimension> sb;
            sb.parse_inputs(argc, argv);
            sb.scan();
            sb.run();
        }
    }
    catch (std::exception &exc)
    {
        std::cerr << std::endl << std::endl
        << "----------------------------------------------------"<< std::endl;
        std::cerr << "Exception on processing: " << std::endl
        << exc.what() << std::endl
        << "Aborting!" << std::endl
        << "----------------------------------------------------"<< std::endl;
        exit_code = -1;
    }
    catch (...)
    {
        std::cerr << std::endl << std::endl
        << "----------------------------------------------------"<< std::endl;
        std::cerr << "Unknown exception!" << std::endl
        << "Aborting!" << std::endl
        << "----------------------------------------------------"<< std::endl;
        exit_code = -1;
    }
    
    return exit_code;  
    
}
