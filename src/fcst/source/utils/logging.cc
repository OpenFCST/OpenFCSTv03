//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Description: Logging objects to be used throughout fcst
//    - Developers: Philip Wardlaw
//
//---------------------------------------------------------------------------

#include "utils/logging.h"

namespace FcstUtilities
{
    // Delcaration of global object log used to ouput data to screen as well as to record data to file.
    FCSTLogStream log(std::cout);
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
FcstUtilities::FCSTLogStream::FCSTLogStream (std::ostream &stream)
:
output_stream (stream),
file(0)
{}

//---------------------------------------------------------------------------
FcstUtilities::FCSTLogStream::~FCSTLogStream()
{        };

//---------------------------------------------------------------------------
void
FcstUtilities::FCSTLogStream::attach (std::ostream &o,
                                      const bool    print_job_id)
{
   // #ifdef OPENFCST_WITH_PETSC
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
   // #endif
    {
        file = &o;            
        o.setf(std::ios::showpoint | std::ios::left);
    }
}

//---------------------------------------------------------------------------
void
FcstUtilities::FCSTLogStream::detach ()
{
    file = 0;    
}

//---------------------------------------------------------------------------
