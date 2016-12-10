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
// - Class: application_base.cc
// - Description: This class implements the interface of applications
// - Developers: Guido Kanschat,     Texas A&M University
//               Valentin N. Zingan, University of Alberta
//               Marc Secanell,      University of Alberta
// - Id: $Id: application_base.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <application_core/application_base.h>

namespace NAME = FuelCell::ApplicationCore;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NAME::ApplicationBase::ApplicationBase(boost::shared_ptr<NAME::ApplicationData> data)
:
Subscriptor(),
data(data)
#ifdef OPENFCST_WITH_PETSC
,mpi_communicator (MPI_COMM_WORLD),
n_mpi_processes (Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
this_mpi_process (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
#endif
{
    FcstUtilities::log << "Application";
    notifications.all();
    
    if(!this->data)
        this->data = boost::shared_ptr<NAME::ApplicationData>(new NAME::ApplicationData);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NAME::ApplicationBase::ApplicationBase(const NAME::ApplicationBase& other)
:
Subscriptor(),
data(other.data)
#ifdef OPENFCST_WITH_PETSC
,mpi_communicator (MPI_COMM_WORLD),
n_mpi_processes (Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
this_mpi_process (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
#endif
{
    FcstUtilities::log << "Base->";
    notifications.all();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NAME::ApplicationBase::~ApplicationBase()
{ }

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

std::string
NAME::ApplicationBase::id() const
{
    return std::string( typeid(*this).name() );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationBase::notify(const Event& reason)
{
    notifications += reason;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationBase::clear()
{
    notifications.all();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationBase::clear_events()
{
    notifications.clear();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned int
NAME::ApplicationBase::get_solution_index()
{
    const std::type_info& info = typeid(*this);
    FcstUtilities::log << "Pure function get_solution_index called in Class " << info.name() << std::endl;
    
    return 100;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

boost::shared_ptr<NAME::ApplicationData>
NAME::ApplicationBase::get_data()
{
    return data;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const boost::shared_ptr<NAME::ApplicationData>
NAME::ApplicationBase::get_data() const
{
    return data;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationBase::print_parameters_to_file(ParameterHandler&                    param,
                                                const std::string&                   filename,
                                                const ParameterHandler::OutputStyle& style)
{
    std::ofstream file_out;
    file_out.open(filename.c_str());
    FcstUtilities::log << "Printing parameter handler" << std::endl;
    param.print_parameters(file_out,
                           style);
    file_out.close();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationBase::print_caller_name(const std::string& caller_name) const
{
    const std::type_info& info = typeid(*this);
    FcstUtilities::log << "Pure function " << caller_name << " called in Class " << info.name() << std::endl;
}