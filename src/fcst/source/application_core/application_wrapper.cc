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
// - Class: application_wrapper.cc
// - Description: This class implements either iterative or time-stepping
//                wrapper of applications
// - Developers: Guido Kanschat,     Texas A&M University
//               Valentin N. Zingan, University of Alberta
//               Marc Secanell,      University of Alberta
//
// ----------------------------------------------------------------------------

#include <application_core/application_wrapper.h>

namespace NAME = FuelCell::ApplicationCore;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NAME::ApplicationWrapper::ApplicationWrapper(NAME::ApplicationBase& app)
:
NAME::ApplicationBase(app.get_data()),
app(&app)
{
  FcstUtilities::log << "->Copy";
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NAME::ApplicationWrapper::~ApplicationWrapper()
{ }

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationWrapper::declare_parameters(ParameterHandler& param)
{
  app->declare_parameters(param);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationWrapper::initialize(ParameterHandler& param)
{
  app->initialize(param);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationWrapper::remesh()
{
  app->remesh();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationWrapper::init_vector(NAME::FEVector& dst) const
{
  app->init_vector(dst);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double
NAME::ApplicationWrapper::residual(NAME::FEVector&        dst,
                                   const NAME::FEVectors& src,
                                   bool                   apply_boundaries)
{
  return app->residual(dst,
                       src,
                       apply_boundaries);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationWrapper::solve(NAME::FEVector&        dst,
                                const NAME::FEVectors& src)
{

  //Prepare the residual
  NAME::FEVector res;
  res.reinit(dst);
  residual(res, src, false);

  NAME::FEVectors srcAndRes;
  srcAndRes.add_vector(res, "residual");
  srcAndRes.merge(src);

  //Notify the application to assemble the matrix
  notify(Event::assign("LinearAssembly"));

  app->solve(dst,
             srcAndRes);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationWrapper::Tsolve(NAME::FEVector&        dst,
                                 const NAME::FEVectors& src)
{
  app->Tsolve(dst,
              src);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double
NAME::ApplicationWrapper::estimate(const NAME::FEVectors& src)
{
  return app->estimate(src);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double
NAME::ApplicationWrapper::evaluate(const NAME::FEVectors& src)
{
  return app->evaluate(src);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationWrapper::grid_out(const std::string& filename) const
{
  app->grid_out(filename);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationWrapper::data_out(const std::string&     filename,
                                   const NAME::FEVectors& src)
{
  app->data_out(filename,
                src);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

std::string
NAME::ApplicationWrapper::id() const
{
  return std::string( typeid(*this).name() );
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationWrapper::notify(const NAME::Event& reason)
{
  NAME::ApplicationBase::notify(reason);
  app->notify(reason);
}
