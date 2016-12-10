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
// - Class: event.cc
// - Description: This class implements notifications of applications
// - Developers: Guido Kanschat,     Texas A&M University
//               Valentin N. Zingan, University of Alberta
//               Marc Secanell,      University of Alberta
// - Id: $Id: event.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <application_core/event.h>

namespace NAME = FuelCell::ApplicationCore;

std::vector<std::string> NAME::Event::names;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NAME::Event::Event()
:
all_true(false),
flags(names.size(), false)
{ }

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NAME::Event
NAME::Event::assign(const char* name)
{
  unsigned int index = names.size();
  names.push_back(name);

  NAME::Event result;
  result.flags[index] = true;

  return result;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::Event::all()
{
  all_true = true;
  std::fill(flags.begin(), flags.end(), true);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::Event::clear()
{
  all_true = false;
  std::fill(flags.begin(), flags.end(), false);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool
NAME::Event::test(const NAME::Event& other) const
{
  // First, test all_true in this object
  if(all_true)
    return true;

  const unsigned int n     = flags.size();
  const unsigned int m     = other.flags.size();
  const unsigned int n_min = std::min(n, m);

  // Now, if all_true set in the
  // other object, then all must be true
  // in this object
  if(other.all_true)
  {
    // Non existing flags are
    // always assumed false
    if(m > n)
      return false;

    // Test all flags separately
    // and return false if one is
    // not set
    for(std::vector<bool>::const_iterator iter  = flags.begin();
                                          iter != flags.end();
                                        ++iter)
      if(!*iter)
        return false;
    // All flags are set
    return true;
  }

  // Finally, compare each flag
  // separately
  for(unsigned int i = 0; i < n_min; ++i)
    if(other.flags[i] && !flags[i])
      return false;
  for(unsigned int i = n_min; i < m; ++i)
    if(other.flags[i])
      return false;
  return true;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool
NAME::Event::any() const
{
  if(all_true)
    return true;
  for(std::vector<bool>::const_iterator iter  = flags.begin();
                                        iter != flags.end();
                                      ++iter)
    if(*iter)
      return true;
  return false;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NAME::Event&
NAME::Event::operator += (const NAME::Event& other)
{
  all_true |= other.all_true;
  if(all_true)
    return *this;

  if(flags.size() < other.flags.size())
    flags.resize(other.flags.size());
  for(unsigned int i = 0; i < flags.size(); ++i)
    flags[i] = flags[i] || other.flags[i];

  return *this;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NAME::Event&
NAME::Event::operator -= (const NAME::Event& other)
{
  if(!other.any())
    return *this;

  all_true = false;
  if(other.all_true)
  {
    for(std::vector<bool>::iterator iter  = flags.begin();
                                    iter != flags.end();
                                  ++iter)
      *iter = false;
    return *this;
  }

  if(flags.size() < other.flags.size())
    flags.resize(other.flags.size());
  for(unsigned int i = 0; i < flags.size(); ++i)
    if(other.flags[i])
      flags[i] = false;

  return *this;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename OS>
void
NAME::Event::print(OS& os) const
{
  if(all_true)
    os << " ALL";

  for(unsigned int i = 0; i < flags.size(); ++i)
    if(flags[i])
      os << ' ' << names[i];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename OS>
void
NAME::Event::print_assigned(OS& os)
{
  for(unsigned int i = 0; i < names.size(); ++i)
    os << i << '\t' << names[i] << std::endl;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename OS>
OS& operator << (OS&                os,
                 const NAME::Event& event)
{
  event.print(os);
  return os;
}