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
// - Class: event.h
// - Description: This class implements notifications of applications
// - Developers: Guido Kanschat,     Texas A&M University
//               Valentin N. Zingan, University of Alberta
//               Marc Secanell,      University of Alberta
//
// ----------------------------------------------------------------------------

#ifndef _FUEL_CELL_APPLICATION_CORE_EVENT_H_
#define _FUEL_CELL_APPLICATION_CORE_EVENT_H_

#include <algorithm>
#include <fstream>

namespace FuelCell
{
namespace ApplicationCore
{

/**
 * Objects of this kind are used to notify interior applications of
 * changes provoked by an outer loop. They are handed to the
 * application through ApplicationBase::notify() and it is up to the
 * actual application how to handle them.
 *
 * Event is organized as an extensible binary enumerator. Every class
 * can add its own events using assign(). A typical code example is
 *
 * <pre>
 * class A
 * {
 *   static Event event;
 * };
 *
 * Event A::event = Event::assign("Event for A");
 * </pre>
 *
 * @author Guido Kanschat
 */

class Event
{
public:

  /**
   * Constructor, generating a
   * clear Event object with all flags set to <tt>false</tt>.
   */
  Event();

  /**
   * This function registers a
   * new event type and assigns a
   * unique identifier to it.
   * All flags are set to <tt>false</tt>.
   * The flag for this new event type
   * is set to <tt>true</tt>.
   * The result of this function
   * should be stored for later
   * use.
   */
  static Event assign(const char* name);

  /**
   * Set all flags to <tt>true</tt>.
   */
  void all();

  /**
   * Set all flags to <tt>false</tt>.
   */
  void clear();

  /**
   * Test whether all the flags
   * set to <tt>true</tt> in the other Event object are
   * also set to <tt>true</tt> in this one.
   */
  bool test(const Event& other) const;

  /**
   * Return <tt>true</tt> if any
   * event is set to <tt>true</tt>.
   */
  bool any() const;

  /**
   * Add the flags of the other Event object.
   */
  Event& operator += (const Event& other);

  /**
   * Clear the flags of the other Event object.
   */
  Event& operator -= (const Event& other);

  /**
   * List the flags to a stream.
   */
  template<typename OS>
  void print(OS& os) const;

  /**
   * List the names to a stream.
   */
  template<typename OS>
  static void print_assigned(OS& os);

private:

  /**
   * Sometimes, actions have to
   * be taken by all
   * means. Therefore, if this
   * value is true, test() always
   * returns true.
   */
  bool all_true;

  /**
   * The actual list of events.
   */
  std::vector<bool> flags;

  /**
   * The actual list of names of events.
   */
  static std::vector<std::string> names;
};

} // ApplicationCore

} // FuelCell

#endif