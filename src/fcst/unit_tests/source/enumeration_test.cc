// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: events_test.cc
// - Description: Test for events class
// - Developers: Marc Secanell Gallart,    University of Alberta
// - $Id: enumeration_test.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include <enumeration_test.h>



 void EnumerationTest::testDefault(){
     //The default "uninitialized" enumerations should be nothing, no reaction

     VariableNames varNothing;
     ReactionNames varNoReaction;

     //Both the following should work equivalently
     TEST_ASSERT(varNothing == VariableNames::nothing);
     TEST_ASSERT(varNothing == nothing);

     //Both the following should work equivalently
     TEST_ASSERT(varNoReaction== ReactionNames::noReaction);
     TEST_ASSERT(varNoReaction == noReaction);


     //Final check
     TEST_ASSERT(varNoReaction != HOR);
     TEST_ASSERT(varNothing != electronic_electrical_potential);

 }



 void EnumerationTest::testComparisons(){

     ReactionNames reaction0 = noReaction;
     ReactionNames reaction1 = HOR;
     ReactionNames reaction2 = ORR;


     TEST_ASSERT(reaction0 != HOR);
     TEST_ASSERT(reaction1 != reaction0);
     TEST_ASSERT(reaction2 != reaction1);
     TEST_ASSERT(reaction2 != reaction0);

     reaction2 = HOR;

     TEST_ASSERT(reaction2 == reaction1);

 }
