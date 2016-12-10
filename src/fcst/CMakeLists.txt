#---------------------------------------------------------------------------
#
#    FCST: Fuel Cell Simulation Toolbox
#
#    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
#
#    This software is distributed under the MIT License.
#    For more information, see the README file in /doc/LICENSE
#
#---------------------------------------------------------------------------

####################### 
# 
# OpenFCST CMake File
#
# This script sets up the necessary variables
# and calls the other folders to do the install
# work.
#
# Developer: Chad Balen, 2014
#
#######################

# NOTE: you may notice unlike the other CMake projects there is no
#PROJECT() command this is because we need to call in the source folder
#AFTER we have set all the compiler options.

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
enable_language(Fortran) #Since COLDAE uses Fortran need to enable this language

#Temporary variable while in transition phase to CMake
set(OPENFCST_INSTALL_CONFIG "" CACHE TYPE PATH) #Dummy variable for setting up the CMake configure file
set(OPENFCST_DIMENSIONS_CONFIG "2" CACHE TYPE STRING) #Dummy variable for setting up the CMake configure file
set(OPENFCST_C11_FLAG "-std=c++11" CACHE TYPE STRING)

# Set variables necessary for installation
# Some are identical to External Super Build CMake since they need to be
# passed to OpenFCST
SET(OPENFCST_DIMENSIONS "2" CACHE TYPE STRING)
SET(PROJECT_NAME "OpenFCST_${OPENFCST_DIMENSIONS}D") #Create Project Name to be used to create project in source folder
SET(FCST_EXE "fuel_cell-${OPENFCST_DIMENSIONS}d.bin") #Name of executable
SET(FCST_GUI_EXE "fcst_gui")
SET(FCST_LIB "fuel_cell-${OPENFCST_DIMENSIONS}d") #Name of library file

set(P4EST_DIR "" CACHE TYPE PATH)
set(METIS_DIR "" CACHE TYPE PATH)
set(PETSC_DIR "" CACHE TYPE PATH)
set(TRILINOS_DIR "" CACHE TYPE PATH)
set(DEALII_DIR "" CACHE TYPE PATH)
set(DAKOTA_DIR "" CACHE TYPE PATH)
set(ALGLIB_DIR "" CACHE TYPE PATH)
set(COLDAE_DIR "" CACHE TYPE PATH)
set(SQLITE_DIR "" CACHE TYPE PATH)
set(CPPTEST_DIR "" CACHE TYPE PATH)

#Variable to determine if CMake will create a library file
#I don't want it to always create library file, b/c it double the time
#required to install this external project
set(OPENFCST_LIBRARY "OFF" CACHE TYPE BOOL)

#Create variables for the path ways to the header files so we can easily
#point to when we are in the source folder
set(HEADER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(UNIT_TESTS_HEADER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/unit_tests/include")

set(OPENFCST_WITH_DAKOTA "OFF" CACHE TYPE BOOL)
set(OPENFCST_WITH_PETSC "OFF" CACHE TYPE BOOL)
set(OPENFCST_WITH_TRILINOS "OFF" CACHE TYPE BOOL)

# NOTE: For some reason when you use an CMake's external project module to install
#OpenFCST it doesn't link the executable and library files to deal.II, Dakota, etc.
#So I copied these lines from http://www.cmake.org/Wiki/CMake_RPATH_handling
#Where they suggest to add lines below to prevent this from occuring
 
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

#Have CMake look into the following directories for CMakeLists.txt
add_subdirectory(source)
add_subdirectory(include)
add_subdirectory(unit_tests)
add_subdirectory(GUI)