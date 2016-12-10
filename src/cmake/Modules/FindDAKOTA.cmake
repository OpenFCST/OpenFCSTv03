## ---------------------------------------------------------------------
## $Id: FindDAKOTA.cmake 32021 2013-12-15 16:51:05Z maier $
##
## Copyright (C) 2012 - 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Try to find the DAKOTA library
#
# This module exports:
#
#   DAKOTA_DIR
#   DAKOTA_INCLUDE_DIR
#   DAKOTA_LIBRARY_DIR
#
#   If Dakota gets their shit together maybe we can set these variables as well:
#   DAKOTA_VERSION
#   DAKOTA_VERSION_MAJOR
#   DAKOTA_VERSION_MINOR
#   DAKOTA_VERSION_SUBMINOR
#   DAKOTA_WITH_MPI
#


#So we can use the FIND_PACKAGE_HANDLE_STANDARD_ARGS later to tell
#CMake if the package was found
INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(DAKOTA_DIR "$ENV{DAKOTA_DIR}")

find_path(DAKOTA_LIBRARY_DIR
  NAMES libdakota_src libdakota_src.so libdakota_src.a
  HINTS
    # petsc is special. Account for that
    ${DAKOTA_DIR}
    ${DAKOTA_LIBRARY_HINT}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib lib64 Install/lib Install/lib64 install/lib install/lib64
  NO_DEFAULT_PATH
)

if(NOT DAKOTA_LIBRARY_DIR STREQUAL "DAKOTA_LIBRARY_DIR-NOTFOUND")
  file(GLOB DAKOTA_SHARED_LIBRARY_FILES "${DAKOTA_LIBRARY_DIR}/*.so")
  file(GLOB DAKOTA_STATIC_LIBRARY_FILES "${DAKOTA_LIBRARY_DIR}/*.a")
endif()
# MESSAGE(STATUS
#   "CMake Cannot Find Dakota library, please confirm you have installed Dakota. If so please provide additional hints through CMake variable(s):\n"
#   "-DDAKOTA_DIR=...\n"
#   "-DDAKOTA_LIBRARY_DIR=...\n")

#Lets go find the path to the folder containing header files
#Note we do not need to find each individual subdirectory folder in include folder
#(Thank goodness) this is because each header file sets the include paths to the main
#folder which we are trying to find.
find_path(DAKOTA_INCLUDE_DIR DakotaStrategy.hpp
  HINTS
    ${DAKOTA_DIR}
    ${DAKOTA_INCLUDE_HINT}
  PATH_SUFFIXES dakota include Install/include install/include
  NO_DEFAULT_PATH
)


# This function is intended to be used in FindXXX.cmake modules files. It handles the 
# REQUIRED, QUIET and version-related arguments to find_package(). It also sets the 
# <packagename>_FOUND variable.
FIND_PACKAGE_HANDLE_STANDARD_ARGS(DAKOTA DEFAULT_MSG
  DAKOTA_LIBRARY_DIR
  DAKOTA_INCLUDE_DIR
  )