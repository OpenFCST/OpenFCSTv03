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
# This is the "Super Build" script that
# installs all external packages including
# OpenFCST
#
# Developer: Chad Balen, 2014
#
#######################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
set(OPENFCST_EXTRA_FLAGS "" CACHE TYPE STRING)
set(CXX_EXTRA_FLAGS "${CXX_EXTRA_FLAGS} ${OPENFCST_EXTRA_FLAGS}")
set(C_EXTRA_FLAGS "${C_EXTRA_FLAGS} ${OPENFCST_EXTRA_FLAGS}")
# WARNING: There's also OpenFCST Extra flags because in the case of the gcc 4.8 bug
#          I cannot send that flag to OpenFCST else segmentation faults occur with Dakota.
#          As a result, it is not recommened that developers use the C++11 thread support in
#          OpenFCST.
set(OPENFCST_CXX_EXTRA_FLAGS "${CXX_EXTRA_FLAGS} ${OPENFCST_EXTRA_FLAGS}")
set(OPENFCST_C_EXTRA_FLAGS "${C_EXTRA_FLAGS} ${OPENFCST_EXTRA_FLAGS}")
MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")

MESSAGE(STATUS
        "\n"
        "=======================================================\n"
        "\n"
        " This script will automatically configure and compile  \n"
        "     the FCST code and its contributing libraries.     \n"
        "\n"
        "     This process will take approx. 1 to 2 hours       \n"
        "\n"
        "     Authors: Chad Balen, Michael Moore,               \n"
        "        Philip Wardlaw, Andreas Putz, Madhur Bhaiya,   \n"
        "        and Marc Secanell                              \n"
        "\n"
        "=======================================================\n"
        "\n"
        )

#Check if user has installed Doxygen for OpenFCST documentation.
#If not give warning and suggest they may want to install Doxygen.
find_package(Doxygen)
if(NOT DOXYGEN_FOUND)
  MESSAGE(STATUS
          "\n"
          "=======================================================\n"
          "\n"
          " WARNING: Doxygen was NOT FOUND as a result OpenFCST   \n"
          "          will be unable to build its documentation.   \n"
          "          If you require OpenFCST documentation,       \n"
          "          please install Doxygen.                      \n"
          "\n"
          "=======================================================\n"
          "\n"
          )
endif()


#----------------------------
#
# Initiate CMake Modules and Macros
#
#----------------------------

#Set CMake Module Path, so CMake knows where to look to find external CMake Modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

#
# Load all macros:
#
#INCLUDE(setup_external_macros)
FILE(GLOB _macro_files "cmake/macros/*.cmake")
FOREACH(_file ${_macro_files})
  INCLUDE(${_file})
ENDFOREACH()


#----------------------------
#
# OpenFCST Find MPI
#
#----------------------------

#Create Hint for where to look incase CMake can't find its path, so user can
#specify its location.
set(MPI_DIR_HINT "" CACHE TYPE FILEPATH)

#Find OpenMPI C and CXX Compilers
find_program(OPENFCST_MPI_C_COMPILER
        NAMES mpicc
        HINTS ${MPI_DIR_HINT} ENV MPIHOME 
        PATH_SUFFIXES bin
        DOC "C Compiler to use with OpenFCST (mpicc preferred)."
        )

find_program(OPENFCST_MPI_CXX_COMPILER
        NAMES mpicxx
        HINTS ${MPI_DIR_HINT} ENV MPIHOME 
        PATH_SUFFIXES bin
        DOC "C++ Compiler to use with OpenFCST (mpicxx preferred)."
        )
        
find_program(OPENFCST_MPI_FORTRAN_COMPILER
        NAMES mpif90
        HINTS ${MPI_DIR_HINT} ENV MPIHOME 
        PATH_SUFFIXES bin
        DOC "Fortran Compiler to use with OpenFCST (mpif90 preferred)."
        )

#Check to see if OpenMPI compilers were found, if not terminate
#CMake and tell them to confirm they have installed OpenMPI and
#possibly give it's path as a hint to CMake
ErrMsgRegular(${OPENFCST_MPI_C_COMPILER} "OPENFCST_MPI_C_COMPILER-NOTFOUND" "mpicc" "MPI_DIR_HINT" "mpi-dir")
ErrMsgRegular(${OPENFCST_MPI_CXX_COMPILER} "OPENFCST_MPI_CXX_COMPILER-NOTFOUND" "mpicxx" "MPI_DIR_HINT" "mpi-dir")
ErrMsgRegular(${OPENFCST_MPI_FORTRAN_COMPILER} "OPENFCST_MPI_FORTRAN_COMPILER-NOTFOUND" "mpif90" "MPI_DIR_HINT" "mpi-dir")

#Hopefully this will fix the Secanell KDevelop Bug
get_filename_component(MPI_PATH ${OPENFCST_MPI_CXX_COMPILER} PATH)
set(LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:${MPI_PATH}")

#Set CMake to use the compilers found
set(CMAKE_C_COMPILER ${OPENFCST_MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${OPENFCST_MPI_CXX_COMPILER})
set(CMAKE_Fortran_COMPILER ${OPENFCST_MPI_FORTRAN_COMPILER})

#----------------------------
#
# Create OpenFCST Project
#
#----------------------------
set(PROJECT_NAME OpenFCST)
PROJECT(${PROJECT_NAME})
enable_language(Fortran)
include(ExternalProject)
include(CheckCXXCompilerFlag)
find_package(MPI) #Must come after PROJECT() else fails


#----------------------------
#
# Check for C++11
#
#----------------------------
#C++11 is mandatory to compile OpenFCST but their are two possible flags
#the compiler may use, so check for which one they have
CHECK_CXX_COMPILER_FLAG("-std=c++11" stdcxx11_flag)
if(stdcxx11_flag) #If stdcxx11_flag returns 1 (ie true)
  set(OPENFCST_C11_FLAG "-std=c++11")
else() #Check for option 2
  CHECK_CXX_COMPILER_FLAG("-std=c++0x" stdcxx0x_flag)
  if(stdcxx0x_flag) #If stdcxx0x_flag returns 1 (ie true)
    set(OPENFCST_C11_FLAG "-std=c++0x")
  else() #Neither compiler flag was found give fatal error
    MESSAGE(FATAL_ERROR 
              "\n"
              "=======================================================\n"
              "\n"
              " ERROR: CMake does not detect that your compiler is    \n"
              "        configured to use the new C++11 standard.      \n"
              "        This is a requirement for parts of OpenFCST    \n"
              "        If your compiler uses a different flag then    \n"
              "        -std=c++0x or -std=c++11 then please specify   \n"
              "        your compiler C++11 flag with the CMake        \n"
              "        variable:                                      \n"
              "            -DOPENFCST_C11_FLAG=...                    \n"
              "        Or if you are using the OpenFCST shell script  \n"
              "        use the flag:                                  \n"
              "            --cxx11-flag=...                           \n"
              "\n"
              "=======================================================\n"
              "\n"
              )
  endif()
endif()

#----------------------------
#
# Check Compiler Type and Version
#
#----------------------------
#I've never tested openFCST requires:
#  i) GNU
#  ii) >= 4.7 (for C++11 support)
#  iii) if = 4.8.X then we need to set some specific compiler flags due to bugs
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.7.0")
    ErrMsgVersionNCCMV("GNU" "Compiler" ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} "4.7" "4.7.2")
  endif()
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9.0" AND 
       CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "4.8.0")
    #In GNU 4.8.X there is a bug when using C++11 threads so need to set this flag

# WARNING: Cannot send flag below to OpenFCST Extra Flags because the gcc 4.8 bug
#          will cause segmentation faults if using Dakota. As a result, it is not 
#          recommened that developers use the C++11 thread support in OpenFCST.
    set(CXX_EXTRA_FLAGS "${CXX_EXTRA_FLAGS} -Wl,--no-as-needed")
  endif()
else()
  ErrMsgVersionNCCMV("GNU" "Compiler" ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} "4.7" "4.7.2")
endif()

#----------------------------
#
# Check OpenMP Compiler Flag
#
#----------------------------
set(OPENFCST_OPENMP "OFF" CACHE TYPE BOOL)
if(OPENFCST_OPENMP) #Only need to check if we want to use it
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(OPENFCST_CXX_EXTRA_FLAGS "${OpenMP_CXX_FLAGS}")
    set(OPENFCST_C_EXTRA_FLAGS "${OpenMP_C_FLAGS}")
  else()
    MESSAGE(FATAL_ERROR 
                "\n"
                "=======================================================\n"
                "\n"
                " ERROR: Some parts of OpenFCST require OpenMP, however \n"
                "        your compiler does not appear to support       \n"
                "        OpenMP. Check that your compiler supports and  \n"
                "        is compiled correctly. You can turn OpenMP off \n"
                "        by setting the CMake variable as follows:      \n"
                "            -DOPENFCST_OPENMP=OFF                      \n"
                "        Or if you are using the OpenFCST shell script  \n"
                "        then DO NOT specify the flag:                  \n"
                "            --with-openmp                              \n"
                "\n"
                "=======================================================\n"
                "\n"

                )
  endif()
endif()

#----------------------------
#
# Find Boost
#
#----------------------------
#CMake's FindBoost package is pretty good at finding the vendor installed
#Boost by machine. If you wish to use your own installation of boost, things
#get tricky. So we use BOOST_DIR as our hint to point to a boost directory
#other than vendor boost. If BOOST_DIR not set (ie "") then find vendor Boost
#else we change BOOST_ROOT and turn off Boost_NO_SYSTEM_PATHS ON so that way
#when it does find our non-vendor boost it DOESN'T USE IT!!!!!
set(BOOST_DIR "" CACHE TYPE PATH)
if(NOT BOOST_DIR STREQUAL "")
  set(BOOST_ROOT ${BOOST_DIR})
  set(Boost_NO_SYSTEM_PATHS ON) #So Boost doesn't use vendor boost instead of another boost
endif()

find_package(Boost
               COMPONENTS iostreams #For deal.ii
                 serialization #For deal.ii & Dakota
                 system #For deal.ii & Dakota
                 thread #For deal.ii
                 filesystem #For Dakota
                 regex #For Dakota
                 signals #For Dakota
                 program_options #For command line argument parsing
                 )

ErrMsgRegular(${Boost_FOUND} "0" "Boost" "BOOST_DIR" "boost-dir")

#We require a Boost version of atleast 1.49.0
if(Boost_MAJOR_VERSION LESS 1) 
  set(BOOST_MAJOR_VERSION_BAD "ON")
  set(BOOST_MINOR_VERSION_BAD "ON")
else()
  set(BOOST_MAJOR_VERSION_BAD "OFF")
  if(Boost_MINOR_VERSION LESS 49)
    set(BOOST_MINOR_VERSION_BAD "ON")
  else()
    set(BOOST_MINOR_VERSION_BAD "OFF")
  endif()
endif()


if(BOOST_MAJOR_VERSION_BAD OR BOOST_MINOR_VERSION_BAD)
  ErrMsgVersionCCMV("Boost" "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}" "1.49.0" "BOOST_DIR" "boost-dir")
endif()

#----------------------------
#
# Find BLAS & LAPACK
#
#----------------------------
find_package(BLAS)
ErrMsgRegularNoSol(${BLAS_FOUND} "TRUE" "BLAS")

find_package(LAPACK)
ErrMsgRegularNoSol(${LAPACK_FOUND} "TRUE" "LAPACK")

#If CMake cannot find it and user has a self-installed version then from the looks of it
#they have to hardcode the paths to CMake based on deal.ii documentation so then they'll
#have to alter our bash script to get it to work. But if they are using a self-compiled
#version of BLAS & LAPACK this should be easy for them.

# NOTE: the created by these packages are:
# i) BLAS_FOUND / LAPACK_FOUND - set to true if a library implementing the BLAS/LAPACK interface is found
# ii) BLAS_LINKER_FLAGS / LAPACK_LINKER_FLAGS - uncached list of required linker flags (excluding -l and -L).
# iii) BLAS_LIBRARIES / LAPACK_LIBRARIES - uncached list of libraries (using full path name) to link against to use BLAS95/LAPACK95
# iv) BLAS95_LIBRARIES / LAPACK95_LIBRARIES - uncached list of libraries (using full path name) to link against to use BLAS95/LAPACK95
# v) BLAS95_FOUND / LAPACK95_FOUND - set to true if a library implementing the BLAS/LAPACK f95 interface is found
# vi) BLA_STATIC - if syoet on this determines what kind of linkage we do (static)
# vii) BLA_VENDOR - if set checks only the specified vendor, if not set checks all the possibilities
# viii) BLA_F95 - if set on tries to find the f95 interfaces for BLAS/LAPACK

#----------------------------
#
# Find FLEX
#
#----------------------------
#So this is required for Dakota, but Dakota does not give any documentation on CMake variables
#to give for FLEX, and after looking through their CMake they made none. Since Dakota has
#it so you can only used a FLEX installed by system (not user installed) in the /usr/ directory
#we cannot set any CMake variables for help with finding FLEX. So we will just try
#to find FLEX if we can't find then neither will Dakota. Then an error is sent to the user
#to tell them to install it before they get to the Dakota installation.
find_package(FLEX)
ErrMsgRegularNoSol(${FLEX_FOUND} "TRUE" "FLEX")

#----------------------------
#
# Find SciPy & NumPy
#
#----------------------------
find_package(SciPy)
ErrMsgRegularNoSol(${SCIPY_FOUND} "TRUE" "SciPy")

find_package(NumPy)
ErrMsgRegularNoSol(${NUMPY_FOUND} "TRUE" "NumPy")

#----------------------------
#
# Find Bison
#
#----------------------------
find_package(BISON)
ErrMsgRegularNoSol(${BISON_FOUND} "TRUE" "BISON")

#----------------------------
#
# Find Qt Packages
#
#----------------------------
#Currently we only the Qt4-designer for GUI
find_package(Qt4 QUIET 
               COMPONENTS QtDesigner)

if(DEFINED Qt4_FOUND) #Most recent variable for checking
  set(qt_found "${Qt4_FOUND}")
elseif(DEFINED QT4_FOUND) #Compatible variable
  set(qt_found "${QT4_FOUND}")
elseif(DEFINED QT_FOUND) #Older Compatible variable
  set(qt_found "${QT_FOUND}")
else()
  MESSAGE(FATAL_ERROR
                  "\n"
                  "=======================================================\n"
                  "\n"
                  " ERROR: The variable used by CMake for determining     \n"
                  "        if Qt4 was found does not appear to exist      \n"
                  "\n"
                  "=======================================================\n"
                  "\n"
                  )
endif()
ErrMsgRegularNoSol(${qt_found} "TRUE" "Qt4-designer")

#----------------------------
#
# Find Sphinx
#
#----------------------------
find_package(Sphinx)
if(NOT SPHINX_FOUND)
  MESSAGE(STATUS
          "\n"
          "=======================================================\n"
          "\n"
          " WARNING: Sphinx was NOT FOUND as a result OpenFCST    \n"
          "          will be unable to build its documentation.   \n"
          "          If you require OpenFCST documentation,       \n"
          "          please install Sphinx.                       \n"
          "\n"
          "=======================================================\n"
          "\n"
          )
endif()

#----------------------------
#
# OpenFCST CMake Variables
#
#----------------------------
#Temporary variable while in transition phase to CMake
set(ONLY_INSTALL_OPENFCST "OFF" CACHE TYPE BOOL)
set(ONLY_INSTALL_P4EST "OFF" CACHE TYPE BOOL)
set(ONLY_INSTALL_METIS "OFF" CACHE TYPE BOOL)
set(ONLY_INSTALL_PETSC "OFF" CACHE TYPE BOOL)
set(ONLY_INSTALL_TRILINOS "OFF" CACHE TYPE BOOL)
set(ONLY_INSTALL_DEALII "OFF" CACHE TYPE BOOL)
# set(ONLY_INSTALL_DAKOTA "OFF" CACHE TYPE BOOL)
set(ONLY_INSTALL_CPPTEST "OFF" CACHE TYPE BOOL)
set(ONLY_INSTALL_PARAVIEW "OFF" CACHE TYPE BOOL)
set(ONLY_INSTALL_SQLITE "OFF" CACHE TYPE BOOL)
set(ONLY_INSTALL_ALGLIB "OFF" CACHE TYPE BOOL)
set(ONLY_INSTALL_COLDAE "OFF" CACHE TYPE BOOL)

#set bool for determining if OpenFCST will create its own library file
#typically set to OFF because this doubles the time required to Build the
#OpenFCST external project
set(OPENFCST_LIBRARY "OFF" CACHE TYPE BOOL)

#Set dimensions for OpenFCST
set(OPENFCST_DIMENSIONS "2" CACHE TYPE STRING)

#Restrict OPENFCST_DIMENSIONS to 2 options to be set in GUI
set_property(CACHE OPENFCST_DIMENSIONS PROPERTY STRINGS 2 3)

#Confirm it is only one of the set options if not then end installation
if(NOT (OPENFCST_DIMENSIONS STREQUAL "1" OR 
        OPENFCST_DIMENSIONS STREQUAL "2" OR 
        OPENFCST_DIMENSIONS STREQUAL "3"))
  MESSAGE(FATAL_ERROR 
            "\n"
            "=======================================================\n"
            "\n"
            " ERROR: The OpenFCST dimension can only be set to:     \n"
            "            1 - for both 2D and 3D simulations\n"
            "            2 - for 2D simulations\n"
            "            3 - for 3D simulations\n"
            "        You have it set to:\n"
            "            -DOPENFCST_DIMENSIONS = ${OPENFCST_DIMENSIONS}\n"
            "        Please specify the appropriate value with the  \n"
            "        CMake variable:"
            "            -DOPENFCST_DIMENSIONS=...                  \n"
            "        Or if you are using the OpenFCST shell script  \n"
            "        use the flag:                                  \n"
            "            --openfcst-dimen=...                       \n"
            "\n"
            "=======================================================\n"
            "\n"
            )
endif()

#For the following programs:
#1. Set Build Type
#2. Restriction options to store in the variable
#3. Check if there is an error

#Trilinos
set(TRILINOS_BUILD_TYPE "Release" CACHE TYPE STRING)
set_property(CACHE TRILINOS_BUILD_TYPE PROPERTY STRINGS Release Debug)
ErrMsgBuildTypeRD(${TRILINOS_BUILD_TYPE} "Trilinos" "TRILINOS_BUILD_TYPE" "trilinos-debug")

#OpenFCST
set(OPENFCST_BUILD_TYPE "Release" CACHE TYPE STRING)
set_property(CACHE OPENFCST_BUILD_TYPE PROPERTY STRINGS Release DebugRelease Debug)
ErrMsgBuildTypeRDB(${OPENFCST_BUILD_TYPE} "OpenFCST" "OPENFCST_BUILD_TYPE" "openfcst-debug")

#PETSc
set(PETSC_BUILD_TYPE "Release" CACHE TYPE STRING)
set_property(CACHE PETSC_BUILD_TYPE PROPERTY STRINGS Release Debug)
ErrMsgBuildTypeRD(${PETSC_BUILD_TYPE} "PETSc" "PETSC_BUILD_TYPE" "petsc-debug")

#Set PETSc Architecture required by deal.II
#Set default values to those used on OpenSUSE linux
if(PETSC_BUILD_TYPE STREQUAL "Release")
  set(PETSC_ARCH arch-linux2-c-opt CACHE TYPE STRING) #May need to change depending on system architecture
else()
  set(PETSC_ARCH arch-linux2-c-debug CACHE TYPE STRING) #May need to change depending on system architecture
endif()


#Turn ON/OFF optional External Project
set(OPENFCST_WITH_DAKOTA "OFF" CACHE TYPE BOOL)
set(OPENFCST_WITH_PETSC "OFF" CACHE TYPE BOOL)
set(OPENFCST_WITH_TRILINOS "OFF" CACHE TYPE BOOL)
set(OPENFCST_WITH_PARAVIEW "OFF" CACHE TYPE BOOL)

#Create versions for external projects, so user can select different versions
#Default values should be the newest working version
set(DEALII_VER "v8.4.1" CACHE TYPE STRING)
# set(DAKOTA_VER "v5.4_stable" CACHE TYPE STRING)
set(METIS_VER "v5.1" CACHE TYPE STRING)
set(P4EST_VER "v3.4.2" CACHE TYPE STRING)
set(PARAVIEW_VER "v4.3.1")
set(PETSC_VER "v3.4.4" CACHE TYPE STRING)
set(TRILINOS_VER "v11.4.1" CACHE TYPE STRING)

#Set path to pre-compiled External Project Path directories
set(ALGLIB_DIR "${CMAKE_INSTALL_PREFIX}/contrib/ALGLIB")
set(COLDAE_DIR "${CMAKE_INSTALL_PREFIX}/contrib/COLDAE")
set(CPPTEST_DIR "${CMAKE_INSTALL_PREFIX}/contrib/cpptest")
# set(DAKOTA_DIR "${CMAKE_INSTALL_PREFIX}/contrib/Dakota/${DAKOTA_VER}" CACHE TYPE PATH)
set(DEALII_DIR "${CMAKE_INSTALL_PREFIX}/contrib/deal.II/${DEALII_VER}" CACHE TYPE PATH)
set(METIS_DIR "${CMAKE_INSTALL_PREFIX}/contrib/Metis/${METIS_VER}" CACHE TYPE PATH)
set(P4EST_DIR "${CMAKE_INSTALL_PREFIX}/contrib/p4est/${P4EST_VER}" CACHE TYPE PATH)
set(PARAVIEW_DIR "${CMAKE_INSTALL_PREFIX}/contrib/ParaView/${PARAVIEW_VER}")
set(PETSC_DIR "${CMAKE_INSTALL_PREFIX}/contrib/PETSc/${PETSC_VER}" CACHE TYPE PATH)
set(SQLITE_DIR "${CMAKE_INSTALL_PREFIX}/contrib/sqlite")
set(TRILINOS_DIR "${CMAKE_INSTALL_PREFIX}/contrib/Trilinos/${TRILINOS_VER}" CACHE TYPE PATH)

#Go install all External projects
if(NOT (ONLY_INSTALL_OPENFCST OR
          ONLY_INSTALL_P4EST OR
          ONLY_INSTALL_METIS OR
          ONLY_INSTALL_PETSC OR
          ONLY_INSTALL_TRILINOS OR
          ONLY_INSTALL_DEALII OR
#           ONLY_INSTALL_DAKOTA OR
          ONLY_INSTALL_CPPTEST OR
          ONLY_INSTALL_PARAVIEW OR
          ONLY_INSTALL_SQLITE OR
          ONLY_INSTALL_ALGLIB OR
          ONLY_INSTALL_COLDAE
          )
  )
  add_subdirectory(contrib)
endif()

set(INSTALL_OPENFCST "TRUE")
# set(PACKAGE_LIST "P4EST" "METIS" "PETSC" "TRILINOS" "DEALII" "DAKOTA" "CPPTEST" "PARAVIEW" "SQLITE" "ALGLIB" "COLDAE")
set(PACKAGE_LIST "P4EST" "METIS" "PETSC" "TRILINOS" "DEALII" "CPPTEST" "PARAVIEW" "SQLITE" "ALGLIB" "COLDAE") 
foreach(ELEMENT ${PACKAGE_LIST})
  if(ONLY_INSTALL_${ELEMENT})
    STRING(TOLOWER ${ELEMENT} TEMP) #Set var TEMP to be lower case version of ELEMENT used in Error Message
    #Check that --library-fcst wasn't also specified at the same time, if so throw error
    ErrMsgConflictVars("${ONLY_INSTALL_OPENFCST}" "${ONLY_INSTALL_${ELEMENT}}" "ON" "ON" "ONLY_INSTALL_OPENFCST" "ONLY_INSTALL_${ELEMENT}" "library-fcst" "library-${TEMP}")
    
    #Since we are not setting dependencies then packages would be installed at the same time, so
    #check that no other packages are also being specified with --libary-XXX command
    foreach(OTHER_PACKAGE ${PACKAGE_LIST})
      if( NOT (OTHER_PACKAGE STREQUAL ${ELEMENT}))
        STRING(TOLOWER ${OTHER_PACKAGE} TEMP2) #Set var TEMP2 to be lower case version of OTHER_PACKAGE used in Error Message
        ErrMsgConflictVars("${ONLY_INSTALL_${OTHER_PACKAGE}}" "${ONLY_INSTALL_${ELEMENT}}" "ON" "ON" "ONLY_INSTALL_${OTHER_PACKAGE}" "ONLY_INSTALL_${ELEMENT}" "library-${TEMP2}" "library-${TEMP}")
      endif()
    endforeach()
    
    #Since we only want to install OpenFCST set bool to turn off its installation later
    set(INSTALL_OPENFCST "FALSE")
    
    if(ELEMENT STREQUAL "P4EST")
      add_subdirectory(contrib/p4est)
    elseif(ELEMENT STREQUAL "METIS")
      add_subdirectory(contrib/Metis)
    elseif(ELEMENT STREQUAL "PETSC")
      add_subdirectory(contrib/PETSc)
    elseif(ELEMENT STREQUAL "TRILINOS")
      add_subdirectory(contrib/Trilinos)
    elseif(ELEMENT STREQUAL "DEALII")
      add_subdirectory(contrib/deal.II)
#     elseif(ELEMENT STREQUAL "DAKOTA")
#       add_subdirectory(contrib/Dakota)
    elseif(ELEMENT STREQUAL "CPPTEST")
      add_subdirectory(contrib/cpptest)
    elseif(ELEMENT STREQUAL "PARAVIEW")
      add_subdirectory(contrib/ParaView)
    elseif(ELEMENT STREQUAL "SQLITE")
      add_subdirectory(contrib/sqlite)
    elseif(ELEMENT STREQUAL "ALGLIB")
      add_subdirectory(contrib/ALGLIB)
    elseif(ELEMENT STREQUAL "COLDAE")
      add_subdirectory(contrib/COLDAE)
    endif()

  endif()
endforeach()


#-----------------------
#
# Copy/Install files into the Install folder
# 
#-----------------------
# Copy files to Install folder
file(COPY README run_tests tests.cmake DESTINATION ${CMAKE_INSTALL_PREFIX})
# Copy folders to Install folder
file(COPY doc examples my_data databases test ../python DESTINATION ${CMAKE_INSTALL_PREFIX})



#The configure_file() command will copy the target file to the Install folder.
#It will also check the file to see if there are any variables of the form @VAR_NAME@
#and replace them with the CMake variable value of the same variable name
configure_file(fcst_env.sh.in
  "${CMAKE_INSTALL_PREFIX}/fcst_env.sh" @ONLY)
configure_file(fcst_run_nohup.sh
  "${CMAKE_INSTALL_PREFIX}/bin/fcst_run_nohup.sh" @ONLY)

  
#-----------------------
#
# Install OpenFCST
# 
#-----------------------

# Make sure that all openFCST dependencies are specified, in particular
# packages needed
# Note: KDevelop tries to install both 2D and 3D at same time leading to 
#       broken mpi compilers, so need to set OpenFCST 2D/3D Depends
#       so they are guarenteed to install separately.
if(ONLY_INSTALL_OPENFCST)
  set(OPENFCST2D_DEPENDS "") #default value
else()
  set(OPENFCST2D_DEPENDS "COLDAE") #default value
endif()

set(OPENFCST3D_DEPENDS "OpenFCST_2D") #default value
if(OPENFCST_DIMENSIONS STREQUAL "3") #if not installing 2D then change 3D default value
  if(ONLY_INSTALL_OPENFCST)
    set(OPENFCST3D_DEPENDS "")
  else()
    set(OPENFCST3D_DEPENDS "COLDAE")
  endif() 
endif()

if(OPENFCST_BUILD_TYPE STREQUAL "Release")
  set(SHOW_LOG "1")
else()
  #So if we are running in debug mode for developers, lets show them the build log
  #as its happening for their convenience
  set(SHOW_LOG "0")
endif()

if(OPENFCST_DIMENSIONS STREQUAL "1")
  set(OPENFCST_INSTALL_DIMS "2" "3")
elseif(OPENFCST_DIMENSIONS STREQUAL "2")
  set(OPENFCST_INSTALL_DIMS "2")
elseif(OPENFCST_DIMENSIONS STREQUAL "3")
  set(OPENFCST_INSTALL_DIMS "3")
endif()

#If we are only install one of the other packages then need to set
#OPENFCST_INSTALL_DIMS to null so then foreach loop below will not run
if(NOT INSTALL_OPENFCST)
  set(OPENFCST_INSTALL_DIMS "")
endif()

foreach(OPENFCST_INSTALL_DIM ${OPENFCST_INSTALL_DIMS})
  
  if(OPENFCST_INSTALL_DIM STREQUAL "2")
    set(OPENFCST_DEPENDS "${OPENFCST2D_DEPENDS}")
  elseif(OPENFCST_INSTALL_DIM STREQUAL "3")
    set(OPENFCST_DEPENDS "${OPENFCST3D_DEPENDS}")
  endif()
  
  set(EXTERNAL_PROJECT_NAME "OpenFCST_${OPENFCST_INSTALL_DIM}D")
  ExternalProject_Add(
    ${EXTERNAL_PROJECT_NAME}
    DEPENDS ${OPENFCST_DEPENDS}
    #URL ${CMAKE_SOURCE_DIR}/fcst/${EXTERNAL_PROJECT_NAME}.tar.gz
    BINARY_DIR fcst/Build
    PREFIX fcst
    UPDATE_COMMAND ""
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/fcst
      -DOPENFCST_INSTALL_CONFIG=${CMAKE_INSTALL_PREFIX}
      -DOPENFCST_DIMENSIONS_CONFIG=${OPENFCST_DIMENSIONS}
      -DOPENFCST_DIMENSIONS=${OPENFCST_INSTALL_DIM}
      -DCMAKE_BUILD_TYPE=${OPENFCST_BUILD_TYPE}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DOPENFCST_C11_FLAG=${OPENFCST_C11_FLAG}
      -DCMAKE_EXE_LINKER_FLAGS=${OPENFCST_CXX_EXTRA_FLAGS}
      -DCMAKE_CXX_FLAGS=${OPENFCST_CXX_EXTRA_FLAGS}
      -DCMAKE_C_FLAGS=${OPENFCST_C_EXTRA_FLAGS}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DOPENFCST_LIBRARY=${OPENFCST_LIBRARY}
      -DOPENFCST_WITH_PETSC=${OPENFCST_WITH_PETSC}
      -DOPENFCST_WITH_TRILINOS=${OPENFCST_WITH_TRILINOS}
      -DDEALII_DIR=${DEALII_DIR}
      -DOPENFCST_WITH_DAKOTA=${OPENFCST_WITH_DAKOTA}
#       -DDAKOTA_DIR=${DAKOTA_DIR}
      -DP4EST_DIR=${P4EST_DIR}
      -DMETIS_DIR=${METIS_DIR}
      -DTRILINOS_DIR=${TRILINOS_DIR}
      -DPETSC_DIR=${PETSC_DIR}
      -DALGLIB_DIR=${ALGLIB_DIR}
      -DCOLDAE_DIR=${COLDAE_DIR}
      -DSQLITE_DIR=${SQLITE_DIR}
      -DCPPTEST_DIR=${CPPTEST_DIR}
      -DBoost_INCLUDE_DIR=${Boost_INCLUDE_DIRS}
      -DBoost_LIBRARY_DIRS=${Boost_LIBRARY_DIRS}
      -DBoost_IOSTREAMS_LIBRARY=${Boost_IOSTREAMS_LIBRARY}
      -DBoost_IOSTREAMS_LIBRARY_DEBUG=${Boost_IOSTREAMS_LIBRARY_DEBUG}
      -DBoost_IOSTREAMS_LIBRARY_RELEASE=${Boost_IOSTREAMS_LIBRARY_RELEASE}
      -DBoost_SERIALIZATION_LIBRARY=${Boost_SERIALIZATION_LIBRARY}
      -DBoost_SERIALIZATION_LIBRARY_DEBUG=${Boost_SERIALIZATION_LIBRARY_DEBUG}
      -DBoost_SERIALIZATION_LIBRARY_RELEASE=${Boost_SERIALIZATION_LIBRARY_RELEASE}
      -DBoost_SYSTEM_LIBRARY=${Boost_SYSTEM_LIBRARY}
      -DBoost_SYSTEM_LIBRARY_DEBUG=${Boost_SYSTEM_LIBRARY_DEBUG}
      -DBoost_SYSTEM_LIBRARY_RELEASE=${Boost_SYSTEM_LIBRARY_RELEASE}
      -DBoost_THREAD_LIBRARY=${Boost_THREAD_LIBRARY}
      -DBoost_THREAD_LIBRARY_DEBUG=${Boost_THREAD_LIBRARY_DEBUG}
      -DBoost_THREAD_LIBRARY_RELEASE=${Boost_THREAD_LIBRARY_RELEASE}
      -DBoost_FILESYSTEM_LIBRARY=${Boost_FILESYSTEM_LIBRARY}
      -DBoost_FILESYSTEM_LIBRARY_DEBUG=${Boost_FILESYSTEM_LIBRARY_DEBUG}
      -DBoost_FILESYSTEM_LIBRARY_RELEASE=${Boost_FILESYSTEM_LIBRARY_RELEASE}
      -DBoost_REGEX_LIBRARY=${Boost_REGEX_LIBRARY}
      -DBoost_REGEX_LIBRARY_DEBUG=${Boost_REGEX_LIBRARY_DEBUG}
      -DBoost_REGEX_LIBRARY_RELEASE=${Boost_REGEX_LIBRARY_RELEASE}
      -DBoost_SIGNALS_LIBRARY=${Boost_SIGNALS_LIBRARY}
      -DBoost_SIGNALS_LIBRARY_DEBUG=${Boost_SIGNALS_LIBRARY_DEBUG}
      -DBoost_SIGNALS_LIBRARY_RELEASE=${Boost_SIGNALS_LIBRARY_RELEASE}
      -DBoost_PROGRAM_OPTIONS_LIBRARY=${Boost_PROGRAM_OPTIONS_LIBRARY}
      -DBoost_PROGRAM_OPTIONS_LIBRARY_DEBUG=${Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG}
      -DBoost_PROGRAM_OPTIONS_LIBRARY_RELEASE=${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE}
      
  
    STAMP_DIR fcst/stamp
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/fcst
    TMP_DIR fcst/tmp
    #BUILD_IN_SOURCE 1
  
    LOG_DOWNLOAD 1 # Wrap download in script to log output
    LOG_UPDATE 1 # Wrap update in script to log output
    LOG_CONFIGURE 1 # Wrap configure in script to log output
    LOG_BUILD ${SHOW_LOG} # Wrap build in script to log output
    LOG_INSTALL 1 # Wrap install in script to log output
  )

  ExternalProject_Add_Step(${EXTERNAL_PROJECT_NAME} BEGIN_INSTALL # Names of project and custom step
    COMMAND echo "-----------------------------------"
    COMMAND echo ""
    COMMAND echo "Installing ${EXTERNAL_PROJECT_NAME}"
    COMMAND echo ""
    COMMAND echo "-----------------------------------"
    DEPENDEES  # Steps on which this step depends
    DEPENDERS mkdir # Steps that depend on this step
    ALWAYS 1
  )
  
endforeach()

#-----------------------
#
# Modfiy files with hardcoded paths dependent on each system at Install
# 
#-----------------------

#Fix CTestTestfile so it has the correct path to necessary config.txt file
configure_file(${CMAKE_SOURCE_DIR}/test/CTestTestfile.cmake
   "${CMAKE_INSTALL_PREFIX}/test/CTestTestfile.cmake" @ONLY)
   
configure_file(${CMAKE_SOURCE_DIR}/tests.cmake
   "${CMAKE_INSTALL_PREFIX}/tests.cmake" @ONLY)

#Modify Sphinx Makefiles to have the correct SPHINXBUILD variable
#when installing
configure_file(${CMAKE_SOURCE_DIR}/examples/Makefile
  "${CMAKE_INSTALL_PREFIX}/examples/Makefile" @ONLY)
    
configure_file(${CMAKE_SOURCE_DIR}/examples/make.bat
  "${CMAKE_INSTALL_PREFIX}/examples/make.bat" @ONLY)
  
configure_file(${CMAKE_SOURCE_DIR}/../python/bin/writeVTK.py
  "${CMAKE_INSTALL_PREFIX}/python/bin/writeVTK.py" @ONLY)
  
configure_file(${CMAKE_SOURCE_DIR}/test/vtkLocationCompare.py
  "${CMAKE_INSTALL_PREFIX}/test/vtkLocationCompare.py" @ONLY)



#-----------------------
#
# Create Doxygen documentation
# 
#-----------------------
#If Doxygen was found in the beginning now we want to create the documentation
if(DOXYGEN_FOUND)
#Create custom phony Make command to produce Doxygen Documentation
  add_custom_command(OUTPUT create_doxygen_doc 
                      COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile &> ${CMAKE_BINARY_DIR}/3a_Doxygen_report.txt
                      WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/doc
                      ) 
  
  add_custom_target(doxygen_doc ALL DEPENDS create_doxygen_doc) # See WARNING below
endif()

# WARNING: There is a bug in KDevelop where if you use the ALL command 
#          in multiple add_custom_target commands then you get errors 
#          when compiling

# #-----------------------
# #
# # Create Sphinx documentation
# # 
# #-----------------------
# If SPHINX was found in the beginning now we want to create the examples
if(SPHINX_FOUND)

    if(NOT DEFINED SPHINX_THEME)
        set(SPHINX_THEME default)
    endif()
 
    if(NOT DEFINED SPHINX_THEME_DIR)
        set(SPHINX_THEME_DIR)
    endif()
 
    # configured documentation tools and intermediate build results
    set(BINARY_BUILD_DIR "${CMAKE_INSTALL_PREFIX}/examples/_build")
    
    # Sphinx cache with pickled ReST documents
    set(SPHINX_CACHE_DIR "${CMAKE_INSTALL_PREFIX}/examples/_build/doctrees")
    
    # HTML output directory
    set(SPHINX_HTML_DIR "${CMAKE_INSTALL_PREFIX}/examples/_build/html")
    
    add_custom_command(OUTPUT create_examples
                    COMMAND ${SPHINX_EXECUTABLE} -q -b html -d ${SPHINX_CACHE_DIR} ${CMAKE_INSTALL_PREFIX}/examples ${SPHINX_HTML_DIR} &> ${CMAKE_BINARY_DIR}/4a_examples.txt
                    WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/examples
                    )
                    
    add_custom_target(examples ALL DEPENDS create_examples)


endif()