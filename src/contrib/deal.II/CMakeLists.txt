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
# This script installs deal.II
#
# Developer: Chad Balen, 2014
#
#######################

if(OPENFCST_WITH_PETSC OR OPENFCST_WITH_TRILINOS)
  set(OPENFCST_WITH_P4EST "ON")
  set(OPENFCST_WITH_METIS "ON")
else()
  set(OPENFCST_WITH_P4EST "OFF")
  set(OPENFCST_WITH_METIS "OFF")
endif()

#-----------------------
#
# Decide How to install
# 
#-----------------------
#Check for *.tar.gz file in respective src contrib folder
#If that fails then go and download it from internet

#Set default variables during installation
set(EXTERNAL_PROJECT_NAME "deal.II")
set(DEALII_TAR_FILE "${EXTERNAL_PROJECT_NAME}_${DEALII_VER}.tar.gz")
  
#Give warning to notify user that a pre-installed package was not found as an FYI
# ErrMsgFindPackagePIDNF(${EXTERNAL_PROJECT_NAME} "DEALII_DIR" "deal-dir" ${DEALII_TAR_FILE})
  
find_path(DEALII_TAR_DIR
            ${DEALII_TAR_FILE}
            HINTS ${CMAKE_CURRENT_SOURCE_DIR}
            )
if(DEALII_TAR_DIR STREQUAL "DEALII_TAR_DIR-NOTFOUND")
  # NOTE: This value should be hardcoded because:
  #         1) Too difficult to also manage variable for selecting URL since naming
  #            convention different than *.tar.gz file
  #         2) If we switch to URL convention this means each Version variable
  #            will have a different convention; super confusing to keep track of
  #         3) URL not guarenteed to work forever and ever so this is more of a
  #            last resort
  #         4) If we get to this point lets just make user use software version
  #            we know works
  
  ErrMsgFindPackageTFNFE(${EXTERNAL_PROJECT_NAME} ${DEALII_TAR_FILE} "DEALII_DIR deal-dir")

else()
  set(DEALII_URL_PATH "${DEALII_TAR_DIR}/${DEALII_TAR_FILE}")
endif()

#If pre-installed software was not found then lets now install
#-----------------------
#
# Install deal.II
# 
#-----------------------
ExternalProject_Add(
 ${EXTERNAL_PROJECT_NAME}
 DEPENDS ${DEALII_DEPENDS}
 URL ${DEALII_URL_PATH}
 BINARY_DIR Build
 PREFIX src
 #UPDATE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/configure_dealii --grid_in-path=${CMAKE_BINARY_DIR}/contrib/deal.II/src/source/grid/grid_in.cc --grid_in_patch-path=${CMAKE_CURRENT_SOURCE_DIR}/grid_in_patch
 CMAKE_ARGS
   -DCMAKE_INSTALL_PREFIX=${DEALII_DIR}
   #-DCMAKE_BUILD_TYPE=DebugRelease
   -DDEAL_II_WITH_CXX11=ON
   -DCMAKE_EXE_LINKER_FLAGS=${CXX_EXTRA_FLAGS}
   -DCMAKE_CXX_FLAGS=${CXX_EXTRA_FLAGS}
   -DDEAL_II_WITH_MPI=ON
   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
   -DMPI_CXX_COMPILER=${CMAKE_CXX_COMPILER}
   -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
   -DDEAL_II_WITH_METIS=${OPENFCST_WITH_METIS}
   -DMETIS_DIR=${METIS_DIR}
   -DDEAL_II_WITH_UMFPACK=ON 
   -DDEAL_II_WITH_MUMPS=OFF 
   -DDEAL_II_WITH_P4EST=${OPENFCST_WITH_P4EST} 
   -DP4EST_DIR=${P4EST_DIR}
   -DDEAL_II_WITH_PETSC=${OPENFCST_WITH_PETSC} 
   -DPETSC_DIR=${PETSC_DIR}
   -DPETSC_ARCH=${PETSC_ARCH}
   -DDEAL_II_WITH_TRILINOS=${OPENFCST_WITH_TRILINOS}
   -DTRILINOS_DIR=${TRILINOS_DIR}
   -DDEAL_II_FORCE_BUNDLED_UMFPACK=ON
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
   -DDEAL_II_COMPONENT_PARAMETER_GUI=${DEAL_II_COMPONENT_PARAMETER_GUI}
   
 #BUILD_COMMAND ""
 #INSTALL_COMMAND ""
 TEST_AFTER_INSTALL 1
   
 STAMP_DIR stamp
 SOURCE_DIR src
 TMP_DIR tmp
 #BUILD_IN_SOURCE 1
 
 LOG_DOWNLOAD 1 # Wrap download in script to log output
 LOG_UPDATE 1 # Wrap update in script to log output
 LOG_CONFIGURE 1 # Wrap configure in script to log output
 LOG_BUILD 1 # Wrap build in script to log output
 LOG_INSTALL 1 # Wrap install in script to log output
 LOG_TEST 1 # Wrap test in script to log output
)

ExternalProject_Add_Step(${EXTERNAL_PROJECT_NAME} BEGIN_INSTALL # Names of project and custom step
 COMMAND echo "-----------------------------------"
 COMMAND echo ""
 COMMAND echo "Installing ${EXTERNAL_PROJECT_NAME} ${DEALII_VER}"
 COMMAND echo ""
 COMMAND echo "-----------------------------------"
 DEPENDEES # Steps on which this step depends
 DEPENDERS mkdir # Steps that depend on this step
 ALWAYS 1
)