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
# This script installs PETSc
#
# Developer: Chad Balen, 2014
#
#######################

#-----------------------
#
# Decide How to install
# 
#-----------------------
#Check for *.tar.gz file in respective src contrib folder
#If that fails then go and download it from internet

#Set default variables during installation
set(EXTERNAL_PROJECT_NAME "PETSc")
set(PETSC_TAR_FILE "${EXTERNAL_PROJECT_NAME}_${PETSC_VER}.tar.gz")

#Give warning to notify user that a pre-installed package was not found as an FYI
# ErrMsgFindPackagePIDNF(${EXTERNAL_PROJECT_NAME} "PETSc_DIR" "petsc-dir" ${PETSC_TAR_FILE})

find_path(PETSC_TAR_DIR
            ${PETSC_TAR_FILE}
            HINTS ${CMAKE_CURRENT_SOURCE_DIR}
            )
if(PETSC_TAR_DIR STREQUAL "PETSC_TAR_DIR-NOTFOUND")
  # NOTE: This value should be hardcoded because:
  #         1) Too difficult to also manage variable for selecting URL since naming
  #            convention different than *.tar.gz file
  #         2) If we switch to URL convention this means each Version variable
  #            will have a different convention; super confusing to keep track of
  #         3) URL not guarenteed to work forever and ever so this is more of a
  #            last resort
  #         4) If we get to this point lets just make user use software version
  #            we know works
  set(PETSC_URL_PATH "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.4.4.tar.gz")
  
#   ErrMsgFindPackageTFNFW(${EXTERNAL_PROJECT_NAME} ${PETSC_TAR_FILE})
  
else()
  set(PETSC_URL_PATH "${PETSC_TAR_DIR}/${PETSC_TAR_FILE}")
endif()

#Get path to BLAS so we can tell PETSc were to look
#get_filename_component(BLAS_PATH ${BLAS_LIBRARIES} PATH)

#If pre-installed software was not found then lets now install
#-----------------------
#
# Install PETSc
# 
#-----------------------
if(PETSC_BUILD_TYPE STREQUAL "Release")
  set(PETSC_DEBUG "0" TYPE STRING)
else()
  set(PETSC_DEBUG "1" TYPE STRING)
endif()

ExternalProject_Add(
  ${EXTERNAL_PROJECT_NAME}
  DEPENDS ${PETSC_DEPENDS}
  URL ${PETSC_URL_PATH}
  PREFIX ${PETSC_DIR}
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ./config/configure.py --with-debugging=${PETSC_DEBUG} CXXOPTFLAGS=${OPENFCST_C11_FLAG} --with-shared-libraries=1 --with-x=0 --with-mpi=1 --download-hypre=1 --download-scalapack=1 --download-mumps=1 --download-ptscotch=1
  BUILD_COMMAND make
  INSTALL_COMMAND make all test

  STAMP_DIR stamp
  SOURCE_DIR ${PETSC_DIR}
  TMP_DIR tmp
  BUILD_IN_SOURCE 1

  LOG_DOWNLOAD 1 # Wrap download in script to log output
  LOG_UPDATE 1 # Wrap update in script to log output
  LOG_CONFIGURE 1 # Wrap configure in script to log output
  LOG_BUILD 1 # Wrap build in script to log output
  LOG_INSTALL 1 # Wrap install in script to log output
)
  
ExternalProject_Add_Step(${EXTERNAL_PROJECT_NAME} BEGIN_INSTALL # Names of project and custom step
  COMMAND echo "-----------------------------------"
  COMMAND echo ""
  COMMAND echo "Installing ${EXTERNAL_PROJECT_NAME} ${PETSC_BUILD_TYPE} ${PETSC_VER}"
  COMMAND echo ""
  COMMAND echo "-----------------------------------"
  DEPENDEES  # Steps on which this step depends
  DEPENDERS mkdir # Steps that depend on this step
  ALWAYS 1
)