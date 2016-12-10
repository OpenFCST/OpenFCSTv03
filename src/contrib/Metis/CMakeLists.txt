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
# This script installs Metis
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
set(EXTERNAL_PROJECT_NAME "Metis")
set(METIS_TAR_FILE "${EXTERNAL_PROJECT_NAME}_${METIS_VER}.tar.gz")

#Give warning to notify user that a pre-installed package was not found as an FYI
# ErrMsgFindPackagePIDNF(${EXTERNAL_PROJECT_NAME} "METIS_DIR" "metis-dir" ${METIS_TAR_FILE})

find_path(METIS_TAR_DIR
            ${METIS_TAR_FILE}
            HINTS ${CMAKE_CURRENT_SOURCE_DIR}
            )
if(METIS_TAR_DIR STREQUAL "METIS_TAR_DIR-NOTFOUND")
  # NOTE: This value should be hardcoded because:
  #         1) Too difficult to also manage variable for selecting URL since naming
  #            convention different than *.tar.gz file
  #         2) If we switch to URL convention this means each Version variable
  #            will have a different convention; super confusing to keep track of
  #         3) URL not guarenteed to work forever and ever so this is more of a
  #            last resort
  #         4) If we get to this point lets just make user use software version
  #            we know works
  set(METIS_URL_PATH "http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz")
  
#   ErrMsgFindPackageTFNFW(${EXTERNAL_PROJECT_NAME} ${METIS_TAR_FILE})

else()
  set(METIS_URL_PATH "${METIS_TAR_DIR}/${METIS_TAR_FILE}")
endif()

#If pre-installed software was not found then lets now install
#-----------------------
#
# Install Metis
# 
#-----------------------
ExternalProject_Add(
  ${EXTERNAL_PROJECT_NAME}
  DEPENDS ${METIS_DEPENDS}
  URL ${METIS_URL_PATH}
  #BINARY_DIR ${EXTERNAL_PROJECT_NAME}/src
  PREFIX src
  UPDATE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/configure_Metis --metis-makefile=${CMAKE_BINARY_DIR}/contrib/Metis/src/Makefile
  CONFIGURE_COMMAND make CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} config shared=1 prefix=${METIS_DIR} 
  #BUILD_COMMAND ""
  #INSTALL_COMMAND ""

  STAMP_DIR stamp
  SOURCE_DIR src
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
  COMMAND echo "Installing ${EXTERNAL_PROJECT_NAME} ${METIS_VER}"
  COMMAND echo ""
  COMMAND echo "-----------------------------------"
  DEPENDEES  # Steps on which this step depends
  DEPENDERS mkdir # Steps that depend on this step
  ALWAYS 1
)