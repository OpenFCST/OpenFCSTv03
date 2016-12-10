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
# This script installs sqlite
#
# Developer: Chad Balen, 2014
#
#######################

#-----------------------
#
# Install sqlite
# 
#-----------------------
set(EXTERNAL_PROJECT_NAME "sqlite")
ExternalProject_Add(
 ${EXTERNAL_PROJECT_NAME}
 DEPENDS ${SQLITE_DEPENDS}
 URL ${CMAKE_CURRENT_SOURCE_DIR}/${EXTERNAL_PROJECT_NAME}.tar.gz
 BINARY_DIR Build
 PREFIX src
 UPDATE_COMMAND ""
 CMAKE_ARGS 
   -DCMAKE_INSTALL_PREFIX=${SQLITE_DIR}
   #-DCMAKE_BUILD_TYPE=${OPENFCST_BUILD_TYPE}
   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
 #BUILD_COMMAND ""
 #INSTALL_COMMAND ""
 
 STAMP_DIR stamp
 SOURCE_DIR src
 TMP_DIR tmp
 #BUILD_IN_SOURCE 1
 
 LOG_DOWNLOAD 1 # Wrap download in script to log output
 LOG_UPDATE 1 # Wrap update in script to log output
 LOG_CONFIGURE 1 # Wrap configure in script to log output
 LOG_BUILD 1 # Wrap build in script to log output
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