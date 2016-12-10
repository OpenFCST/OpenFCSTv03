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
# This script installs ParaView
#
# Developer: Chad Balen, 2015
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
set(EXTERNAL_PROJECT_NAME "ParaView")
set(PARAVIEW_URL_PATH "http://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v4.3&type=source&os=all&downloadFile=ParaView-v4.3.1-source.tar.gz")

#If pre-installed software was not found then lets now install
#-----------------------
#
# Install ParaView
# 
#-----------------------

ExternalProject_Add(
  ${EXTERNAL_PROJECT_NAME}
  DEPENDS ${PARAVIEW_DEPENDS}
  URL ${PARAVIEW_URL_PATH}
  BINARY_DIR ${PARAVIEW_DIR}/Build
  PREFIX src
  UPDATE_COMMAND ""
  CMAKE_ARGS 
    -DCMAKE_INSTALL_PREFIX=${PARAVIEW_DIR}
    -DCMAKE_BUILD_TYPE=Release
    -DPARAVIEW_USE_MPI=ON
    -DPARAVIEW_ENABLE_PYTHON=ON
    -DBUILD_SHARED_LIBS=ON
    -DPARAVIEW_BUILD_QT_GUI=OFF
    -DQT_QMAKE_EXECUTABLE=OFF
    -DPARAVIEW_ENABLE_FFMPEG=OFF
    -DPARAVIEW_USE_VISITBRIDGE=OFF
   
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
  LOG_TEST 1 # Wrap test in script to log output
)
  
ExternalProject_Add_Step(${EXTERNAL_PROJECT_NAME} BEGIN_INSTALL # Names of project and custom step
  COMMAND echo "-----------------------------------"
  COMMAND echo ""
  COMMAND echo "Installing ${EXTERNAL_PROJECT_NAME} ${PARAVIEW_VER}"
  COMMAND echo ""
  COMMAND echo "-----------------------------------"
  DEPENDEES  # Steps on which this step depends
  DEPENDERS mkdir # Steps that depend on this step
  ALWAYS 1
)