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
# This script installs the Unit Test header files into the
# installation directory.
#
# Developer: Chad Balen, 2014
#
#######################

# NOTE: If you wish to add more header folders to unit_tests to be installed, please check
#       the OpenFCST include folder CMakeLists.txt for the instructions on how to add to this list
#       will have to change parts to match the variables below. I'm sure you will be able to figure
#       out the problem. If not God help you.

#Create variable SRCFiles containing a list of all header files in this directory
file(GLOB HFiles ${UNIT_TESTS_HEADER_DIR}/*.h)
file(GLOB include_HFiles ${UNIT_TESTS_HEADER_DIR}/include/*.h)

#Installation of the header files
INSTALL(FILES ${HFiles}
         DESTINATION unit_tests
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(FILES ${include_HFiles}
         DESTINATION unit_tests/include
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)