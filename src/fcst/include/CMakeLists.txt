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
# This script installs the header files into the
# installation directory.
#
# Developer: Chad Balen, 2014
#
#######################

# NOTE: In the future if you wish to create NEW FOLDERS then please append
#the two sections below (preferably in alphabetical order like they are shown
#when you open the includes folder in Dolphin. The steps are:
#
# 1. This is done by copying the following line into Section A:
#
# file(GLOB XXX_HFiles XXX/*.h)
#
#to where the new folder would be placed alphabetically in the list.
#Then replace the XXX's with the name of the folder
#
# 2. Copy the following lines into Section B (again alphabetically):
#
# INSTALL(FILES ${XXX_HFiles}
#          DESTINATION include/XXX
#          PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
#
# Then replace the XXX's with the name of the folder again.

#----------------------------
#
# Section A
#
#----------------------------

#Create variable SRCFiles containing a list of all header files in their 
#respective directories
file(GLOB HFiles ${HEADER_DIR}/*.h)
file(GLOB AFCC_HFiles ${HEADER_DIR}/AFCC/*.h)
file(GLOB application_core_HFiles ${HEADER_DIR}/application_core/*.h)
file(GLOB applications_HFiles ${HEADER_DIR}/applications/*.h)
file(GLOB contribs_HFiles ${HEADER_DIR}/contribs/*.h)
file(GLOB equations_HFiles ${HEADER_DIR}/equations/*.h)
file(GLOB grid_HFiles ${HEADER_DIR}/grid/*.h)
file(GLOB layers_HFiles ${HEADER_DIR}/layers/*.h)
file(GLOB materials_HFiles ${HEADER_DIR}/materials/*.h)
file(GLOB microscale_HFiles ${HEADER_DIR}/microscale/*.h)
file(GLOB postprocessing_HFiles ${HEADER_DIR}/postprocessing/*.h)
file(GLOB reactions_HFiles ${HEADER_DIR}/reactions/*.h)
file(GLOB solvers_HFiles ${HEADER_DIR}/solvers/*.h)
file(GLOB utils_HFiles ${HEADER_DIR}/utils/*.h)


#----------------------------
#
# Section B
#
#----------------------------

#Installation of the header files from above
INSTALL(FILES ${HFiles}
         DESTINATION include
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(FILES ${AFCC_HFiles}
         DESTINATION include/AFCC
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(FILES ${application_core_HFiles}
         DESTINATION include/application_core
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(FILES ${applications_HFiles}
         DESTINATION include/applications
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(FILES ${contribs_HFiles}
         DESTINATION include/contribs
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(FILES ${equations_HFiles}
         DESTINATION include/equations
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(FILES ${grid_HFiles}
         DESTINATION include/grid
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(FILES ${layers_HFiles}
         DESTINATION include/layers
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(FILES ${materials_HFiles}
         DESTINATION include/materials
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(FILES ${microscale_HFiles}
         DESTINATION include/microscale
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(FILES ${postprocessing_HFiles}
         DESTINATION include/postprocessing
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(FILES ${reactions_HFiles}
         DESTINATION include/reactions
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(FILES ${solvers_HFiles}
         DESTINATION include/solvers
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(FILES ${utils_HFiles}
         DESTINATION include/utils
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)