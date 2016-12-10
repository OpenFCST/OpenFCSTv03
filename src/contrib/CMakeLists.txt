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
# This script starts the individual scripts for installing
# the various external projects. 
#
# Developer: Chad Balen, 2014
#
#######################

# NOTE: This script is broken into 4 sections:
#         a) Use find_package and determine if we are going to throw fatal error at user
#         b) Check to see what packages need to be installed
#         c) Setup the dependency tree for installation
#         d) Install packages


#-----------------------
#
# SECTION 1: Perform find_package Checking
# 
#-----------------------

#Check to see if programs or found on computer
find_package(P4EST QUIET) #Keep it Quiet we'll give them the warning ourselves
find_package(METIS QUIET) #Keep it Quiet we'll give them the warning ourselves
find_package(PETSC QUIET) #Keep it Quiet we'll give them the warning ourselves
find_package(TRILINOS QUIET) #Keep it Quiet we'll give them the warning ourselves
find_package(deal.II QUIET HINTS ${DEALII_DIR}) #Keep it Quiet we'll give them the warning ourselves
# find_package(DAKOTA QUIET) #Keep it Quiet we'll give them the warning ourselves

#If they specified a different directory and the package was not found throw
#Fatal Error. Else we have issues because when it tries to install OpenFCST with
#the wrong locations for files will be given.
ErrMsgFindPackagePID("p4est" ${P4EST_FOUND} ${P4EST_DIR} "${CMAKE_INSTALL_PREFIX}/contrib/p4est/${P4EST_VER}" "P4EST_DIR" "p4est-dir")
ErrMsgFindPackagePID("Metis" ${METIS_FOUND} ${METIS_DIR} "${CMAKE_INSTALL_PREFIX}/contrib/Metis/${METIS_VER}" "METIS_DIR" "metis-dir")
ErrMsgFindPackagePID("PETSc" ${PETSC_FOUND} ${PETSC_DIR} "${CMAKE_INSTALL_PREFIX}/contrib/PETSc/${PETSC_VER}" "PETSC_DIR" "petsc-dir")
ErrMsgFindPackagePID("Trilinos" ${TRILINOS_FOUND} ${TRILINOS_DIR} "${CMAKE_INSTALL_PREFIX}/contrib/Trilinos/${TRILINOS_VER}" "TRILINOS_DIR" "trilinos-dir")
ErrMsgFindPackagePID("deal.II" ${deal.II_FOUND} ${DEALII_DIR} "${CMAKE_INSTALL_PREFIX}/contrib/deal.II/${DEALII_VER}" "DEALII_DIR" "deal-dir")
# ErrMsgFindPackagePID("Dakota" ${DAKOTA_FOUND} ${DAKOTA_DIR} "${CMAKE_INSTALL_PREFIX}/contrib/Dakota/${DAKOTA_VER}" "DAKOTA_DIR" "dakota-dir")

#Check to see OpenFCST has installed Paraview previously
find_path(PARAVIEW_BIN_DIR
  NAMES pvserver
  HINTS
    # petsc is special. Account for that
    ${PARAVIEW_DIR}
  PATH_SUFFIXES bin
  NO_DEFAULT_PATH
)

if(PARAVIEW_BIN_DIR STREQUAL "PARAVIEW_BIN_DIR-NOTFOUND")
  if(OPENFCST_WITH_PARAVIEW STREQUAL "ON")
    set(PARAVIEW_FOUND "FALSE")
    set(PARAVIEW_INSTALL "ON")
  else()
    set(PARAVIEW_FOUND "FALSE")
    set(PARAVIEW_INSTALL "OFF")
  endif()
else()
  set(PARAVIEW_FOUND "TRUE")
  set(PARAVIEW_INSTALL "OFF")
endif()

#Create Configuration strings for user output
if(deal.II_FOUND STREQUAL "0")
  set(DEALII_FOUND "FALSE")
else()
  set(DEALII_FOUND "TRUE")
endif()

#-----------------------
#
# SECTION 2: Check what needs to be installed
# 
#-----------------------
# NOTE: If statement structure is as follows:
#         Step 1: Check to see if user selected to install package (i.e. set OPENFCST_WITH_PETSC to ON, etc.)
#         Step 2: If true, then check to see if package was found and whether it is in the Install/contrib/... folder
#         Step 3: If Step 2 is true then user has specified their own pre-installed package so set the 
#                   that package installation to OFF
#         Step 4: If Step 2 is false then user has not specified their own pre-installed package so set the
#                   that package installation to ON
#         Step 5: If Step 1 was false then user does not want to install package so go to the last else 
#                   statement and set package to OFF

# if(OPENFCST_WITH_PETSC)
if(OPENFCST_WITH_PETSC OR OPENFCST_WITH_TRILINOS)
  #Check p4est installation
  if(NOT P4EST_DIR STREQUAL "${CMAKE_INSTALL_PREFIX}/contrib/p4est/${P4EST_VER}")
    set(P4EST_INSTALL "OFF")
  else()
    set(P4EST_INSTALL "ON")
  endif()

  #Check Metis installation
  if(NOT METIS_DIR STREQUAL "${CMAKE_INSTALL_PREFIX}/contrib/Metis/${METIS_VER}")
    set(METIS_INSTALL "OFF")
  else()
    set(METIS_INSTALL "ON")
  endif()
else()
  set(P4EST_INSTALL "OFF")
  set(METIS_INSTALL "OFF")
endif() 
  
#Check PETSc installation
if(OPENFCST_WITH_PETSC)
  if(NOT PETSC_DIR STREQUAL "${CMAKE_INSTALL_PREFIX}/contrib/PETSc/${PETSC_VER}")
    set(PETSC_INSTALL "OFF")
  else()
    set(PETSC_INSTALL "ON")
  endif()
else()
  set(PETSC_INSTALL "OFF")
endif()
  
# Check Trilinos installation
if(OPENFCST_WITH_TRILINOS)
  if(NOT TRILINOS_DIR STREQUAL "${CMAKE_INSTALL_PREFIX}/contrib/Trilinos/${TRILINOS_VER}")
    set(TRILINOS_INSTALL "OFF")
  else()
    set(TRILINOS_INSTALL "ON")
  endif()
else()
  set(TRILINOS_INSTALL "OFF")
endif()

#Check deal.II installation
if(NOT DEALII_DIR STREQUAL "${CMAKE_INSTALL_PREFIX}/contrib/deal.II/${DEALII_VER}")
  set(DEALII_INSTALL "OFF")
else()
  set(DEALII_INSTALL "ON")
endif()


#Check Dakota installation
# if(OPENFCST_WITH_DAKOTA)
#   if(NOT DAKOTA_DIR STREQUAL "${CMAKE_INSTALL_PREFIX}/contrib/Dakota/${DAKOTA_VER}")
#     set(DAKOTA_INSTALL "OFF")
#   else()
#     set(DAKOTA_INSTALL "ON")
#   endif()
# else()
#   set(DAKOTA_INSTALL "OFF")
# endif()


#Initialize config strings for each optional package & deal.ii
set(P4EST_CONFIG_STRING "#     p4est   ")
set(METIS_CONFIG_STRING "#     Metis   ")
set(PETSC_CONFIG_STRING "#     PETSc   ")
set(TRILINOS_CONFIG_STRING "#     Trilinos ")
set(DEALII_CONFIG_STRING "#     deal.II ")
set(PARAVIEW_CONFIG_STRING "#     ParaView")
# set(DAKOTA_CONFIG_STRING "#     Dakota  ")
#Create list (ie Array) of all packages to be shown in configuration output
# set(PACKAGE_LIST "P4EST" "METIS" "PETSC" "TRILINOS" "DEALII" "DAKOTA" "PARAVIEW")
set(PACKAGE_LIST "P4EST" "METIS" "PETSC" "TRILINOS" "DEALII" "PARAVIEW") 

#Loop through each element of PACKAGE_LIST
foreach(ELEMENT ${PACKAGE_LIST})
  # Check bool for XXX_FOUND to see if package was found.
  # If it was found then add YES to XXX_CONFIG_STRING else add NO
  if(${ELEMENT}_FOUND STREQUAL "TRUE")
    set(${ELEMENT}_CONFIG_STRING "${${ELEMENT}_CONFIG_STRING}              YES")
  else()
    set(${ELEMENT}_CONFIG_STRING "${${ELEMENT}_CONFIG_STRING}              NO")
  endif()
  
  # Check to see if XXX_INSTALL is ON (ie do we want to install the package)
  # If so then add YES and package version, else add NO and N/A
  if(${ELEMENT}_INSTALL STREQUAL "ON")
    set(${ELEMENT}_CONFIG_STRING "${${ELEMENT}_CONFIG_STRING}                   YES          ${${ELEMENT}_VER}")
  else()
    set(${ELEMENT}_CONFIG_STRING "${${ELEMENT}_CONFIG_STRING}                    NO          N/A")
  endif()  
endforeach()

MESSAGE(STATUS "
###
#
#  OpenFCST configuration:
#    OpenFCST Build Type:       ${OPENFCST_BUILD_TYPE}
#    Build OpenFCST Shared Lib: ${OPENFCST_LIBRARY}      
#    CMAKE_INSTALL_PREFIX:      ${CMAKE_INSTALL_PREFIX}
#    CMAKE_SOURCE_DIR:          ${CMAKE_SOURCE_DIR}
#    CMAKE_BINARY_DIR:          ${CMAKE_BINARY_DIR}
#    CXX11 Flag:                ${OPENFCST_C11_FLAG}
#    OpenFCST With OpenMP:      ${OPENFCST_OPENMP}
#    OpenFCST Extra Flags:      ${OPENFCST_CXX_EXTRA_FLAGS}
#
#  CMAKE settings:
#    Compiler:                  ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}
#    CXX Compiler:              ${CMAKE_CXX_COMPILER}
#    C Compiler:                ${CMAKE_C_COMPILER}
#    Fortran Compiler:          ${CMAKE_Fortran_COMPILER}
#
#  Configured Features:
#
#    PACKAGE       PRE-INSTALLATION FOUND      INSTALLING      VERSION
#   -------------------------------------------------------------------
${P4EST_CONFIG_STRING}
${METIS_CONFIG_STRING}
${PETSC_CONFIG_STRING}
${TRILINOS_CONFIG_STRING}
${DEALII_CONFIG_STRING}
${PARAVIEW_CONFIG_STRING}
#
###
"
)
#${DAKOTA_CONFIG_STRING} -- this line was added between ${DEALII_CONFIG_STRING} and ${PARAVIEW_CONFIG_STRING}

#-----------------------
#
# SECTION 3: Setup the dependency tree for installation
# 
#-----------------------
# NOTE: Now that we know what needs to be installed lets figure out how the dependency tree should look.
#       Special notes are:
#         i)   Do not need to set p4est dependency because if it is to be installed its the first package
#              so it has no dependency. So we start with Metis.
#         ii)  Packages like deal.ii that are required we skip Step 1, if statement tree structure below:
#         iii) If all packages were to be installed the order they are installed in is:
#                a) p4est
#                b) Metis
#                c) PETSc
#                d) Trilinos
#                e) deal.ii
#                f) Dakota
#                g) appframe
#                h) cpptest
#                i) sqlite
#                j) ALGLIB
#                k) COLDAE
#                l) OpenFCST 2D
#                m) OpenFCST 3D
#       The if statement tree structure is as follows:
#         Step 1: Check to see if we want to install package. If not we do not need to bother setting its dependency
#         Step 2: Check to see if the previous package is going to be installed; i.e if we are installing PETSc
#                   check if we are also installing Metis; see Special NOTE 3 for package installation order
#         Step 3: If Step 2 is true then stop if statements and set package dependency to this package
#         Step 4: If Step 2 is false check next previous package to see if its installed, if true stop and set
#                   dependency to this package
#         Step 5: Repeat Steps 2-4 till all packages are checked and then the case must be to set the package to 
#                   have no depency since it will be the first package to be installed

#Set Metis dependency
if(METIS_INSTALL)
  if(P4EST_INSTALL)
    set(METIS_DEPENDS "p4est")
  else()
    set(METIS_DEPENDS "")
  endif()
endif()

#Set PETSc dependency
if(PETSC_INSTALL)
  if(METIS_INSTALL)
    set(PETSC_DEPENDS "Metis")
  elseif(P4EST_INSTALL)
    set(PETSC_DEPENDS "p4est")
  else()
    set(PETSC_DEPENDS "")
  endif()
endif()

#Set Trilinos dependency
if(TRILINOS_INSTALL)
  if(PETSC_INSTALL)
    set(TRILINOS_DEPENDS "PETSc")
  elseif(METIS_INSTALL)
    set(TRILINOS_DEPENDS "Metis")
  elseif(P4EST_INSTALL)
    set(TRILINOS_DEPENDS "p4est")
  else()
    set(TRILINOS_DEPENDS "")
  endif()
endif()

#Set deal.ii dependency
if(DEALII_INSTALL)
  if(TRILINOS_INSTALL)
    set(DEALII_DEPENDS "Trilinos")
  elseif(PETSC_INSTALL)
    set(DEALII_DEPENDS "PETSc")
  elseif(METIS_INSTALL)
    set(DEALII_DEPENDS "Metis")
  elseif(P4EST_INSTALL)
    set(DEALII_DEPENDS "p4est")
  else()
    set(DEALII_DEPENDS "")
  endif()
endif()

#Set Dakota dependency
# if(DAKOTA_INSTALL)
#   if(DEALII_INSTALL)
#     set(DAKOTA_DEPENDS "deal.II")
#   elseif(TRILINOS_INSTALL)
#     set(DAKOTA_DEPENDS "Trilinos")
#   elseif(PETSC_INSTALL)
#     set(DAKOTA_DEPENDS "PETSc")
#   elseif(METIS_INSTALL)
#     set(DAKOTA_DEPENDS "Metis")
#   elseif(P4EST_INSTALL)
#     set(DAKOTA_DEPENDS "p4est")
#   else()
#     set(DAKOTA_DEPENDS "")
#   endif()
# endif()

#Set cpptest dependency
# if(DAKOTA_INSTALL)
#   set(CPPTEST_DEPENDS "Dakota")
if(DEALII_INSTALL)
  set(CPPTEST_DEPENDS "deal.II")
elseif(TRILINOS_INSTALL)
  set(CPPTEST_DEPENDS "Trilinos")
elseif(PETSC_INSTALL)
  set(CPPTEST_DEPENDS "PETSc")
elseif(METIS_INSTALL)
  set(CPPTEST_DEPENDS "Metis")
elseif(P4EST_INSTALL)
  set(CPPTEST_DEPENDS "p4est")
else()
  set(CPPTEST_DEPENDS "")
endif()

#set ParaView dependency
if(PARAVIEW_INSTALL)
  set(PARAVIEW_DEPENDS "cpptest")
  set(SQLITE_DEPENDS   "ParaView")
else()
  set(SQLITE_DEPENDS   "cpptest")
endif()

set(ALGLIB_DEPENDS   "sqlite")
set(COLDAE_DEPENDS   "ALGLIB")

#-----------------------
#
# SECTION 4: Install packages
# 
#-----------------------
# Install p4est
if(P4EST_INSTALL)
  add_subdirectory(p4est)
endif()

# Install Metis
if(METIS_INSTALL)
  add_subdirectory(Metis)
endif()
  
# Install PETSc
if(PETSC_INSTALL)
  add_subdirectory(PETSc)
endif()

# Install Trilinos
if(TRILINOS_INSTALL)
  add_subdirectory(Trilinos)
endif()

# Install deal.II
if(DEALII_INSTALL)
  add_subdirectory(deal.II)
endif()

# Install Dakota
if(DAKOTA_INSTALL)
  add_subdirectory(Dakota)
endif()

#Install cpptest
add_subdirectory(cpptest)

#Install ParaView
if(PARAVIEW_INSTALL)
  add_subdirectory(ParaView)
endif()

# Install sqlite
add_subdirectory(sqlite)

# Install ALGLIB
add_subdirectory(ALGLIB)

# Install COLDAE
add_subdirectory(COLDAE)