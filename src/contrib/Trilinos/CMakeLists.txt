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
# This script installs Trilinos
#
# Developer: Chad Balen, 2014
#
#######################

#-----------------------
#
# Decide How to install
# 
#-----------------------
#First we want to check if we can detect if user has already installed software
#If that fails then check for *.tar.gz file in respective src contrib folder
#If that fails then go and download it from internet


#Want to find packages quietly; i.e. if can't be found don't tell user during 
#configure, we got two other methods for installing
# NOTE: If user wished to use pre-installed software and CMake can't find it
#       they can provide a hint to it by setting the <SOFTWARE_NAME>_DIR
#       installation
#Set default variables during installation
set(EXTERNAL_PROJECT_NAME "Trilinos")
set(TRILINOS_TAR_FILE "${EXTERNAL_PROJECT_NAME}_${TRILINOS_VER}.tar.gz")

#Give warning to notify user that a pre-installed package was not found as an FYI
# ErrMsgFindPackagePIDNF(${EXTERNAL_PROJECT_NAME} "TRILINOS_DIR" "trilinos-dir" ${TRILINOS_TAR_FILE})

find_path(TRILINOS_TAR_DIR
            ${TRILINOS_TAR_FILE}
            HINTS ${CMAKE_CURRENT_SOURCE_DIR}
            )
if(TRILINOS_TAR_DIR STREQUAL "TRILINOS_TAR_DIR-NOTFOUND")
  # NOTE: This value should be hardcoded because:
  #         1) Too difficult to also manage variable for selecting URL since naming
  #            convention different than *.tar.gz file
  #         2) If we switch to URL convention this means each Version variable
  #            will have a different convention; super confusing to keep track of
  #         3) URL not guarenteed to work forever and ever so this is more of a
  #            last resort
  #         4) If we get to this point lets just make user use software version
  #            we know works
  set(TRILINOS_URL_PATH "http://trilinos.sandia.gov/download/files/trilinos-11.4.1-Source.tar.gz")
  
#   ErrMsgFindPackageTFNFW(${EXTERNAL_PROJECT_NAME} ${TRILINOS_TAR_FILE})
  
else()
  set(TRILINOS_URL_PATH "${TRILINOS_TAR_DIR}/${TRILINOS_TAR_FILE}")
endif()

if(NOT METIS_FOUND)
  set(METIS_INCLUDE_DIRS ${METIS_DIR}/include)
  set(METIS_LIBRARY_DIRS ${METIS_DIR}/lib)
endif()

#If pre-installed software was not found then lets now install
#-----------------------
#
# Install Trilinos
# 
#-----------------------
ExternalProject_Add(
  ${EXTERNAL_PROJECT_NAME}
  DEPENDS ${TRILINOS_DEPENDS}
  URL ${TRILINOS_URL_PATH}
  BINARY_DIR Build
  PREFIX src
  UPDATE_COMMAND ""
  CMAKE_ARGS 
    -D CMAKE_BUILD_TYPE:STRING=${TRILINOS_BUILD_TYPE}
    -D Trilinos_ENABLE_Fortran:BOOL=ON
    -D TPL_ENABLE_MPI=ON
    -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -D CMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
    -D CMAKE_CXX_FLAGS=-g -O3 ${OPENFCST_C11_FLAG}
    -D CMAKE_C_FLAGS=-g -O3
    -D CMAKE_FORTRAN_FLAGS=-g -O5
    -D BUILD_SHARED_LIBS=ON
    -D Trilinos_VERBOSE_CONFIGURE=FALSE
    -D CMAKE_VERBOSE_MAKEFILE=FALSE
    -D Trilinos_ENABLE_TESTS:BOOL=ON
    -D Trilinos_ENABLE_Sacado:BOOL=ON
    -D Trilinos_ENABLE_Stratimikos:BOOL=ON
    -D Trilinos_ENABLE_Teuchos:BOOL=ON
    -D Trilinos_ENABLE_EpetraExt:BOOL=ON
    -D Trilinos_ENABLE_Amesos:BOOL=ON
    -D Trilinos_ENABLE_Zoltan:BOOL=OFF
    -D Amesos_ENABLE_PARAKLETE:BOOL=OFF #Paraklete is broken don't bother trying to turn on
    -D Amesos_ENABLE_KLU:BOOL=ON
    -D TPL_ENABLE_SuperLU:BOOL=OFF
    -D TPL_ENABLE_METIS:BOOL=ON
    -D METIS_LIBRARY_DIRS=${METIS_LIBRARY_DIRS}
    -D METIS_INCLUDE_DIRS=${METIS_INCLUDE_DIRS}
    #Everything pertaining to SuperLU and SuperLUDist should be turned commented out till we can test their ability to work with deal.ii
#     -D SuperLU_INCLUDE_DIRS:FILEPATH="/usr/include"
#     -D SuperLU_LIBRARY_DIRS:PATH="/usr/lib64"
#     -D TPL_SuperLU_LIBRARIES:FILEPATH="/usr/lib64/libsuperlu.so"
#     -D TPL_ENABLE_SuperLUDist:BOOL=ON
#     -D SuperLUDist_INCLUDE_DIRS:PATH="/home/cbalen/Installation/SuperLU_DIST/src/SRC"
#     -D SuperLUDist_LIBRARY_DIRS:PATH="/home/cbalen/Installation/SuperLU_DIST/src/lib"
#     -D TPL_SuperLUDist_LIBRARIES:FILEPATH="/home/cbalen/Installation/SuperLU_DIST/src/lib/libsuperlu_dist_2.5.a"
#     -D TPL_ENABLE_ParMETIS:BOOL=ON #ParMETIS is required for SuperLUDist
#     -D ParMETIS_LIBRARY_DIRS="/home/cbalen/Installation/ParMetis/Install/lib"
#     -D ParMETIS_INCLUDE_DIRS="/home/cbalen/Installation/ParMetis/Install/include"
#     -D TPL_ENABLE_METIS:BOOL=ON
#     -D METIS_LIBRARY_DIRS="/home/cbalen/CMake_Tutorial/Mine/ExternalProjects/Install/contrib/Metis/lib"
#     -D METIS_INCLUDE_DIRS="/home/cbalen/CMake_Tutorial/Mine/ExternalProjects/Install/contrib/Metis/include"
#     -D BLAS_LIBRARY_NAMES:STRING=libblas.so
#     -D BLAS_LIBRARY_DIRS:STRING=/usr/lib64
#     -D LAPACK_LIBRARY_NAMES:STRING=liblapack.so
#     -D LAPACK_LIBRARY_DIRS:STRING=/usr/lib64
    -DCMAKE_INSTALL_PREFIX=${TRILINOS_DIR}
  #BUILD_COMMAND ""
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
  COMMAND echo "Installing ${EXTERNAL_PROJECT_NAME} ${TRILINOS_BUILD_TYPE} ${TRILINOS_VER}"
  COMMAND echo ""
  COMMAND echo "-----------------------------------"
  DEPENDEES # Steps on which this step depends
  DEPENDERS mkdir # Steps that depend on this step
  ALWAYS 1
)