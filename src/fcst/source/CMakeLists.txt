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
# This script is the meat of the OpenFCST
# installation and installs OpenFCST. 
#
# Developer: Chad Balen, 2014
#
#######################

# Collect all source files the target consists of:
file(GLOB AFCC_CCFiles "${CMAKE_CURRENT_SOURCE_DIR}/AFCC/*.cc")
file(GLOB application_core_CCFiles "${CMAKE_CURRENT_SOURCE_DIR}/application_core/*.cc")
file(GLOB applications_CCFiles "${CMAKE_CURRENT_SOURCE_DIR}/applications/*.cc")
file(GLOB contribs_CCFiles "${CMAKE_CURRENT_SOURCE_DIR}/contribs/*.cc")
file(GLOB equations_CCFiles "${CMAKE_CURRENT_SOURCE_DIR}/equations/*.cc")
file(GLOB grid_CCFiles "${CMAKE_CURRENT_SOURCE_DIR}/grid/*.cc")
file(GLOB layers_CCFiles "${CMAKE_CURRENT_SOURCE_DIR}/layers/*.cc")
file(GLOB materials_CCFiles "${CMAKE_CURRENT_SOURCE_DIR}/materials/*.cc")
file(GLOB microscale_CCFiles "${CMAKE_CURRENT_SOURCE_DIR}/microscale/*.cc")
file(GLOB postprocessing_CCFiles "${CMAKE_CURRENT_SOURCE_DIR}/postprocessing/*.cc")
file(GLOB reactions_CCFiles "${CMAKE_CURRENT_SOURCE_DIR}/reactions/*.cc")
file(GLOB solvers_CCFiles "${CMAKE_CURRENT_SOURCE_DIR}/solvers/*.cc")
file(GLOB utils_CCFiles "${CMAKE_CURRENT_SOURCE_DIR}/utils/*.cc")
file(GLOB UnitTest_CCFiles "${CMAKE_SOURCE_DIR}/unit_tests/source/*.cc")

#Combine all the source files into a single variable
SET(SRCFiles
  main.cc
  ${AFCC_CCFiles}
  ${application_core_CCFiles}
  ${applications_CCFiles}
  ${contribs_CCFiles}
  ${equations_CCFiles}
  ${grid_CCFiles}
  ${layers_CCFiles}
  ${materials_CCFiles}
  ${microscale_CCFiles}
  ${postprocessing_CCFiles}
  ${reactions_CCFiles}
  ${solvers_CCFiles}
  ${utils_CCFiles}
  ${UnitTest_CCFiles}
  # You can specify additional files or variables of lists of files here!
  )

#Tell CMake where to find the header files the SRCFiles and all other
#external projects. NOTE: deal.II not included in this list because deal.II
#has macros to do all this work for us that are used later.
include_directories(${HEADER_DIR}
                      ${HEADER_DIR}/AFCC
                      ${HEADER_DIR}/application_core 
                      ${HEADER_DIR}/applications
                      ${HEADER_DIR}/contribs
                      ${HEADER_DIR}/equations
                      ${HEADER_DIR}/grid
                      ${HEADER_DIR}/layers
                      ${HEADER_DIR}/materials
                      ${HEADER_DIR}/microscale
                      ${HEADER_DIR}/postprocessing
                      ${HEADER_DIR}/reactions
                      ${HEADER_DIR}/solvers
                      ${HEADER_DIR}/utils
                      ${UNIT_TESTS_HEADER_DIR}
                      ${ALGLIB_DIR}/include
                      ${COLDAE_DIR}/include
                      ${CPPTEST_DIR}/include
                      ${SQLITE_DIR}/include
                      ${DAKOTA_DIR}/include
                      )
    
#Tell CMake to find the deal.II project
FIND_PACKAGE(deal.II REQUIRED
  HINTS
    ${DEALII_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR 
            "\n"
            "=======================================================\n"
            "\n"
            " ERROR: Could not locate deal.II. Please specifying    \n"
            "        the appropriate value with the CMake variable: \n"
            "            -DDEALII_DIR=...                           \n"
            "        Or if you are using the OpenFCST shell script  \n"
            "        use the flag:                                  \n"
            "            --deal-dir=...                             \n"
            "\n"
            "=======================================================\n"
            "\n"
            )
ENDIF ()

#
# Are all dependencies fullfilled?
#
if(OPENFCST_WITH_PETSC) #This if statement required b/c user may not install with PETSc
  add_definitions(-DOPENFCST_WITH_PETSC)
  IF(NOT DEAL_II_WITH_PETSC)
    MESSAGE(FATAL_ERROR 
            "\n"
            "=======================================================\n"
            "\n"
            " ERROR: The deal.II library found at:                  \n"
            "            ${DEALII_DIR}                              \n"
            "        was not configured with PETSc. Since you       \n"
            "        used your own installation of deal.ii I leave\n"
            "        it to you to fix the problem."
            "\n"
            "=======================================================\n"
            "\n"
            )
  ENDIF()
endif()

if(OPENFCST_WITH_TRLINOS) #This if statement required b/c user may not install with PETSc
  IF(NOT DEAL_II_WITH_TRILINOS)
    MESSAGE(FATAL_ERROR 
            "\n"
            "=======================================================\n"
            "\n"
            " ERROR: The deal.II library found at:                  \n"
            "            ${DEALII_DIR}                              \n"
            "        was not configured with Trilinos. Since you    \n"
            "        used your own installation of deal.ii I leave  \n"
            "        it to you to fix the problem."
            "\n"
            "=======================================================\n"
            "\n"
            )
  ENDIF()
endif()


#Add our own compiler flags, required by OpenFCST                      
add_definitions(-Ddeal_II_dimension=${OPENFCST_DIMENSIONS} -DGLIBCXX_DEBUG -D_GLIBCXX_USE_NANOSLEEP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENFCST_C11_FLAG}") #Add C++11 compiler flag

# This sets some cached variables to the values used for compiling the
# deal.II library.
DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(${PROJECT_NAME}) #Now we can create the Project

#----------------------------
#
# Find Boost
#
#----------------------------
#Don't need to do as much work since all we want to do is pass all the info from 
#Super Build CMake script. We just want to try and find boost, but we will
#force all the values from the Super Build CMake script.
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

include_directories(${Boost_INCLUDE_DIRS}
                      ${DEAL_II_INCLUDE_DIRS}
                      )

#Set the OpenFCST Version Number
set(OPENFCST_MAJOR_VERSION 0)
set(OPENFCST_MINOR_VERSION 2)
set(OPENFCST_VERSION ${OPENFCST_MAJOR_VERSION}.${OPENFCST_MINOR_VERSION})

#Create the executable
ADD_EXECUTABLE(${FCST_EXE} ${SRCFiles})

# DEAL_II_SETUP_TARGET() appends necessary include directories, linker flags, compile
# definitions and the deal.II library link interface to the given target.
DEAL_II_SETUP_TARGET(${FCST_EXE})


# Make sure the linker can find the external project libraries once it is built. 
LINK_DIRECTORIES(${ALGLIB_DIR}/lib
                  ${COLDAE_DIR}/lib
                  ${CPPTEST_DIR}/lib
                  ${CPPTEST_DIR}/lib64
                  ${DAKOTA_DIR}/lib
                  ${SQLITE_DIR}/lib
                  ${Boost_LIBRARY_DIRS}
                  ${DEALII_DIR}/lib
                  )

#Collect all library files to be used in TARGET_LINK_LIBRARIES() command             
file(GLOB ALGLIBLIB "${ALGLIB_DIR}/lib/*.so")
file(GLOB COLDAELIB "${COLDAE_DIR}/lib/*.so")
file(GLOB cpptestLIB "${CPPTEST_DIR}/lib/*.so" "${CPPTEST_DIR}/lib64/*.so") #DO NOT USE THE .a library else errors occur
file(GLOB sqliteLIB "${SQLITE_DIR}/lib/*.so")

#Specify librariess used to link to executable
TARGET_LINK_LIBRARIES(${FCST_EXE}
                         ${ALGLIBLIB} 
                         ${COLDAELIB} 
                         ${cpptestLIB} 
                         ${sqliteLIB}
                         ) 
                         
IF(CMAKE_BUILD_TYPE STREQUAL Debug)
  TARGET_LINK_LIBRARIES(${FCST_EXE}
                         ${DEAL_II_LIBRARIES_DEBUG}) 
ELSE()
  TARGET_LINK_LIBRARIES(${FCST_EXE}
                         ${DEAL_II_LIBRARIES_RELEASE}) 
ENDIF()

#Check to see if you are including Dakota. If so also link Dakota libraries
if(OPENFCST_WITH_DAKOTA)
  add_definitions(-D_WITH_DAKOTA) #Add this flag to tell OpenFCST we are using Dakota
  file(GLOB DakotaLIB_SHARED "${DAKOTA_DIR}/lib/*.so")
  file(GLOB DakotaLIB_STATIC "${DAKOTA_DIR}/lib/*.a")
  TARGET_LINK_LIBRARIES(${FCST_EXE}
                         ${DakotaLIB_SHARED}
                         ${DakotaLIB_STATIC}
                         )
endif()

#If OpenFCST not compiled in Release mode then some units tests may fail
#So we add the following flag if the build type is not Release
#As well, lets link to deal.II at the same time since it requires the same if
#statement condition
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  TARGET_LINK_LIBRARIES(${FCST_EXE}
                          ${DEAL_II_LIBRARIES_RELEASE}
                          ${Boost_IOSTREAMS_LIBRARY_RELEASE}
                          ${Boost_SERIALIZATION_LIBRARY_RELEASE}
                          ${Boost_SYSTEM_LIBRARY_RELEASE}
                          ${Boost_THREAD_LIBRARY_RELEASE}
                          ${Boost_FILESYSTEM_LIBRARY_RELEASE}
                          ${Boost_REGEX_LIBRARY_RELEASE}
                          ${Boost_SIGNALS_LIBRARY_RELEASE}
                          ${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE}
                          )
else()
  add_definitions(-O)
  TARGET_LINK_LIBRARIES(${FCST_EXE}
                          ${DEAL_II_LIBRARIES_DEBUG}
                          ${Boost_IOSTREAMS_LIBRARY_DEBUG}
                          ${Boost_SERIALIZATION_LIBRARY_DEBUG}
                          ${Boost_SYSTEM_LIBRARY_DEBUG}
                          ${Boost_THREAD_LIBRARY_DEBUG}
                          ${Boost_FILESYSTEM_LIBRARY_DEBUG}
                          ${Boost_REGEX_LIBRARY_DEBUG}
                          ${Boost_SIGNALS_LIBRARY_DEBUG}
                          ${Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG}
                          )
endif()

#If user select to have an OpenFCST library file be created then next few lines
#create a library file in the same way as you create an executable
if(OPENFCST_LIBRARY)
  ADD_LIBRARY(${FCST_LIB} SHARED ${SRCFiles})
  DEAL_II_SETUP_TARGET(${FCST_LIB})
  TARGET_LINK_LIBRARIES(${FCST_LIB} 
                         ${ALGLIBLIB} 
                         ${COLDAELIB} 
                         ${cpptestLIB} 
                         ${sqliteLIB}
                         )
                         
  if(OPENFCST_WITH_DAKOTA)
    TARGET_LINK_LIBRARIES(${FCST_LIB} 
                           ${DakotaLIB})
  endif()
  IF(CMAKE_BUILD_TYPE STREQUAL Debug)
    TARGET_LINK_LIBRARIES(${FCST_LIB}
                           ${DEAL_II_LIBRARIES_DEBUG}
                           ${Boost_IOSTREAMS_LIBRARY_DEBUG}
                           ${Boost_SERIALIZATION_LIBRARY_DEBUG}
                           ${Boost_SYSTEM_LIBRARY_DEBUG}
                           ${Boost_THREAD_LIBRARY_DEBUG}
                           ${Boost_FILESYSTEM_LIBRARY_DEBUG}
                           ${Boost_REGEX_LIBRARY_DEBUG}
                           ${Boost_SIGNALS_LIBRARY_DEBUG}
                           ${Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG}
                           ) 
  ELSE()
    TARGET_LINK_LIBRARIES(${FCST_LIB}
                           ${DEAL_II_LIBRARIES_RELEASE}
                           ${Boost_IOSTREAMS_LIBRARY_RELEASE}
                           ${Boost_SERIALIZATION_LIBRARY_RELEASE}
                           ${Boost_SYSTEM_LIBRARY_RELEASE}
                           ${Boost_THREAD_LIBRARY_RELEASE}
                           ${Boost_FILESYSTEM_LIBRARY_RELEASE}
                           ${Boost_REGEX_LIBRARY_RELEASE}
                           ${Boost_SIGNALS_LIBRARY_RELEASE}
                           ${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE}                           
                           ) 
  ENDIF()
  INSTALL(TARGETS ${FCST_LIB}
         DESTINATION lib #this location is Instal_Folder/fcst/lib
         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_WRITE WORLD_EXECUTE)
endif()

#Actually Install executable now
INSTALL(TARGETS ${FCST_EXE}
         RUNTIME
         DESTINATION ../bin
         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_WRITE WORLD_EXECUTE)
 
#-----------------------
#
# Create Config file for other CMake projects to find
# 
#-----------------------
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE OpenFCST)

# Copy fixed template files to the Build folder
configure_file(${CMAKE_SOURCE_DIR}/cmake/OpenFCSTConfig.cmake.in
  "${CMAKE_BINARY_DIR}/lib/OpenFCSTConfig.cmake" @ONLY)

configure_file(${CMAKE_SOURCE_DIR}/cmake/OpenFCSTConfigVersion.cmake.in
  "${CMAKE_BINARY_DIR}/lib/OpenFCSTConfigVersion.cmake" @ONLY)
 
# Now Install these files in the Install/fcst/lib folder
install(FILES
  "${CMAKE_BINARY_DIR}/lib/OpenFCSTConfig.cmake"
  "${CMAKE_BINARY_DIR}/lib/OpenFCSTConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" COMPONENT dev)
