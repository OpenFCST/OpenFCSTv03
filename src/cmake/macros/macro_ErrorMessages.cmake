## ---------------------------------------------------------------------
## $Id: macro_set_if_empty.cmake 31527 2013-11-03 09:58:45Z maier $
##
## Copyright (C) 2012 - 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

# NOTE: When using these functions the first Variable is the function name
#       followed by its inputs. It takes inputs in the form of strings so if
#       value is stored as a variabled remember to send to function 
#       ${FUNC_NAME_HERE}. A function cannot have an input missing else error
#       occurs. So if you wish to leave an input empty the following shows
#       what is acceptable:
#
#       Acceptable: " "
#       NOT Acceptable: ""
#
#       So using the ErrMsgRegular as an example:
#       Ex. 1 set(OPENFCST_MPI_C_COMPILER "OPENFCST_MPI_C_COMPILER-NOTFOUND")
#             ErrMsgRegular(${OPENFCST_MPI_C_COMPILER} OPENFCST_MPI_C_COMPILER-NOTFOUND mpicc MPI_DIR_HINT mpi-dir)
#
#           Output: =======================================================
#
#                    ERROR: CMake CANNOT find mpicc, please
#                           confirm you have installed mpicc. If so please
#                           specify its location with the CMake variable:
#                               -DMPI_DIR_HINT=...
#                           Or if you are using the OpenFCST shell script
#                           use the flag:
#                               --mpi-dir=...
#                   
#                   =======================================================
#
#       Ex. 2 set(OPENFCST_MPI_C_COMPILER " ")
#             ErrMsgRegular(${OPENFCST_MPI_C_COMPILER} OPENFCST_MPI_C_COMPILER-NOTFOUND mpicc MPI_DIR_HINT mpi-dir)
#
#           Output: #Nothing, condition was not met to produce output
#
# NOTE: It is recommended that all literal strings that are sent to a function to be encapsulated in quotes.
#       This is not necessary but highly recommended for readability.

#-----------------------
#
# Regular Error Message
# 
#-----------------------
# CONDITION: variable value that tells cmake if program was found or not
# COMPARISON: string that CONDITION equals when program is NOT found
# PROGRAM_NAME: string containing the name of program to be displayed in error message
# OPENFCST_VAR: string containing what the CMake variable is they can use to resolve problem
# BASH_VAR: string containing what the bash script flag is they can use to resolve problem
function(ErrMsgRegular CONDITION COMPARISON PROGRAM_NAME OPENFCST_VAR BASH_VAR)
    if(${CONDITION} STREQUAL ${COMPARISON})
        MESSAGE(FATAL_ERROR 
                  "\n"
                  "=======================================================\n"
                  "\n"
                  " ERROR: CMake CANNOT find ${PROGRAM_NAME}, please       \n"
                  "        confirm you have installed ${PROGRAM_NAME}. If so please \n"
                  "        specify its location with the CMake variable:  \n"
                  "            -D${OPENFCST_VAR}=...                         \n"
                  "        Or if you are using the OpenFCST shell script  \n"
                  "        use the flag:                                  \n"
                  "            --${BASH_VAR}=...                              \n"
                  "\n"
                  "=======================================================\n"
                  "\n"
                  )
    endif()
endfunction()

# NOTE: So this one is a bit different than above function. In the case where we don't have
#       CMake variables to specify paths when we don't want user to use a self-compiled
#       version as is the case with BLAS, LAPACK, etc. we will tell them exactly that
#       have a program like YaST install it for them so CMake has the best chances of
#       finding it.
# CONDITION: variable value that tells cmake if program was found or not

# COMPARISON: string that CONDITION equals when program IS found
# NOTE: The comparsion works a bit differently here, b/c most of the programs that require
#       this function for testing I do not want to have to deal with the wrath of uninstalling
#       them from my computer. So since I know exactly what the value is if it is found, then
#       anything else must be that it wasn't found. So COMPARISON must be what the value
#       is if the program IS found.

# PROGRAM_NAME: string containing the name of program to be displayed in error message
function(ErrMsgRegularNoSol CONDITION COMPARISON PROGRAM_NAME)
    if(NOT ${CONDITION} STREQUAL ${COMPARISON})
        MESSAGE(FATAL_ERROR 
                  "\n"
                  "=======================================================\n"
                  "\n"
                  " ERROR: CMake CANNOT find ${PROGRAM_NAME}, please      \n"
                  "        confirm you have installed ${PROGRAM_NAME}     \n"
                  "        with your operating system setup and           \n"
                  "        configuration tool (i.e. OpenSUSE's YaST,      \n"
                  "        Ubuntu's SoftwareCentre, etc.)                 \n"
                  "\n"
                  "=======================================================\n"
                  "\n"
                  )
    endif()
endfunction()


#-----------------------
#
# Build Type Error Message, when condition is EITHER Release or Debug
# 
#-----------------------
# CONDITION: value of CMake variable, used to check that user typed Release and Debug
#            the way we want
# PROGRAM_NAME: string containing the name of program to be displayed in error message
# OPENFCST_VAR: string containing what the CMake variable is they can use to resolve problem
# BASH_VAR: string containing what the bash script flag is they can use to resolve problem
function(ErrMsgBuildTypeRD CONDITION PROGRAM_NAME OPENFCST_VAR BASH_VAR)
    if(NOT (${CONDITION} STREQUAL "Release" OR 
        ${CONDITION} STREQUAL "Debug"))
      MESSAGE(FATAL_ERROR 
                "\n"
                "=======================================================\n"
                "\n"
                " ERROR: The ${PROGRAM_NAME} build type can only be set to:       \n"
                "            Release or Debug                           \n"
                "        You have it set to:\n"
                "            -D${OPENFCST_VAR} = ${CONDITION}   \n"
                "        Please specify the appropriate value with the  \n"
                "        CMake variable:                                \n"
                "            -D${OPENFCST_VAR}=...                     \n"
                "        Or if you are using the OpenFCST shell script  \n"
                "        use the flag:                                  \n"
                "            --${BASH_VAR} to use Debug mode else       \n"
                "            PETSc is automatically set to Release mode \n"
                "\n"
                "=======================================================\n"
                "\n"
                )
    endif()
endfunction()

#-----------------------
#
# Build Type Error Message, when condition is EITHER Release, Debug, or DebugRelease
# 
#-----------------------
# CONDITION: value of CMake variable, used to check that user typed Release, Debug, & DebugRelease
#            the way we want
# PROGRAM_NAME: string containing the name of program to be displayed in error message
# OPENFCST_VAR: string containing what the CMake variable is they can use to resolve problem
# BASH_VAR: string containing what the bash script flag is they can use to resolve problem
function(ErrMsgBuildTypeRDB CONDITION PROGRAM_NAME OPENFCST_VAR BASH_VAR)
    if(NOT (${CONDITION} STREQUAL "Release" OR 
        ${CONDITION} STREQUAL "DebugRelease" OR
        ${CONDITION} STREQUAL "Debug"))
      MESSAGE(FATAL_ERROR 
                "\n"
                "=======================================================\n"
                "\n"
                " ERROR: The ${PROGRAM_NAME} build type can only be set to:       \n"
                "            Release, DebugRelease, or Debug            \n"
                "        You have it set to:\n"
                "            -D${OPENFCST_VAR} = ${CONDITION}   \n"
                "        Please specify the appropriate value with the  \n"
                "        CMake variable:                                \n"
                "            -D${OPENFCST_VAR}=...                     \n"
                "        Or if you are using the OpenFCST shell script  \n"
                "        use the flag:                                  \n"
                "            --${BASH_VAR} to use Debug mode else       \n"
                "            PETSc is automatically set to Release mode \n"
                "\n"
                "=======================================================\n"
                "\n"
                )
    endif()
endfunction()

#-----------------------
#
# Version Error Message, when there is No Changing of CMake Variable (NCCMV) to solve problem
# 
#-----------------------
# NOTE: To give some more information when we have a case like say the compiler, if the compiler does not meet OpenFCST
#       requirements we don't have any built variables to change the gcc compiler that CMake finds so we cannot give
#       a quick solution to move on. Instead they will have to spend some time changing their system settings.

# PROGRAM_NAME: string containing the prefered program for OpenFCST, ex gcc
# TYPE_OF_PROGRAM: string containing the type of program to be displayed in error message, ex compiler
# WHAT_WAS_FOUND: what program was found by CMake, ex gcc or clang for compiler type
# FOUND_VERSION: version of program found by CMake, ex v4.8 for gcc compiler or v2.1 for clang compiler
# VERSION_MIN: minimum version that OpenFCST requires
# VERSION_PREFERRED: the prefered version for OpenFCST, ie what do us developers use?
# OPENFCST_VAR: string containing what the CMake variable is they can use to resolve problem
# BASH_VAR: string containing what the bash script flag is they can use to resolve problem
function(ErrMsgVersionNCCMV PROGRAM_NAME TYPE_OF_PROGRAM WHAT_WAS_FOUND FOUND_VERSION VERSION_MIN VERSION_PREFERRED)
    MESSAGE(FATAL_ERROR 
          "\n"
          "=======================================================\n"
          "\n"
          " ERROR: OpenFCST requires the ${PROGRAM_NAME} ${TYPE_OF_PROGRAM}. You have: \n"
          "                  ${WHAT_WAS_FOUND} v${FOUND_VERSION} \n"
          "        Please install a ${PROGRAM_NAME} ${TYPE_OF_PROGRAM} >= ${VERSION_MIN} \n"
          "        (${VERSION_PREFERRED} is preferred) and try again. \n"
          "\n"
          "=======================================================\n"
          "\n"
          )
endfunction()

#-----------------------
#
# Version Error Message, when user can Change CMake Variable (CCMV) to solve problem
# 
#-----------------------
# NOTE: To give some more information when we have a case like say Boost. Since boost is a program that the user
#       may have built by his system or even have self-compiled versions as well we have a variable to easily
#       change the directory for which version of the program that CMake will use. In this case all they need to
#       do is tell OpenFCST what the different directory is through a CMake variable change.

# PROGRAM_NAME: string containing the name of program to be displayed in error message
# VERSION_MIN: minimum version that OpenFCST requires
# FOUND_VERSION: version of program found by CMake, ex v4.8 for gcc compiler
# OPENFCST_VAR: string containing what the CMake variable is they can use to resolve problem
# BASH_VAR: string containing what the bash script flag is they can use to resolve problem
function(ErrMsgVersionCCMV PROGRAM_NAME FOUND_VERSION VERSION_MIN CMAKE_VAR BASH_VAR)
    MESSAGE(FATAL_ERROR 
          "\n"
          "=======================================================\n"
          "\n"
          " ERROR: OpenFCST requires a ${PROGRAM_NAME} version of atleast   \n"
          "        ${VERSION_MIN}, the ${PROGRAM_NAME} version found is: ${FOUND_VERSION}. \n"
          "        Please update your ${PROGRAM_NAME} library or specify a  \n"
          "        newer version's location with the CMake        \n"
          "        variable:                                      \n"
          "            -D${CMAKE_VAR}=...                            \n"
          "        Or if you are using the OpenFCST shell script  \n"
          "        use the flag:                                  \n"
          "            --${BASH_VAR}=...                          \n"
          "\n"
          "=======================================================\n"
          "\n"
          )
endfunction()

#-----------------------
#
# Find Package Error Message, when user sets a Pre-Installed Directory (PID) to use
# 
#-----------------------
# PROGRAM_NAME: string containing the name of program to be displayed in error message
# PROGRAM_FOUND: True/False bool saying if program was found or not
# PROGRAM_DIR: set by user to where to look for software
# TYPPATH: typical/default path where program is installed

# NOTE: PROGRAM_DIR & TYPPATH are used to check if user has specified another location for
#       where to look for software. If they are the same, not need to worry let CMake install software

# OPENFCST_VAR: string containing what the CMake variable is they can use to resolve problem
# BASH_VAR: string containing what the bash script flag is they can use to resolve problem
function(ErrMsgFindPackagePID PROGRAM_NAME PROGRAM_FOUND PROGRAM_DIR TYPPATH OPENFCST_VAR BASH_VAR)
    if(NOT ${PROGRAM_FOUND} AND 
          NOT ${PROGRAM_DIR} STREQUAL ${TYPPATH})
        MESSAGE(FATAL_ERROR
                  "\n"
                  "=======================================================\n"
                  "\n"
                  " ERROR: You have specified a different location for a  \n"
                  "        pre-installed ${PROGRAM_NAME} package but ${PROGRAM_NAME} was NOT  \n"
                  "        found. Please try giving a different hint.     \n"
                  "        This can be done by specifying its location    \n"
                  "        with the CMake variable:                       \n"
                  "            -D${OPENFCST_VAR}=...                            \n"
                  "        Or if you are using the OpenFCST shell script  \n"
                  "        use the flag:                                  \n"
                  "            --${BASH_VAR}=...                            \n"
                  "\n"
                  "=======================================================\n"
                  "\n"
                  )
    endif()
endfunction()

#-----------------------
#
# Find Package Error Message, when Pre-Installed Directory is Not Found (PIDNF)
# 
#-----------------------
# PROGRAM_NAME: string containing the name of program to be displayed in error message
# OPENFCST_VAR: string containing what the CMake variable is they can use to resolve problem
# BASH_VAR: string containing what the bash script flag is they can use to resolve problem
# TAR_FILE: tar file name that CMake will look for; let's user know if maybe they named the tar file wrong
function(ErrMsgFindPackagePIDNF PROGRAM_NAME OPENFCST_VAR BASH_VAR TAR_FILE)
      MESSAGE(STATUS
          "\n"
          "=======================================================\n"
          "\n"
          " WARNING: A pre-installed version of ${PROGRAM_NAME} was not     \n"
          "   found, CMake will now look for ${TAR_FILE}    \n"
          "   in the source's contrib/${PROGRAM_NAME} directory. \n"
          "   If you intended to use a pre-installed ${PROGRAM_NAME} \n"
          "   package please specify its location with the CMake  \n"
          "   variable:                                           \n"
          "       -D${OPENFCST_VAR}=...                                 \n"
          "   Or if you are using the OpenFCST shell script use   \n"
          "   the flag:                                           \n"
          "       --${BASH_VAR}=...\n"
          "\n"
          "=======================================================\n"
          "\n"
          )
endfunction()

#-----------------------
#
# Find Package Error Message, when Tar File is Not Found Warning (TFNFW)
# 
#-----------------------
# NOTE: This is different from one below because use this function when we 
#       can download from internet, in this case there still is an option
#       to install package

# PROGRAM_NAME: string containing the name of program to be displayed in error message
# TAR_FILE: tar file name that CMake will look for; let's user know if maybe they named the tar file wrong
function(ErrMsgFindPackageTFNFW PROGRAM_NAME TAR_FILE)
      MESSAGE(STATUS
            "\n"
            "=======================================================\n"
            "\n"
            " WARNING: ${TAR_FILE} was not found in the       \n"
            "   source's contrib/${PROGRAM_NAME} directory so we will now \n"
            "   download it from the internet.                      \n"
            "\n"
            "=======================================================\n"
            "\n"
            )
endfunction()

#-----------------------
#
# Find Package Error Message, when Tar File is Not Found Error (TFNFE)
# 
#-----------------------
# NOTE: This is different from one above because if we cannot download it from
#       internet (like with deal.II) then throw and error so they have to 
#       give the tar file OR use the CMake/Bash variables to give a pre-installed
#       package.

# PROGRAM_NAME: string containing the name of program to be displayed in error message
# TAR_FILE: tar file name that CMake will look for; let's user know if maybe they named the tar file wrong
# OPENFCST_VAR: string containing what the CMake variable is they can use to resolve problem
# BASH_VAR: string containing what the bash script flag is they can use to resolve problem
function(ErrMsgFindPackageTFNFE PROGRAM_NAME TAR_FILE OPENFCST_VAR BASH_VAR)
      MESSAGE(FATAL_ERROR
                  "\n"
                  "=======================================================\n"
                  "\n"
                  " ERROR: ${TAR_FILE} was not found, cannot install      \n"
                  "        ${PROGRAM_NAME}. Please install ${PROGRAM_NAME}\n"
                  "        or try giving a different hint. This can be    \n"
                  "        done by specifying its location with the CMake \n"
                  "        variable:                                      \n"
                  "            -D${OPENFCST_VAR}=...                      \n"
                  "        Or if you are using the OpenFCST shell script  \n"
                  "        use the flag:                                  \n"
                  "            --${BASH_VAR}=...                          \n"
                  "\n"
                  "=======================================================\n"
                  "\n"
                  )
endfunction()

#-----------------------
#
# User Passes Flag They Shouldnt Error Message
# 
#-----------------------
# NOTE: This is to check two CMake Variables that the user sets, if there is a conflict return
#       error message.

# CONDITION: value of CMake variable, used to check that user typed Release, Debug, & DebugRelease
#            the way we want
# PROGRAM_NAME: string containing the name of program to be displayed in error message
# OPENFCST_VAR: string containing what the CMake variable is they can use to resolve problem
# BASH_VAR: string containing what the bash script flag is they can use to resolve problem
function(ErrMsgConflictVars VAR1 VAR2 CONDITION1 CONDITION2 OPENFCST_VAR1 OPENFCST_VAR2 BASH_VAR1 BASH_VAR2)
    if(${VAR1} STREQUAL ${CONDITION1} AND
         ${VAR2} STREQUAL ${CONDITION2})
      MESSAGE(FATAL_ERROR
                "\n"
                "=======================================================\n"
                "\n"
                " ERROR: You have specified:                            \n"
                "           ${OPENFCST_VAR1} = ${VAR1}                  \n"
                "           ${OPENFCST_VAR2} = ${VAR2}                  \n"
                "        These variables cannot be set at the same time.\n"
                "        Please turn off one of the two variables and   \n"
                "        re-install. If using the OpenFCST shell script \n"
                "        they are controlled with the respective flags: \n"
                "            --${BASH_VAR1}                             \n"
                "            --${BASH_VAR2}                             \n"
                "\n"
                "=======================================================\n"
                "\n"
                )
    endif()
endfunction()