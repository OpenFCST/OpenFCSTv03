# -----------------------------------------------------------  
# -- Get environment
# -----------------------------------------------------------  

## -- Set hostname
## --------------------------
find_program(HOSTNAME_CMD NAMES hostname)
exec_program(${HOSTNAME_CMD} ARGS OUTPUT_VARIABLE HOSTNAME)

set(CTEST_SITE                          "${HOSTNAME}")

## -- Set site / build name
## --------------------------

find_program(UNAME NAMES uname)
macro(getuname name flag)
  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(getuname)

getuname(osname -s)
getuname(osrel  -r)
getuname(cpu    -m)

set(CTEST_BUILD_NAME                    "${osname}-${cpu}-prod")

## -- SVN command
## ----------------
find_program(CTEST_SVN_COMMAND NAMES svn)

## -- make command
## -----------------
find_program(MAKE NAMES make)

# -----------------------------------------------------------  
# -- build specific
# -----------------------------------------------------------  

set(MODEL                               "Continuous")

## -- DashBoard Root
set(CTEST_DASHBOARD_ROOT                "$ENV{PWD}")

## -- SRC Dir
set(CTEST_SOURCE_DIRECTORY              "${CTEST_DASHBOARD_ROOT}/test")

## -- BIN Dir                                            
set(CTEST_BINARY_DIRECTORY              "${CTEST_SOURCE_DIRECTORY}/build-${CTEST_BUILD_NAME}") 

## -- Build options
##set(OPTION_BUILD                        "-j16")

# -----------------------------------------------------------  
# -- commands
# -----------------------------------------------------------  

## -- Update Command
set(CTEST_UPDATE_COMMAND               "${CTEST_SVN_COMMAND}")

## -- Configure Command
set(CTEST_CONFIGURE_COMMAND            "./fcst_install --with-configure=off") 

## -- Build Command
set(CTEST_BUILD_COMMAND                "./fcst_install --library=fcst")

# -----------------------------------------------------------  
# -- Configure CTest
# -----------------------------------------------------------  

# -- Write flags to /test/config.txt for CTestTestfile.cmake to read them
file(APPEND ${CTEST_DASHBOARD_ROOT}/test/config.txt "PETSc@OPENFCST_WITH_PETSC@" "\n")  # Turn PETSc option on
file(APPEND ${CTEST_DASHBOARD_ROOT}/test/config.txt "DAKOTA@OPENFCST_WITH_DAKOTA@" "\n") # Turn Dakota option on

## -- read CTestCustom.cmake file
ctest_read_custom_files("${CTEST_SOURCE_DIRECTORY}")

# -----------------------------------------------------------  
# -- Settings
# -----------------------------------------------------------  

## -- Process timeout in seconds
set(CTEST_TEST_TIMEOUT           "14400")

## -- Set output to english
set( $ENV{LC_MESSAGES}      "en_EN" )

# -----------------------------------------------------------  
# -- Run CTest
# -----------------------------------------------------------  

##----Remove comments made using a single hash key in order to
##----get CTest to test that component.

## -- Start
message(" -- Start dashboard ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_start(${MODEL} TRACK ${MODEL})

## -- Update
#message(" -- Update ${MODEL} - ${CTEST_BUILD_NAME} --")
#ctest_update(           SOURCE "${CTEST_DASHBOARD_ROOT}" RETURN_VALUE res)

## -- Configure, this will actually configure and compile the entire code
#message(" -- Configure ${MODEL} - ${CTEST_BUILD_NAME} --")
#ctest_configure(BUILD  "${CTEST_DASHBOARD_ROOT}" RETURN_VALUE res)

## -- BUILD, this will actually compile fcst only
#message(" -- Build ${MODEL} - ${CTEST_BUILD_NAME} --")
#ctest_build(    BUILD  "${CTEST_DASHBOARD_ROOT}" RETURN_VALUE res)



## -- TEST
message(" -- Test ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_test(     BUILD  "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)

## -- SUBMIT
##message(" -- Submit ${MODEL} - ${CTEST_BUILD_NAME} --")
##ctest_submit(                                              RETURN_VALUE res)

message(" -- Finished ${MODEL}  - ${CTEST_BUILD_NAME} --")
