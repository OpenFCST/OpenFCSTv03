#!/bin/bash

############################################################################################################
# This file will run the app_diffusion for 2D  test that is contained in the
# data/app_diffusion/app_diffusion/testing/2D folder.
# The script is run by ctest, keep in mind that the results that would normally be printed to screen will be
# suppressed by ctest # and printed to its own file.

# The script will (by line number):
# - first navigate to the folder where the test is
# - will run the code with the correct data files.
# - The result from the simulation is queried to see if it ran without error by checking ${PIPESTATUS[0]}.
#   A zero means the test ran without error.
# - If a non-zero result is returned by the simulation, the code will print out a message saying that there
#   was an error. As this will also be suppressed by the code, the message is also printed to the
#   tests_summary.log file. The first line containing a 'tee' command will create the file, subsequent
#   calls will append their output to the end of the file, so as not to overwrite it.
# - Before exiting, the test_summary.log file is copied to the fcst main folder where it will be opened by the
#   run_tests script.
# - If the simulation ran correctly, then the results from the simulation are compared to expected results.
#   Both sets of results are stored in a texts files, appended with .dat. The simpleCompare.py script is
#   a python file that will read in the two files and compare them. If they are not within a reasonable
#   agreement the python script will return a non-zero and an error is printed. If they are in reasonable
#   agreement, a zero is returned indicated that all is well and a message is printed. Again the message
#   is captured by ctest, so it is also appended to the tests_summary.log file.
############################################################################################################
cd ../data/app_diffusion/testing/2D

rm solution.dat

../../../../bin/fuel_cell-2d.bin main.prm
echo  "called!"
if [ "${PIPESTATUS[0]}" != "0" ]; then
    echo  2>&1 | tee tests_summary.log
    echo  "Results summary from the app_diffusion_2D test:" 2>&1 | tee --append tests_summary.log
    echo  2>&1 | tee --append tests_summary.log
    echo "-------------------------------------------------------------------------------" 2>&1 | tee --append tests_summary.log
    echo "The simulation did not run correctly. Please review the tests_output.log file  "  2>&1 | tee --append tests_summary.log
    echo "-------------------------------------------------------------------------------" 2>&1 | tee --append tests_summary.log
    echo  2>&1 | tee --append tests_summary.log
    cat tests_summary.log >> ../../../../tests_summary.log
    exit 2
else
    python ../../../../test/simpleCompare.py solution_ref.dat solution.dat
    if [ "${PIPESTATUS[0]}" != "0" ]; then
        echo  2>&1 | tee tests_summary.log
        echo  "Results summary from the app_diffusion_2D test:" 2>&1 | tee --append tests_summary.log
        echo 2>&1 | tee --append tests_summary.log
        echo "-------------------------------------------------------------------------------" 2>&1 | tee --append tests_summary.log
        echo "Results from the test do not match expected results" 2>&1 | tee --append tests_summary.log
        echo "Please check the results in solution_ref.dat against that of solution.dat" 2>&1 | tee --append tests_summary.log
        echo "-------------------------------------------------------------------------------" 2>&1 | tee --append tests_summary.log
        echo 2>&1 | tee --append tests_summary.log
        cat tests_summary.log >> ../../../../tests_summary.log
        exit 2
    else
        echo  2>&1 | tee tests_summary.log
        echo  "Results summary from the app_diffusion_2D test:" 2>&1 | tee --append tests_summary.log
        echo 2>&1 | tee --append tests_summary.log
        echo "--------------------------------------------" 2>&1 | tee --append tests_summary.log
        echo "Results from the test match expected results" 2>&1 | tee --append tests_summary.log
        echo  "--------------------------------------------" 2>&1 | tee --append tests_summary.log
        echo 2>&1 | tee --append tests_summary.log
        cat tests_summary.log >> ../../../../tests_summary.log
        exit 0
    fi
fi