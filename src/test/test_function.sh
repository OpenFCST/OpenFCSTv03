###########################################################################################
#
# Function used by test applications that compare different solution files. The solution is
# stored in two files:
#  - solution_ref.dat : Reference solution for regression test
#  - solution.dat : Actual solution from running the program.
#
#
# The function takes two arguments:
# - Argument 1: Name of the test you are running. It must be a single string
# - Argument 2: Relative or absolute path to Install/ folder
#
# Usage:
#
# In your test script first load this script using
#   source ./test_function.sh
# 
# Then, create the variables and call the function:
#   test_name="Backward_step"
#   path="../../../../.."
#   test_function $test_name $path
#
# Author: M. Secanell, 2015
#
###########################################################################################

test_function () {

test_name=$1
test_script="$2/test/simpleCompare.py"
log_name="$2/tests_summary.log"
ref_data_name="$3/solution_ref.dat"

if [ "${PIPESTATUS[0]}" != "0" ]; then
    echo                                                                                   2>&1 | tee          tests_summary.log
    echo "Results summary from the $test_name test:"                                2>&1 | tee --append tests_summary.log
    echo                                                                                   2>&1 | tee --append tests_summary.log
    echo "-------------------------------------------------------------------------------" 2>&1 | tee --append tests_summary.log
    echo "The simulation did not run correctly. Please review the tests_output.log file  " 2>&1 | tee --append tests_summary.log
    echo "-------------------------------------------------------------------------------" 2>&1 | tee --append tests_summary.log
    echo                                                                                   2>&1 | tee --append tests_summary.log
    cat tests_summary.log >> $log_name
    exit 2
else
    python $test_script $ref_data_name solution.dat
    error_status=${PIPESTATUS[0]}
    if [ "$error_status" != "0" ] && [ "$4" == "NoPETSc" ]; then
        echo                                                                                   2>&1 | tee          tests_summary.log
        echo "Results summary from the $test_name test:"                                2>&1 | tee --append tests_summary.log
        echo                                                                                   2>&1 | tee --append tests_summary.log
        echo "-------------------------------------------------------------------------------" 2>&1 | tee --append tests_summary.log
        echo "Results from the test do not match expected results"                             2>&1 | tee --append tests_summary.log
        echo "Please check the results in solution_ref.dat against that of solution.dat"       2>&1 | tee --append tests_summary.log
        echo "NOTE: This test cannot pass when using the solver MUMPS with adaptive"           2>&1 | tee --append tests_summary.log
        echo "refinement; i.e. for this test to pass (or to use this application) do NOT use"  2>&1 | tee --append tests_summary.log
        echo "the flag --with-petsc."                                                          2>&1 | tee --append tests_summary.log
        echo "-------------------------------------------------------------------------------" 2>&1 | tee --append tests_summary.log
        echo                                                                                   2>&1 | tee --append tests_summary.log
        cat tests_summary.log >> $log_name
        exit 2
    elif [ "${error_status}" != "0" ]; then
        echo                                                                                   2>&1 | tee          tests_summary.log
        echo "Results summary from the $test_name test:"                                2>&1 | tee --append tests_summary.log
        echo                                                                                   2>&1 | tee --append tests_summary.log
        echo "-------------------------------------------------------------------------------" 2>&1 | tee --append tests_summary.log
        echo "Results from the test do not match expected results"                             2>&1 | tee --append tests_summary.log
        echo "Please check the results in solution_ref.dat against that of solution.dat"       2>&1 | tee --append tests_summary.log
        echo "-------------------------------------------------------------------------------" 2>&1 | tee --append tests_summary.log
        echo                                                                                   2>&1 | tee --append tests_summary.log
        cat tests_summary.log >> $log_name
    exit 2
    else
        echo                                                    2>&1 | tee          tests_summary.log
        echo "Results summary from the $test_name test:" 2>&1 | tee --append tests_summary.log
        echo                                                    2>&1 | tee --append tests_summary.log
        echo "--------------------------------------------"     2>&1 | tee --append tests_summary.log
        echo "Results from the test match expected results"     2>&1 | tee --append tests_summary.log
        echo "--------------------------------------------"     2>&1 | tee --append tests_summary.log
        echo                                                    2>&1 | tee --append tests_summary.log
        cat tests_summary.log >> $log_name
        exit 0
    fi
fi

}