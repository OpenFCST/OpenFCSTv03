cd ./unit_tests/
rm logfile.log
../../bin/fuel_cell-2d.bin main_unit_tests.prm

if grep -q "Unit Tests Successfully Completed" "logfile.log"; then
    echo  2>&1 | tee tests_summary.log
    echo  "Results summary from the units test:" 2>&1 | tee --append tests_summary.log
    echo 2>&1 | tee --append tests_summary.log
    echo "--------------------------------------------" 2>&1 | tee --append tests_summary.log
    echo  "Unit Tests Successfully Completed" 2>&1 | tee --append tests_summary.log
    echo "--------------------------------------------" 2>&1 | tee --append tests_summary.log
    echo 2>&1 | tee --append tests_summary.log
    cp tests_summary.log ../../tests_summary.log 
    exit 0
else
    echo  2>&1 | tee tests_summary.log
    echo  "Results summary from the units test:" 2>&1 | tee --append tests_summary.log
    echo 2>&1 | tee --append tests_summary.log
    echo "--------------------------------------------" 2>&1 | tee --append tests_summary.log
    echo  "Unit Tests Completed with Errors!" 2>&1 | tee --append tests_summary.log
    echo "--------------------------------------------" 2>&1 | tee --append tests_summary.log
    echo 2>&1 | tee --append tests_summary.log
    cp tests_summary.log ../../tests_summary.log 
    exit 2
fi
