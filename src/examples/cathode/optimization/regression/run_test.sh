#!/bin/bash

##############################
#
# Note: This script is expected to be run from /test
# folder, not its current location.
#
# 
# If you want to check current directory: echo ${PWD}
#
##############################

# Define name of test:
test_name="AppCathode>>Optimization"

# Define relative (or absolute path to Install) from testing folder above:
path=$FCST_DIR

# Enter testing folder
cd $path/examples/cathode/optimization
rm dakota_tabular.dat

#Get number of cores to use
argnumcores=`echo "$*" | perl -n -e 'm/--cores=(\S+)/; print $1'`
if [ "$argnumcores" == "" ];then
  argnumcores=1
fi

# Load the test function:
. $path/test/test_function_optimization.sh

# Call application:
mpirun -np $argnumcores $path/bin/fuel_cell-2d.bin main.prm

# Run test function:
test_function_optimization $test_name $path


