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
test_name="AppDiffusion3D>>Regression"

# Define relative (or absolute path to Install) from testing folder above:
path=$FCST_DIR

# Enter testing folder
cd $path/examples/microscale/test
rm solution.dat

# Load the test function:
. $path/test/test_function.sh

# Call application:
$path/bin/fuel_cell-3d.bin main.prm

# Run test function:
test_function $test_name $path ./regression