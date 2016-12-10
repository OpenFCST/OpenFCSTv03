#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

This scripts prepares configures a polarization curve and runs it.
OPTIONS:
"-h"     - Prints the help for this scripts
"-r"     - executes the regression test
"-c"     - Configure only. Generate parameter scripts and exits.
"-i 0.1" - Maximal increment on the polcurve
"-u 0.35"- Lower limit of the cell voltage
"-a 2"    - Number of adaptive refinement steps

Example; polcurve with high resolution and mass transport knee:
$0 -i 0.05 -u 0.1
EOF
}

# Parse commandline

flag_regression=0
flag_cleanup=0

while getopts "hrci:u:a:" OPTION
do
  case $OPTION in
    h)
      usage
      exit 1
      ;;
    c)
      flag_cleanup=1
      ;;
    r)
      flag_regression=1
      ;;
    i)
      size_step=$OPTARG
      echo: "Incement in the polcurve $size_step"
      ;;
    u)
      u_l=$OPTARG
      echo: "Lower potential point in the polcurve set to $u_l [V]"
      ;;
    a)
      num_adapt=$OPTARG
      echo: "Number of adaption steps set to $num_adapt"
      ;;
    ?)
      usage
      exit 1
      ;;
  esac
done

# Preset defaults:
if [ -z "${$size_step}" ]; then
    $size_step=0.1
fi
if [ -z "${num_adapt}" ]; then
    num_adapt=2
fi
if [ -z "${u_l}" ]; then
    u_l=0.4
fi

# Define name of test:
test_name="PemfcNIThermal>>PolarizationCurve"

# Retrieve test directory (directory of file location):
test_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Enter testing folder
cd $test_dir

# Cleanup
rm polarization_curve.dat
if [ "$flag_cleanup" = 1 ]; then
  rm *.log
  exit 1
fi
  
# Directory configuration
cp ../template/data.prm data.prm
# I need to figure out the slahs problem
#sed -i "s/include ..//template//data.prm/include/.//data_param.prm/g" data.prm
sed -i "s/  set Number of Refinements                = 3/  set Number of Refinements                = $num_adapt/g" data.prm
sed -i "s/    set Final voltage [V] = 0.35/    set Final voltage [V] = $size_step/g" main.prm
# Load the test function:
source "$FCST_DIR/test/test_function_polarization.sh"

# Call application:
$FCST_DIR/bin/fuel_cell-2d.bin main.prm

# Run test function:
test_function_polarization $test_name $FCST_DIR