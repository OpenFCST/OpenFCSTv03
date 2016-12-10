#!/bin/bash

######################################################################
# Run a simulation with a specified file out of any directory
# Add --  alias fcst_run='sh /home/andput/AFCCEXT/fcst/trunk/fcst/run_simulation'  -- to .alias or .bashrc
#
# This script will take a simulation input file in the command line 
#	or prompt for one if none is entered.
######################################################################

source fcst_env.sh

run_directory=`pwd`

sim_file=`echo "$1"`

if [ "$sim_file" = "" ]; then
read -p 'Enter simulation parameter file:' sim_file
fi

echo -e "Running FCST in $run_directory with $sim_file \n" 
nohup fuel_cell-2d.bin $sim_file &> fcst.log &
PID=$!

echo -e "Running FCST with PID $PID"
echo $PID > fcst.pid
