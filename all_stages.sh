#!/bin/bash

# Command line arguments
FRB=$1      # FRB name
a_or_m=$2   # AIPS or MIRIAD solutions for calibration
# Both polarisations are executed
n=$3

# Execute stage123.sh for x polarisation, save jobid3 as jobid3_x, then execute for y.
# This allows us to provide jobid3_x and jobid3_y to stage456.sh and have all stages
# submitted at once.
source stage123.sh $FRB $a_or_m x $n
jobid3_x=$jobid3

source stage123.sh $FRB $a_or_m y $n
jobid3_y=$jobid3

# Launch stage456.sh, providing the optional jobid3_x and jobid3_y arguments
./stage456.sh $FRB $n $jobid3_x $jobid3_y
