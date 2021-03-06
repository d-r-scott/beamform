#!/bin/bash

# Command line arguments
FRB=$1
# Following two are optional, only used when all_stages.sh is run.
# They are used to allow all stages to be submitted at once.
jobid3_x=$3
jobid3_y=$4

# NOTE: BOTH POLARISATIONS (X AND Y) WILL BE PROCESSED

logdir=./log          # Directory for log files
logpre=${logdir}/${FRB}

source FRBdata.sh $FRB $a_or_m $pol
echo "FRB$FRB"
echo "DM= $DM"
echo "f0= $f0"
n=$2
# Stage 4: Dedispersion
out4=${logpre}_stage4.out
args4_x="$FRB x $DM $f0 $n"
args4_y="$FRB y $DM $f0 $n"

# Check if we are provided with jobid3_x and jobid4_x, and if so, set them as dependencies for stage4.
# If they're not provided, stage4 has no dependencies
if [ "$jobid3_x" != "" ] && [ "$jobid3_y" != "" ]; then
  dependency4_x="--dependency=afterok:$jobid3_x"
  dependency4_y="--dependency=afterok:$jobid3_y"
else
  dependency4_x=""
  dependency4_y=""
fi

echo "sbatch --output=$out4 --error=$out4 $dependency4_x stage4_dedispersion.sh $args4_x"
jobid4_x=$(sbatch --output=$out4 --error=$out4 $dependency4_x stage4_dedispersion.sh $args4_x | cut -d " " -f 4)

echo "sbatch --output=$out4 --error=$out4 $dependency4_y stage4_dedispersion.sh $args4_y"
jobid4_y=$(sbatch --output=$out4 --error=$out4 $dependency4_y stage4_dedispersion.sh $args4_y | cut -d " " -f 4)

# Stage 5: IFFT
out5=${logpre}_stage5.out
args5_x="$FRB x $DM $n"
args5_y="$FRB y $DM $n"

echo "sbatch --output=$out5 --error=$out5 --dependency=afterok:$jobid4_x stage5_ifft.sh $args5_x"
jobid5_x=$(sbatch --output=$out5 --error=$out5 --dependency=afterok:$jobid4_x stage5_ifft.sh $args5_x | cut -d " " -f 4)

echo "sbatch --output=$out5 --error=$out5 --dependency=afterok:$jobid4_y stage5_ifft.sh $args5_y"
jobid5_y=$(sbatch --output=$out5 --error=$out5 --dependency=afterok:$jobid4_y stage5_ifft.sh $args5_y | cut -d " " -f 4)

# Stage 6: Generate dynamic spectra and calculate Stokes parameters
out6=${logpre}_stage6.out
args6="$FRB $DM $n"

echo "sbatch --output=$out6 --error=$out6 --dependency=afterok:$jobid5_x:$jobid5_y stage6_dynspecs.sh $args6"
sbatch --output=$out6 --error=$out6 --dependency=afterok:$jobid5_x:$jobid5_y stage6_dynspecs.sh $args6
