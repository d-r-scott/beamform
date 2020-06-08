#!/bin/bash

# Command line arguments
FRB=$1

# NOTE: BOTH POLARISATIONS (X AND Y) WILL BE PROCESSED

logdir=./log          # Directory for log files
logpre=${logdir}/${FRB}

source FRBdata.sh $FRB $a_or_m $pol
echo "FRB$FRB"
echo "DM=		$DM"
echo "f0=		$f0"

# Stage 4: Dedispersion
out4=${logpre}_stage4.out
args4_x="$FRB x $DM $f0"
args4_y="$FRB y $DM $f0"

echo "sbatch --output=$out4 --error=$out4 stage4_dedispersion.sh $args4_x"
jobid4_x=$(sbatch --output=$out4 --error=$out4 stage4_dedispersion.sh $args4_x | cut -d " " -f 4)

echo "sbatch --output=$out4 --error=$out4 stage4_dedispersion.sh $args4_y"
jobid4_y=$(sbatch --output=$out4 --error=$out4 stage4_dedispersion.sh $args4_y | cut -d " " -f 4)

# Stage 5: IFFT
out5=${logpre}_stage5.out
args5_x="$FRB x $DM"
args5_y="$FRB y $DM"

echo "sbatch --output=$out5 --error=$out5 --dependency=afterok:$jobid4_x stage5_ifft.sh $args5_x"
jobid5_x=$(sbatch --output=$out5 --error=$out5 --dependency=afterok:$jobid4_x stage5_ifft.sh $args5_x | cut -d " " -f 4)

echo "sbatch --output=$out5 --error=$out5 --dependency=afterok:$jobid4_y stage5_ifft.sh $args5_y"
jobid5_y=$(sbatch --output=$out5 --error=$out5 --dependency=afterok:$jobid4_y stage5_ifft.sh $args5_y | cut -d " " -f 4)

# Stage 6: Generate dynamic spectra and calculate Stokes parameters
out5=${logpre}_stage6.out
args6="$FRB $DM"

echo "sbatch --output=$out6 --error=$out6 --dependency=afterok:$jobidy5_x,$jobid5_y stage6_dynspecs.sh $args6"
sbatch --output=$out6 --error=$out6 --dependency=afterok:$jobidy5_x,$jobid5_y stage6_dynspecs.sh $args6
