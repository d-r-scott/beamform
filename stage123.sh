#!/bin/bash

# Command line arguments
FRB=$1      # FRB name
a_or_m=$2   # AIPS or MIRIAD solutions for correlation
pol=$3      # Polarisation (x or y)
n_ant=$4    # Number of antennae there is for this FRB

# KEY DIFFERENCE BETWEEN THIS VERSION AND HYERIN'S ORIGINAL VERSION:
# Antenna number is not specified, all antennas are always processed

logdir=./log          # Directory for log files
logpre=${logdir}/${FRB}_${pol}

source FRBdata.sh $FRB $a_or_m $pol
echo "FRB$FRB"
echo "offset=		$offset"
echo "DM=		$DM"
echo "f0=		$f0"
echo "calcfile=	$calcfile"
echo "fcm=		$fcm"
echo "hwfile=	$hwfile"
echo "aips=		$aips"
echo "mir=	$mir"
echo "f_vcraft=	$f_vcraft"

# Stage 1: Per-antenna correlation
out1=${logpre}_stage1.out

args1="$FRB $a_or_m $pol $offset $calcfile $fcm $f_vcraft"
if [ "$a_or_m" == "AIPS" ]; then
  args1="$args1 $aips"
elif [ "$a_or_m" == "MIRIAD" ]; then
  args1="$args1 $mir"
else
  echo "ERROR: Must provide AIPS or MIRIAD exactly!"
  echo "Exiting..."
  exit
fi
# We'll add hwfile in inside stage1_correlation.sh

max_ant=$(( n_ant - 1 ))
jobid1=""
# TODO: Need to put this in a for loop for all antennae
for ant in `seq 0 $max_ant`; do
  echo "sbatch --output=$out1 --error=$out1 stage1_correlation.sh $args1 $ant $hwfile"
  new_jobid=$(sbatch --output=$out1 --error=$out1 stage1_correlation.sh $args1 $ant $hwfile | cut -d " " -f 4)
  jobid1="$jobid1:$new_jobid"
done

# Stage 2: Summing fine channel spectra
out2=${logpre}_stage2.out
args2="$FRB $pol"

echo "sbatch --output=$out2 --error=$out2 --dependency=afterok$jobid1 stage2_summing.sh $args2"
jobid2=$(sbatch --output=$out2 --error=$out2 --dependency=afterok$jobid1 stage2_summing.sh $args2 | cut -d " " -f 4)

# Stage 3: Derippling
out3=${logpre}_stage3.out
args3="$FRB $pol $fftlen"   # fftlen was exported by stage1_correlation.sh

echo "sbatch --output=$out3 --error=$out3 --dependency=afterok:$jobid2 stage3_derippling.sh $args3"
sbatch --output=$out3 --error=$out3 --dependency=afterok:$jobid2 stage3_derippling.sh $args3
