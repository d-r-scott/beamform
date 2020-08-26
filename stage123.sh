#!/bin/bash

# Command line arguments
FRB=$1      # FRB name
a_or_m=$2   # AIPS or MIRIAD solutions for beamforming
pol=$3      # Polarisation (x or y)

# KEY DIFFERENCE BETWEEN THIS VERSION AND HYERIN'S ORIGINAL VERSION:
# Antenna number is not specified, all antennas are always processed

logdir=./log          # Directory for log files
logpre=${logdir}/${FRB}_${pol}

# Create data directories
basedir=./output
outdir=${basedir}/${FRB}
f_outdir=${outdir}/f

if [ ! -d $outdir ]; then
  mkdir $outdir
fi

if [ ! -d $f_outdir ]; then
  mkdir $f_outdir
fi

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
echo "i=  $i"
echo "n=  $n"
echo "n_ant=  $n_ant"

# Stage 1: Per-antenna beamforming

args1="$FRB $a_or_m $pol $offset $calcfile $fcm"
if [ "$a_or_m" == "AIPS" ]; then
  args1="$args1 $aips"
elif [ "$a_or_m" == "MIRIAD" ]; then
  args1="$args1 $mir"
else
  echo "ERROR: Must provide AIPS or MIRIAD exactly!"
  echo "Exiting..."
  exit
fi
# We'll add hwfile in inside stage1_beamforming.sh

# Processing parameters
fftlen=$(( $n * 64 ))

args1="$args1 $i $n"

max_ant=$(( n_ant - 1 ))
jobid1=""
# TODO: Need to put this in a for loop for all antennae
for ant in `seq 0 $max_ant`; do
  out1=${logpre}_stage1_${ant}.out
  echo "$args1"
  echo "sbatch --output=$out1 --error=$out1 stage1_beamforming.sh $args1 $f_vcraft $ant $hwfile"
  new_jobid=$(sbatch --output=$out1 --error=$out1 stage1_beamforming.sh $args1 "$f_vcraft" $ant $hwfile | cut -d " " -f 4)
  jobid1="$jobid1:$new_jobid"
done

# Stage 2: Summing fine channel spectra
out2=${logpre}_stage2.out
args2="$FRB $pol"

echo "sbatch --output=$out2 --error=$out2 --dependency=afterok$jobid1 stage2_summing.sh $args2"
jobid2=$(sbatch --output=$out2 --error=$out2 --dependency=afterok$jobid1 stage2_summing.sh $args2 | cut -d " " -f 4)

# Stage 3: Derippling
out3=${logpre}_stage3.out
args3="$FRB $pol $fftlen"   # fftlen was exported by stage1_beamforming.sh

echo "sbatch --output=$out3 --error=$out3 --dependency=afterok:$jobid2 stage3_derippling.sh $args3"
jobid3=$(sbatch --output=$out3 --error=$out3 --dependency=afterok:$jobid2 stage3_derippling.sh $args3 | cut -d " " -f 4)

# We set jobid3 above so that it can be used by all_stages.sh to chain stages 1, 2, 3 with 4, 5, 6
