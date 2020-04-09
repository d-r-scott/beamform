#!/bin/bash
# ---------------------------------------------------------------------------------
# Original Author: Hyerin Cho
# Last Updated by: David Scott	[david.r.scott@postgrad.curtin.edu.au]
# ---------------------------------------------------------------------------------
# Shell script to reconstruct a localized ASKAP FRB to its full time resolution (~3 ns)
# This script consists of 2 parts:
#   PART 1) Create high resolution 1D array in frequency domain with delay corrections
#   PART 2) Remove ripples generated at PFB stage, coherently dedisperse, and ifft back to time series

# NOTE : before doing all the reconstruction, it is useful to check the calibration solutions first with 0407 data
# NOTE : i=1, n=16384 takes a long time, so it is recommended to run parallel for each antenna and then sum up the time series at the very last step
# (for each antenna: max. time ~ 9 hrs, max. memory ~ 30g)
# TODO: maybe consider overlap-save method in the future
# Please refer to the google docs README for more information: https://docs.google.com/document/d/1uSbsudW5tkjy-GbJaCduERS36AgCEgWox-V8vY-YdA8/edit?usp=sharing
# ---------------------------------------------------------------------------------

# USAGE: bash fulltimeres.sh [antenna] [FRB ID] [AIPS/MIRIAD] [Polarisation (x/y)]
if (( $# != 4 )); then
	echo "USAGE: bash fulltimeres.sh [antenna] [FRB ID] [AIPS/MIRIAD] [Polarisation (x/y)]"
	exit
fi

## Specify one antenna number, or a - if you want all antennas
# TODO: Implement the above
an=$1
FRB=$2		# FRB to do
a_or_m=$3	# AIPS or MIRIAD
pol=$4		# polarisation (x or y)

args=()
## set array size and offset
i=1
n=16384 # $n * 54 * 336 is the total length of the output array

source FRBdata.sh $FRB $a_or_m $pol

echo $FRB
echo $offset
echo $DM
echo $f0
echo $calcfile
echo $fcm
echo $hwfile
echo $aips
echo $mir
echo $f_vcraft
echo $f_outfile
echo $t_outfile

exit


args+=("-i $i")
args+=("-n $n")
args+=("--offset $offset")
if [ "$an" != "-" ]; then
	args+=("--an ${an}")
fi
args+=("--calcfile $calcfile")
args+=("--parset $fcm")
args+=("--hwfile $hwfile")
if [ "$a_or_m" == "AIPS" ]; then
	args+=("--aips_c $aips") # optional, comment out if you want to use MIRIAD bandpass corrections
elif [ "$a_or_m" == "MIRIAD" ]; then
	args+=("--mirsolutions $mir") # optional, comment out if you want to use AIPS phase corrections
else
	echo "ERROR: Must provide AIPS or MIRIAD exactly!"
	echo "Exiting..."
	exit
fi
args+=("-o $f_outfile")

# PART 1
python craftcor.py ${args[@]} --tab $f_vcraft


# PART 2
fftlen=$(( $n*64 ))
python freq2time.py -f $f_outfile -d $DM --f0 $f0 -o $t_outfile -l $fftlen -q
