#!/bin/bash
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

# FRB details (for other FRBs, refer to the README file)

#echo "FRB181112"
#offset=2106003
#DM=589.25
#f0=1297.5
## Calibration solution directories
#basedr=/fred/oz002/users/hcho/craft/
#calcfile=${basedr}Calibration/aipscal/frb181112/c1_f0/craftfrb.im # geometric delays
#fcm=${basedr}Calibration/aipscal/frb181112/fcm.txt # clock delays
#hwfile=${basedr}Calibration/mircal/frb181112/hwdelays_0407_SB7031_20181112222901_round-8.txt # hardware delays
#mir=${basedr}Calibration/mircal/frb181112/20181112222901_call_beam01_i4096_f9.uvlin # MIRIAD gain, bandpass
#aips=${basedr}Calibration/aipscal/frb181112/bandpass.bp.txt # AIPS gain, bandpass
## VCRAFT file directory - change
#f_vcraft=${basedr}python/voltages/FRB181112/ak**/beam01/*.vcraft
## output file directories
#f_outfile="./test_output/f_y_an${an}.npy" # part 1 output (frequency domain)
#t_outfile="./test_output/t_y_an${an}.npy" # part 2 output (time doamin)

# Above is made obsolete by FRBdata.sh
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
if [ $# != 0 ]
then
args+=("--an ${an}")
fi
args+=("--calcfile $calcfile")
args+=("--parset $fcm")
args+=("--hwfile $hwfile")
#args+=("--mirsolutions $mir") # optional, comment out if you want to use AIPS phase corrections
args+=("--aips_c $aips") # optional, comment out if you want to use MIRIAD bandpass corrections
args+=("-o $f_outfile")

# PART 1
python craftcor.py ${args[@]} --tab $f_vcraft


# PART 2
fftlen=$(( $n*64 ))
python freq2time.py -f $f_outfile -d $DM --f0 $f0 -o $t_outfile -l $fftlen -q
