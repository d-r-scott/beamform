#!/bin/bash --login

#SBATCH --job-name=beamforming
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=0:30:00
#SBATCH --mem=128g

# Command line arguments
FRB=$1      # FRB name
a_or_m=$2   # AIPS or MIRIAD solutions for correlation
pol=$3      # Polarisation (x or y)
offset=$4   # Data offset (how many microseconds at the start to skip)
calcfile=$5
fcm=$6
a_m_file=$7 # AIPS/MIRIAD file (both are in this variable, the contents depends on $a_or_m
i=$8
n=$9
f_vcraft=${10}
ant=${11}      # Antenna number
hwfile=${12}   # Hardware delays. Probably not there for newer FRBs.

# Set data directories - stage123.sh has already ensured they exist
basedir=./output
outdir=${basedir}/${FRB}_n${n}
f_outdir=${outdir}/f

# Get modules to load and load them
source modules.sh
module load $modules_1

outfile=${f_outdir}/${FRB}_${ant}_${pol}_f.npy
echo $outfile
echo "$offset"

args="-i $i -n $n --offset $offset --calcfile $calcfile --parset $fcm"
echo "$args"
if [ "$a_or_m" == "AIPS" ]; then
  args="$args --aips_c $a_m_file"
else
  args="$args --mirsolutions $a_m_file"
fi
if [ "$hwfile" != "" ]; then
  args="$args --hwfile $hwfile"
fi
args="$args --an $ant -o $outfile"
echo "$args"

echo "python craftcor_tab.py $args --tab $f_vcraft"
python craftcor_tab.py $args --tab $f_vcraft
