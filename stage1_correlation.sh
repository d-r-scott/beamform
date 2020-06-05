#!/bin/bash --login

#SBATCH --job-name=correlation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=32g

# Command line arguments
FRB=$1      # FRB name
a_or_m=$2   # AIPS or MIRIAD solutions for correlation
pol=$3      # Polarisation (x or y)
offset=$4   # Data offset (how many microseconds at the start to skip)
calcfile=$5
fcm=$6
f_vcraft=$7
a_m_file=$8 # AIPS/MIRIAD file (both are in this variable, the contents depends on $a_or_m
hwfile=$9   # Hardware delays. Probably not there for newer FRBs.

# Processing parameters
i=1
n=40960   # $n * 54 * 336 is the total length of the output array. Try to make n % 32 == 0.

# Get data directories
source dir_vars.sh $FRB

# Get modules to load and load them
source modules.sh
module load $modules_1

args="-i $i -n $n --offset $offset --calcfile $calcfile --parset $fcm --tab $f_vcraft"
if [ "$a_or_m" == "AIPS" ]; then
  args="$args --aips_c $a_m_file"
else
  args="$args --mirsolutions $a_m_file"
fi
if [ "$hwfile" != "" ]; then
  args="$args --hwfile $hwfile"
fi
args="$args -o $f_outdir"

echo "python craftcor_tab.py $args"
python craftcor_tab.py $args