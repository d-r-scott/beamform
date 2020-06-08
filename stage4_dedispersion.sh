#!/bin/bash

#SBATCH --job-name=dedispersion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=32g

# Command line arguments
FRB=$1  # FRB name
pol=$2  # Polarisation (x or y)
DM=$3   # Dispersion measure in pc/cm3
f0=$4   # Central frequency in MHz

# Set data directories
basedir=./output
outdir=${basedir}/${FRB}
f_outdir=${outdir}/f

# Get modules to load and load them
source modules.sh
module load $modules_4

args="-f ${outdir}/${FRB}_sum_${pol}_f_derippled.npy --DM $DM --f0 $f0 --bw 336 -o ${outdir}/${FRB}_sum_${pol}_f_dedispersed_${DM}.npy"

echo "python3 dedisperse.py $args"
python3 dedisperse.py $args
