#!/bin/bash

#SBATCH --job-name=ifft
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0:10:00
#SBATCH --mem=32g

# Command line arguments
FRB=$1  # FRB name
pol=$2  # Polarisation (x or y)
DM=$3   # Dispersion measure in pc/cm3

# Get data directories
source dir_vars.sh $FRB

# Get modules to load and load them
source modules.sh
module load $modules_5

args="-f ${outdir}/${FRB}_sum_${pol}_f_dedispersed_${DM}.npy -o ${outdir}/${FRB}_sum_${pol}_t_${DM}.npy"

echo "python3 ifft.py $args"
python3 ifft.py $args