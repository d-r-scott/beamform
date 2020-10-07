#!/bin/bash

#SBATCH --job-name=dynspecs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=64g

# Command line arguments
FRB=$1  # FRB name
DM=$2   # Dispersion measure in pc/cm3
n=$3

# Set data directories
basedir=./output
outdir=${basedir}/${FRB}_n${n}
f_outdir=${outdir}/f

# Get modules to load and load them
source modules.sh
module load $modules_6

args="-x ${outdir}/${FRB}_sum_x_t_${DM}.npy -y ${outdir}/${FRB}_sum_y_t_${DM}.npy -o ${outdir}/${FRB}_sum_!_@_${DM}.npy"

echo "python3 dynspecs.py $args"
python3 dynspecs.py $args
