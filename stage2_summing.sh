#!/bin/bash

#SBATCH --job-name=summing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0:15:00
#SBATCH --mem=16g

# Command line arguments
FRB=$1  # FRB name
pol=$2  # Polarisation (x or y)

# Set data directories
basedir=./output
outdir=${basedir}/${FRB}
f_outdir=${outdir}/f

# Get modules to load and load them
source modules.sh
module load $modules_2

args="--f_dir $f_outdir -f ${FRB} -p ${pol} -o $outdir/${FRB}_sum_${pol}_f.npy"

echo "python3 sum.py $args"
python3 sum.py $args
