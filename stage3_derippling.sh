#!/bin/bash

#SBATCH --job-name=derippling
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=32g

# Command line arguments
FRB=$1  # FRB name
pol=$2  # Polarisation (x or y)
fftlen=$3

# Set data directories
basedir=./output
outdir=${basedir}/${FRB}
f_outdir=${outdir}/f

# Get modules to load and load them
source modules.sh
module load $modules_3

args="-f ${outdir}/${FRB}_sum_${pol}_f.npy -l $fftlen -o ${outdir}/${FRB}_sum_${pol}_f_derippled.npy"

echo "python deripple.py $args"
python deripple.py $args
