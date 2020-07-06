#!/bin/bash

#SBATCH --job-name=fine_DM
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0:30:00
#SBATCH --mem=128g

# Command line arguments
x=$1        # Dedispersed 3ns x time series
y=$2        # Dedispersed 3ns y time series
DM_min=$3   # Minimum Delta DM
DM_max=$4   # Maximum Delta DM
DM_res=$5   # DM resolution
t_min=$6    # Minimum time to crop to (s)
t_max=$7    # Maximum time to crop to (s)
t_res=$8    # Time resolution to use when calculating S/N (us)
w=$9        # Burst width to use when calculating S/N (us)
b=${10}     # Bandwidth (MHz)
f=${11}     # Central frequency (MHz)

module load gcc/7.3.0 openmpi/3.0.0 python/3.7.4 numpy/1.18.2-python-3.7.4 scipy/1.4.1-python-3.7.4 matplotlib/3.2.1-python-3.7.4 astropy/4.0.1-python-3.7.4

echo "python3 fine_DM.py -x $x -y $y --DM_min $DM_min --DM_max $DM_max --DM_res $DM_res --t_min $t_min --t_max $t_max --t_res $t_res -w $w -b $b -f $f"
python3 fine_DM.py -x $x -y $y --DM_min $DM_min --DM_max $DM_max --DM_res $DM_res --t_min $t_min --t_max $t_max --t_res $t_res -w $w -b $b -f $f

