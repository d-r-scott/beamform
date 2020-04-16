#!/bin/bash

# Crops all files for a given FRB

module load gcc/6.4.0
module load python/3.6.4
module load numpy/1.16.3-python-3.6.4

FRB=$1
max_ant=$2	# highest antenna number (don't forget they're zero-indexed!)

# Use ../FRBdata.sh just to get the offset for this FRB - only the FRB argument matters
source ../FRBdata.sh $FRB AIPS x 0

if [ ! -d $FRB/cropped ]; then
	mkdir $FRB/cropped
fi

for pol in x y; do
	for ant in $(seq 0 $max_ant); do
		echo "./crop_t.py -t $offset -o ${FRB}/cropped/${FRB}_${pol}_${ant}_crop.npy ${FRB}/${FRB}_${pol}_${ant}_t.npy"
		python crop_t.py -t $offset -o ${FRB}/cropped/${FRB}_${pol}_${ant}_crop.npy ${FRB}/${FRB}_${pol}_${ant}_t.npy
		echo ""
	done
done
