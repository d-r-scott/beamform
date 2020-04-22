#!/bin/bash

FRB=$1
pol=$2
out=$3

module load gcc/6.4.0 
module load python/3.7.4
module load numpy/1.18.2-python-3.7.4

python sum.py -p $pol -o $out $FRB
