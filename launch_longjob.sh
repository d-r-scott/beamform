#!/bin/bash

#module load cython
mode=recon #sumant #cal #dmtest #sumant #get_snr #tab_each #testoffset #recon #dd #
here=`pwd`

if [ $mode = recon ]; then
  script=fulltimeres.sh
  echo "Reconstructing to its full time resolution"

  frb=$1
  a_or_m=$2
  pol=$3

  if [ $frb = 190102 ]; then
    nant=23
  elif [ $frb = 181112 ]; then
    nant=12
  elif [ $frb = 180924 ]; then
    nant=24
  elif [ $frb = 190608 ]; then
    nant=21 #25
  elif [ $frb = 190611 ]; then
    nant=25
  elif [ $frb = 190711 ]; then
    nant=28
  elif [ $frb = 190714 ]; then
    nant=25
  elif [ $frb = 191001 ]; then
    nant=30
  elif [ $frb = 191228 ]; then
    nant=21
  elif [ $frb = 200430 ]; then
    nant=26
  else
    nant=0
    echo "Please specify total # of antenna for this FRB${frb}"
  fi

  for (( iant=0; iant<$nant; iant++ )); do
    fn=${here}/log/recon${iant}${pol}_${frb}.log
    t=1:30:00
    mem=64g
    jobname=${iant}${pol}recon

    echo "sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script $iant $frb $a_or_m $pol"
    sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script $iant $frb $a_or_m $pol
  done
fi
