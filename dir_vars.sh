#!/bin/bash
# Script to export consistent data directories for a given FRB
# Also ensures the directories exist

FRB=$1

basedir=./output
outdir=${basedir}/${FRB}
f_outdir=${outdir}/f

if [ ! -d $outdir ]; then
  mkdir $outdir
fi

if [ ! -d $f_outdir ]; then
  mkdir $f_outdir
fi