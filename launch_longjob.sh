#!/bin/bash

#module load cython
mode=recon #sumant #cal #dmtest #sumant #get_snr #tab_each #testoffset #recon #dd #
here=`pwd`

if [ $mode = testoffset ]
then
script=test_fftoffset.sh
echo "testing fft offset"
for offset in 0 #{26..32} # {0..32} 
do
for an in {0..23}
do
    fn=log/test_offset${offset}_an${an}.log #${script}
    t=6:00:00
    jobname=an${an}
    mem=7g
    echo "sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script $offset"
    sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script $offset $an
done
done
fi

if [ $mode = tab_each ]
then
script=do_corr_cho.sh
echo "TAB-ing for each antenna"
for an in 0 #{0..23}
do
    fn=${here}/log/sf_in_offset1.log #${an}.log #${script}
    t=8:00:00 #10:00:00
    mem=14g
    jobname=sf_in #${an}

    echo "sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script $an"
    sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script $an
done
fi

if [ $mode = sumant ]
then
script=sum.sh
echo "TAB-ing for each antenna"
for offset in 0 #{26..32}
do
    fn=${here}/log/sumant.log #${script}
    t=6:00:00
    mem=20g
    jobname=recon_sumant

    echo "sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script $offset"
    sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script #$offset
done
fi

if [ $mode = get_snr ]
then
script=get_snr.sh
echo "TAB-ing for each antenna"
for an in 0 #{0..31} #{0..23}
do
    fn=${here}/log/snrref_an${an}.log #${script}
    t=6:00:00
    mem=30g # 15g 40g
    jobname=${an}ansnr

    echo "sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script $offset"
    sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script #$an
done
fi

if [ $mode = dmtest ]
then
echo "DM test"
script=do_dmtest.sh
t=6:00:00
for decimal in {0..9}
do
DM=$(bc <<< 364.0+0.1*$decimal)
fn=${here}/log/dm_test_${DM}.log
jobname=${DM}dm_test
mem=70g
echo "sbatch --job-name=${jobname} --output=${fn} --error=${fn} --mem=${mem} --time=${t} $script $DM"
sbatch --job-name=${jobname} --output=${fn} --error=${fn} --mem=${mem} --time=${t} $script $DM
done
fi

if [ $mode = cal ]
then
script=calcorr.sh
echo "Correlation for calibration data"
refant=0
frb=191001
pol=y
mode=nobp #aips_maginv #aips_noxpol, nobp, mir, aips
fn=${here}/log/cal_${mode}_${pol}_${frb}_refant${refant}.log
t=5:00:00
mem=15g #15g #7g
jobname=${mode}${pol}cal

echo "sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script $frb $pol $refant $mode"
sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script $frb $pol $refant $mode

fi

if [ $mode = dd ]
then
script=dd_timeseries.sh
echo "FFFF to dedispersed time series"

fn=${here}/log/dd_181112.log
t=5:00:00
mem=30g
jobname=dd

echo "sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script"
sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script
fi

if [ $mode = polcal ]
then
script=do_polcal.sh
echo "Applying polarization calibration to VELA"
pol=y

fn=${here}/log/polcal_${pol}.log
t=5:00:00
mem=15g #30g # 
jobname=${pol}polcal

echo "sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script"
sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script $pol
fi

if [ $mode = recon ]
then
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
echo "$frb"
fi
for (( iant=0; iant<$nant; iant++  ))
#for iant in {0..22} #23} #11} #
do
fn=${here}/log/recon${iant}${pol}_${frb}.log
t=3:00:00
mem=100g
jobname=${iant}recon

echo "sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script $iant $frb $a_or_m $pol"
sbatch --job-name=${jobname} --output=${fn} --error=${fn} --time=${t} --mem=${mem} $script $iant $frb $a_or_m $pol
done
fi
