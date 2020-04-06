#!/bin/bash
# Script to export FRB metadata
# Much of this is sourced from Hyerin Cho's README: https://docs.google.com/document/d/1uSbsudW5tkjy-GbJaCduERS36AgCEgWox-V8vY-YdA8/edit

basedr=/fred/oz002/users/hcho/craft/

# Should be two arguments:
# source FRBdata.sh [FRB name] [AIPS or MIRIAD] [Polarisation (x/y)] [antenna]
if (( $# != 3 )); then
	echo "ERROR: Usage: source FRBdata.sh [FRB name] [AIPS or MIRIAD] [Polarisation (x/y)] [antenna]"
else
	FRB=$1
	a_or_m=$2
	pol=$3
	an=$4

	if [ "$FRB" == "180924" ]; then
		offset=1874193
		DM=362.193
		f0=1320.5
		calcfile="${basedr}Calibration/aipscal/frb180924/c1_f0/craftfrb.im"
		hwfile="${basedr}Calibration/SB6635_b18_neweop.hwdelays"
	
		if [ "$a_or_m" == "AIPS" ]; then
			aips="${basedr}Calibration/aipscal/frb180924/noxpol/bandpasses.bp.txt"
			fcm="${basedr}Calibration/aipscal/frb180924/fcm.txt"
		elif [ "$a_or_m" == "MIRIAD" ]; then
			mir="${basedr}Calibration/0407_allfreq/20180924212734_call_beam37_i1024_f9.uvaver"
			fcm="${basedr}Calibration/FRB180924-calcfiles/fcm_release2_ak25mod.txt"
		else
			echo "ERROR: Must provide AIPS or MIRIAD exactly!!!"
		fi

		if [ "$pol" == "x" ]; then
			f_vcraft="${basedr}python/voltages/FRB180924/ak**/beam36/*.vcraft"
			f_outfile="./f_x_an${an}.npy"
			t_outfile="./t_x_an${an}.npy"
		elif [ "$pol" == "y" ]; then
			f_vcraft="${basedr}python/voltages/FRB180924/ak**/beam37/*.vcraft"
			f_outfile="./f_y_an${an}.npy"
			t_outfile="./t_y_an${an}.npy"
		else
			echo "ERROR: Must provide polarisation as x or y!"
		fi

	elif [ "$FRB" == "181112" ]; then
		offset=2106003
		DM=589.265
		f0=1297.5
		calcfile="${basedr}Calibration/aipscal/frb181112/c1_f0/craftfrb.im"
		fcm="${basedr}Calibration/aipscal/frb181112/fcm.txt"
		hwfile="${basedr}Calibration/mircal/frb181112/hwdelays_0407_SB7031_20181112222901_round-8.txt"
		aips="${basedr}Calibration/aipscal/frb181112/bandpass.bp.txt"

		if [ "$pol" == "x" ]; then
			mir="${basedr}Calibration/mircal/frb181112/20181112222901_call_beam00_i4096_f9.uvlin"
			f_vcraft="${basedr}python/voltages/FRB181112/ak**/beam00/*.vcraft"
			f_outfile="./f_x_an${an}.npy"
			t_outfile="./t_x_an${an}.npy"
		elif [ "$pol" = "y" ]; then
			mir="${basedr}Calibration/mircal/frb181112/20181112222901_call_beam01_i4096_f9.uvlin"
			f_vcraft="${basedr}python/voltages/FRB181112/ak**/beam01/*.vcraft"
			f_outfile="./f_y_an${an}.npy"
			t_outfile="./t_y_an${an}.npy"
		else
			echo "ERROR: Must provide polarisation as x or y!"
		fi
	
	elif [ "$FRB" == "190102" ]; then
		offset=1880759
		DM=364.544
		f0=1271.5
		calcfile="${basedr}Calibration/aipscal/frb190102/craftfrb.im"
		fcm="${basedr}Calibration/aipscal/frb190102/fcm.txt"
		hwfile="${basedr}Calibration/aipscal/frb190102/hwdelays_converted.txt"
		aips="${basedr}Calibration/aipscal/frb190102/noxpol/bandpasses.bp.txt"

		if [ "$pol" == "x" ]; then
			f_vcraft="${basedr}python/voltages/FRB190102/ak**/beam16/*.vcraft"
			f_outfile="./f_x_an${an}.npy"
			t_outfile="./t_x_an${an}.npy"
		elif [ "$pol" == "y" ]; then
			f_vcraft="${basedr}python/voltages/FRB190102/ak**/beam17/*.vcraft"
			f_outfile="./f_y_an${an}.npy"
			t_outfile="./t_y_an${an}.npy"
		else
			echo "ERROR: Must provide polarisation as x or y!"
		fi

	else
		echo "ERROR: FRB not recognised!"
	fi
fi
