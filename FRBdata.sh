#!/bin/bash
# Script to export FRB metadata
# Much of this is sourced from Hyerin Cho's README: https://docs.google.com/document/d/1uSbsudW5tkjy-GbJaCduERS36AgCEgWox-V8vY-YdA8/edit

basedr=/fred/oz002/users/dscott/
basedr2=/fred/oz002/users/hcho/craft/

# Should be two arguments:
# source FRBdata.sh [FRB name] [AIPS or MIRIAD] [Polarisation (x/y)] [antenna]
if (( $# != 4 )); then
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
			f_vcraft="${basedr2}python/voltages/FRB180924/ak**/beam36/*.vcraft"
		elif [ "$pol" == "y" ]; then
			f_vcraft="${basedr2}python/voltages/FRB180924/ak**/beam37/*.vcraft"
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
			f_vcraft="${basedr2}python/voltages/FRB181112/ak**/beam00/*.vcraft"
		elif [ "$pol" = "y" ]; then
			mir="${basedr}Calibration/mircal/frb181112/20181112222901_call_beam01_i4096_f9.uvlin"
			f_vcraft="${basedr2}python/voltages/FRB181112/ak**/beam01/*.vcraft"
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
			f_vcraft="${basedr2}python/voltages/FRB190102/ak**/beam16/*.vcraft"
		elif [ "$pol" == "y" ]; then
			f_vcraft="${basedr2}python/voltages/FRB190102/ak**/beam17/*.vcraft"
		else
			echo "ERROR: Must provide polarisation as x or y!"
		fi


	elif [ "$FRB" == "190608" ]; then
		offset=2018819
		DM=338.8
		f0=1271.5
		calcfile="${basedr}Calibration/aipscal/frb${FRB}/craftfrb.im" # geometric delays
		fcm="${basedr}Calibration/aipscal/frb${FRB}/fcm.txt" # clock delays
		hwfile= # No hardware delays
		#mir= # MIR gain, bandpass
		aips="${basedr}Calibration/aipscal/frb${FRB}/noxpol/200225/bandpasses.bp.txt" # AIPS gain, bandpass

		if [ "$pol" == "x" ]; then
			f_vcraft=$(find ${basedr2}python/voltages/FRB${FRB}/ak**/beam36/*.vcraft -not -name "*ak13*" -not -name "*ak19*" -not -name "*ak20*" -not -name "*ak28*" )
		elif [ "$pol" == "y" ]; then
			f_vcraft=$(find ${basedr2}python/voltages/FRB${FRB}/ak**/beam37/*.vcraft -not -name "*ak13*" -not -name "*ak19*" -not -name "*ak20*" -not -name "*ak28*" )
		else
			echo "ERROR: Must provide polarisation as x or y!"
		fi

	elif [ "$FRB" == "190611" ]; then
		offset=2018819
		DM=338.8
		f0=1271.5
		calcfile="${basedr}Calibration/aipscal/frb${FRB}/craftfrb.im"
		fcm="${basedr}Calibration/aipscal/frb${FRB}/fcm.txt"
		hwfile= # No hardware delays
		#mir= # MIR gain, bandpass
		aips="${basedr}Calibration/aipscal/frb${FRB}/noxpol/200225/bandpasses.bp.txt" # AIPS gain, bandpass

		if [ "$pol" == "x" ]; then
			f_vcraft=$(find ${basedr2}python/voltages/FRB${FRB}/ak**/beam36/*.vcraft -not -name "*ak13*" -not -name "*ak19*" -not -name "*ak20*" -not -name "*ak28*" )
		elif [ "$pol" == "y" ]; then
			f_vcraft=$(find ${basedr2}python/voltages/FRB${FRB}/ak**/beam37/*.vcraft -not -name "*ak13*" -not -name "*ak19*" -not -name "*ak20*" -not -name "*ak28*" )
		else
			echo "ERROR: Must provide polarisation as x or y!"
		fi

	elif [ "$FRB" == "190711" ]; then
		offset=302100
		DM=587.5
		f0=1271.5
		calcfile="${basedr}Calibration/aipscal/frb${FRB}/craftfrb.im" # geometric delays
		fcm="${basedr}Calibration/aipscal/frb${FRB}/fcm.txt" # clock delays
		hwfile= # No hardware delays
		#mir= # MIR gain, bandpass
		aips="${basedr}Calibration/aipscal/frb${FRB}/noxpol/bandpasses_noxpol.bp.txt" # AIPS gain, bandpass

		if [ "$pol" == "x" ]; then
			f_vcraft="${basedr2}python/voltages/FRB${FRB}/ak**/beam30/*.vcraft"
		elif [ "$pol" == "y" ]; then
			f_vcraft="${basedr2}python/voltages/FRB${FRB}/ak**/beam31/*.vcraft"
		else
			echo "ERROR: Must provide polarisation as x or y!"
		fi

	elif [ "$FRB" == "190714" ]; then
		offset=1887164
		DM=504.13		# Uncertain
		f0=1271.5
		calcfile="${basedr}Calibration/aipscal/frb${FRB}/craftfrb.im" # geometric delays
		fcm="${basedr}Calibration/aipscal/frb${FRB}/fcm.txt" # clock delays
		hwfile= # No hardware delays
		#mir= # MIR gain, bandpass
		aips="/fred/oz002/users/dscott/Calibration/aipscal/frb190714/bandpasses_noxpol_FRB.bp.txt"

		if [ "$pol" == "x" ]; then
			f_vcraft="${basedr2}python/voltages/FRB${FRB}/ak**/beam56/*.vcraft"
		elif [ "$pol" == "y" ]; then
			echo "y pol corrupted"
		    #f_vcraft="${basedr2}python/voltages/FRB${FRB}/ak**/beam27/*.vcraft"
		else
			echo "ERROR: Must provide polarisation as x or y!"
		fi

	elif [ "$FRB" == "191228" ]; then
		offset=1576224
		DM=296.5
		f0=1271.5
		calcfile="${basedr}Calibration/aipscal/frb${FRB}/craftfrb.im" # geometric delays
		fcm="${basedr}Calibration/aipscal/frb${FRB}/fcm.txt" # clock delays
		hwfile= # No hardware delays
		#mir= # MIR gain, bandpass
		aips="${basedr}Calibration/aipscal/frb${FRB}/noxpol/bandpasses_noxpol.bp.txt" # AIPS gain, bandpass

		if [ "$pol" == "x" ]; then
			f_vcraft="${basedr2}python/voltages/FRB${FRB}/ak**/beam42/*.vcraft"
		elif [ "$pol" == "y" ]; then
			f_vcraft="${basedr2}python/voltages/FRB${FRB}/ak**/beam43/*.vcraft"
		else
			echo "ERROR: Must provide polarisation as x or y!"
		fi

	else
		echo "ERROR: FRB not recognised!"
	fi

	f_outfile="./output/${FRB}_${pol}_${an}_f.npy"
	t_outfile="./output/${FRB}_${pol}_${an}_t.npy"

fi
