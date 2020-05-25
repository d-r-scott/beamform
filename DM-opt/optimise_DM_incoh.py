#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy import fft

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(description='Optimises the DM of a burst by S/N via coherent dedispersion', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', help='Input dispersed dynamic Stokes I spectrum')
parser.add_argument('-m', '--DM_min', help='Minimum of DM range', default=0)
parser.add_argument('-M', '--DM_max', help='Maximum of DM range', default=4000)
parser.add_argument('-p', '--prec', help='Precision of DM to optimise to', default='1')
parser.add_argument('-f', '--f_mid', help='Central frequency of band (MHz)')
parser.add_argument('-b', '--bw', help='Bandwidth (MHz)')
parser.add_argument('-o', help='Output dynamic Stokes I spectrum')
values = parser.parse_args()

dt		= 1e-6							# Time resolution in s
df		= 1								# Time resolution in MHz
n_chans = int(values.bw/df)				# Number of frequency channels
f_min	= values.f_mid - values.bw/2	# Bottom of band
f_max	= values.f_mid + values.bw/2	# Top of band
k_dm    = 4.148808e3    				# MHz^2 pc^-1 cm^3 s

# Functions for consistent calculation of idt from DM and vice versa
idt_f   = lambda DM, fi : k_dm * (fi**(-2) - f_max**(-2)) * DM / dt
DM_f    = lambda idt, fi : idt * dt / (k_dm * (fi**(-2) - f_max**(-2)))



def dedisperse(dynspec, DM):
	F, T = dynspec.shape

	idt_max = int(idt_f(DM, f_min))

	dd_dynspec = np.zeros((F, T - idt_max))

	for f in range(F):
		fi = f_min + f*df
		shift = int(idt_f(DM, fi))
		for t in range(T - idt_max):
			dd_dynspec[-(f+1), t] = dynspec[-(f+1), t + shift]