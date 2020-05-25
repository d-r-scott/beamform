#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy import fft

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(description='Optimises the DM of a burst by S/N via coherent dedispersion', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', help='Input dispersed fourier spectrum')
parser.add_argument('-m', '--DM_min', help='Minimum of DM range', default=0)
parser.add_argument('-M', '--DM_max', help='Maximum of DM range', default=4000)
parser.add_argument('-p', '--prec', help='Precision of DM to optimise to', default='1')
parser.add_argument('-f', '--f_mid', help='Central frequency of band (MHz)')
parser.add_argument('-b', '--bw', help='Bandwidth (MHz)')
parser.add_argument('-o', help='Output time series file')
values = parser.parse_args()

spec = np.load(values.i, mmap_mode='r')
n_sam = spec.shape[0]

f_start = values.f_mid - values.bw/2
f_stop = values.f_mid + values.bw/2

dedisperse = lambda DM: np.exp(2j*np.pi*DM/2.41e-4*np.array([(f-f_mid)**2/f_mid**2/f*1e6 for f in np.linspace(f_stop, f_start, n_sam)]))

def reduce(A, n):
	"""
	Reduces the time resolution of a given array by a factor n.
	Reduces along the 0th axis, make sure this one is the time axis!

	:param A: Input array to be reduced
	:param n: Factor to reduce by
	:return: None
	"""
	if n > 1:	# if n == 1 (or below, and is therefore invalid) don't bother reducing and just return the input
		A_red = []
		with progressbar.ProgressBar(max_value=int(A.shape[0]/n)) as bar:
			for i in range(int(A.shape[0]/n)):
				A_red.append(np.sum(A[i*n:(i+1)*n], axis=0))
				bar.update(i)
		return np.array(A_red)
	else:
		return(A)

def opt_wrapper(DM):
	dd_arr = dedisperse(DM)
	dd_spec = spec*dd_arr
	t_series = fft.ifft(dd_spec)
	red_t_series = reduce(t_series, 336000)	# Reduce time resolution to 1ms
	peak = np.max(red_t_ser)
	return -1*peak

x0 = [ np.mean((values.DM_min, values.DM_max)) ]
res = minimize(opt_wrapper, x0, bounds=(values.DM_min, values.DM_max), tol=values.prec)
print(res)