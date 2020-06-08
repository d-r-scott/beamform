#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import Bounds
import sys
import progressbar as pb
from joblib import Parallel, delayed
import multiprocessing
sys.path.append('../output')
import dynspec as ds

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(description='Optimises the DM of a burst by S/N via coherent dedispersion', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', help='Input dispersed dynamic Stokes I spectrum')
parser.add_argument('-m', '--DM_min', type=float, help='Minimum of DM range', default=0)
parser.add_argument('-M', '--DM_max', type=float, help='Maximum of DM range', default=4000)
parser.add_argument('-p', '--prec', type=float, help='Precision of DM to optimise to', default='1')
parser.add_argument('-w', type=int, help='Width of FRB in 1us time samples', default=1000)
parser.add_argument('-f', '--f_mid', type=float, help='Central frequency of band (MHz)')
parser.add_argument('-b', '--bw', type=int, help='Bandwidth (MHz)')
parser.add_argument('--dt', type=float, help='Time resolution (us)', default=1)
parser.add_argument('-o', help='Output dynamic Stokes I spectrum')
values = parser.parse_args()

dt		= values.dt*1e-6				# Time resolution in s
df		= 1								# Time resolution in MHz
n_chans = int(values.bw/df)				# Number of frequency channels
f_min	= values.f_mid - values.bw/2	# Bottom of band
f_max	= values.f_mid + values.bw/2	# Top of band
k_dm    = 4.148808e3    				# MHz^2 pc^-1 cm^3 s

# Functions for consistent calculation of idt from DM and vice versa
idt_f   = lambda DM, fi : k_dm * (fi**(-2) - f_max**(-2)) * DM / dt
DM_f    = lambda idt, fi : idt * dt / (k_dm * (fi**(-2) - f_max**(-2)))

in_dynspec = np.load(values.i, mmap_mode='r')


def dedisperse(dynspec, DM):
	F, T = dynspec.shape

	idt_max = int(idt_f(DM, f_min))

	dd_dynspec = np.zeros((F, T - idt_max))

	for f in range(F):
		fi = f_min + f*df
		shift = int(idt_f(DM, fi))
		for t in range(T - idt_max):
			#prof[t] += dynspec[-(f+1), t + shift]
			dd_dynspec[-(f+1), t] = dynspec[-(f+1), t + shift]

	return dd_dynspec


def get_sn(dynspec, w):
	# Sum dynspec across frequencies to get time series
	# Dividing  my sqrt(num freq. chann.) keeps mean=0, std=1
	prof = np.sum(dynspec, axis=0)/np.sqrt(dynspec.shape[0])

	# Sum over the width of the burst
	#print("Rolling prof")
	#roll_prof = np.sum([np.roll(prof, i) for i in range(w)], axis=0)/np.sqrt(w)
	sum_prof = np.array([ np.sum(prof[t:t+w]) for t in range(len(prof)-w) ])/np.sqrt(w)

	# Get the maximum value of the profile and return as S/N
	return max(sum_prof)



def opt_f(x):
	DM = x
	w = int(values.w/values.dt)

	sn = get_sn(dedisperse(in_dynspec, DM), w)
	return sn


def do_opt():
	DM_min = values.DM_min
	DM_max = values.DM_max

	# This ain't working
	"""
	x0 = np.array([(DM_min + DM_max)/2])
	print(f'Staring optimisation with DM={x0[0]}')
	#res = minimize(opt_f, x0, bounds=Bounds([DM_min], [DM_max]), tol=values.prec)
	res = minimize(opt_f, x0, tol=values.prec)

	print('Optimisation finished!')
	print(res.x)
	print(res.fun)
	print(res.success)
	"""
	DMs = np.arange(DM_min, DM_max, values.prec)

	num_cores = multiprocessing.cpu_count()

	sns = Parallel(n_jobs=num_cores-1)(delayed(opt_f)(DM) for DM in DMs)

	"""""
	with pb.ProgressBar(max_value=len(sns)) as bar:
		bar.update(0)
		for i, DM in enumerate(DMs):
			sns[i] = opt_f([DM])
			bar.update(i)
	"""

	opt_DM = DMs[np.argmax(sns)]
	print(opt_DM)
	plt.plot(DMs, sns)
	plt.show()
	opt_dynspec = dedisperse(in_dynspec, opt_DM)
	opt_dynspec_red = ds.reduce(opt_dynspec, int(500/values.dt), transpose=True)
	ds.plot_dynspec(opt_dynspec_red, freq_range=(862.5-336/2, 862.5+336/2), tick_space=(20, 50), axis_res=(0.5, 1), peak_idx=np.argmax(np.sum(opt_dynspec_red, axis=0)))
	np.save(values.o, opt_dynspec)

if __name__ == '__main__':
	do_opt()