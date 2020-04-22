#!/usr/bin/env python3

"""
Given an FRB's summed series for x and/or y, reduce the time resolution by a given factor and plot the series
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from scipy import signal


def _main():
	args = get_args()

	if not ( args.x or args.y or args.i or args.q or args.u or args.v):
		print("Need at least one of x, y, i, q, u, or v!")
	elif args.x or args.y:
		X, Y = load_xy(args)
		stokes = get_stokes(X, Y, args)
		peak = np.argmax(stokes[0])	# get the peak as the highest point in Stokes I
		reduced_stokes = reduce_all(stokes, args.n, args.c, peak=peak)
		plot_stokes(reduced_stokes, smooth=args.s)
	else:
		stokes = load_stokes(args)
		peak = np.argmax(stokes[0])
		reduced_stokes = reduce_all(stokes, args.n, args.c, peak=peak)
		plot_stokes(reduced_stokes, smooth=args.s)


def get_args():
	"""
	Get command line arguments
	"""
	parser = ArgumentParser(description="Given an FRB's summed series for x and/or y, reduce the time resolution by a given factor and plot the Stokes parameters", formatter_class=ArgumentDefaultsHelpFormatter)

	parser.add_argument('-x', type=str, help='X polarisation filename')
	parser.add_argument('-y', type=str, help='Y polarisation filename')
	parser.add_argument('-i', type=str, help='Stokes I filename. If -x and/or -y are provided, this is the file the newly-calculated Stokes I is saved to')
	parser.add_argument('-q', type=str, help='Stokes Q filename. If -x and -y are provided, this is the file the newly-calculated Stokes Q is saved to')
	parser.add_argument('-u', type=str, help='Stokes U filename. If -x and -y are provided, this is the file the newly-calculated Stokes U is saved to')
	parser.add_argument('-v', type=str, help='Stokes V filename. If -x and -y are provided, this is the file the newly-calculated Stokes V is saved to')
	parser.add_argument('-n', type=int, help='Factor to reduce time resolution by [default=1000000]', default=1000000)
	parser.add_argument('-s', type=int, help='Smooth to apply to the plotted series [default=None]', default=None)
	parser.add_argument('-c', action='store_true', help='Crop plotted series to the locale of the peak', default=False)

	return parser.parse_args()


def load_xy(args):
	"""
	Load X and Y polarisation series. Numpy's mmap_mode argument lets the actual data stay on disk instead of loading it
	into memory.
	"""

	if args.x:
		print(f"Loading {args.x}...")
		X = np.load(args.x, mmap_mode='r')
	else:
		X = None

	if args.y:
		print(f"Loading {args.y}...")
		Y = np.load(args.y, mmap_mode='r')
	else:
		Y = None

	return X, Y


def load_stokes(args):
	"""
	Loads the stokes parameters from the given filenames in args.
	If the file doesn't exist, that parameter will be None.
	"""
	ret = []

	for fname in (args.i, args.q, args.u, args.v):
		if fname:
			try:
				print(f"Loading {fname}...")
				par = np.load(fname, mmap_mode='r')
			except FileNotFoundError:
				par = None
		else:
			par = None
		ret.append(par)

	return ret


def get_stokes(X, Y, args):
	"""
	Calculates, saves, and returns the Stokes parameters I, Q, U, and V from provided X and/or Y polarisations.
	If only one of X or Y is not None, only I is calculated and saved. Q, U, and V are returned as None
	"""

	# The Stokes arrays are generally large. This function saves a given array and loads it again in mmap_mode to save on memory
	def save_load(fname, A):
		np.save(fname, A)
		del A
		return np.load(fname, mmap_mode='r')

	if X is not None and Y is not None:
		print("Calculating I...")
		I = np.abs(X)**2 + np.abs(Y)**2
		I = save_load(args.i, I)

		print("Calculating Q...")
		Q = np.abs(X)**2 - np.abs(Y)**2
		Q = save_load(args.q, Q)

		print("Calculating U...")
		U = 2*np.real(np.conj(X) * Y)
		U = save_load(args.u, U)

		print("Calculating V...")
		V = 2*np.imag(np.conj(X) * Y)
		V = save_load(args.v, V)
	
	elif Y is None:
		print("Calculating I...")
		I = np.abs(X)**2
		I = save_load(args.i, I)

		Q = U = V = None

	else:
		print("Calculating I...")
		I = np.abs(Y)**2
		I = save_load(args.i, I)

		Q = U = V = None

	return I, Q, U, V


def reduce_all(series_list, n, crop, peak=None):
	"""
	Reduces the time resolution of all the arrays in the given list by a factor of n.
	If crop is True, then only the 2*crop_range values surrounding the peak of the first array in the list is returned.
	"""
	dt = ((336*u.MHz)**(-1)).to(u.ns)	# Base time resolution
	ret = []
	crop_range = 1000000	# Space around the peak to crop if c is True
	for series in series_list:
		if series is not None:
			if crop:
				ret.append(reduce(series[peak-crop_range:peak+crop_range], n))
			else:
				ret.append(reduce(series, n))

	print("Time resolution reduced to " + str(dt*n))

	return ret


def reduce(A, n):
	"""
	Reduces the time resolution of a given array by a factor n.
	"""
	A_red = []
	for i in range(int(len(A)/n)):
		A_red.append(np.sum(A[i*n:(i+1)*n]))

	return np.array(A_red)


def plot_stokes(stokes, smooth=None):
	"""
	Plots the given stokes parameters. If smooth is not none, they are smooth with a Gaussian.
	"""
	# Little function to generate a smooth series by convolving with a Gaussian
	# https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
	def smooth_f(y, box_pts):
		box = signal.gaussian(box_pts, box_pts/4)	# The "box" is a Gaussian that fits out to two standard deviations
		y_smooth = np.convolve(y, box, mode='same')
		return y_smooth

	if smooth:
		for par in stokes:
			smooth_par = smooth_f(par, smooth)
			plt.plot(smooth_par - np.mean(par))		# Always zero-meaning plotted lines
	else:
		for par in stokes:
			plt.plot(par - np.mean(par))

	plt.legend(['I', 'Q', 'U', 'V'])
	plt.show()


if __name__ == "__main__":
	_main()

