#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
import astropy.units as u
from scipy.fft import fft, ifft


def _main():
	args = get_args()

	x_3ns_t = load(args.x)
	y_3ns_t = load(args.y)
	i_3ns_t = load(args.i)

	x_3ns_t_crop, y_3ns_t_crop, t_res_us = crop(x_3ns_t, y_3ns_t, i_3ns_t)

	Delta_DMs = np.arange(args.DM_min, args.DM_max, args.DM_res)

	H_array_fname = create_H_array(x_3ns_t_crop.shape[0], 1297.5, 336, Delta_DMs, args.DM_min, args.DM_max, args.DM_res)

	dedisperse_many(x_3ns_t_crop, y_3ns_t_crop, H_array_fname, Delta_DMs, t_res_us)

	os.remove(H_array_fname)


def get_args():
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description='Calculates S/N(DM) for a given FRB', formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument('-x', help='Dedispersed 3ns x time series')
	parser.add_argument('-y', help='Dedispersed 3ns y time series')
	parser.add_argument('-i', help='Dedispersed 3ns i time series')
	parser.add_argument('--DM_min', help='Minimum Delta DM', type=float)
	parser.add_argument('--DM_max', help='Maximum Delta DM', type=float)
	parser.add_argument('--DM_res', help='DM resolution', type=float)
	return parser.parse_args()


def load(fname):
	return np.load(fname, mmap_mode='r')


def crop(x_3ns_t, y_3ns_t, i_3ns_t):
	assert len(x_3ns_t) == len(y_3ns_t)

	while True:
		t_res_us = int(input("Time resolution to plot with in us: "))

		# Original time resolution is (336 MHz)^-1, which is exactly 1/336 of 1 us
		reduction_factor = 336 * t_res_us

		i_red_t = reduce(i_3ns_t, 336*t_res_us)

		t_red = (np.arange(0, len(i_red_t)) * (t_res_us * u.us)).to(u.s)

		plt.plot(t_red, i_red_t)
		plt.show()

		in_str = input("Input t_min, t_max (both in s) or X for new time resolution: ")

		if in_str != "X":
			t_min, t_max = map(float, in_str.split(", "))

			t_3ns = (np.arange(0, len(x_3ns_t)) / (336*u.MHz)).to(u.s)

			t_mask = (t_3ns.value >= t_min) & (t_3ns.value < t_max)
			x_3ns_t_crop = x_3ns_t[t_mask]
			y_3ns_t_crop = y_3ns_t[t_mask]

			return x_3ns_t_crop, y_3ns_t_crop, t_res_us


def create_H_array(n_sam, f0, bw, Delta_DMs, Delta_DM_min, Delta_DM_max, DM_res):
	f_min = f0 - float(bw)/2
	f_max = f0 + float(bw)/2
	freqs = np.linspace(f_max, f_min, n_sam)

	k_DM = 2.41e-4

	H_array = np.zeros((Delta_DMs.shape[0], n_sam), dtype=np.complex64)

	for i, DM in enumerate(Delta_DMs):
		if i % 100 == 0:
			print(DM)
		H_array[i, :] = np.exp(2j*np.pi*DM/k_DM*((freqs-f0)**2/f0**2/freqs*1e6))

	H_array_fname = f'H_{f0}_{n_sam}_{Delta_DM_min}-{Delta_DM_max}-{DM_res}.npy'

	np.save(H_array_fname, H_array)

	del H_array
	return H_array_fname


def dedisperse_many(x_3ns_t_crop, y_3ns_t_crop, H_array_fname, Delta_DMs, t_res_us):
	H_array = np.load(H_array_fname, mmap_mode='r')

	x_3ns_f_crop = fft(x_3ns_t_crop)
	y_3ns_f_crop = fft(y_3ns_t_crop)

	peak_sns = np.zeros(Delta_DMs.shape)

	for i, DM in enumerate(Delta_DMs):
		x_3ns_t_dd = dedisperse_ifft(x_3ns_f_crop, H_array[i,:])
		y_3ns_t_dd = dedisperse_ifft(y_3ns_f_crop, H_array[i,:])

		i_3ns_t_dd = np.abs(x_3ns_t_dd)**2 + np.abs(y_3ns_t_dd)**2

		reduction_factor = 336*t_res_us
		i_red_t_dd = reduce(i_3ns_t_dd, reduction_factor)
		i_norm = normalise(i_red_t_dd)
		peak_sns[i] = np.max(i_norm)

	plt.plot(Delta_DMs, peak_sns)
	plt.xlabel(r'$\Delta$ DM')
	plt.ylabel('S/N')
	plt.show()
	np.save('SN_DM.npy', peak_sns)


def dedisperse_ifft(f, H):
	return ifft(f*H)


def normalise(a):
	return (a - np.mean(a))/np.std(a)


def reduce(A, n):
	"""
	Reduces the time resolution of a given array by a factor n.
	"""
	A_red = []
	for i in range(int(A.shape[0]/n)):
		A_red.append(np.sum(A[i*n:(i+1)*n], axis=0))

	return np.array(A_red)


if __name__ == '__main__':
	_main()
