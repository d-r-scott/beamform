#!/usr/bin/env python3

import numpy as np
import dynspec as ds
import sys
from scipy.fft import ifft
import gc


def old_coh_dedisp(FFFF, DM, f_mid=864.5, bw=336, quiet=False):
	nSam = len(FFFF)

	# ASKAP Parameters
	f_start = f_mid - float(bw) / 2  # 1153.
	f_stop = f_mid + float(bw) / 2  # 1488.

	if not quiet:
		print('dedispersing....')
	dedisperse = np.exp(2j * np.pi * DM / 2.41e-4 * np.array([(f - f_mid) ** 2 / f_mid ** 2 / f * 1e6 for f in np.linspace(f_stop, f_start, nSam)]))

	FFFF *= dedisperse

	return FFFF


def new_coh_dedisp(FFFF, DM, f_mid=864.5, bw=336, quiet=False):
	nSam = len(FFFF)
	k_DM = 4148.808

	f_start = f_mid - float(bw) / 2
	f_stop = f_mid + float(bw) / 2

	dedisperse = np.exp(2j * np.pi * DM * k_DM * np.array([(f - f_mid) ** 2 / f_mid ** 2 / f * 1e6 for f in np.linspace(f_stop, f_start, nSam)]))

	FFFF *= dedisperse

	return FFFF


def do_dedisp(fname, dedisp_func):
	f = np.load(fname)

	print('dedispersing')
	f_dd = dedisp_func(f, 380)
	del f

	print('ifft')
	t_dd = ifft(f_dd)[0,:,0]
	del f_dd

	print('generating dynspec')
	dynspec = ds.generate_dynspec(t_dd)
	del t_dd
	gc.collect()
	return dynspec


def do_dedisp_wrapper(old_new, x_fname, y_fname, dedisp_func):
	print(old_new)
	print(x_fname)
	print(y_fname)
	print('x')
	x_ds = do_dedisp(x_fname, dedisp_func)
	np.save(f'x_ds_{old_new}.npy', x_ds)
	del x_ds
	gc.collect()
	x_ds = np.load(f'x_ds_{old_new}.npy', mmap_mode='r')

	print('y')
	y_ds = do_dedisp(y_fname, dedisp_func)
	np.save(f'y_ds_{old_new}.npy', y_ds)
	del y_ds
	gc.collect()
	y_ds = np.load(f'y_ds_{old_new}.npy', mmap_mode='r')

	ds.save_stokes_dynspec(x_ds, y_ds, f'compare_dedisp_{old_new}', dir='.')


old_new = sys.argv[1]
x_fname = sys.argv[2]
y_fname = sys.argv[3]

dedisp_func = old_coh_dedisp if old_new == 'old' else new_coh_dedisp

do_dedisp_wrapper(old_new, x_fname, y_fname, dedisp_func)