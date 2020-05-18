"""
A collection of functions for creating, manipulating, and plotting dynamic spectra from polarimetric data
"""

import numpy as np
import matplotlib.pyplot as plt

def generate_dynspec(t_ser):
	n = 336
	dynspec = np.zeros((int(t_ser.shape[0]/n), n), dtype=np.complex64)
	for i in range(int(t_ser.shape[0]/n)):
		dynspec[i, :] = np.fft.fft(t_ser[i*n:(i+1)*n])
	return dynspec

def IQUV(x, y):
	i = np.abs(x)**2 + np.abs(y)**2
	q = np.abs(x)**2 - np.abs(y)**2
	u = 2*np.real(np.conj(x) * y)
	v = 2*np.imag(np.conj(x) * y)
	return i, q, u, v

def save_stokes_dynspec(x, y, frb):
	print("Generating x dynspec")
	x_ds = generate_dynspec(x)
	#if y != 0:
	print("Generating y dynspec")
	y_ds = generate_dynspec(y)
	#else:
	#	y_ds = 0
	print("Calculating Stokes parameters")
	i_ds, q_ds, u_ds, v_ds = IQUV(x_ds, y_ds)
	stk_str = ['i', 'q', 'u', 'v']
	for j, par in enumerate((i_ds, q_ds, u_ds, v_ds)):
		print("Saving {}_ds_{}.npy".format(stk_str[j], frb))
		np.save('{}_ds_{}.npy'.format(stk_str[j], frb), par)
	del i_ds, q_ds, u_ds, v_ds

def load_stokes_dynspec(frb):
	stk_str = ['i', 'q', 'u', 'v']
    ret = ( np.load('{}/{}_ds_{}.npy'.format(frb, stk, frb) ) for stk in stk_str )
    return ret
