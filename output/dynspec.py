"""
A collection of functions for creating, manipulating, and plotting dynamic spectra from polarimetric data
"""

import numpy as np
import matplotlib.pyplot as plt

def generate_dynspec(t_ser):
	"""
	Creates a dynamic spectrum at the highest time resolution from the given time series.
	Assumes 336 frequency channels.

	:param t_ser: input time series of voltages
	:return: dynamic spectrum of voltages
	"""
	n = 336
	dynspec = np.zeros((int(t_ser.shape[0]/n), n), dtype=np.complex64)
	for i in range(int(t_ser.shape[0]/n)):
		dynspec[i, :] = np.fft.fft(t_ser[i*n:(i+1)*n])
	return dynspec

def IQUV(x, y):
	"""
	Calculates and returns the Stokes polarisation parameters I, Q, U, and V from input X and Y polarisations.

	:param x: Input X polarisation
	:param y: Input Y polarisation
	:return: A tuple (I, Q, U, V)
	"""
	i = np.abs(x)**2 + np.abs(y)**2
	q = np.abs(x)**2 - np.abs(y)**2
	u = 2*np.real(np.conj(x) * y)
	v = 2*np.imag(np.conj(x) * y)
	return i, q, u, v

def save_stokes_dynspec(x, y, frb):
	"""
	Generates X and Y dynamic spectra, calculates Stokes IQUV dynamic spectra, and saves IQUV to save on memory.
	If one of the polarisations is not available, import y as None.

	:param x: Input x time series
	:param y: Input y time series
	:param frb: String with name of FRB, for the filenames
	:return: None
	"""
	print("Generating x dynspec")
	x_ds = generate_dynspec(x)

	if y is None:
		print("Y not given, setting to 0")
		y_ds = 0
	else:
		print("Generating y dynspec")
		y_ds = generate_dynspec(y)

	print("Calculating Stokes parameters")
	i_ds, q_ds, u_ds, v_ds = IQUV(x_ds, y_ds)

	stk_str = ['i', 'q', 'u', 'v']

	for j, par in enumerate((i_ds, q_ds, u_ds, v_ds)):
		print("Saving {}_ds_{}.npy".format(stk_str[j], frb))
		np.save('{}_ds_{}.npy'.format(stk_str[j], frb), par)
		# Delete once we're done with the parameter so it doesn't waste memory
		del par

def load_stokes_dynspec(frb):
	"""
	Loads stokes parameters from provided FRB directory:
		{frb}/[iquv]_ds_{frb}.npy

	:param frb: String with name of FRB, for the filenames
	:return: A tuple (I, Q, U, V)
	"""
	stk_str = ['i', 'q', 'u', 'v']
	ret = ( np.load('{}/{}_ds_{}.npy'.format(frb, stk, frb), mmap_mode='r' ) for stk in stk_str )
	return ret


def reduce(A, n):
	"""
	Reduces the time resolution of a given array by a factor n.
	Reduces along the 0th axis, make sure this one is the time axis!

	:param A: Input array to be reduced
	:param n: Factor to reduce by
	:return: None
	"""
	A_red = []
	for i in range(int(A.shape[0]/n)):
		A_red.append(np.sum(A[i*n:(i+1)*n], axis=0))

	return np.array(A_red)

def plot_dynspec(dynspec,
				 xmin=None, xmax=None, ymin=None, ymax=None,
				 xlabel=None, ylabel=None,
				 vmin=None, vmax=None,
				 cmap='inferno'
):
	"""
	Plots a dynamic spectrum in a nice way with summed axes on the sides.

	:param dynspec: Dynamic spectrum to plot, with shape (F, T) where F and T are the lengths of the frequency and time
	                axes respectively
	:param xmin: Min time sample to crop to
	:param xmax: Max time sample to crop to
	:param ymin: Min freq channel to crop to
	:param ymax: Max freq channel to crop to
	:param xlabel: Label for x (time) axis
	:param ylabel: Label for y (freq) axis
	:param vmin: Minimum value for colormap
	:param vmax: Maximum value for colormap
	:param cmap: Colormap to plot with [default: inferno]
	:return: None
	"""
	rows = 3
	cols = 5
	ax_shape = (rows, cols)

	F, T = dynspec.shape

	xmin = xmin if xmin else 0
	xmax = xmax if xmax else T
	ymin = ymin if ymin else 0
	ymax = ymax if ymax else F

	ds_ax = plt.subplot2grid(ax_shape, (1, 0), colspan = cols - 1, rowspan = rows - 1)  #dyanmic spectrum
	t_ax = plt.subplot2grid(ax_shape, (0, 0), colspan = cols - 1, sharex = ds_ax)   #time series
	f_ax = plt.subplot2grid(ax_shape, (1, cols - 1), rowspan = rows - 1, sharey = ds_ax)	#freq series

	ds_ax.imshow(dynspec, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax, origin='lower', interpolation='none')
	ds_ax.set_xlim(left=xmin, right=xmax)
	ds_ax.set_ylim(bottom=ymax, top=ymin)
	ds_ax.set_xlabel(xlabel)
	ds_ax.set_ylabel(ylabel)

	t_ser = [ np.sum(dynspec[:, t]) for t in range(T) ]
	t_ax.plot(range(T), t_ser, c='k', lw=0.5)

	f_ser = [ np.sum(dynspec[f, :]) for f in range(F) ]
	f_ax.plot(f_ser, range(F), c='k', lw=0.5)

	plt.tight_layout()
	plt.show()
