"""
A collection of functions for creating, manipulating, and plotting dynamic spectra from polarimetric data
"""

import numpy as np
import matplotlib.pyplot as plt
import progressbar
import gc

def generate_dynspec(t_ser):
	"""
	Creates a dynamic spectrum at the highest time resolution from the given time series.
	Assumes 336 frequency channels.

	:param t_ser: input time series of voltages
	:return: dynamic spectrum of voltages
	"""
	n = 336
	dynspec = np.zeros((int(t_ser.shape[0]/n), n), dtype=np.complex64)
	with progressbar.ProgressBar(max_value=int(t_ser.shape[0]/n)) as bar:
		for i in range(int(t_ser.shape[0]/n)):
			dynspec[i, :] = np.fft.fft(t_ser[i*n:(i+1)*n])
			bar.update(i)

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

def save_stokes_dynspec(x_ds, y_ds, frb, dir=None):
	"""
	Generates X and Y dynamic spectra, calculates Stokes IQUV dynamic spectra, and saves IQUV to save on memory.
	If one of the polarisations is not available, import y as None.
	Also normalises and transposes the dynspec so that it's usable.

	:param x_ds: Input x dynspec
	:param y_ds: Input y dynspec
	:param frb: String with name of FRB, for the filenames
	:return: None
	"""

	if dir is None:
		dir=frb

	print("Calculating Stokes parameters")

	# lambda functions for each of the Stokes parameters
	stokes = {
		"i" : lambda x, y: np.abs(x)**2 + np.abs(y)**2,
		"q" : lambda x, y: np.abs(x)**2 - np.abs(y)**2,
		"u" : lambda x, y: 2*np.real(np.conj(x) * y),
		"v" : lambda x, y: 2*np.imag(np.conj(x) * y)
	}

	stk_str = ['i', 'q', 'u', 'v']

	for stk in stk_str:
		print(f'Calculating {stk}')
		par = stokes[stk](x_ds, y_ds)
		par_norm = normalise(par)
		del par
		par_tran = par_norm.transpose()
		del par_norm

		print(f'Saving {stk}_ds_{frb}.npy')
		np.save(f'{dir}/{stk}_ds_{frb}.npy', par_tran)
		del par_tran
		gc.collect()

def load_stokes_dynspec(frb, dir=None):
	"""
	Loads stokes parameters from provided FRB directory:
		{frb}/[iquv]_ds_{frb}.npy

	:param frb: String with name of FRB, for the filenames
	:param dir: If the FRB directory is different from the FRB name, provide dir with the directory name
	:return: A tuple (I, Q, U, V)
	"""
	if dir is None:
		dir=frb

	stk_str = ['i', 'q', 'u', 'v']
	ret = ( np.load('{}/{}_ds_{}.npy'.format(dir, stk, frb), mmap_mode='r' ) for stk in stk_str )
	return ret


def reduce(A, n, transpose=False):
	"""
	Reduces the time resolution of a given array by a factor n.
	Reduces along the 0th axis, make sure this one is the time axis!

	:param A: Input array to be reduced
	:param n: Factor to reduce by
	:param transpose: Does the function need to transpose the array before reducing? Set to True iff axis 0 is not time
	:return: None
	"""

	if n > 1:	# if n == 1 (or below, and is therefore invalid) don't bother reducing and just return the input
		A = A.transpose() if transpose else A
		A_red = []
		with progressbar.ProgressBar(max_value=int(A.shape[0]/n)) as bar:
			for i in range(int(A.shape[0]/n)):
				A_red.append(np.sum(A[i*n:(i+1)*n], axis=0))
				bar.update(i)
		return np.array(A_red).transpose() if transpose else np.array(A_red)
	else:
		return(A)

def plot_dynspec(dynspec,
				 xmin=None, xmax=None, ymin=None, ymax=None,
				 xlabel=None, ylabel=None,
				 vmin=None, vmax=None,
				 freq_range=None, tick_space=None, axis_res=None, peak_idx=None,
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
	:param freq_range: 2-tuple (f_min, f_max) of ranges for nicer axis ticks
	:param tick_space: 2-tuple (t_space, f_space) of tick frequency in pixels for nicer axis ticks
	:param axis_res: 2-tuple (t_res, f_res) of axis resolutions (in ms and MHz respectively) for nicer axis ticks
	:param peak_idx: Index of burst peak
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

	t_ser = [ np.sum(dynspec[:, t]) for t in range(T) ]
	t_ax.plot(range(T), t_ser, c='k', lw=0.5)

	f_ser = [ np.sum(dynspec[f, :]) for f in range(F) ]
	f_ax.plot(f_ser, range(F), c='k', lw=0.5)

	# If freq_range, tick_space, axis_res, and peak_idx are given, make the tick marks nicer
	if freq_range is not None and tick_space is not None and axis_res is not None and peak_idx is not None:
		f_min, f_max = freq_range
		t_space, f_space = tick_space
		t_res, f_res = axis_res

		# Time axis
		#	Axis positions - 10 pre- and post- peak of time series (peak is at 0), t_space apart
		first_xtick = peak_idx - 10*t_space
		last_xtick = peak_idx + 10*t_space

		xticks = np.arange(first_xtick, last_xtick+1, t_space)
		ds_ax.set_xticks(xticks)

		#	Tick labels
		xtick_labels = [ '%.0f'%(lbl) for lbl in np.arange(-10*t_space*t_res, 10*t_space*t_res+1, t_space*t_res) ]
		ds_ax.set_xticklabels(xtick_labels)

		# Frequency axis
		#	Axis positions - One tick every f_space, on values where freq % f_res == 0
		#					 i.e. if f_res == 25, put ticks at x00, x25, x50, x75...
		#					 Once upon a time you wrote down why these are the way they are. I hope you kept your notes
		#					 Hint: the terms with int in them are effectively floor and ceil
		first_ytick = f_max - int(f_max/f_res)*f_res
		last_ytick = (F-1) - ( (int(f_min/f_res + 1) * f_res) - f_min )

		yticks = np.arange(last_ytick, first_ytick, -1*f_space)	# this is backwards because imshow be weird
		ds_ax.set_yticks(yticks)

		#	Tick labels
		ytick_labels = [ '%.0f'%(lbl) for lbl in np.arange((int(f_min/f_space)+1) * f_space, int(f_max/f_space)*f_space+1, f_space*f_res) ]
		ds_ax.set_yticklabels(ytick_labels)
		
	elif freq_range is not None or tick_space is not None or axis_res is not None or peak_idx is not None:
		print("ERROR: Need all of freq_range, tick_space, axis_res, and peak_idx to make axes nice!")

	ds_ax.set_xlim(left=xmin+0.5, right=xmax-0.5)
	ds_ax.set_ylim(bottom=ymax-0.5, top=ymin+0.5)
	ds_ax.set_xlabel(xlabel)
	ds_ax.set_ylabel(ylabel)

	plt.tight_layout()
	plt.show()

def normalise(ds):
	"""
	Normalises a dynamic spectrum along frequency channels to cancel out RFI

	:param ds: Input dynspec to normalise
	"""
	T, F = ds.shape
	out_ds = ds.copy()
	norm = lambda x: (x - np.mean(x))/np.std(x)
	with progressbar.ProgressBar(max_value = F) as bar:
		for f in range(F):
			out_ds[:,f] = norm(ds[:,f])
			bar.update(f)
	del ds
	return out_ds
