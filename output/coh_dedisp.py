#!/usr/bin/env python3

import numpy as np
import dynspec as ds
import sys
import gc
sys.path.append('..')
import freq2time as ftt


def _main():
	args = get_args()

	process(args.x, 'x', args)
	process(args.y, 'y', args)


def get_args():
	"""
	Get arguments from command line

	:returns: args object with arguments as fields
	"""
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description='Coherently dedisperses polarimetric fourier spectra, iffts and saves result.\nOptionally converts to dynamic spectra.', formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument('-x', help='Input x fourier spectrum')
	parser.add_argument('-y', help='Input y fourier spectrum')
	parser.add_argument('-d', '--DM', type=float, help='DM to coherently dedisperse to')
	parser.add_argument('-f', type=float, help='Central frequency of band (MHz)')
	parser.add_argument('-b', type=float, help='Bandwidth (MHz)')
	parser.add_argument('--dynspec', help='Flag to enable conversion to dynamic spectra', action='store_true', default=False)
	parser.add_argument('-o', help='Output filename. Will be prefixed with [xy]_[ft]_, and appended with _cohDM{DM}(_ds).npy')
	args = parser.parse_args()

	return args


def save(arr, o, xy, ft, DM, ds):
	"""
	Saves a given arr with specified filename parameters.
	Filename format:
		if not ds:
			{xy}_{ft}_{o}_cohDM{DM}.npy
		else:
			{xy}_{ft}_{o}_cohDM{DM}_ds.npy

	:param arr: Array to be saved
	:param o: Base of filename
	:param xy: One of 'x' or 'y'
	:param ft: One of 'f' or 't'
	:param DM: DM dedispersed to
	:param ds: Boolean indicating if _ds suffix should be included
	"""
	fname = f'{xy}_{ft}_{o}_cohDM{DM}_ds.npy' if ds else f'{xy}_{ft}_{o}_cohDM{DM}.npy'
	print(f'Saving {fname}...')
	np.save(fname, arr)


def process(fname, xy, args):
	"""
	Processes the given fourier spectrum

	:param fname: filename for fourier spectrum
	:param xy: 'x' or 'y' indicating which polarisation is being processed
	:param args: Command line arguments
	"""

	print(xy)
	f = np.load(fname)
	f_dd = ftt.coh_dedisp(f, args.DM, f_mid=args.f, bw=args.b)
	del f
	t_dd = ftt.ifft_long(f_dd)[0, :, 0]
	if args.dynspec:
		t_dd_ds = ds.generate_dynspec(t_dd)
	save(f_dd, args.o, xy, 'f', args.DM, False)
	del f_dd
	save(t_dd, args.o, xy, 't', args.DM, False)
	del t_dd
	if args.dynspec:
		save(t_dd_ds, args.o, xy, 't', args.DM, True)
	del t_dd_ds
	gc.collect()


if __name__ == '__main__':
	_main()