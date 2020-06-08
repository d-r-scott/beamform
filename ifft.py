#!/usr/bin/env python3

import numpy as np
from scipy.fft import ifft

def _main():
	args = get_args()
	f = load(args.f)
	t = do_ifft(f)
	save(t, args.o)


def get_args():
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description='Performs ifft on given spectrum to obtain time series')
	parser.add_argument('-f', help='Spectrum file to ifft')
	parser.add_argument('-o', help='Output file to save time series to')
	return parser.parse_args()


def load(fname):
	print(f'Loading {fname}')
	return np.load(fname)


def do_ifft(f):
	print('IFFTing')
	return ifft(f)


def save(t, fname):
	print(f'Saving {fname}')
	np.save(fname, t)


if __name__ == '__main__':
	_main()