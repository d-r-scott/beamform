#!/usr/bin/env python3

import numpy as np

def _main():
	args = get_args()
	f = load(args.f)
	f_dd = dedisperse(f, args.DM, args.f0, args.bw)
	save(f_dd, args.o)


def get_args():
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description='Coherently dedisperses given fine spectrum', formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument('-f', help='Fine spectrum file to dedisperse')
	parser.add_argument('--DM', type=float, help='Dispersion measure to dedisperse to in pc/cm3')
	parser.add_argument('--f0', type=float, help='Central frequency of observation in MHz')
	parser.add_argument('--bw', type=float, help='Bandwidth of observation in MHz', default=336)
	return parser.parse_args()


def load(fname):
	return np.load(fname)


def dedisperse(f, DM, f0, bw):
	"""
	Takes heavy inspiration from Hyerin Cho's coh_dedisp function from freq2time.py
	"""
	n_sam = len(f)
	k_DM = 4818.808


if __name__ == '__main__':
	_main()