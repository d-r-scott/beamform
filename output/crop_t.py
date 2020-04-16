"""
Crops a given time series (assumed to be very large, ~GB) to the region around a specified time index.
Intended for use on full time resolution voltage dumps of FRB detections, post-coherent dedispersion.
The time index is generally the offset of the burst.
"""

import numpy as np
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def _main():
	args = get_args()
	full_t_ser = load(args.fname)
	crop_t_ser = crop(full_t_ser, args.t, args.b, args.a)
	save(crop_t_ser, args.o)

def get_args():
	parser = ArgumentParser(description='Crops a given time series around the provided time index', formatter_class=ArgumentDefaultsHelpFormatter)

	parser.add_argument(dest='fname', type=str, help='full time series .npy file')
	parser.add_argument('--t0', '-t', dest='t', type=int, help='time index to crop about [default=0]', default=0)
	parser.add_argument('--bef', '-b', dest='b', type=int, help='number of samples to include before t0 [default=100]', default=100)
	parser.add_argument('--aft', '-a', dest='a', type=int, help='number of samples to include after t0 [default=1000]', default=1000)
	parser.add_argument('-o', dest='o', type=str, help='output filename [default=./crop_t_ser.npy]', default='./crop_t_ser.npy')

	return parser.parse_args()

def load(fname):
	print("Loading {}...".format(fname))
	return np.load(fname, mmap_mode='r')

def crop(full_t_ser, t, b, a):
	print("Cropping...")
	crop_t_ser = full_t_ser[t-b:t+a].copy()
	del full_t_ser	# These usually are huuuge, delete to be safe
	return crop_t_ser
	
def save(crop_t_ser, dest_fname):
	print("Saving {}...".format(dest_fname))
	np.save(dest_fname, crop_t_ser)
	print("Saved!")

if __name__ == "__main__":
	_main()
