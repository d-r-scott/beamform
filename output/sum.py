"""
Sums time series together for a given FRB
"""

import numpy as np
from argpase import ArgumentParser, ArgumentDefaultsHelpFormatter
import glob

def _main():
	args = get_args()
	
	# Check we have a valid polarisation
	if args.p in ('x', 'y', 'both'):
		x_fnames, y_fnames = find_files(args.FRB, args.p)

		# the above will either be lists of filenames or None, depending on if we're doing that polarisation

		if x_fnames is not None:
			x_sum = do_sum(x_fnames)
			save(x_sum, args.o, 'x_')

		if y_fnames is not None:
			y_sum = do_sum(y_fnames)
			save(y_sum, args.o, 'y_')

	else:
		print("Invalid polarisation")

def get_args()
	parser = ArgumentParser(description="Sums time series together for a given FRB", formatter_class=ArgumentDefaultsHelpFormatter)

	parser.add_argument(dest='FRB', type=str, help='FRB to sum')
	parser.add_argument('-c', type='store_true', help='Use cropped time series instead of full', default=False)
	# TODO: Implement -c flag
	parser.add_argument('-o', type=str, help='Output filename (polarisation will be appended as prefix)', default='sum_t.npy')
	parser.add_argument('-p', type=str, help='Polarisation to do (default is both)', default='both')

	return parser.parse_args()

def find_files(FRB):
	glob_str_x = '{}/*x*t.npy'.format(FRB)
	glob_str_y = '{}/*y*t.npy'.format(FRB)

	x_fnames = glob.glob(glob_str_x)
	y_fnames = glob.glob(glob_str_y)

	if args.p == 'both':
		# Return lists for x and y
		return x_fnames, y_fnames
	elif args.p == 'x':
		return x_fnames, None
	else:	# must be 'y'
		return None, y_fnames

def do_sum(fnames):
	# Sum one file at a time, since they're very large
	# Numpy's mmap_mode helps by allowing us to do this without moving them into memory

	# Initialise summation array with the first file
	sum_arr = np.load(fnames[0]).copy()

	# Go over the rest of the files, summing their contents into sum_arr
	for fname in fnames[1:]:
		# mmap_mode='r' keeps the loaded array on disk, instead of loading it into memory
		sum_arr += np.load(fname, mmap_mode='r')
	
	return sum_arr

def save(arr, outfile, prefix):
	dest_fname = prefix + outfile
	np.save(darr,dest_fname)

if __name__ == "__main__":
	_main()
