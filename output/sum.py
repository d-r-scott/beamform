"""
Sums time series together for a given FRB
"""

import numpy as np
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import glob

def _main():
	args = get_args()
	
	# Check we have a valid polarisation
	if args.p in ('x', 'y', 'both'):
		x_f_fnames, y_f_fnames, x_t_fnames, y_t_fnames = find_files(args.FRB, args.p)

		# the above will either be lists of filenames or None, depending on if we're doing that polarisation

		if x_f_fnames is not None:
			#x_f_sum = do_sum(x_f_fnames)
			#save(x_f_sum, args.o, 'x_f_')
			x_t_sum = do_sum(x_t_fnames)
			save(x_t_sum, args.o, 'x_t_')

		if y_f_fnames is not None:
			#y_f_sum = do_sum(y_f_fnames)
			#save(y_f_sum, args.o, 'y_f_')
			y_t_sum = do_sum(y_t_fnames)
			save(y_t_sum, args.o, 'y_t_')

	else:
		print("Invalid polarisation")

def get_args():
	parser = ArgumentParser(description="Sums time series together for a given FRB", formatter_class=ArgumentDefaultsHelpFormatter)

	parser.add_argument(dest='FRB', type=str, help='FRB to sum')
	#parser.add_argument('-c', type='store_true', help='Use cropped time series instead of full', default=False)
	# TODO: Implement -c flag
	parser.add_argument('-o', type=str, help='Output filename (polarisation and t/f will be appended as prefix)', default='sum_t.npy')
	parser.add_argument('-p', type=str, help='Polarisation to do (default is both)', default='both')

	return parser.parse_args()

def find_files(FRB, p):
	glob_str_x_f = '{}/*x*f.npy'.format(FRB)
	glob_str_y_f = '{}/*y*f.npy'.format(FRB)
	glob_str_x_t = '{}/*x*t.npy'.format(FRB)
	glob_str_y_t = '{}/*y*t.npy'.format(FRB)

	x_f_fnames = glob.glob(glob_str_x_f)
	y_f_fnames = glob.glob(glob_str_y_f)
	x_t_fnames = glob.glob(glob_str_x_t)
	y_t_fnames = glob.glob(glob_str_y_t)

	if p == 'both':
		# Return lists for x and y
		return x_f_fnames, y_f_fnames, x_t_fnames, y_t_fnames
	elif p == 'x':
		return x_f_fnames, None, x_t_fnames, None
	else:	# must be 'y'
		return None, y_f_fnames, None, y_t_fnames

def do_sum(fnames):
	# Sum one file at a time, since they're very large
	# Numpy's mmap_mode helps by allowing us to do this without moving them into memory

	print(fnames)

	# Initialise summation array with the first file
	print("Initialising sum array...")
	sum_arr = np.load(fnames[0]).copy()

	# Go over the rest of the files, summing their contents into sum_arr
	for fname in fnames[1:]:
		print(fname)
		# mmap_mode='r' keeps the loaded array on disk, instead of loading it into memory
		sum_arr += np.load(fname, mmap_mode='r')
	
	return sum_arr

def save(arr, outfile, prefix):
	print("Saving...")
	dest_fname = prefix + outfile
	np.save(dest_fname, arr)

if __name__ == "__main__":
	_main()
