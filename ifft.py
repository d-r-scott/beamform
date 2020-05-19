import numpy as np
from scipy.fft import ifft
import sys

FFFF = np.load(sys.argv[1])
print('iffting...')
t_series = ifft(FFFF)
np.save(sys.argv[2], t_series)
