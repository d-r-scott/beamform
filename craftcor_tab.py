#!/usr/bin/env python
"""
Tied-array beamforming vcraft files, based on "craftcor.py".

Copyright (C) CSIRO 2017
"""

__author__ = ['Keith Bannister, CSIRO <keith.bannister@csiro.au>',
              'Hyerin Cho, Curtin University/Swinburne University '
              + '<chyerin1996@gmail.com>',
              'David Scott, Curtin University '
              + '<david.r.scott@postgrad.curtin.edu.au']

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
import multiprocessing
import numpy as np
import os
import signal
import time

from astropy.coordinates import SkyCoord

from calc11 import ResultsFile
from miriad import MiriadGainSolutions
import vcraft


# constants
C_LIGHT = 299792458.0   # Speed of light (m/s)
CHAN_BWIDTH = 27.0      # Channel bandwidth
OS_NYQ_BWIDTH = 32.0    # Oversampled Nyquist bandwidth
F_OS = OS_NYQ_BWIDTH / CHAN_BWIDTH  # Oversampling factor
NUM_GUARD_CHAN = OS_NYQ_BWIDTH - CHAN_BWIDTH    # Number of guard channels

# TODO: Cleanup (not in any particular order)
#   IN GENERAL: MAKE THINGS CONSISTENT AND NICE!
#   (1) PEP-008ify everything
#         - Probably do first, easiest to do the others with readable code!
#   (2) Document all functions and classes
#         - Includes docstrings and line comments
#   (3) Ensure only necessary functions remain
#   (4) Remove magic numbers!
#   (5) Make as much as possible compatible with Python 3


class AntennaSource(object):
    # TODO: (1, 2, 4, 5)
    def __init__(self, vfile):
        # TODO: (2, 4, 5)
        self.vfile = vfile
        self.ant_name = self.vfile.hdr['ANT'][0].lower()
        self.antno = int(self.vfile.hdr['ANTENNA_NO'][0])
        self.mjd_start = self.vfile.start_mjd
        self.trigger_frame = self.vfile.start_frameid
        self.hdr = self.vfile.hdr
        self.init_geom_delay_us = None
        self.all_geom_delays = []
        self.all_mjds = []
        self.pol = self.vfile.pol.lower()
        self.fringe_rot_params = None
        print('antenna {} {}'.format(self.ant_name, self.vfile.freqconfig))

    def do_f_tab(self, corr, i_ant):
        """Perform the tied-array beamforming for this antenna.

        :param corr: Correlator object containing constants and
                     functions required for beamforming
        :param i_ant: Incremental antenna number
        :return: beamformed complex voltage time series
        """
        # TODO: (2, 4, 5)
        self.fringe_rot_params = FringeRotParams(corr, self)

        # calculate sample start
        framediff_samp = corr.ref_ant.trigger_frame - self.trigger_frame
        framediff_us = framediff_samp / corr.fs

        geom_delay_us, geom_delay_rate_us = \
            corr.get_geometric_delay_delayrate_us(self)
        self.all_geom_delays.append(geom_delay_us)
        self.all_mjds.append(corr.curr_mjd_mid)

        geom_delay_samp = geom_delay_us * corr.fs
        fixed_delay_us = corr.get_fixed_delay_usec(self.antno)
        fixed_delay_samp = fixed_delay_us * corr.fs
        total_delay_samp = framediff_samp
        whole_delay = int(np.round(total_delay_samp))
        total_delay_us = total_delay_samp / corr.fs
        whole_delay_us = whole_delay / corr.fs

        frac_delay_samp = total_delay_samp - whole_delay
        frac_delay_us = frac_delay_samp * corr.fs

        data_out = np.zeros((corr.n_int, corr.n_fine_chan, corr.n_pol_in),
                            dtype=np.complex64)
        n_fine = corr.n_fft - 2 * corr.nguard_chan

        n_samp = corr.n_int * corr.n_fft

        # time-dependent geometric delays
        # np.linspace(0, 1, n_samp) == time in units of integrations
        geom_delays_us = (geom_delay_us
                         + geom_delay_rate_us
                         * np.linspace(0, 1, n_samp)
                         - fixed_delay_us)

        np.save('delays/geom_delays_us_{}'.format(i_ant), geom_delays_us)

        print('')
        print_var('framediff_samp', framediff_samp)
        print_var('framediff_us', framediff_us)
        print_var('geom_delay_us', geom_delay_us)
        print_var('geom_delay_rate_us', geom_delay_rate_us)
        print_var('geom_delay_samp', geom_delay_samp)
        print_var('np.mean(geom_delays_us)', np.mean(geom_delays_us))
        print_var('fixed_delay_us', fixed_delay_us)
        print_var('fixed_delay_samp', fixed_delay_samp)
        print_var('total_delay_samp', total_delay_samp)
        print_var('whole_delay', whole_delay)
        print_var('total_delay_us', total_delay_us)
        print_var('whole_delay_us', whole_delay_us)
        print_var('frac_delay_samp', frac_delay_samp)
        print_var('frac_delay_us', frac_delay_us)
        print('')

        print("antenna #: ", i_ant, self.ant_name)
        sample_offset = whole_delay + corr.abs_delay
        frameid = self.vfile.start_frameid + sample_offset
        print('FRAMEID: ' + str(frameid) + ', remainder from 32: '
              + str(frameid % 32))

        '''
        To avoid iPFB fractional delay, set FRAMEID such that the 
        remainder is 0
        '''

        raw_data = self.vfile.read(sample_offset, n_samp)

        assert raw_data.shape == (n_samp, corr.ncoarse_chan),\
            'Unexpected shape from vfile: {} expected ({},{})'.format(
                raw_data.shape, n_samp, corr.ncoarse_chan)

        turn_fracs = np.zeros((336, n_fine+1))

        for i, chan in enumerate(xrange(corr.ncoarse_chan)):
            # Channel frequency
            centre_freq = corr.freqs[chan]
            turn_fracs[i, 0] = centre_freq

            '''
            Array of fine channel frequencies relative to centre frequency
            More appropriately named Delta_f?
            Goes from -0.5 to 0.5
            '''
            delta_freq = (corr.fine_chan_bwidth
                          * (np.arange(n_fine, dtype=np.float)
                          - float(n_fine)/2.0))

            # TODO: what is sideband?
            if corr.sideband == -1:
                delta_freq = -delta_freq

            '''
            raw_data's shape: (n_samp, corr.ncoarse_chan)
            n_samp = input.i * (64 * input.n)
            '''
            x1 = raw_data[:, chan].reshape(-1, corr.n_fft)

            '''
            fixed_delay_us = corr.get_fixed_delay_usec(self.antno)
            fixed_delay_us is contained in fcm.txt for each antenna
            geom_delay_us = corr.get_geometric_delay_delayrate_us(self)[0]
            geom_delay_us accounts for Earth's rotation
            '''
            # Fringe rotation for Earth's rotation
            turn_fringe = centre_freq * geom_delays_us
            phasor_fringe = np.exp(2j * np.pi * turn_fringe,
                                   dtype=np.complex64)
            x1 *= phasor_fringe

            '''
            corr.nguard_chan = NUM_GUARD_CHAN * input.n 
                             = number of fine channels on each side to cut off
            xfguard is xf1 with the ends trimmed off
            '''
            xf1 = np.fft.fft(x1, axis=1)
            xf1 = np.fft.fftshift(xf1, axes=1)
            # scale because otherwise it overflows
            xfguard_f = xf1[:, corr.nguard_chan:corr.nguard_chan+n_fine:]

            # Fractional sample phases
            turn_frac = delta_freq * np.mean(geom_delays_us)
            turn_fracs[i, 1:] = turn_frac

            # phasors to rotate the data with constant amplitude = 1
            phasor = np.exp(np.pi * 2j * turn_frac, dtype=np.complex64)

            # get absolute frequencies in gigahertz
            freq_ghz = (centre_freq+delta_freq) / 1e3

            # get the calibration solutions and apply them to the phasors
            mir_cor = corr.mir.get_solution(i_ant, 0, freq_ghz)

            if mir_cor[0] == 0: # if correction is 0, flag data
                phasor *= 0
            else:
                phasor /= mir_cor

            xfguard_f *= phasor

            # select the channels for this coarse channel
            fine_chan_start = chan * n_fine
            fine_chan_end = (chan+1) * n_fine

            '''
            RECAP
            xfguard is a "dynamic" spectrum of the current coarse 
            channel in only the fine channels not trimmed 
            (oversampled PFB).
            xfguard.shape == (input.i, 64 * input.n)
            '''

            '''
            Slot xfguard (a trimmed spectrum for this coarse channel) 
            into the corresponding slice of the fine channels
            '''
            data_out[:, fine_chan_start:fine_chan_end, 0] = xfguard_f

        np.save('delays/turn_fracs_{}'.format(i_ant), turn_fracs)

        return data_out


class FringeRotParams(object):
    # TODO: (1, 2, 4, 5)
    cols = ('U (m)', 'V (m)', 'W (m)', 'DELAY (us)')

    def __init__(self, corr, ant):
        # TODO: (1, 2, 4, 5)
        mid_data = corr.fringe_rot_data_mid[ant.ant_name]
        self.u,self.v,self.w,self.delay = map(float, [mid_data[c] for c in FringeRotParams.cols])
        self.delay_start = float(corr.fringe_rot_data_start[ant.ant_name]['DELAY (us)'])
        self.delay_end = float(corr.fringe_rot_data_end[ant.ant_name]['DELAY (us)'])
        self.delay_rate = (self.delay_end - self.delay_start)/float(corr.n_int)
        self.ant = ant
        self.corr = corr

    def __str__(self):
        # TODO: (1, 2, 4, 5)
        s = 'FR {} uvw=({},{},{}) m = {} us'.format(self.ant.ant_name, self.u, self.v, self.w, self.delay)
        return s

    __repr__ = __str__


class Correlator(object):
    # TODO: (1, 2, 4, 5)
    def __init__(self, ants, sources, values, abs_delay=0):
        # TODO: (1, 2, 4, 5)
        self.running = True
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)
        self.ants = ants
        self.values = values
        self.pool = None
        if self.values.num_threads > 1:
            self.pool = multiprocessing.Pool(processes=values.num_threads)

        self.parse_par_set()

        for ia, a in enumerate(self.ants):
            a.ia = ia
            a.antpos = self.get_ant_location(a.antno)

        self.abs_delay = abs_delay
        self.ref_ant = ants[0]
        self.calc_results = ResultsFile(values.calcfile)
        self.d_utc = 0
        self.mjd0 = self.ref_ant.mjd_start + self.d_utc / 86400.0
        self.frame0 = self.ref_ant.trigger_frame
        self.n_int = values.n_int
        self.n_fft = 64 * values.fft_size
        self.nguard_chan = NUM_GUARD_CHAN * values.fft_size
        self.oversamp = F_OS
        self.fs = self.oversamp     # samples per microsecond
        self.ncoarse_chan = len(self.ref_ant.vfile.freqs)
        self.sideband = -1
        self.coarse_chan_bwidth = 1.0
        self.n_fine_per_coarse = self.n_fft - 2 * self.nguard_chan
        self.n_fine_chan = self.ncoarse_chan * self.n_fine_per_coarse
        self.fine_chan_bwidth = (self.coarse_chan_bwidth
                                 / float(self.n_fine_per_coarse))
        self.full_bw = self.fine_chan_bwidth * self.n_fine_chan
        self.freq_scrunch = values.freq_scrunch
        assert self.freq_scrunch >= 1
        assert self.n_fine_per_coarse % self.freq_scrunch == 0, \
            'Freq scrunch must yield an integer number of fine channels ' \
            'per coarse channel'
        self.n_fine_out_per_coarse = self.n_fine_per_coarse / self.freq_scrunch
        self.n_fine_out_chan = self.n_fine_out_per_coarse * self.ncoarse_chan
        self.out_chan_bwidth = (self.coarse_chan_bwidth
                                / float(self.n_fine_out_per_coarse))
        self.n_pol_in = 1
        self.n_pol_out = 1
        self.f0 = self.ants[0].vfile.freqs[0]
        self.freqs = self.ants[0].vfile.freqs
        self.freq_mid = self.freqs.mean()
        self.int_time_secs = float(self.n_int * self.n_fft) / (self.fs * 1e6)
        self.int_time_days = self.int_time_secs / 86400.
        self.curr_int_no = 0
        self.curr_samp = self.curr_int_no * self.n_int + 1000
        self.calc_mjd()
        self.get_fringe_rot_data()
        self.pol = self.ants[0].pol
        self.parse_mir()
        self.par_set = None
        self.mir = None
        self.curr_mjd_start = None
        self.curr_mjd_mid = None
        self.curr_mjd_end = None
        self.fringe_rot_data_start = None
        self.fringe_rot_data_mid = None
        self.fringe_rot_data_end = None

        logging.debug('F0 %f FINE CHANNEL %f kHz num=%d freqs=%s', self.f0,
                      self.fine_chan_bwidth * 1e3, self.n_fine_chan,
                      self.freqs)

    def exit_gracefully(self, signum, frame):
        # TODO: (1, 2, 4, 5)
        self.running = False

    def parse_par_set(self):
        # TODO: (1, 2, 4, 5)
        self.par_set = {}

        # open the fcm file
        with open(self.values.par_set, 'rU') as f:
            for line in f:
                if '=' not in line or line.startswith('#'):
                    continue

                name, value = line.strip().split('=')
                name = name.strip()
                value = value.strip()
                self.par_set[name] = value

    def parse_mir(self):
        # TODO: (1, 2, 4, 5)
        self.mir = MiriadGainSolutions(self.values.mirsolutions,
                                       self.values.aips_c, self.pol,
                                       self.freqs)

    def get_ant_location(self, antno):
        # TODO: (1, 2, 4, 5)
        key = 'common.antenna.ant{}.location.itrf'.format(antno)
        value = self.par_set[key]
        location = map(float,
                       value.replace('[', '').replace(']', '').split(','))
        return location

    def get_fixed_delay_usec(self, antno):
        # TODO: (1, 2, 4, 5)
        key = 'common.antenna.ant{}.delay'.format(antno)
        value = self.par_set[key]
        delay_ns = float(value.replace('ns', ''))
        delay_us = delay_ns / 1e3

        return delay_us

    def get_geometric_delay_delayrate_us(self, ant):
        # TODO: (1, 2, 4, 5)
        fr1 = FringeRotParams(self, ant)
        fr2 = FringeRotParams(self, self.ref_ant)

        # TODO: There is a discrepancy here, below comment says fr1 is
        #       ref ant, but above suggests fr2 is?
        # fr1: reference antenna
        # Account for effects of Earth's rotation
        delay = fr1.delay_start - fr2.delay_start
        delayrate = fr1.delay_rate - fr2.delay_rate

        with open('delays/{}_ant_delays.dat'.format(ant.antno), 'w') as f:
            f.write('#field fr1({}) fr2({})\n'.format(ant, self.ref_ant))
            f.write('delay_start {} {}\n'.format(fr1.delay_start,
                                                 fr2.delay_start))
            f.write('delay {} {}\n'.format(fr1.delay, fr2.delay))
            f.write('delay_end {} {}\n'.format(fr1.delay_end, fr2.delay_end))
            f.write('delay_rate {} {}\n'.format(fr1.delay_rate,
                                                fr2.delay_rate))

        return delay, delayrate

    def calc_mjd(self):
        # TODO: (1, 2, 4, 5)
        i = float(self.curr_int_no)
        abs_delay_days = float(self.abs_delay)/86400./(self.fs*1e6)
        self.curr_mjd_start = self.mjd0 + self.int_time_days * (i + 0.0) + abs_delay_days
        self.curr_mjd_mid = self.mjd0 + self.int_time_days * (i + 0.5) + abs_delay_days
        self.curr_mjd_end = self.mjd0 + self.int_time_days * (i + 1.0) + abs_delay_days

    def get_calc_results(self, mjd):
        # TODO: (1, 2, 4, 5)
        res = self.calc_results.scans[0].eval_src0_poly(mjd)

        return res

    def get_fringe_rot_data(self):
        # TODO: (1, 2, 4, 5)
        self.fringe_rot_data_start = self.get_calc_results(self.curr_mjd_start)
        self.fringe_rot_data_mid = self.get_calc_results(self.curr_mjd_mid)
        self.fringe_rot_data_end = self.get_calc_results(self.curr_mjd_end)

    def do_tab(self, an=None):
        # TODO: (1, 2, 4, 5)
        # Tied-array beamforming
        
        n_samp = self.n_int
        n_chan = self.ncoarse_chan*self.n_fine_per_coarse
            
        sum_aligned = np.zeros((n_samp, n_chan, self.n_pol_in),
                               dtype=np.complex64)
        
        if an == None: # add all antennas
            print('## Summing up all ' + str(len(self.ants)) + ' antennas')
            if self.values.tab:
                print('coherent sum')
            else:
                print('incoherent sum')
            for iant, ant in enumerate(self.ants):
                if not self.running:
                    raise KeyboardInterrupt()
            
                temp = ant.do_f_tab(self, iant)
                
                if self.values.tab:
                    sum_aligned += temp
                else:
                    sum_aligned += np.power(abs(temp), 2)

            return sum_aligned
        else:
            print('## Operate on only antenna #: ' + str(an))
            ant = self.ants[an]
            iant = an
            temp = ant.do_f_tab(self, iant)
            return temp
    

def parse_delays(values):
    # TODO: (2, 4, 5)
    delay_file = values.calcfile.replace('.im', '.hwdelays')

    if os.path.exists(delay_file) is False:
        delay_file = values.hwfile

    delays = {}
    if delay_file is not None and os.path.exists(delay_file):
        with open(delay_file, 'rU') as dfile:
            for line in dfile:
                bits = line.split()
                if not line.startswith('#') and len(bits) == 2:
                    raw = -int(bits[1])
                    if raw % 8 != 0:    # if it is not a multiple of 8, round
                        new = int(8 * round(float(raw)/8))
                        print('hwdelay ', raw, ' rounded to ', new)
                        delays[bits[0].strip()] = new
                    else:
                        delays[bits[0].strip()] = raw

        logging.info('Loaded %s delays from %s', len(delays), delay_file)
    else:
        logging.info('No delays loaded. %s does not exist', delay_file)

    return delays


def load_sources(calc_file):
    """Loads the source in the given calc file

    :param calc_file: File containing source information
    :return: List containing one source dict.
             Source dict contains source name, RA, and DEC (in radians)
    """
    # TODO: (2, 5)
    calc_input = calc_file.replace('.im', '.calc')

    d = {}
    for line in open(calc_input, 'rU'):
        if len(line) == 0 or line.startswith('#'):
            continue

        bits = line.split(':')
        if len(bits) != 2:
            continue
        
        k, v = bits

        d[k.strip()] = v.strip()

    assert d['NUM SOURCES'] == '1'
    name = d['SOURCE 0 NAME']
    # ra/dec in radians
    ra = float(d['SOURCE 0 RA'])
    dec = float(d['SOURCE 0 DEC'])
    pos = SkyCoord(ra, dec, unit=('rad', 'rad'), frame='icrs')
    sources = [{'name': name, 'ra': pos.ra.deg, 'dec': pos.dec.deg}]

    return sources


def get_antennas(values):
    """Calculate delays and return antennas

    :param values: Command line input values
    :return: list of AntennaSource objects representing the list of antennas
    """

    delay_map = parse_delays(values)
    antennas = [AntennaSource(mux) for mux in
                vcraft.mux_by_antenna(values.files, delay_map)]

    print('NUMBER OF ANTENNAS', len(antennas))

    return antennas


def print_var(name, value):
    """Print a variable name and value in a consistent way

    :param name: Name of variable to be printed
    :param value: Value of variable to be printed
    """
    base_str = '{} : {}'
    print(base_str.format(name, value))


def _main():
    # TODO: (2, 5)
    parser = ArgumentParser(description='Script description',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='Be verbose', default=True)
    parser.add_argument('-o', '--outfile', help='Output fits/.npy file',
                        default='corr.fits')
    parser.add_argument('-c', '--channel', type=int, help='Channel to plot',
                        default=0)
    parser.add_argument('-n', '--fft-size', type=int,
                        help='Multiple of 64 channels to make channels - '
                             + 'default=1', default=1)
    parser.add_argument('-t', '--num-threads', type=int,
                        help='Number of threads to run with', default=1)
    parser.add_argument('--calc_file', help='Calc file for fringe rotation')
    parser.add_argument('-w', '--hwfile', help='Hw delay file')
    parser.add_argument('-p', '--par_set', help='Parset for delays')
    parser.add_argument('--show', help='Show plot', action='store_true',
                        default=False)
    parser.add_argument('-i', '--n_int',
                        help='Number of fine spectra to average', type=int,
                        default=128)
    parser.add_argument('-f', '--freq_scrunch',
                        help='Frequency average by this factor', default=1,
                        type=int)
    parser.add_argument('--rfidelay', type=int,
                        help='Delay in fine samples to add to second component'
                             + ' to make an RFI data set', default=0)
    parser.add_argument('--mirsolutions',
                        help='Root file name for miriad gain solutions')
    parser.add_argument('--aips_c', help='AIPS bandpass polynomial fit coeffs',
                        default=None)
    parser.add_argument('--an', type=int, help='Specific antenna',
                        default=None)
    parser.add_argument('--offset', type=int, help='FFT offset to add',
                        default=0)
    parser.add_argument(dest='files', nargs='+')
    values = parser.parse_args()

    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    sources = load_sources(values.calcfile)

    antennas = get_antennas(values)

    given_offset = values.offset
    corr = Correlator(antennas, sources, values, abs_delay=given_offset)

    t0 = time.time()
    try:
        print('PERFORMING TIED-ARRAY BEAMFORMING')

        temp = corr.do_tab(values.an)
        fn = values.outfile

        print('saving output to ' + fn)

        np.save(fn, temp)
    except Exception as e:
        print('ERROR OCCURRED')
        print(e)
    finally:
        print('craftcor_tab.py running time: ' + str(time.time() - t0))
        print('done')


if __name__ == '__main__':
    _main()
