#!/usr/bin/env python
"""
Tied-array beamforming vcraft files, based on "craftcor.py".

Copyright (C) CSIRO 2017
"""
#import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import vcraft
from calc11 import ResultsFile
#from corruvfits import CorrUvFitsFile
from astropy.coordinates import SkyCoord
import multiprocessing
import signal
import time
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from miriad import MiriadGainSolutions

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

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


def print_delay(xx):
    # TODO: (2, 3, 5)
    xxang = np.angle(xx)
    punwrap = np.unwrap(xxang)
    f = np.arange(len(punwrap)) - len(punwrap)/2
    gradient, phase = np.polyfit(f, punwrap, 1)

    '''
    pylab.figure(20)
    pylab.plot(f, xxang)
    pylab.plot(f, punwrap, 'x')
    pylab.plot(f, np.polyval((gradient, phase), f))
    pylab.show()
    '''

    delay = gradient / 2. / np.pi * len(punwrap)
    delayns = delay / F_OS * 1e3 * (54./len(punwrap))

    print('Unwrapped phase = {},\n'.format(phase)
          + 'rad = {} deg,\n'.format(np.degrees(phase))
          + 'gradient = {} rad per channel,\n'.format(gradient)
          + 'delay = {},\n'.format(delay)
          + 'samples = {} ns,\n'.format(delay / F_OS * 1e3)
          + 'nsamp = {}'.format(len(punwrap)))

    return delay, np.degrees(phase)


class PlotOut(object):
    # TODO: (1, 2, 3, 4, 5)
    def __init__(self, corr):
        # TODO: (1, 2, 3, 4, 5)
        self.stuff = []
        self.corr = corr
        self.delayout = open(corr.values.outfile.replace('.fits','')+'.delays.parset', 'w')

    def put_product(self, a1, a2, xxp):
        # TODO: (1, 2, 3, 4, 5)
        if self.corr.values.show:
            if a1 != a2:
                xx= xxp[:,0]
                self.stuff.append(xx)
                self.last_xx = (a1, a2, xx)
                if len(self.stuff) >0:
                    self.plot_stuff(a1, a2)

    def plot_stuff(self, a1, a2):
        # TODO: (1, 2, 3, 4, 5)
        '''(a1, a2, xx) = self.last_xx
        fig, (ax1, ax2, ax3, ax4, ax5,ax6) = pylab.subplots(6,1)
        xxang = np.angle(xx)
        print_delay(xx)
        ax1.plot(abs(xx))
        ax2.plot(np.degrees(xxang), 'o')
        nf = self.corr.nfine_out_per_coarse
        for i in xrange(self.corr.ncoarse_chan):
            ax2.axvline(nf*i, c='r')
            ax1.axvline(nf*i, c='r')
            print 'PLOTSTUFF', nf, i, nf*i, nf*(i+1), self.corr.ncoarse_chan, xx.shape
            print_delay(xx[nf*i:nf*(i+1)])

        lag = np.fft.fftshift(abs(np.fft.fft(xx)))
        lagidx = np.argmax(lag)
        lagoff = lagidx - len(xx)/2.0
        lagns= lagoff/self.corr.full_bw*1e3
        lagsn = lag[lagidx]/np.std(lag[0:lagidx-100])
        lagx = np.arange(len(lag)) - len(lag)/2.0
        ax3.plot(lagx/self.corr.full_bw*1e3, lag,label='lag')
        ax3.axvline(0, c='r')
        ax3.set_xlabel('Delay (ns)')
        angxx = np.angle(np.array(self.stuff))
        ax4.imshow(angxx, aspect='auto')
        phase_vt = np.degrees(angxx[:, 10:50].mean(axis=1))
        t = np.arange(len(phase_vt))
        rate, offset = np.polyfit(t, phase_vt, 1)
        fit_std = (np.polyval((rate, offset), t) - phase_vt).std()
        ax5.plot(phase_vt)
        ax6.imshow(abs(np.array(self.stuff)), aspect='auto')
        print 'Phase rate={} offset={} std={} deg lagoff={}samples = {}ns lag S/N={}'.format(rate, offset, fit_std, lagoff, lagns, lagsn)
        fig.suptitle('a1 {} a2 {}'.format(a1.ant_name, a2.ant_name))
        curr_delay = self.corr.get_fixed_delay_usec(a2.antno)*1e3
        if a1 != a2:
            self.delayout.write('#{}={} S/N={} delay={}ns\n'.format(a1.ant_name, a2.ant_name, lagsn, lagns))
            self.delayout.write('common.antenna.ant{}.delay={}ns\n'.format(a2.antno, curr_delay+lagns))
            self.delayout.flush()

        pylab.show()'''
        return # added by hyerin 

    def finish(self):
        # TODO: (1, 2, 3, 4, 5)
        stuff = np.array(self.stuff)
        #pylab.imshow(np.angle(stuff))
        #pylab.show()


class AntennaSource(object):
    # TODO: (1, 2, 3, 4, 5)
    def __init__(self, vfile):
        # TODO: (1, 2, 3, 4, 5)
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

    def do_f(self, corr):
        # TODO: (1, 2, 3, 4, 5)
        self.fringe_rot_params = FringeRotParams(corr, self)
        # calculate sample start
        framediff_samp = corr.refant.trigger_frame - self.trigger_frame
        framediff_us = framediff_samp / corr.fs
        (geom_delay_us, geom_delay_rate_us) = corr.get_geometric_delay_delayrate_us(self)
        self.all_geom_delays.append(geom_delay_us)
        self.all_mjds.append(corr.curr_mjd_mid)

        #logging.debug'ALL GEOM DELAYS', self.all_geom_delays, type(geom_delay_us)
        #print 'ALL MJDs', self.all_mjds
        # test with 0 delay rate
        #geom_delay_us = self.all_geom_delays[0]
        
        geom_delay_samp = geom_delay_us * corr.fs
        fixed_delay_us = corr.get_fixed_delay_usec(self.antno)
        fixed_delay_samp = fixed_delay_us*corr.fs
        total_delay_samp = framediff_samp
        whole_delay = int(np.round(total_delay_samp))
        total_delay_us = total_delay_samp / corr.fs
        whole_delay_us = whole_delay / corr.fs

        frac_delay_samp = (total_delay_samp - whole_delay)
        frac_delay_us = frac_delay_samp *corr.fs

        # get data
        nsamp = corr.nint*corr.nfft
        logging.debug('F %s sample delays: frame: %f geo %f fixed %f total %f whole: %d frac: %f nsamp: %d phase=%f deg',
                      self.ant_name, framediff_samp, geom_delay_samp, fixed_delay_samp,
                      total_delay_samp, whole_delay, frac_delay_samp, nsamp,
                      360. * frac_delay_us * corr.f0)
        logging.debug('F %s us delays: frame: %f geo %f fixed %f total %f whole: %d frac: %f nsamp: %d phase=%f deg rate=%e',
                      self.ant_name, framediff_us, geom_delay_us, fixed_delay_us,
                      total_delay_us, whole_delay_us, frac_delay_us, nsamp,
                      360. * frac_delay_us * corr.f0, geom_delay_rate_us)

        sampoff = corr.curr_samp*corr.nfft + whole_delay
        rawd = self.vfile.read(sampoff, nsamp)
        assert rawd.shape == (nsamp, corr.ncoarse_chan), 'Unexpected shape from vfile: {} expected ({},{})'.format(rawd.shape, nsamp, corr.ncoarse_chan)
        self.data = np.zeros((corr.nint, corr.nfine_chan, corr.npol_in), dtype=np.complex64)
        d1 = self.data
        nfine = corr.nfft - 2*corr.nguard_chan

        for c in xrange(corr.ncoarse_chan):
            cfreq = corr.freqs[c]
            freqs = (np.arange(nfine, dtype=np.float) - float(nfine)/2.0)*corr.fine_chanbw
            if corr.sideband == -1:
                freqs = -freqs
                
            x1 = rawd[:, c].reshape(-1, corr.nfft)
            xf1 = np.fft.fft(x1, axis=1)
            xf1 = np.fft.fftshift(xf1, axes=1)
            xfguard = xf1[:, corr.nguard_chan:corr.nguard_chan+nfine:] # scale because oterhwise it overflows
            delta_t = -fixed_delay_us + geom_delay_us
            phases = delta_t * (freqs + cfreq)
            logging.debug('PHASOR %s[%s] chan=%s freq=%sfixed=%f us geom=%f us delta_t %s us coff*fixed = %f deg coff*geom = %f deg',
                          self.ant_name, self.ia, c, cfreq, fixed_delay_us, geom_delay_us, delta_t, cfreq * fixed_delay_us * 360., cfreq * geom_delay_us * 360.)


            # If you plot the phases you're about to correct, after adding a artificial
            # 1 sample delay ad tryig to get rid of it with a phase ramp, it becaomes
            # blatetly clear what you should do
            phasor = np.exp(np.pi*2j*phases, dtype=np.complex64)
            '''
            pylab.figure(10)
            pylab.plot(np.angle(phasor[0, :]))
            pylab.plot(np.angle(phasor[-1:, :]))
            '''


            xfguard *= phasor
            # slice out only useful channels
            fcstart = c*nfine
            fcend = (c+1)*nfine
            self.data[:, fcstart:fcend, 0] = xfguard

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
        framediff_samp = corr.refant.trigger_frame - self.trigger_frame
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

        data_out = np.zeros((corr.nint, corr.nfine_chan, corr.npol_in),
                            dtype=np.complex64)
        n_fine = corr.nfft - 2*corr.nguard_chan

        n_samp = corr.nint * corr.nfft

        # time-dependent geometric delays
        # np.linspace(0, 1, n_samp) == time in units of integrations
        geom_delays_us = (geom_delay_us
                         + geom_delay_rate_us
                         * np.linspace(0, 1, n_samp)
                         - fixed_delay_us)

        np.save('delays/geom_delays_us_{}'.format(i_ant), geom_delays_us)

        print('')
        base_str = '{} : {}'
        print(base_str.format('framediff_samp', framediff_samp))
        print(base_str.format('framediff_us', framediff_us))
        print(base_str.format('geom_delay_us', geom_delay_us))
        print(base_str.format('geom_delay_rate_us', geom_delay_rate_us))
        print(base_str.format('geom_delay_samp', geom_delay_samp))
        print(base_str.format('np.mean(geom_delays_us)',
                              np.mean(geom_delays_us)))
        print(base_str.format('fixed_delay_us', fixed_delay_us))
        print(base_str.format('fixed_delay_samp', fixed_delay_samp))
        print(base_str.format('total_delay_samp', total_delay_samp))
        print(base_str.format('whole_delay', whole_delay))
        print(base_str.format('total_delay_us', total_delay_us))
        print(base_str.format('whole_delay_us', whole_delay_us))
        print(base_str.format('frac_delay_samp', frac_delay_samp))
        print(base_str.format('frac_delay_us', frac_delay_us))
        print('')

        print("antenna #: ", i_ant, self.ant_name)
        sample_offset = whole_delay + corr.abs_delay
        frameid = self.vfile.start_frameid + sample_offset
        print('FRAMEID: '
              + str(frameid)
              + ', remainder from 32: '
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
            delta_freq = (corr.fine_chanbw
                          * (np.arange(n_fine, dtype=np.float)
                          - float(n_fine)/2.0))

            # TODO: what is sideband?
            if corr.sideband == -1:
                delta_freq = -delta_freq

            '''
            raw_data's shape: (n_samp, corr.ncoarse_chan)
            n_samp = input.i * (64 * input.n)
            '''
            x1 = raw_data[:, chan].reshape(-1, corr.nfft)

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
            
    def get_data(self, chan, pol):
        # TODO: (1, 2, 3, 4, 5)
        return self.data[:, chan, pol]


class FringeRotParams(object):
    # TODO: (1, 2, 3, 4, 5)
    cols = ('U (m)', 'V (m)', 'W (m)', 'DELAY (us)')

    def __init__(self, corr, ant):
        # TODO: (1, 2, 3, 4, 5)
        mid_data = corr.frdata_mid[ant.ant_name]
        self.u,self.v,self.w,self.delay = map(float, [mid_data[c] for c in FringeRotParams.cols])
        self.delay_start = float(corr.frdata_start[ant.ant_name]['DELAY (us)'])
        self.delay_end = float(corr.frdata_end[ant.ant_name]['DELAY (us)'])
        self.delay_rate = (self.delay_end - self.delay_start)/float(corr.nint)
        self.ant = ant
        self.corr = corr

    def __str__(self):
        # TODO: (1, 2, 3, 4, 5)
        s = 'FR {} uvw=({},{},{}) m = {} us'.format(self.ant.ant_name, self.u, self.v, self.w, self.delay)
        return s

    __repr__ = __str__


class Correlator(object):
    # TODO: (1, 2, 3, 4, 5)
    def __init__(self, ants, sources, values, abs_delay=0):
        # TODO: (1, 2, 3, 4, 5)
        self.running = True
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)
        self.ants = ants
        self.values = values
        self.pool = None
        if self.values.num_threads > 1:
            self.pool = multiprocessing.Pool(processes=values.num_threads)

        self.parse_parset()

        for ia, a in enumerate(self.ants):
            a.ia = ia
            a.antpos = self.get_ant_location(a.antno)

        refantname = self.parset['cp.ingest.tasks.FringeRotationTask.params.refant'].lower()
        self.abs_delay = abs_delay
        #self.refant = filter(lambda a:a.ant_name == refantname, ants)[0]
        self.refant = ants[0]
        self.calcresults = ResultsFile(values.calcfile)
        self.dutc = 0
        self.mjd0 = self.refant.mjd_start + self.dutc / 86400.0
        self.frame0 = self.refant.trigger_frame
        self.nint = values.nint
        self.nfft = 64*values.fft_size
        self.nguard_chan = 5*values.fft_size
        self.oversamp = 32./27.
        self.fs = self.oversamp # samples per microsecnd
        self.ncoarse_chan = len(self.refant.vfile.freqs)
        self.sideband = -1
        self.coarse_chanbw = 1.0
        self.nfine_per_coarse = self.nfft - 2*self.nguard_chan
        self.nfine_chan = self.ncoarse_chan*self.nfine_per_coarse
        self.fine_chanbw = self.coarse_chanbw / float(self.nfine_per_coarse)
        self.full_bw = self.fine_chanbw * self.nfine_chan
        self.fscrunch = values.fscrunch
        assert self.fscrunch >= 1
        assert self.nfine_per_coarse % self.fscrunch == 0, 'Fsrunch must yield an integer number of fine channels per coarse channel'
        self.nfine_out_per_coarse = self.nfine_per_coarse / self.fscrunch
        self.nfine_out_chan = self.nfine_out_per_coarse*self.ncoarse_chan
        self.out_chanbw = self.coarse_chanbw / float(self.nfine_out_per_coarse)
        self.npol_in = 1
        self.npol_out = 1
        #self.f0 = self.ants[0].vfile.freqs.mean() # centre frequency for fringe rotation
        self.f0 = self.ants[0].vfile.freqs[0]
        self.freqs = self.ants[0].vfile.freqs
        self.fmid = self.freqs.mean()
        self.inttime_secs = float(self.nint*self.nfft)/(self.fs*1e6)
        self.inttime_days = self.inttime_secs/86400.
        self.curr_intno = 0
        self.curr_samp = self.curr_intno*self.nint + 1000
        self.prodout = PlotOut(self)
        self.calcmjd()
        self.get_fr_data()
        self.pol = self.ants[0].pol
        self.parse_mir()
        #self.fileout = CorrUvFitsFile(values.outfile, self.fmid, self.sideband*self.out_chanbw, \
                                      #self.nfine_out_chan, self.npol_out, self.mjd0, sources, ants, self.sideband)

        logging.debug('F0 %f FINE CHANNEL %f kHz num=%d freqs=%s', self.f0, self.fine_chanbw*1e3, self.nfine_chan, self.freqs)

    def exit_gracefully(self, signum, frame):
        # TODO: (1, 2, 3, 4, 5)
        self.running = False

    def parse_parset(self):
        # TODO: (1, 2, 3, 4, 5)
        self.parset = {}

        # open the fcm file
        with open(self.values.parset, 'rU') as f:
            for line in f:
                if '=' not in line or line.startswith('#'):
                    continue

                name, value = line.strip().split('=')
                name = name.strip()
                value = value.strip()
                self.parset[name] = value

    def parse_mir(self):
        # TODO: (1, 2, 3, 4, 5)
        #self.mir = None
        #if self.values.mirsolutions is not None or self.values.aips_c is not None:
        self.mir = MiriadGainSolutions(self.values.mirsolutions,self.values.aips_c, self.pol, self.freqs)

    def get_ant_location(self, antno):
        # TODO: (1, 2, 3, 4, 5)
        key = 'common.antenna.ant{}.location.itrf'.format(antno)
        value = self.parset[key]
        location = map(float, value.replace('[','').replace(']','').split(','))
        return location

    def get_fixed_delay_usec(self, antno):
        # TODO: (1, 2, 3, 4, 5)
        key = 'common.antenna.ant{}.delay'.format(antno)
        value = self.parset[key]
        delayns =  float(value.replace('ns',''))
        delayus = delayns/1e3

        return delayus

    def get_uvw(self, ant1, ant2):
        # TODO: (1, 2, 3, 4, 5)
        fr1 = FringeRotParams(self, ant1)
        fr2 = FringeRotParams(self, ant2)
        uvw = np.array([fr1.u - fr2.u, fr1.v - fr2.v, fr1.w - fr2.w]) / C_LIGHT

        return uvw

    def get_geometric_delay_delayrate_us(self, ant):
        # TODO: (1, 2, 3, 4, 5)
        fr1 = FringeRotParams(self, ant)
        fr2 = FringeRotParams(self, self.refant)

        # TODO: There is a discrepancy here, below comment says fr1 is ref ant, but above suggests fr2 is?
        # fr1: reference antenna
        # Account for effects of Earth's rotation
        #delay = fr1.delay - fr2.delay
        #delayrate = fr1.delay_rate - fr2.delay_rate
        delay = fr1.delay_start - fr2.delay_start
        delayrate = fr1.delay_rate - fr2.delay_rate

        with open('delays/{}_ant_delays.dat'.format(ant.antno), 'w') as f:
            f.write('#field fr1({}) fr2({})\n'.format(ant, self.refant))
            f.write('delay_start {} {}\n'.format(fr1.delay_start, fr2.delay_start))
            f.write('delay {} {}\n'.format(fr1.delay, fr2.delay))
            f.write('delay_end {} {}\n'.format(fr1.delay_end, fr2.delay_end))
            f.write('delay_rate {} {}\n'.format(fr1.delay_rate, fr2.delay_rate))

        return (delay, delayrate)

    def calcmjd(self):
        # TODO: (1, 2, 3, 4, 5)
        i = float(self.curr_intno)
        abs_delay_days = float(self.abs_delay)/86400./(self.fs*1e6)
        self.curr_mjd_start = self.mjd0 + self.inttime_days*(i + 0.0) + abs_delay_days
        self.curr_mjd_mid = self.mjd0 + self.inttime_days*(i + 0.5) + abs_delay_days
        self.curr_mjd_end = self.mjd0 + self.inttime_days*(i + 1.0) + abs_delay_days

    def next_integration(self):
        # TODO: (1, 2, 3, 4, 5)
        #self.curr_intno +=
        self.curr_intno += 1
        self.curr_samp += self.nint
        self.calcmjd()
        self.get_fr_data()

    def get_calc_results(self, mjd):
        # TODO: (1, 2, 3, 4, 5)
        #res = self.calcresults.scans[0].eval_src0_poly_delta(mjd, self.refant.ant_name.lower())
        res = self.calcresults.scans[0].eval_src0_poly(mjd)  # Calcresults is a ResultsFile

        return res

    def get_fr_data(self):
        # TODO: (1, 2, 3, 4, 5)
        self.frdata_start = self.get_calc_results(self.curr_mjd_start)
        self.frdata_mid = self.get_calc_results(self.curr_mjd_mid)
        self.frdata_end = self.get_calc_results(self.curr_mjd_end)


    def do_f(self):
        # TODO: (1, 2, 3, 4, 5)
        for iant, ant in enumerate(self.ants):
            if not self.running:
                raise KeyboardInterrupt()

            ant.do_f(self)

    def do_x_corr(self, a1, a2):
        # TODO: (1, 2, 3, 4, 5)
        npolout = self.npol_out
        xx = np.empty([self.nfine_chan, npolout], dtype=np.complex64)
        #np.seterr(all='raise')
        rfidelay = self.values.rfidelay
        for p1 in xrange(self.npol_in):
            for p2 in xrange(self.npol_in):
                d1 = a1.data[:, :, p1]
                d2 = a2.data[:, :, p2]
                pout = p2 + p1*self.npol_in
                ntimesi = d1.shape[0] - rfidelay
                ntimes = float(ntimesi)

                try:
                    for c in xrange(self.nfine_chan):
                        #xx[:,pout] = (d1 * np.conj(d2)).mean(axis=0)
                        # vdot conjugates the first argument
                        # this is equivalent to (d1 * conj(d2)).mean(axis=0)
                        if self.sideband == -1:
                            xx[c, pout] = np.vdot(d2[:ntimesi, c], d1[rfidelay:, c])/ntimes
                        else:# take complex conjugate if inverted
                            xx[c, pout] = np.vdot(d1[:ntimesi, c], d2[rfidelay:, c])/ntimes
                            
                except Exception, e:
                    print 'Error', e
                    import ipdb
                    ipdb.set_trace()

        if self.fscrunch > 1:
            chans = np.arange(self.nfine_chan, step=self.fscrunch)
            xx = np.add.reduceat(xx, chans, axis=0)/float(self.fscrunch)

        assert xx.shape == (self.nfine_out_chan, self.npol_out)
            
        self.put_product(a1, a2, xx)

    def do_tab(self, an=None):
        # TODO: (1, 2, 3, 4, 5)
        # Tied-array beamforming
        
        nsamp = self.nint
        nchan = self.ncoarse_chan*self.nfine_per_coarse
            
        sum_aligned = np.zeros((nsamp, nchan, self.npol_in), dtype=np.complex64)
        
        if an == None: # add all antennas
            print('## Summing up all '+str(len(self.ants))+' antennas')
            if self.values.tab:
                print('coherent sum')
            else:
                print('incoherent sum')
            for iant, ant in enumerate(self.ants):
                if not self.running:
                    raise KeyboardInterrupt()
            
                temp = ant.do_f_tab(self,iant)
                
                if self.values.tab:
                    sum_aligned += temp
                else:
                    sum_aligned += np.power(abs(temp),2)

            return sum_aligned
        else:
            print('## Operate on only antenna #: '+str(an))
            ant = self.ants[an]
            iant = an
            temp = ant.do_f_tab(self,iant)
            return temp
    
    def put_product(self, a1, a2, xx):
        # TODO: (1, 2, 3, 4, 5)
        self.prodout.put_product(a1, a2, xx)
        uvw = self.get_uvw(a1, a2)
        #self.fileout.put_data(uvw, self.curr_mjd_mid, a1.ia, a2.ia,
            #self.inttime_secs, xx)


def parse_delays(values):
    # TODO: (1, 2, 3, 4, 5)
    delayfile = values.calcfile.replace('.im','.hwdelays')
    if os.path.exists(delayfile)==False:
        delayfile = values.hwfile
    #print(delayfile)
    delays = {}
    if delayfile is not None and os.path.exists(delayfile):
        with open(delayfile, 'rU') as dfile:
            for line in dfile:
                bits = line.split()
                if not line.startswith('#') and len(bits) == 2:
                    raw = -int(bits[1])
                    if raw % 8 !=0: # if it is not a multiple of 8, round
                        new = int(8 * round(float(raw)/8))
                        print('hwdelay ',raw,' rounded to ',new)
                        delays[bits[0].strip()] = new
                    else:
                        delays[bits[0].strip()] = raw

        logging.info('Loaded %s delays from %s', len(delays), delayfile)
    else:
        logging.info('No delays loaded. %s does not exist', delayfile)
    

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
    sources = [{'name' : name, 'ra' : pos.ra.deg, 'dec' : pos.dec.deg}]

    return sources


def get_antennas(values):
    """Calculate delays and return antennas

    :param values: Command line input values
    :return: list of AntennaSource objects representing the list of antennas
    """

    delay_map = parse_delays(values)
    antennas = [AntennaSource(mux) for mux in
                vcraft.mux_by_antenna(values.files, delay_map)]

    print('NUMBER OF ANTENNAS',len(antennas))

    return antennas


def _main():
    # TODO: (2, 5)
    parser = ArgumentParser(description='Script description',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='Be verbose', default=True)
    parser.add_argument('-o','--outfile', help='Output fits/.npy file',
                        default='corr.fits')
    parser.add_argument('-c','--channel', type=int, help='Channel to plot',
                        default=0)
    parser.add_argument('-n','--fft-size', type=int,
                        help='Multiple of 64 channels to make channels - '
                             + 'default=1', default=1)
    parser.add_argument('-t','--num-threads', type=int,
                        help='Number of threads to run with', default=1)
    parser.add_argument('--calc_file', help='Calc file for fringe rotation')
    parser.add_argument('-w','--hwfile', help='Hw delay file')
    parser.add_argument('-p','--parset', help='Parset for delays')
    parser.add_argument('--show', help='Show plot', action='store_true',
                        default=False)
    parser.add_argument('-i','--nint',
                        help='Number of fine spectra to average', type=int,
                        default=128)
    parser.add_argument('-f','--fscrunch',
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
