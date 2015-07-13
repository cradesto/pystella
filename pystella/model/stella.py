import os
import numpy as np
import math
from pystella.rf.spectrum import SeriesSpectrum, Spectrum

__author__ = 'bakl'


class Stella:
    def __init__(self, name, path='./', info=False):
        """Creates a Stella model instance.  Required parameters:  name."""
        self.name = name
        self.path = path  # path to model files
        if info:
            self.show_info()

    def __str__(self):
        return "%s, path: %s" % (self.name, self.path)

    def __repr__(self):
        return "%s, path: %s" % (self.name, self.path)
        # return "%s" % self.name

    def show_info(self):
        ext = ('tt', 'ph', 'res', 'swd')
        for e in ext:
            fname = os.path.join(self.path, self.name + '.' + e)
            if os.path.isfile(fname):
                print "Exist %s-file: %s" % (e, fname)
            else:
                print "No %s-file: %s" % (e, fname)

    @property
    def is_any_data(self):
        ext = ('tt', 'ph', 'res', 'swd')
        return any(map(os.path.isfile, [os.path.join(self.path, self.name + '.' + e) for e in ext]))

    def read_serial_spectrum(self, t_diff=1.05):
        fname = os.path.join(self.path, self.name + '.ph')
        serial = SeriesSpectrum(self.name)

        # read first line with frequencies
        f = open(fname, 'r')
        # read NX,NY,NZ
        header1 = f.readline()
        f.close()

        freqs = [map(float, header1.split())]
        freqs = np.array(freqs).reshape(-1)
        freqs = np.exp(math.log(10) * freqs)
        serial.set_freq(freqs)

        # read time timeph, nfrus, dum, ttt(Nfreq)
        data = np.loadtxt(fname, comments='!', skiprows=1)

        times = np.array(data[:, 0])
        is_times = np.zeros(len(times), dtype=bool)
        k = 1
        for i in range(len(times)):
            if np.abs(times[i] / times[k]) > t_diff:
                is_times[k] = True
                k = i
        is_times[0] = True
        is_times[-1] = True
        times_thin = times[is_times]

        sdata = []
        for i in range(len(times)):
            if is_times[i]:
                fl = np.array(data[i, 3:])
                fl[fl < 0] = 0.
                fl = np.exp(math.log(10) * fl)
                s = Spectrum(self.name, freq=serial.freq, flux=fl, is_sort_wl=True)
                sdata.append(s)

        serial.set_data(sdata)
        serial.set_times(times_thin)
        return serial
