import os
from rf.spectrum import SeriesSpectrum, Spectrum
import numpy as np

__author__ = 'bakl'


class Stella:
    def __init__(self, name, path='./', info=False):
        """Creates a Stella model instance.  Required parameters:  name."""
        self.name = name
        self.path = path  # path to model files
        if info:
            self.show_info()

    def show_info(self):
        ext = ('tt', 'ph', 'res')
        for e in ext:
            fname = os.path.join(self.path, self.name+'.'+e)
            if os.path.isfile(fname):
                print "Exist %s-file: %s" % (e, fname)

    def read_serial_spectrum(self):
        fname = os.path.join(self.path, self.name+'.ph')
        serial = SeriesSpectrum(self.name)

        # TODO read first line with frequencies
        f = open(fname, 'r')
        # read NX,NY,NZ
        header1 = f.readline()
        f.close()

        freqs = [map(float, header1.split())]
        freqs = np.log10(freqs)
        serial.set_freq(freqs)

        # TODO read time timeph, nfrus, dum, ttt(Nfreq)
        data = np.loadtxt(fname, comments='!', skiprows=1)
        times = data[:, 1]
        serial.set_times(times)

        sdata = []
        for i in range(len(times)):
            s = Spectrum(self.name, wl=serial.wl, flux=data[i, 3:])
            sdata.append(s)
        serial.set_data(sdata)
        return serial

