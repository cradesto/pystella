import os
import numpy as np
import math
from pystella.model.sn_res import StellaRes
from pystella.rf.spectrum import SeriesSpectrum, Spectrum

__author__ = 'bakl'


class Stella:
    def __init__(self, name, path='./', info=False):
        """Creates a Stella model instance.  Required parameters:  name."""
        self.name = name
        self.path = path  # path to model files
        if info:
            self.show_info()
        self._serial_spec = None

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

    @property
    def is_ph_data(self):
        fname = os.path.join(self.path, self.name + '.ph')
        return os.path.isfile(fname)

    @property
    def is_res_data(self):
        ext = ['res']
        return any(map(os.path.isfile, [os.path.join(self.path, self.name + '.' + e) for e in ext]))

    @property
    def is_tt_data(self):
        ext = ['tt']
        return any(map(os.path.isfile, [os.path.join(self.path, self.name + '.' + e) for e in ext]))

    @property
    def series_spec(self):
        return self._serial_spec

    def set_series_spec(self, ss):
        self._serial_spec = ss

    def get_res(self):
        return StellaRes(self.name, self.path)

    def read_series_spectrum(self, t_diff=1.05, t_beg=0.0, is_nfrus=True):

        # read first line with frequencies
        fname = os.path.join(self.path, self.name + '.ph')
        f = open(fname, 'r')
        try:
            header1 = f.readline()
        finally:
            f.close()

        freqs = [map(float, header1.split())]
        freqs = np.array(freqs).reshape(-1)
        freqs = np.exp(math.log(10) * freqs)

        data = np.loadtxt(fname, comments='!', skiprows=1)

        times = np.array(data[:, 0])
        is_times = np.zeros(len(times), dtype=bool)
        k = 1
        for i in range(len(times)):
            if times[i] < t_beg:
                k = i
                continue
            if np.abs(times[i] / times[k]) > t_diff:
                is_times[k] = True
                k = i
        is_times[0] = True  # times[0] > 0.
        is_times[-1] = True
        times_thin = times[is_times]

        series = SeriesSpectrum(self.name)
        for i in range(len(times)):
            if is_times[i]:
                t = times[i]
                if is_nfrus:
                    nfrus = data[i, 1]  # exact number of used (saved) freqs
                    freqs = freqs[:nfrus]
                    fl = np.array(data[i, 3:nfrus+3])
                else:
                    fl = np.array(data[i, 3:])
                fl[fl < 0] = 0.
                fl = np.exp(math.log(10) * fl)
                s = Spectrum(self.name, freq=freqs, flux=fl, is_sort_wl=True)
                series.add(t, s)

        series.set_freq(freqs)
        # series.set_times(times_thin)

        self.set_series_spec(series)
        return series

    def read_tt_data(self):
        num_line_header = 87
        fname = os.path.join(self.path, self.name + '.tt')
        header = ''
        i = 0
        with open(fname, "r") as f:
            for line in f:
                i += 1
                if i < num_line_header:
                    continue
                if line.lstrip().startswith("time"):
                    header = line
                    num_line_header = i
                    break
        # time Tbb rbb Teff Rlast_sc R(tau2/3) Mbol MU MB MV MI MR Mbolavg  gdepos
        # time Tbb rbb Teff Rlast_sc R(tau2/3) Mbol MU MB MV MI MR   Mbolavg  gdepos
        if header != '':
            names = map(str.strip, header.split())
            names = [w.replace('R(tau2/3)', 'Rph') for w in names]
            dtype = np.dtype({'names': names, 'formats': [np.float64] * len(names)})
            block = np.loadtxt(fname, skiprows=num_line_header + 1, dtype=dtype)
            return block
        else:
            return None
