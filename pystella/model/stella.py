import os
import numpy as np
import math
from pystella.model.sn_res import StellaRes
from pystella.model.sn_swd import StellaShockWaveDetail
from pystella.model.sn_tt import StellaTt
from pystella.rf.spectrum import SeriesSpectrum, Spectrum
from pystella.model import sn_eve

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
                print("Exist %s-file: %s" % (e, fname))
            else:
                print("No %s-file: %s" % (e, fname))

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
        fname = os.path.join(self.path, self.name + '.res')
        return os.path.isfile(fname)

    @property
    def is_tt_data(self):
        fname = os.path.join(self.path, self.name + '.tt')
        return os.path.isfile(fname)

    @property
    def series_spec(self):
        return self._serial_spec

    def set_series_spec(self, ss):
        self._serial_spec = ss

    def get_eve(self, name, path=None):
        if path is None:
            path = self.path
        eve = sn_eve.load_rho(os.path.join(path, name+'.rho'))
        return eve

    def get_res(self):
        return StellaRes(self.name, self.path)

    def get_tt(self):
        return StellaTt(self.name, self.path)

    def get_swd(self):
        swd = StellaShockWaveDetail(self.name, self.path)
        return swd

    def read_series_spectrum(self, t_diff=1.005, t_beg=0.0, t_end=None, is_nfrus=True):
        if t_end is None:
            t_end = float('inf')
        # read first line with frequencies
        fname = os.path.join(self.path, self.name + '.ph')
        f = open(fname, 'r')
        try:
            header1 = f.readline()
        finally:
            f.close()

        freqs = [float(x) for x in header1.split()]
        freqs = np.array(freqs).reshape(-1)
        freqs = [10.**nu for nu in freqs]
        # freqs = np.exp(math.log(10) * freqs)

        data = np.loadtxt(fname, comments='!', skiprows=1)

        times = np.array(data[:, 0])
        is_times = np.zeros(len(times), dtype=bool)
        k = 1
        for i in range(len(times)):
            if times[i] < t_beg:
                k = i
                continue
            if times[i] > t_end:
                break
            if np.abs(times[i] / times[k]) > t_diff:
                is_times[k] = True
                k = i
        is_times[0] = True  # times[0] > 0.
        is_times[-1] = True

        series = SeriesSpectrum(self.name)
        for i in range(len(times)):
            if is_times[i]:
                t = times[i]
                if is_nfrus:
                    nfrus = int(data[i, 1])  # exact number of used (saved) freqs
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

    # old
    def read_tt_data(self):
        return self.get_tt().read()
        # num_line_header = 87
        # fname = os.path.join(self.path, self.name + '.tt')
        # header = ''
        # i = 0
        # with open(fname, "r") as f:
        #     for line in f:
        #         i += 1
        #         if i < num_line_header:
        #             continue
        #         if line.lstrip().startswith("time"):
        #             header = line
        #             num_line_header = i
        #             break
        # # time Tbb rbb Teff Rlast_sc R(tau2/3) Mbol MU MB MV MI MR Mbolavg  gdepos
        # # time Tbb rbb Teff Rlast_sc R(tau2/3) Mbol MU MB MV MI MR   Mbolavg  gdepos
        # if header != '':
        #     names = map(str.strip, header.split())
        #     names = [w.replace('R(tau2/3)', 'Rph') for w in names]
        #     dtype = np.dtype({'names': names, 'formats': [np.float64] * len(names)})
        #     block = np.loadtxt(fname, skiprows=num_line_header + 1, dtype=dtype)
        #     return block
        # else:
        #     return None


def show_info(path, cond=lambda i: True):
    """Print information list about models in the path
    :param path: working directory
    :param cond: condition function, like lambda i: 30 < i.M < 1000 and i.R > 100
    :return: None
   """
    from os import listdir
    from os.path import isfile, join

    files = [f for f in listdir(path) if isfile(join(path, f)) and f.endswith('.tt')]
    for f in files:
        # print 'Read: %s' % f
        name, ext = os.path.splitext(f)
        stella = Stella(name, path=path)
        info = stella.get_tt().Info
        #         print(info.Data)
        if cond(info):
            # #         if 30 < info.R < 40:
            info.show()
