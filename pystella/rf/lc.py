import csv
import numpy as np
# import scipy.optimize as opt
# from scipy import interpolate, integrate

from pystella.rf import band

__author__ = 'bakl'


class LightCurve:
    def __init__(self, b, time, mags, tshift=0):
        """Creates a Light Curve instance.  Required parameters:  b (band), time, mags."""
        if isinstance(b, basestring):  # convert band-name to band instance
            self._b = band.band_by_name(b)
        else:
            self._b = b
        self._t = np.array(time, copy=True)  # [days]
        self._m = np.array(mags, copy=True)  # magnitudes
        self._tshift = tshift
        self._mshift = 0

    @property
    def Time(self):
        return self._t + self._tshift

    @property
    def Length(self):
        return len(self._t)

    @property
    def Mag(self):
        return self._m + self._mshift

    @property
    def Band(self):
        return self._b

    @property
    def tshift(self):
        return self._tshift

    @tshift.setter
    def tshift(self, shift):
        self._tshift = shift

    @property
    def mshift(self):
        return self._mshift

    @mshift.setter
    def mshift(self, shift):
        self._mshift = shift


class SetLightCurve:
    def __init__(self, name=''):
        """Creates a Set of Light Curves."""
        self._name = name
        self._set = {}

    @property
    def Name(self):
        return self._name

    @property
    def Set(self):
        return self._set

    @property
    def Length(self):
        return len(self._set)

    def __getitem__(self, index):
        return self.Set[index]

    @property
    def Bands(self):
        if len(self.Set) == 0:
            raise ValueError('There are no bands in SetLightCurve.')
        res = [lc.Band for name, lc in self.Set.items()]
        return res

    @property
    def BandNames(self):
        res = [b.Name for b in self.Bands]
        return res

    @property
    def TimeDef(self):
        b = self.BandNames[0]
        return self.Set[b].Time

    def add(self, lc):
        self._set[lc.Band.Name] = lc

    def save(self, fname):
        with open(fname, 'wb') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['{:^8s}'.format(x) for x in ['time'] + self.BandNames])
            for i, (row) in enumerate(zip(self.TimeDef, *[self.Set[b] for b in self.BandNames])):
                # row = row[-1:] + row[:-1]  # make time first column
                writer.writerow(['{:8.3f}'.format(x) for x in row])
                # writer.writerow(['{:3.4e}'.format(x) for x in row])
