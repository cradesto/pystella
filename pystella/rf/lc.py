import numpy as np
from pystella.rf import band

__author__ = 'bakl'


class LightCurve(object):
    def __init__(self, b, time, mags, errs=None, tshift=0):
        """Creates a Light Curve instance.  Required parameters:  b (band), time, mags."""
        if isinstance(b, str):  # convert band-name to band instance
            if band.band_is_exist(b):
                self._b = band.band_by_name(b)
            else:
                raise ValueError("No such band: {}".format(b))
        else:
            self._b = b
        self._t = np.array(time, copy=True)  # [days]
        self._m = np.array(mags, copy=True)  # magnitudes

        self._e = None
        if errs is not None:
            self._e = np.array(errs, copy=True)  # magnitudes
        self._tshift = tshift
        self._mshift = 0

    @property
    def Time(self):
        return self.T + self.tshift

    @property
    def T(self):
        return self._t

    @property
    def Length(self):
        return len(self._t)

    @property
    def Mag(self):
        return self.M + self.mshift

    @property
    def M(self):
        return self._m

    @property
    def IsErr(self):
        return self._e is not None

    @property
    def MagErr(self):
        if self.IsErr:
            return self._e
        else:
            return np.zeros(self.Length)

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

    @property
    def tmin(self):
        return np.min(self._t)

    @property
    def TimeMin(self):
        return np.min(self.Time)

    @property
    def TimeMax(self):
        return np.max(self.Time)

    @property
    def TimeLcMax(self):
        idx = np.argmin(self.Mag)
        return self.Time[idx]

    def copy(self, tlim=None):
        if tlim is not None:
            is_good = np.where((self.Time >= tlim[0]) & (self.Time <= tlim[1]))
            time = self.T[is_good]
            mags = self.M[is_good]
            errs = None
            if self.IsErr:
                errs = self.MagErr[is_good]
        else:
            time = self.T
            mags = self.M
            errs = None
            if self.IsErr:
                errs = self.MagErr

        lc = LightCurve(self.Band, time, mags, errs)
        lc.tshift = self.tshift
        lc.mshift = self.mshift


class SetLightCurve(object):
    """Set of the Light Curves"""
    def __init__(self, name=''):
        """Creates a Set of Light Curves."""
        self._name = name
        self._set = {}
        self._loop = 0

    @property
    def Name(self):
        return self._name

    @Name.setter
    def Name(self, v):
        self._name = v

    @property
    def Set(self):
        return self._set

    @property
    def Length(self):
        return len(self._set)

    @property
    def Bands(self):
        if len(self.Set) == 0:
            raise ValueError('There are no bands in SetLightCurve.')
        # for name, lc in self.Set.items():
        #     yield lc.Band
        res = (lc.Band for name, lc in self.Set.items())
        return res

    # @property
    # def Bands(self):
    #     if len(self.Set) == 0:
    #         raise ValueError('There are no bands in SetLightCurve.')
    #     res = [lc.Band for name, lc in self.Set.items()]
    #     return res
    #
    # @property
    # def BandNames(self):
    #     for b in self.Bands:
    #         yield b.Name
    @property
    def BandNames(self):
        res = [b.Name for b in self.Bands]
        return res

    @property
    def TimeDef(self):
        # b = next(self.BandNames)
        b = self.BandNames[0]
        return self.Set[b].Time

    @property
    def tmin(self):
        res = [lc.tmin for name, lc in self.Set.items()]
        return min(res)

    # for cycle
    def __getitem__(self, index):
        return self.Set[index]

    # def __iter__(self):
    #     self._loop = 0
    #     return self
    def __iter__(self):
        for k, v in self.Set.items():
            yield v
    #
    # def next(self):
    #     idx = self._loop
    #     if idx >= self.Length:
    #         raise StopIteration
    #     # b = self.BandNames[idx]
    #     self._loop += 1
    #     # return self.Set[b]
    #     return next(self.Set)

    def __len__(self):
        return len(self.Set)

    def add(self, lc):
        self._set[lc.Band.Name] = lc

    def rm(self, bname):
        return self._set.pop(bname, None)

    def get(self, bn):
        for n, lc in self.Set.items():
            if lc.Band.Name == bn:
                return lc
        return None

    def set_tshift(self, tshift):
        for n, lc in self.Set.items():
            lc.tshift = tshift

    def set_mshift(self, mshift):
        for n, lc in self.Set.items():
            lc.mshift = mshift
