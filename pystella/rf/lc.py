import numpy as np

from pystella.rf import band
from pystella.rf.ts import TimeSeries

__author__ = 'bakl'


class LightCurve(TimeSeries):
    def __init__(self, b, time, mags, errs=None, tshift=0., mshift=0.):
        """Creates a Light Curve instance.  Required parameters:  b (band), time, mags."""
        if isinstance(b, str):  # convert band-name to band instance
            if band.band_is_exist(b):
                self._b = band.band_by_name(b)
            else:
                raise ValueError("No such band: {}".format(b))
        else:
            self._b = b

        super().__init__(time, mags, errs, name=self._b.Name, tshift=tshift)
        self._mshift = mshift

    @property
    def Mag(self):
        return self.V + self.mshift

    @property
    def M(self):
        return self.V

    @property
    def MagErr(self):
        return self.Err

    @property
    def Band(self):
        return self._b

    @property
    def mshift(self):
        return self._mshift

    @mshift.setter
    def mshift(self, shift):
        self._mshift = shift

    @property
    def TimeLcMax(self):
        idx = np.argmin(self.Mag)
        return self.Time[idx]

    def copy(self, tlim=None):
        errs = None
        if tlim is not None:
            is_good = np.where((self.Time >= tlim[0]) & (self.Time <= tlim[1]))
            time = self.T[is_good]
            mags = self.V[is_good]
            if self.IsErr:
                errs = self.MagErr[is_good]
        else:
            time = self.T
            mags = self.V
            if self.IsErr:
                errs = self.MagErr

        lc = LightCurve(self.Band, time, mags, errs)
        lc.tshift = self.tshift
        lc.mshift = self.mshift

    def shifted(self):
        errs = None
        if self.IsErr:
            errs = self.MagErr

        return  LightCurve(self.Band, self.Time, self.Mag, errs)

    @classmethod
    def Merge(cls, lc1, lc2):
        if lc1.Band.Name != lc2.Band.Name:
            raise ValueError("Merging is possible only for the same filters: {} VS {}".
                             format(lc1.Band.Name, lc2.Band.Name))
        bname = lc1.Band.Name
        t = np.concatenate((lc1.Time, lc2.Time))
        m = np.concatenate((lc1.Mag, lc2.Mag))

        sorti = np.argsort(t)
        time = t[sorti]
        mags = m[sorti]

        errs = None
        if lc1.IsErr and lc2.IsErr:
            e = np.concatenate((lc1.MagErr, lc2.MagErr))
            errs = e[sorti]
        res = LightCurve(bname, time, mags, errs=errs)
        return res


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

    def IsBand(self, bname):
        return bname in self.BandNames

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

    def is_band(self, bn):
        return bn in self.BandNames

    def set_tshift(self, tshift):
        for n, lc in self.Set.items():
            lc.tshift = tshift

    def set_mshift(self, mshift):
        for n, lc in self.Set.items():
            lc.mshift = mshift

    @classmethod
    def Merge(cls, curves1, curves2):
        if curves1 is None:
            return curves2
        if curves2 is None:
            return curves1

        res = SetLightCurve("{}+{}".format(curves1.Name, curves2.Name))
        # Add Light Curves from the first set
        for lc1 in curves1:
            lc2 = curves2.get(lc1.Band.Name)
            if lc2 is None:
                res.add(lc1)
            else:
                lc = LightCurve.Merge(lc1, lc2)
                res.add(lc)
        # Add remaining Light Curves from the second set
        for lc in curves2:
            if not res.IsBand(lc.Band.Name):
                res.add(lc)

        return res
