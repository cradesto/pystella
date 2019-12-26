import numpy as np

from pystella.rf import band
from pystella.rf.ts import TimeSeries, SetTimeSeries

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

        super().__init__(self._b.Name, time, mags, errs, tshift=tshift)
        self._mshift = mshift
        self._attrs = {}

    @property
    def Mag(self):
        return self.V + self.mshift

    @property
    def M(self):
        return self.V

    @M.setter
    def M(self, v):
        self.V = v

    @property
    def MagErr(self):
        return self.Err

    @property
    def Band(self):
        return self._b

    @property
    def BName(self):
        return self.Band.Name

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

    def attrs(self, nm, *val):
        if not val:
            return self._attrs[nm]
        else:
            self._attrs[nm] = val

    def toarray(self, is_err=True):
        if is_err and self.IsErr:
            res = np.array([self.Time, self.Mag, self.MagErr])
        else:
            res = np.array([self.Time, self.Mag])
        return res.T

    def copy_tlim(self, tlim=None):
        errs = None

        if tlim is not None:
            is_good = np.where((self.Time >= tlim[0]) & (self.Time <= tlim[1]))
            time = self.T[is_good]
            mags = self.V[is_good]
            if self.IsErr:
                errs = self.Err[is_good]
        else:
            time = self.T
            mags = self.V
            if self.IsErr:
                errs = self.Err

        lc = LightCurve(self.Band, time, mags, errs)
        lc.tshift = self.tshift
        lc.mshift = self.mshift
        return lc

    def copy(self, name=None, f=None):
        lc = super(type(self), self).copy(name=name, f=f)
        lc.mshift = self.mshift
        return lc

    def clone(self):
        errs = None

        if self.IsErr:
            errs = self.Err

        return LightCurve(self.Band, self.Time, self.Mag, errs), self.tshift, self.mshift

    def sorted_time(self, order=None):
        ind = np.argsort(self.Time, order=order)
        time = self.T[ind]
        mags = self.V[ind]
        errs = None
        if self.IsErr:
            errs = self.Err[ind]
        lc = LightCurve(self.Band, time, mags, errs)
        lc.tshift = self.tshift
        lc.mshift = self.mshift
        return lc

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
            e = np.concatenate((lc1.Err, lc2.Err))
            errs = e[sorti]
        res = LightCurve(bname, time, mags, errs=errs)
        return res


def LC_interp(orig, time, is_spline=True):
    if is_spline:
        from scipy.interpolate import InterpolatedUnivariateSpline
        s = InterpolatedUnivariateSpline(orig.Time, orig.Mag, k=1)
        mags = s(time)
    else:
        mags = np.interp(time, orig.Time, orig.Mag)
    if orig.IsErr:
        if is_spline:
            from scipy.interpolate import InterpolatedUnivariateSpline
            s = InterpolatedUnivariateSpline(orig.Time, orig.MagErr, k=1)
            errs = s(time)
        else:
            errs = np.interp(time, orig.Time, orig.MagErr)
        lc = LightCurve(orig.Band, time, mags, errs)
    else:
        lc = LightCurve(orig.Band, time, mags)

    # lc.tshift = orig.tshift
    # lc.mshift = orig.mshift
    return lc


class SetLightCurve(SetTimeSeries):
    """Set of the Light Curves"""

    def __init__(self, name=''):
        """Creates a Set of Light Curves."""
        super().__init__(name)
        # self._loop = 0

    @property
    def Bands(self):
        if len(self.Set) == 0:
            raise ValueError('There are no bands in SetLightCurve.')
        # for name, lc in self.Set.items():
        #     yield lc.Band
        res = (lc.Band for name, lc in self.Set.items())
        return res

    @property
    def BandNames(self):
        res = [b.Name for b in self.Bands]
        return res

    def IsBand(self, bname):
        return bname in self.BandNames

    def add(self, lc):
        self._set[lc.Band.Name] = lc

    def get(self, bn, default=None):
        for n, lc in self.Set.items():
            if lc.Band.Name == bn:
                return lc
        return default

    # def __getattr__(self, attr):
    #     lc = self.get(attr, None)
    #     if lc is None:
    #         raise AttributeError(attr)
    #     return lc
    #
    def is_band(self, bn):
        return bn in self.BandNames

    def set_mshift(self, mshift):
        for n, lc in self.Set.items():
            lc.mshift = mshift

    def clone(self):
        res = SetLightCurve(self.Name)
        for lc in self:
            clone = lc.clone()
            res.add(clone)
        return res

    def sorted_time(self, order=None):
        res = SetLightCurve(self.Name)
        for lc in self:
            clone = lc.sorted_time(order=order)
            res.add(clone)
        return res

    def copy(self, name=None, f=None):
        if name is None:
            name = self.Name
        res = SetLightCurve(name)
        for lc in self:
            cp = lc.copy(f=f)
            res.add(cp)
        return res

    def copy_tmlim(self, tlim=None, mlim=None):
        """
        Copy SetLightCurve to other SetLightCurve
        :param mlim: time limits, default None
        :param tlim: magnitude limits, default None
        :return:
        """
        if tlim is not None:
            res = self.copy(f=lambda x: (tlim[0] <= x.Time) & (x.Time <= tlim[1]))
        else:
            res = self.copy()

        if mlim is not None:
            res = res.copy(f=lambda x: (mlim[0] >= x.Mag) & (x.Mag >= mlim[1]))

        return res

    def merge(self, curves2):
        res = SetLightCurve.Merge(self, curves2)
        return res

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
