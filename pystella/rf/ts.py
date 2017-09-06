import numpy as np


class TimeSeries(object):
    def __init__(self, time, values, errs=None, name=None, tshift=0.):
        self._name = name
        self._t = np.array(time, copy=True)
        self._v = np.array(values, copy=True)
        self._e = None
        self._tshift = tshift
        if errs is not None:
            self._e = np.array(errs, copy=True)

    @property
    def Name(self):
        return self._name

    @property
    def tshift(self):
        return self._tshift

    @tshift.setter
    def tshift(self, shift):
        self._tshift = shift


    @property
    def T(self):
        return self._t

    @property
    def Time(self):
        return self.T + self.tshift

    @property
    def TimeMin(self):
        return np.min(self.Time)

    @property
    def TimeMax(self):
        return np.max(self.Time)

    @property
    def Length(self):
        return len(self._t)

    @property
    def V(self):
        return self._v

    @property
    def IsErr(self):
        return self._e is not None

    @property
    def Err(self):
        if self.IsErr:
            return self._e
        else:
            return np.zeros(self.Length)

    @property
    def Tmin(self):
        return np.min(self.T)

    @property
    def Tmax(self):
        return np.max(self.T)

    @property
    def TVmax(self):
        idx = np.argmin(self.V)
        return self.T[idx]