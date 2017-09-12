import numpy as np

from pystella.util.arr_dict import first


class TimeSeries(object):
    """Data of f(t)"""
    def __init__(self, name, time, values, errs=None, tshift=0.):
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

    @V.setter
    def V(self, value):
        self._v = value

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

    @classmethod
    def Merge(cls, ts1, ts2):
        if ts1.Name != ts2.Name:
            raise ValueError("Merging is possible only for the same filters: {} VS {}".
                             format(ts1.Name, ts2.Name))
        nm = ts1.Name
        t = np.concatenate((ts1.Time, ts2.Time))
        v = np.concatenate((ts1.V, ts2.V))

        sorti = np.argsort(t)
        time = t[sorti]
        values = v[sorti]

        errs = None
        if ts1.IsErr and ts2.IsErr:
            e = np.concatenate((ts1.Err, ts2.Err))
            errs = e[sorti]
        res = TimeSeries(nm, time, values, errs=errs)
        return res


class SetTimeSeries(object):
    """Set of the TimeSeries"""
    def __init__(self, name=''):
        """Creates a Set of TimeSeries."""
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
    def Names(self):
        if len(self.Set) == 0:
            raise ValueError('There are no bands in SetLightCurve.')
        # for name, lc in self.Set.items():
        #     yield lc.Band
        res = (name for name, ts in self.Set.items())
        return res

    @property
    def TimeDef(self):
        f = first(self.Set)
        return f.Time

    @property
    def tmin(self):
        res = [ts.tmin for name, ts in self.Set.items()]
        return min(res)

    def IsName(self, name):
        return name in self.Names

    # for cycle
    def __getitem__(self, index):
        return self.Set[index]

    # def __iter__(self):
    #     self._loop = 0
    #     return self
    def __iter__(self):
        for k, ts in self.Set.items():
            yield ts

    def __len__(self):
        return len(self.Set)

    def add(self, ts):
        self._set[ts.Name] = ts

    def rm(self, name):
        return self._set.pop(name, None)

    def get(self, bn):
        for n, ts in self.Set.items():
            if n == bn:
                return ts
        return None

    def set_tshift(self, tshift):
        for n, ts in self.Set.items():
            ts.tshift = tshift
