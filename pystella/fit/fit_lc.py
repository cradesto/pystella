class FitLc:
    def __init__(self, name):
        self._name = name
        self._par = {'is_info': False, 'is_debug': False}

    def fit_lc(self, lc_o, lc_m):
        pass

    def fit_curves(self, curves_m, curves_o):
        pass

    def fit_tss(self, tss_m, tss_o):
        pass

    @property
    def Name(self):
        return self._name

    @property
    def is_info(self):
        return self.get('is_info')

    @is_info.setter
    def is_info(self, v):
        self._par['is_info'] = v

    @property
    def is_debug(self):
        return self.get('is_debug')

    @is_debug.setter
    def is_debug(self, v):
        self._par['is_debug'] = v

    def get(self, k, default=None):
        if k in self._par:
            return self._par[k]
        return default

    def set_param(self, k, v):
        self._par[k] = v


class FitLcResult:
    def __init__(self):
        self._mshift = 0.
        self._msigma = 0.
        self._tshift = 0.
        self._tsigma = 0.
        self._measure = None
        self._comm = None

    @property
    def tshift(self):
        return self._tshift

    @tshift.setter
    def tshift(self, shift):
        self._tshift = shift

    @property
    def tsigma(self):
        return self._tsigma

    @tsigma.setter
    def tsigma(self, v):
        self._tsigma = v

    @property
    def mshift(self):
        return self._mshift

    @mshift.setter
    def mshift(self, v):
        self._mshift = v

    @property
    def msigma(self):
        return self._msigma

    @msigma.setter
    def msigma(self, v):
        self._msigma = v

    @property
    def measure(self):
        return self._measure

    @measure.setter
    def measure(self, v):
        self._measure = v

    @property
    def comm(self):
        return self._comm

    @comm.setter
    def comm(self, v):
        self._comm = v
