class FitLc:
    def __init__(self, is_debug=False):
        self._is_debug = is_debug
        self._mshift = 0.
        self._tshift = 0.
        self._tsigma = 0.

    def fit(self, lc_m, lc_o):
        pass

    @property
    def is_debug(self):
        return self._is_debug

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
