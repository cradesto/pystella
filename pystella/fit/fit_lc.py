import numpy as np


class FitLc:
    def __init__(self, name):
        self._name = name
        self._par = {'is_info': False, 'is_debug': False, 'is_quiet': True}

    def print_parameters(self):
        print(f'Parmeters of {self.Name}')
        for k, v in self._par.items():
            print(f'{k:20s}: {v}')

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
    def is_quiet(self):
        return self.get('is_quiet')

    @is_quiet.setter
    def is_quiet(self, v):
        self._par['is_quiet'] = v

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

    @staticmethod
    def theta2arr(theta, nb):
        if len(theta) == 1 + nb:
            return FitLc.thetaDt2arr(theta, nb)
        elif len(theta) == 2 + nb:
            return FitLc.thetaDtDm2arr(theta, nb)

        raise (ValueError("Len(theta) is not in [{},{}].".format(nb+1, nb+2), theta))

    @staticmethod
    def thetaDt2arr(theta, nb):
        if len(theta) != 1+nb:
            raise (ValueError("Len(theta) is not {}.".format(1+nb), theta))
        return theta[0], theta[1:]

    @staticmethod
    def thetaDtDm2arr(theta, nb):
        if len(theta) != 2+nb:
            raise (ValueError("Len(theta) is not {}.".format(2 + nb), theta))
        return theta[0], theta[1], theta[2:]

    @staticmethod
    def thetaDtNames(bnames):
        res = ['dt']
        for j, bn in enumerate(bnames):
            res.append('sig_{}'.format(bn))
        return tuple(res)

    @staticmethod
    def thetaDtDmNames(bnames):
        res = ['dt', 'dm0']
        for j, bn in enumerate(bnames):
            res.append('sig_{}'.format(bn))
        return tuple(res)

    @staticmethod
    def print_thetaDt(theta, bnames):
        dt, sigs = FitLc.thetaDt2arr(theta, len(bnames))
        txt = 'dt= {:4.1f} '.format(dt)
        txt += ' '.join(['sig_{}= {:5.2f}'.format(bn, sigs[j]) for j, bn in enumerate(bnames)])
        print(txt)

    @staticmethod
    def print_thetaDtDm(theta, bnames):
        dt, dm0, sigs = FitLc.thetaDtDm2arr(theta, len(bnames))
        txt = 'dt= {:4.1f} m_0= {:4.1f} '.format(dt, dm0)
        txt += ' '.join(['sig_{}= {:5.2f}'.format(bn, sigs[j]) for j, bn in enumerate(bnames)])
        print(txt)

    @staticmethod
    def log_priorDt(theta, nb, tlim, siglim):
        dt, sigs = FitLc.thetaDt2arr(theta, nb)
        if (dt < tlim[0]) or (dt > tlim[1]):
            return -np.inf
        for x in sigs:
            if (x < -0.) or (x > siglim):
                return -np.inf
        return 0  # flat prior

    @staticmethod
    def log_priorDtDm(theta, tlim, maglim, siglim):
        dt, dm, *sigs = theta
        if (dt < tlim[0]) or (dt > tlim[1]):
            return -np.inf
        if (dm < maglim[0]) or (dm > maglim[1]):
            return -np.inf
        for x in sigs:
            if (x < -0.) or (x > siglim):
                return -np.inf
        return 0  # flat prior


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
