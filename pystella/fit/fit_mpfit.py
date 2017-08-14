import numpy as np
from sklearn.metrics import mean_squared_error

from pystella.fit import mpfit

from pystella.fit.fit_lc import FitLc


class FitMPFit(FitLc):
    def __init__(self, is_debug=False):
        super().__init__(is_debug=is_debug)

    def fit(self, lc_o, lc_m):
        t, tsigma = self.run(lc_o, lc_m, is_debug=super().is_debug)
        self.tshift = t
        self.tsigma = tsigma
        return self

    def run(self, lc_o, lc_m, dt0=0., dm0=None, is_debug=True, xtol=1e-10, ftol=1e-10, gtol=1e-10):
        def least_sq(p, fjac):
            lc_o.tshift = p[0]
            lc_o.mshift = p[1]
            # model interpolation
            m = np.interp(lc_o.Time, lc_m.Time, lc_m.Mag)  # One-dimensional linear interpolation.
            w = np.ones(len(m))
            # w = np.abs(1. - lc_o.Time / max(lc_o.Time))  # weight
            # w = np.abs(lc_o.Time / max(lc_o.Time))  # weight
            res = np.abs((lc_o.Mag - m) ** 2 / m) * w
            # res = np.abs((lc_o.Mag - m) ** 2 / m) * w
            # res = np.sqrt(mean_squared_error(lc_o.Mag, m))
            # w = np.exp(-(max(abs(lc_o.Mag)) - abs(lc.Mag)) * 2)  # weight
            # w = w / max(w)
            # if lc_o.MagErr is not None:
            #     res = res * w / lc_o.MagErr
            return 0, res

        parinfo = [{'value': dt0, 'limited': [1, 1], 'limits': [0., 250.]}]
        if dm0 is not None:
            parinfo.append({'value': dm0, 'limited': [1, 1], 'limits': [-50., 50.]})
        else:
            parinfo.append({'value': lc_o.mshift, 'fixed': 1})

        result = mpfit.mpfit(least_sq, parinfo=parinfo, quiet=not is_debug, maxiter=200)
        # result = mpfit.mpfit(least_sq, parinfo=parinfo, quiet=not is_debug, maxiter=200,
        #                      ftol=ftol, gtol=gtol, xtol=xtol)
        if result.status == 5:
            print('Maximum number of iterations exceeded in mangle_spectrum')

        tshift = result.params[0]
        dof = len(lc_o.Time) - len(result.params)  # deg of freedom
        # scaled uncertainties
        pcerror = result.perror * np.sqrt(result.fnorm / dof)
        tsigma = pcerror[0]  # todo tsigma check

        mshift = result.params[1]

        if is_debug:
            print("The final params are: tshift=%f tsigma=%f mshift=%f" % (tshift, tsigma, mshift))
        return tshift, tsigma

