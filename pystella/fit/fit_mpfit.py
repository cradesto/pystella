import numpy as np

from pystella.fit import mpfit
from pystella.fit.fit_lc import FitLc


class FitMPFit(FitLc):
    def __init__(self, is_debug=False):
        super().__init__(is_debug=is_debug)

    def fit_lc(self, lc_o, lc_m):
        dt0 = lc_o.tshift
        t, tsigma = self.best_lc(lc_o, lc_m, dt0=dt0, is_debug=super().is_debug)
        self.tshift = t
        self.tsigma = tsigma
        return self

    def fit_curves(self, curves_o, curves_m):
        dt0 = curves_o.get(curves_o.BandNames[0]).tshift
        t, tsigma = self.best_curves(curves_o, curves_m, dt0=dt0, is_debug=super().is_debug)
        self.tshift = t
        self.tsigma = tsigma
        return self

    def best_lc(self, lc_o, lc_m, dt0=0., dm0=None, is_debug=True, xtol=1e-10, ftol=1e-10, gtol=1e-10):
        dm = 0.
        if dm0 is not None:
            dm = dm0

        def least_sq(p, fjac):
            lc_o.tshift = dt0 + p[0]
            lc_o.mshift = dm0 + p[1]
            # model interpolation
            m = np.interp(lc_o.Time, lc_m.Time, lc_m.Mag)  # One-dimensional linear interpolation.
            w = np.ones(len(m))
            w = np.abs(1. - (lc_o.Time - lc_o.TimeMin) / (lc_o.TimeMax - lc_o.TimeMin))  # weight
            # w = np.abs(lc_o.Time / max(lc_o.Time))  # weight
            if lc_o.IsErr:
                res = np.abs((lc_o.Mag - m) / (0.1 + lc_o.MagErr)) * w
                # res = np.abs((lc_o.Mag - m) / lc_o.MagErr) * w
            else:
                res = np.abs(lc_o.Mag - m) * w
            # res = np.abs((lc_o.Mag - m) ** 2 / m) * w
            # res = np.abs((lc_o.Mag - m) ** 2 / m) * w
            # res = np.sqrt(mean_squared_error(lc_o.Mag, m))
            # w = np.exp(-(max(abs(lc_o.Mag)) - abs(lc.Mag)) * 2)  # weight
            # w = w / max(w)
            # if lc_o.MagErr is not None:
            #     res = res * w / lc_o.MagErr
            return 0, res

        parinfo = [{'value': 0., 'limited': [1, 1], 'limits': [-250., 250.]}]
        if dm0 is not None:
            parinfo.append({'value': dm, 'limited': [1, 1], 'limits': [-50., 50.]})
        else:
            parinfo.append({'value': lc_o.mshift, 'fixed': 1})

        # result = mpfit.mpfit(least_sq, parinfo=parinfo, quiet=not is_debug, maxiter=200)
        result = mpfit.mpfit(least_sq, parinfo=parinfo, quiet=not is_debug, maxiter=200,
                             ftol=ftol, gtol=gtol, xtol=xtol)
        if result.status == 5:
            print('Maximum number of iterations exceeded in mangle_spectrum')

        tshift = dt0 + result.params[0]
        dof = len(lc_o.Time) - len(result.params)  # deg of freedom
        # scaled uncertainties
        pcerror = result.perror * np.sqrt(result.fnorm / dof)
        tsigma = pcerror[0]  # todo tsigma check

        mshift = dm + result.params[1]

        if is_debug:
            print("The final params are: tshift=%f tsigma=%f mshift=%f" % (tshift, tsigma, mshift))
        return tshift, tsigma

    def best_curves(self, curves_o, curves_m, dt0=0., dm0=None, is_debug=True, xtol=1e-10, ftol=1e-10, gtol=1e-10):
        dm = 0.
        if dm0 is not None:
            dm = dm0

        def least_sq(p, fjac):
            total = []
            for lc_m in curves_m:
                lc_o = curves_o.get(lc_m.Band.Name)
                lc_o.tshift = dt0 + p[0]
                lc_o.mshift = dm + p[1]
                # model interpolation
                m = np.interp(lc_o.Time, lc_m.Time, lc_m.Mag)  # One-dimensional linear interpolation.
                w = np.ones(len(m))
                w = np.abs(1. - (lc_o.Time - lc_o.TimeMin) / (lc_o.TimeMax - lc_o.TimeMin))  # weight
                if lc_o.IsErr:
                    res = np.abs((lc_o.Mag - m) / (0.1 + lc_o.MagErr)) * w
                else:
                    res = np.abs(lc_o.Mag - m) * w
                total = np.append(total, res)
            return 0, total

        parinfo = [{'value': 0., 'limited': [1, 1], 'limits': [-250., 250.]}]
        if dm0 is not None:
            parinfo.append({'value': dm0, 'limited': [1, 1], 'limits': [-50., 50.]})
        else:
            parinfo.append({'value': curves_o.get(curves_o.BandNames[0]).mshift, 'fixed': 1})

        result = mpfit.mpfit(least_sq, parinfo=parinfo, quiet=not is_debug, maxiter=200,
                             ftol=ftol, gtol=gtol, xtol=xtol)
        if result.status == 5:
            print('Maximum number of iterations exceeded in mangle_spectrum')

        tshift = dt0 + result.params[0]
        # scaled uncertainties
        pcerror = result.perror * np.sqrt(result.fnorm / result.dof)
        tsigma = pcerror[0]  # todo tsigma check

        mshift = dm + result.params[1]

        if is_debug:
            print("The final params are: tshift=%f tsigma=%f mshift=%f" % (tshift, tsigma, mshift))
        return tshift, tsigma
