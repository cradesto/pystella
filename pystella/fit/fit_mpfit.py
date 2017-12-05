import numpy as np

from pystella.fit import mpfit
from pystella.fit.fit_lc import FitLc, FitLcResult


class FitMPFit(FitLc):
    def __init__(self):
        super().__init__("MPFit")

    def fit_lc(self, lc_o, lc_m):
        dt0 = lc_o.tshift
        t, tsigma = self.best_lc(lc_o, lc_m, dt0=dt0, is_debug=self.is_debug)
        fit_result = FitLcResult()
        fit_result.tshift = t
        fit_result.tsigma = tsigma
        return fit_result

    def fit_curves(self, curves_o, curves_m):
        dt0 = curves_o.get(curves_o.BandNames[0]).tshift
        dm0 = curves_o.get(curves_o.BandNames[0]).mshift

        result = self.best_curves(curves_o, curves_m, dt0=dt0, dm0=dm0,
                                  is_debug=self.is_debug, is_info=self.is_info)

        tshift = dt0 + result.params[0]
        mshift = dm0 + result.params[1]
        pcerror = result.perror
        # scaled uncertainties
        # pcerror = result.perror * np.sqrt(result.fnorm / result.dof)
        tsigma = pcerror[0]

        if self.is_info:
            print("The final params are: tshift=%f tsigma=%f mshift=%f" % (tshift, tsigma, mshift))

        fit_result = FitLcResult()
        fit_result.tshift = tshift
        fit_result.tsigma = tsigma
        fit_result.measure = result.fnorm / result.dof
        fit_result.comm = 'mpfit: status={:2d} niter={:3d} fnorm={:.2f}'.format(result.status, result.niter, result.fnorm)
        #        fit_result.comm = 'The value of the summed squared residuals for the returned parameter values.'
        return fit_result

    def fit_tss(self, tss_o, tss_m):

        result = self.best_time_series(tss_o, tss_m, is_debug=self.is_debug, is_info=self.is_info)

        tshift = result.params[0]
        pcerror = result.perror
        # scaled uncertainties
        # pcerror = result.perror * np.sqrt(result.fnorm / result.dof)
        tsigma = pcerror[0]

        if self.is_info:
            print("The final params are: tshift=%f tsigma=%f" % (tshift, tsigma))

        fit_result = FitLcResult()
        fit_result.tshift = tshift
        fit_result.tsigma = tsigma
        fit_result.measure = result.fnorm / result.dof
        fit_result.comm = 'mpfit: status={:2d} niter={:3d} fnorm={:.2f}'.format(result.status, result.niter, result.fnorm)
        #        fit_result.comm = 'The value of the summed squared residuals for the returned parameter values.'
        return fit_result

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
                res = np.abs((lc_o.Mag - m) / (0.1 + lc_o.Err)) * w
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

    @staticmethod
    def best_time_series(tss_o, tss_m, is_debug=False, is_info=True, xtol=1e-10, ftol=1e-10, gtol=1e-10):
        def least_sq(p, fjac):
            A = 0.5
            E = 0.002
            total = []
            for ts_m in tss_m:
                ts_o = tss_o[ts_m.Name]
                ts_o.tshift = p[0]
                time_o = ts_o.Time
                # model interpolation
                m = np.interp(time_o, ts_m.Time, ts_m.V)  # One-dimensional linear interpolation.
                w = np.ones(len(m))
                w = np.abs(1. - A * (time_o - min(time_o)) / (max(time_o) - min(time_o)))  # weight
                if ts_o.IsErr:
                    res = np.abs((ts_o.V - m) / (abs(m)*E + ts_o.Err)) * w
                else:
                    res = np.abs(ts_o.V - m) * w
                total = np.append(total, res)
            return 0, total

        parinfo = [{'value': 0., 'limited': [1, 1], 'limits': [-250., 250.]}]
        result = mpfit.mpfit(least_sq, parinfo=parinfo, quiet=not is_debug, maxiter=200,
                             ftol=ftol, gtol=gtol, xtol=xtol)
        if is_info:
            print("status: ", result.status)
            if result.status <= 0:
                print('error message = ', result.errmsg)
            elif result.status == 5:
                print('Maximum number of iterations exceeded in mangle_spectrum')
            else:
                print("Iterations: ", result.niter)
                print("Fitted pars: ", result.params)
                print("Uncertainties: ", result.perror)

        return result

    @staticmethod
    def best_curves(curves_o, curves_m, dt0=0., dm0=0., is_dm=False,
                    is_debug=False, is_info=True, xtol=1e-10, ftol=1e-10, gtol=1e-10):
        dm = 0.
        if dm0 is not None:
            dm = dm0

        def least_sq(p, fjac):
            A = 0.5
            total = []
            for lc_m in curves_m:
                lc_o = curves_o.get(lc_m.Band.Name)
                lc_o.tshift = dt0 + p[0]
                lc_o.mshift = dm + p[1]
                # model interpolation
                m = np.interp(lc_o.Time, lc_m.Time, lc_m.Mag)  # One-dimensional linear interpolation.
                w = np.ones(len(m))
                w = np.abs(1. - A * (lc_o.Time - lc_o.TimeMin) / (lc_o.TimeMax - lc_o.TimeMin))  # weight
                if lc_o.IsErr:
                    res = np.abs((lc_o.Mag - m) / (0.03 + lc_o.Err)) * w
                else:
                    res = np.abs(lc_o.Mag - m) * w
                total = np.append(total, res)
            return 0, total

        parinfo = [{'value': 0., 'limited': [1, 1], 'limits': [-250., 250.]}]
        if is_dm:
            parinfo.append({'value': 0., 'limited': [1, 1], 'limits': [-50., 50.]})
        else:
            parinfo.append({'value': 0., 'fixed': 1})

        result = mpfit.mpfit(least_sq, parinfo=parinfo, quiet=not is_debug, maxiter=200,
                             ftol=ftol, gtol=gtol, xtol=xtol)
        if is_info:

            print("status: ", result.status)
            if result.status <= 0:
                print('error message = ', result.errmsg)
            elif result.status == 5:
                print('Maximum number of iterations exceeded in mangle_spectrum')
            else:
                print("Iterations: ", result.niter)
                print("Fitted pars: ", result.params)
                print("Uncertainties: ", result.perror)

        return result
