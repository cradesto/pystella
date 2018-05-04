import numpy as np

from pystella.fit import mpfit
from pystella.fit.fit_lc import FitLc, FitLcResult
from scipy.interpolate import InterpolatedUnivariateSpline

class FitMPFit(FitLc):
    def __init__(self):
        super().__init__("MPFit")

    def fit_lc(self, lc_o, lc_m, dt0=0., dm0=None):
        #   res = {'dt': tshift, 'dtsig': tsigma, 'dm': dmshift, 'dmsig': dmsigma}
        res = self.best_lc(lc_o, lc_m, dt0=dt0, dm0=dm0, is_debug=self.is_debug)

        fit_result = FitLcResult()
        fit_result.tshift = res['dt']
        fit_result.tsigma = res['dtsig']
        fit_result.mshift = res['dm']
        fit_result.msigma = res['dmsig']
        return fit_result

    def fit_curves(self, curves_o, curves_m):
        dt0 = curves_o.get(curves_o.BandNames[0]).tshift
        dm0 = curves_o.get(curves_o.BandNames[0]).mshift

        result = self.best_curves(curves_o, curves_m, dt0=dt0, dm0=dm0)

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

    def best_lc(self, lc_o, lc_m, dt0=None, dm0=None, At=0., xtol=1e-10, ftol=1e-10, gtol=1e-10, is_spline=True):
        dt_init = lc_o.tshift
        dm_init = lc_o.mshift

        def least_sq(p, fjac):
            lc_o.tshift = dt_init + p[0]
            lc_o.mshift = dm_init + p[1]
            # model interpolation
            if is_spline:
                s = InterpolatedUnivariateSpline(lc_m.Time, lc_m.Mag, k=1)
                m = s(lc_o.Time )
            else:
                m = np.interp(lc_o.Time, lc_m.Time, lc_m.Mag)  # One-dimensional linear interpolation.
            w = np.ones(len(m))
            w = np.abs(1. - At*(lc_o.Time - lc_o.TimeMin) / (lc_o.TimeMax - lc_o.TimeMin))  # weight
            # w = np.abs(lc_o.Time / max(lc_o.Time))  # weight
            diff = lc_o.Mag - m
            if lc_o.IsErr:
                res = diff**2 / lc_o.MagErr**2  * w
                # res = np.abs((lc_o.Mag - m) / lc_o.MagErr) * w
            else:
                res = diff**2 * w
            return 0, res

        parinfo = []
        if dt0 is not None:
            parinfo.append({'value': dt0, 'limited': [1, 1], 'limits': [-250., 250.]})
        else:
            parinfo.append({'value': 0., 'fixed': 1})
            
        if dm0 is not None:
            parinfo.append({'value': dm0, 'limited': [1, 1], 'limits': [-5., 5.]})
        else:
            parinfo.append({'value': 0., 'fixed': 1})

#        print(parinfo)
        result = mpfit.mpfit(least_sq, parinfo=parinfo, quiet=not self.is_info, maxiter=200,
                             ftol=ftol, gtol=gtol, xtol=xtol)
        if result.status == 5:
            print('Maximum number of iterations exceeded in mangle_spectrum')
        elif result.status <= 0:
            print('error message = ', result.errmsg)
            print('parameters = ', result.params)
            raise ValueError(result.errmsg)

        tshift = result.params[0]
        
        #        dof = len(lc_o.Time) - len(result.params)  # deg of freedom
        dof = len(lc_o.Time)
        if dt0 is not None:
            dof -= 1
        if dm0 is not None:
            dof -= 1
            
        # scaled uncertainties
        pcerror = result.perror * np.sqrt(result.fnorm / dof)
        tsigma = pcerror[0]  # todo tsigma check

        dmshift = result.params[1]
        dmsigma = pcerror[1]

        lc_o.tshift = dt_init
        lc_o.mshift = dm_init
        
        if self.is_info:
            print("The final params are: tshift=%f+-%f mshift=%f+-%f  chi2: %e" % (tshift, tsigma, dmshift, dmsigma, result.fnorm))
        res = {'dt': tshift, 'dtsig': tsigma, 'dm': dmshift, 'dmsig': dmsigma, 'chi2': result.fnorm, 'dof': result.dof}
        return res

    @staticmethod
    def best_time_series(tss_o, tss_m, is_debug=False, is_info=True, xtol=1e-10, ftol=1e-10, gtol=1e-10):
        def least_sq(p, fjac):
            A = 0.5
            E = 0.002
            total = []
            for ts_m in tss_m:
                ts_o = tss_o[ts_m.Name]
                ts_o.tshift = p[0]
                # model interpolation
                # check = np.where((ts_o.Time > min(ts_m.Time)) & (ts_o.Time < max(ts_m.Time)))
                # check = range(1, len(ts_o.Time))
                # time_o = ts_o.Time[check]
                # V_o = ts_o.V[check]
                time_o = ts_o.Time
                V_o = ts_o.V

                m = np.interp(time_o, ts_m.Time, ts_m.V)  # One-dimensional linear interpolation.
                # w = np.ones(len(m))
                # w = np.exp(-2*len(check)/len(ts_o.Time)) *7 np.abs(1. - A * (time_o - min(time_o)) / (max(time_o) - min(time_o)))  # weight
                w = np.abs(1. - A * (time_o - min(time_o)) / (max(time_o) - min(time_o)))  # weight
                if ts_o.IsErr:
                    # err_o = ts_o.Err[check]
                    # res = np.abs((V_o - m) / (abs(m)*E + err_o)) * w
                    res = np.abs((V_o - m) / (abs(m)*E + ts_o.Err)) * w
                else:
                    res = np.abs(V_o - m) * w
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

    def best_curves(self, curves_o, curves_m, dt0=None, dm0=None, is_spline=True,
                    At=0., err_mdl=0.,  xtol=1e-10, ftol=1e-10, gtol=1e-10):
        dt_init = {lc.Band.Name: lc.tshift for lc in curves_o}
        dm_init = {lc.Band.Name: lc.mshift for lc in curves_o}

        def least_sq(p, fjac):
            total = []
            for lc_m in curves_m:
                lc_o = curves_o.get(lc_m.Band.Name)
                lc_o.tshift = dt_init[lc_m.Band.Name] + p[0]
                lc_o.mshift = dm_init[lc_m.Band.Name] + p[1]
                # model interpolation
                if is_spline:
                    s = InterpolatedUnivariateSpline(lc_m.Time, lc_m.Mag, k=1)
                    m = s(lc_o.Time)
                else:
                    m = np.interp(lc_o.Time, lc_m.Time, lc_m.Mag)  # One-dimensional linear interpolation.
                #                w = np.ones(len(m))
                w = np.abs(1. - At * (lc_o.Time - lc_o.TimeMin) / (lc_o.TimeMax - lc_o.TimeMin))  # weight
                diff = lc_o.Mag - m
                if lc_o.IsErr:
                    res = diff/(err_mdl + lc_o.MagErr) * w
#                    res = (diff/(err_mdl + lc_o.MagErr))**2 * w
                else:
                    res = diff * w
#                    res = diff**2 * w
                total = np.append(total, res)
            return 0, total

        parinfo = []
        if dt0 is not None:
            parinfo.append({'value': dt0, 'limited': [1, 1], 'limits': [-250., 250.]})
        else:
            parinfo.append({'value': 0., 'fixed': 1})
            
        if dm0 is not None:
            parinfo.append({'value': dm0, 'limited': [1, 1], 'limits': [-5., 5.]})
        else:
            parinfo.append({'value': 0., 'fixed': 1})

        if dt0 is not None or dm0 is not None:
            result = mpfit.mpfit(least_sq, parinfo=parinfo, quiet=not self.is_info, maxiter=200,
                                 ftol=ftol, gtol=gtol, xtol=xtol)
        else:
            dum, r = least_sq([x['value'] for x in parinfo], None)  # just return the weight diff
            r = np.array(r)
            chi2 = np.sum(r**2)
            if self.is_info:
                print("No fit: only run least_sq: len(r)={}  chi2={}".format(len(r), chi2))
                      
            res = {'dt': None, 'dtsig': None, 'dm': None, 'dmsig': None, 'chi2': chi2, 'dof': len(r)}
            return res

        # return initial states
        for lc in curves_o:
            lc.tshift = dt_init[lc.Band.Name]
            lc.mshift = dm_init[lc.Band.Name]

        if result.status <= 0:
            print('error message = ', result.errmsg)
            print('parameters = ', result.params)
            raise ValueError(result.errmsg)
            
        if self.is_info:
            print("status: ", result.status)
            if result.status == 5:
                print('Maximum number of iterations exceeded in mangle_spectrum')
            else:
                print("Iterations: ", result.niter)
                print("Fitted pars: ", result.params)
                print("Uncertainties: ", result.perror)

        # time and magnitude shifts
        tshift = result.params[0]
        dof = np.sum([lc.Length for lc in curves_o])
        if dt0 is not None:
            dof -= 1
        if dm0 is not None:
            dof -= 1
            
        # scaled uncertainties
        pcerror = result.perror * np.sqrt(result.fnorm / dof)
        tsigma = pcerror[0]  # todo tsigma check

        dmshift = result.params[1]
        dmsigma = pcerror[1]

        if self.is_info:
            print("The final params are: tshift=%f+-%f mshift=%f+-%f  chi2: %e" % (tshift, tsigma, dmshift, dmsigma, result.fnorm))
        res = {'dt': tshift, 'dtsig': tsigma, 'dm': dmshift, 'dmsig': dmsigma, 'chi2': result.fnorm, 'dof': result.dof}
        return res
#        return result
