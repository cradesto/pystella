import numpy as np

from pystella.fit import mpfit
from pystella.fit.fit_lc import FitLc, FitLcResult
from scipy.interpolate import InterpolatedUnivariateSpline


class FitMPFit(FitLc):
    def __init__(self):
        super().__init__("MPFit")

    @property
    def xtol(self):
        return self.get('xtol', 1e-10)

    @xtol.setter
    def xtol(self, v):
        self._par['xtol'] = v

    @property
    def ftol(self):
        return self.get('ftol', 1e-10)

    @ftol.setter
    def ftol(self, v):
        self._par['ftol'] = v

    @property
    def gtol(self):
        return self.get('gtol', 1e-10)

    @gtol.setter
    def gtol(self, v):
        self._par['gtol'] = v

    @property
    def maxiter(self):
        return self.get('maxiter', 200)

    @maxiter.setter
    def maxiter(self, v):
        self._par['maxiter'] = v

    def fit_lc(self, lc_m, lc_o, dt0=0., dm0=None):
        #   res = {'dt': tshift, 'dtsig': tsigma, 'dm': dmshift, 'dmsig': dmsigma}
        res = self.best_lc(lc_o, lc_m, dt0=dt0, dm0=dm0)

        fit_result = FitLcResult()
        fit_result.tshift = res['dt']
        fit_result.tsigma = res['dtsig']
        fit_result.mshift = res['dm']
        fit_result.msigma = res['dmsig']
        return fit_result

    def fit_curves(self, curves_m, curves_o):
        dt0 = curves_o.get(curves_o.BandNames[0]).tshift
        dm0 = curves_o.get(curves_o.BandNames[0]).mshift

        result = self.best_curves(curves_o, curves_m, dt0=dt0, dm0=dm0)
        # res = {'dt': tshift, 'dtsig': tsigma, 'dm': dmshift, 'dmsig': dmsigma,
        # 'chi2': result.fnorm, 'dof': result.dof}

        tshift = dt0 + result['dt']  # result.params[0]
        mshift = dm0 + result['dm']  # result.params[1]
        # pcerror = result.perror
        # scaled uncertainties
        # pcerror = result.perror * np.sqrt(result.fnorm / result.dof)
        tsigma = result['dtsig']  # pcerror[0]

        if self.is_info:
            print("The final params are: tshift=%f tsigma=%f mshift=%f" % (tshift, tsigma, mshift))

        fit_result = FitLcResult()
        fit_result.tshift = tshift
        fit_result.tsigma = tsigma
        fit_result.measure = result['chi2'] / result['dof']
        # fit_result.measure = result.fnorm / result.dof
        fit_result.comm = 'mpfit:  fnorm={:.2f} dof={:d}'.format(result['chi2'], result['dof'])
        # fit_result.comm = 'mpfit: status={:2d} niter={:3d} fnorm={:.2f}'.format(result.status, result.niter,
        #                                                                         result.fnorm)
        #        fit_result.comm = 'The value of the summed squared residuals for the returned parameter values.'
        return fit_result

    def fit_tss(self, tss_m, tss_o, A=0.):

        result = self.best_time_series(tss_m, tss_o, At=A, is_debug=self.is_debug, is_info=self.is_info)

        tshift = result.params[0]
        pcerror = result.perror
        # scaled uncertainties
        # pcerror = result.perror * np.sqrt(result.fnorm / result.dof)
        tsigma = pcerror[0]

        err, errsig = result.params[1], pcerror[0]

        if self.is_info:
            print("The final params are: tshift=%f tsigma=%f" % (tshift, tsigma))

        fit_result = FitLcResult()
        fit_result.tshift = tshift
        fit_result.tsigma = tsigma
        fit_result.measure = result.fnorm / result.dof
        fit_result.comm = 'mpfit: status={:2d} niter={:3d} fnorm={:.2f} err= {:.3f}+-{:.4f}'. \
            format(result.status, result.niter, result.fnorm, err, errsig)
        #        fit_result.comm = 'The value of the summed squared residuals for the returned parameter values.'
        return fit_result

    def best_lc(self, lc_m, lc_o, dt0=None, dm0=None, At=0., is_spline=True):
        dt_init = lc_o.tshift
        dm_init = lc_o.mshift

        def least_sq(p, fjac):
            lc_o.tshift = dt_init + p[0]
            lc_o.mshift = dm_init + p[1]
            # model interpolation
            if is_spline:
                s = InterpolatedUnivariateSpline(lc_m.Time, lc_m.Mag, k=1)
                m = s(lc_o.Time)
            else:
                m = np.interp(lc_o.Time, lc_m.Time, lc_m.Mag)  # One-dimensional linear interpolation.
            w = np.ones(len(m))
            w = np.abs(1. - At * (lc_o.Time - lc_o.TimeMin) / (lc_o.TimeMax - lc_o.TimeMin))  # weight
            # w = np.abs(lc_o.Time / max(lc_o.Time))  # weight
            diff = lc_o.Mag - m
            if lc_o.IsErr:
                res = diff ** 2 / lc_o.MagErr ** 2 * w
                # res = np.abs((lc_o.Mag - m) / lc_o.MagErr) * w
            else:
                res = diff ** 2 * w
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
        result = mpfit.mpfit(least_sq, parinfo=parinfo, quiet=self.is_quiet, maxiter=self.maxiter,
                             ftol=self.ftol, gtol=self.gtol, xtol=self.xtol)
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
            print("The final params are: tshift=%f+-%f mshift=%f+-%f  chi2: %e" % (
                tshift, tsigma, dmshift, dmsigma, result.fnorm))
        res = {'dt': tshift, 'dtsig': tsigma, 'dm': dmshift, 'dmsig': dmsigma, 'chi2': result.fnorm, 'dof': result.dof}
        return res

    @staticmethod
    def best_time_series(tss_m, tss_o, At=0., is_debug=False, is_info=True, xtol=1e-10, ftol=1e-10, gtol=1e-10):
        def least_sq(p, fjac):
            E = 0.01
            total = []
            for ts_m in tss_m:
                ts_o = tss_o[ts_m.Name]
                ts_o.tshift = p[0]
                sigma = p[1]
                # model interpolation
                # check = np.where((ts_o.Time > min(ts_m.Time)) & (ts_o.Time < max(ts_m.Time)))
                # check = range(1, len(ts_o.Time))
                # time_o = ts_o.Time[check]
                # V_o = ts_o.V[check]
                time_o = ts_o.Time
                V_o = ts_o.V

                m = np.interp(time_o, ts_m.Time, ts_m.V)  # One-dimensional linear interpolation.
                # sigma = sigma
                w = np.abs(1. - At * (time_o - min(time_o)) / (max(time_o) - min(time_o)))  # weight
                err_m = np.sqrt(np.sum((V_o - m) ** 2) / len(m))
                # err_m = abs(V_o - m)
                if ts_o.IsErr:
                    # err_o = ts_o.Err[check]
                    # res = np.abs((V_o - m) / (abs(m)*E + err_o)) * w
                    res = np.abs((V_o - m) / (err_m + ts_o.Err)) * w
                    # res = np.abs((V_o - m) / (err_m + ts_o.Err+sigma)) * w
                    # res = np.abs(V_o - m) * w
                else:
                    res = np.abs(V_o - m) * w
                total = np.append(total, res)
            return 0, total

        # parinfo = [{'value': 0., 'limited': [1, 1], 'limits': [-250., 250.]}]
        parinfo = [{'value': 0., 'limited': [1, 1], 'limits': [-250., 250.]},
                   {'value': 0.001, 'limited': [1, 1], 'limits': [0.001, 3.]}]
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

    def best_curves(self, curves_m, curves_o, dt0=None, dm0=None, is_spline=True, At=0., **kwargs):
        """
        Find the values of time shift and magnitude shift minimizing the distance between the observational
        and modal light curves
        :param curves_m: modal LCs
        :param curves_o: observational LCs
        :param dt0: initial time shift
        :param dm0: magnitude shift, default: 0.
        :param is_spline:
        :param At:
        :param kwargs:
                dt_limits = kwargs.get("dt_limits", [-250., 250.])
                dm_limits = kwargs.get("dm_limits", [-5., 5])
                band_limits = kwargs.get("band_limits", [0.001, 4.])
                is_fit_sigmas = kwargs.get("is_fit_sigmas", False)
        :return:
        """
        dt_init = {lc.Band.Name: lc.tshift for lc in curves_o}
        dm_init = {lc.Band.Name: lc.mshift for lc in curves_o}
        bnames = curves_o.BandNames

        dt_limits = kwargs.get("dt_limits", [-250., 250.])
        dm_limits = kwargs.get("dm_limits", [-5., 5])
        band_limits = kwargs.get("band_limits", [0.001, 4.])
        is_fit_sigmas = kwargs.get("is_fit_sigmas", False)
        # заменить модели их сплайном
        curves_m_spline = {}
        if is_spline:
            for lc_m in curves_m:
                curves_m_spline[lc_m.Band.Name] = InterpolatedUnivariateSpline(lc_m.Time, lc_m.Mag, k=1)

        def least_sq_simga(p, fjac):
            dt, dm, sigs = FitMPFit.thetaDtDm2arr(p, len(bnames))
            # dt, dm = p
            # dt, dm, sigma = p
            total = []
            for i, bname in enumerate(bnames):
                lc_o = curves_o[bname]
                to, mo, eo = lc_o.Time, lc_o.Mag, lc_o.MagErr
                # lc_o = curves_o.get(lc_m.Band.Name)
                # lc_o.tshift = dt_init[lc_m.Band.Name] + dt
                # lc_o.mshift = dm_init[lc_m.Band.Name] + dm
                # model interpolation
                lc_m = curves_m[bname]
                tm, mm = lc_m.Time, lc_m.Mag
                mm += dm
                tm += dt
                em = sigs[i]
                m = np.interp(to, tm, mm)  # One-dimensional linear interpolation.
                #                w = np.ones(len(m))
                # err_m = abs(min(m)-m) * err_mdl
                # w = np.abs(1. - At * (lc_o.Time - lc_o.TimeMin) / (lc_o.TimeMax - lc_o.TimeMin))  # weight
                diff = lc_o.Mag - m
                sig = np.ones(lc_o.Length) * em
                if lc_o.IsErr:
                    sig = np.sqrt(em**2 + lc_o.MagErr**2)

                    # err_m = np.sqrt(np.sum(diff ** 2) / len(m))
                inv_sigma2 = 1. / sig**2
                chi = np.sqrt(np.abs(diff**2 * inv_sigma2 - np.log(inv_sigma2) + np.log(2 * np.pi)))
                # chi = np.sum(diff**2 * inv_sigma2 - np.log(inv_sigma2)) + np.log(2 * np.pi) * len(diff)
                # chi = np.sum(diff**2 * inv_sigma2 - np.log(inv_sigma2) + np.log(2 * np.pi))
                # chi = diff / sig #+ np.log(sig*np.pi)/lc_o.Length
                # chi += fjac[0]
                # chi = np.ones(lc_o.Length) * chi
                total = np.append(total, chi)
            return 0, total

        def least_sq(p, fjac):
            dt, dm, sigs = FitMPFit.thetaDtDm2arr(p, len(bnames))
            total = []
            for lc_m in curves_m:
                lc_o = curves_o.get(lc_m.Band.Name)
                to, mo = lc_o.Time, lc_o.Mag
                to -= dt
                mo -= dm
                # model interpolation
                if is_spline:
                    s = curves_m_spline[lc_m.Band.Name]
                    mm = s(to)
                else:
                    mm = np.interp(to, lc_m.Time, lc_m.Mag)  # One-dimensional linear interpolation.
                #                w = np.ones(len(m))
                # err_m = abs(min(m)-m) * err_mdl
                w = np.abs(1. - At * (lc_o.Time - lc_o.TimeMin) / (lc_o.TimeMax - lc_o.TimeMin))  # weight
                diff = mo - mm
                err_m = np.sqrt(np.sum(diff ** 2) / len(mm))
                if lc_o.IsErr:
                    chi = diff * w / (err_m + lc_o.MagErr)
                else:
                    chi = diff * w
                #                    res = diff**2 * w
                total = np.append(total, chi)
            return 0, total

        parinfo = []  # [{'value': 0.0, 'limited': [1, 1], 'limits': [0., 3.]}]  # todo changeable parameters
        if dt0 is not None:
            dt_limits = [x+dt0 for x in dt_limits]
            parinfo.append({'value': dt0, 'limited': [1, 1], 'limits': dt_limits})
        else:
            parinfo.append({'value': 0., 'fixed': 1})

        if dm0 is not None:
            dm_limits = [x + dm0 for x in dm_limits]
            parinfo.append({'value': dm0, 'limited': [1, 1], 'limits': dm_limits})
        else:
            parinfo.append({'value': 0., 'fixed': 1})

        for i, bname in enumerate(bnames):
            parinfo.append({'value': 0.1, 'limited': [1, 1], 'limits': band_limits})  # todo changeable parameters

        if dt0 is not None or dm0 is not None:
            func = least_sq
            if is_fit_sigmas:
                func = least_sq_simga
            result = mpfit.mpfit(func, parinfo=parinfo, quiet=self.is_quiet, maxiter=self.maxiter,
                                 ftol=self.ftol, gtol=self.gtol, xtol=self.xtol)
        else:
            dum, r = least_sq([x['value'] for x in parinfo], None)  # just return the weight diff
            r = np.array(r)
            chi2 = np.sum(r ** 2)
            if self.is_info:
                print("No fit: only run least_sq: len(r)={}  chi2={}".format(len(r), chi2))

            res = {'dt': None, 'dtsig': None, 'dm': None, 'dmsig': None, 'chi2': chi2, 'dof': len(r)}
            return res

        # return initial states
        for lc in curves_o:
            lc.tshift = dt_init[lc.Band.Name]
            lc.mshift = dm_init[lc.Band.Name]

        if result.status <= 0:
            print('status = ', result.status)
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

        dof = np.sum([lc.Length for lc in curves_o])
        if dt0 is not None:
            dof -= 1
        if dm0 is not None:
            dof -= 1

        # pcerror = result.perror  # * np.sqrt(result.fnorm / dof)

        # time and magnitude shifts
        dt, dtsig = result.params[0], result.perror[0]
        dm, dmsig = result.params[1], result.perror[1]

        if self.is_info:
            print("The final params are: tshift=%f+-%f mshift=%f+-%f  chi2: %e" % (
                dt, dtsig, dm, dmsig, result.fnorm))

        res = {'dt': dt, 'dtsig': dtsig, 'dm': dm, 'dmsig': dmsig, 'chi2': result.fnorm, 'dof': result.dof}
        # scaled uncertainties
        if len(result.params) > 2:
            err, errsig = result.params[2:], result.perror[2:]
            res['err'] = err
            res['errsig'] = errsig
            if self.is_info:
                for i, bname in enumerate(bnames):
                    print("The sigma_{}= {} +- {}".format(bname, err[i], errsig[i]))

        fit_result = FitLcResult()
        fit_result.tshift = res['dt']
        fit_result.tsigma = res['dtsig']
        fit_result.mshift = dm
        fit_result.msigma = res['dmsig']
        fit_result.measure = res['chi2']
        fit_result.comm = 'result MCMC:dof: {}'.format(res['dof'])

        return fit_result, res, None  # as Fit_MCMC

    #        return result

    def best_curves_gp(self, curves_m, curves_o, dt0=None, dm0=None, Npoints=None,
                       At=0., err_mdl=0.):
        from pystella.fit.fit_gp import FitGP

        dt_init = {lc.Band.Name: lc.tshift for lc in curves_o}
        dm_init = {lc.Band.Name: lc.mshift for lc in curves_o}

        # заменить модель на  сплайны
        curves_m_spline = {}
        for lc_m in curves_m:
            curves_m_spline[lc_m.Band.Name] = InterpolatedUnivariateSpline(lc_m.Time, lc_m.Mag, k=1)

        # Вычислить chi
        def least_sq(p, fjac):
            total = []
            for lc_m in curves_m:
                lc_o = curves_o.get(lc_m.Band.Name)
                lc_o.tshift = dt_init[lc_m.Band.Name] + p[0]
                lc_o.mshift = dm_init[lc_m.Band.Name] + p[1]
                N = Npoints if Npoints is not None else lc_o.Length
                # to = lc_o.Time
                # Chick time ranges
                tmin = max(lc_o.TimeMin, lc_m.TimeMin)
                tmax = min(lc_o.TimeMax, lc_m.TimeMax)

                # заменить наблюдения их гауссовым процессом
                if True and lc_o.IsErr:
                    to_gp = np.linspace(tmin, tmax, N)
                    gp_o = FitGP.lc2gp(lc_o)
                    mo, err_o = gp_o.predict(to_gp.reshape(-1, 1), return_std=True)
                else:
                    to_gp = lc_o.Time
                    mo = lc_o.Mag
                    err_o = lc_o.MagErr

                # model interpolation
                s = curves_m_spline[lc_m.Band.Name]
                m = s(to_gp)
                #                w = np.ones(len(m))
                w = np.abs(1. - At * (to_gp - tmin) / (tmax - tmin))  # weight
                diff = mo - m
                if lc_o.IsErr:
                    deviates = diff * w / (err_mdl + err_o)
                else:
                    deviates = diff * w
                total = np.append(total, deviates)
            return 0, total

        parinfo = []
        if dt0 is not None:
            parinfo.append({'value': dt0, 'limited': [1, 1], 'limits': [-250. + dt0, 250. + dt0]})
        else:
            parinfo.append({'value': 0., 'fixed': 1})

        if dm0 is not None:
            parinfo.append({'value': dm0, 'limited': [1, 1], 'limits': [-5. + dm0, 5. + dm0]})
        else:
            parinfo.append({'value': 0., 'fixed': 1})

        if dt0 is not None or dm0 is not None:
            result = mpfit.mpfit(least_sq, parinfo=parinfo, quiet=self.is_quiet, maxiter=self.maxiter,
                                 ftol=self.ftol, gtol=self.gtol, xtol=self.xtol)
        else:
            dum, r = least_sq([x['value'] for x in parinfo], None)  # just return the weight diff
            r = np.array(r)
            chi2 = np.sum(r ** 2)
            if self.is_info:
                print("No fit: only run least_sq: len(r)={}  chi2={}".format(len(r), chi2))

            res = {'dt': None, 'dtsig': None, 'dm': None, 'dmsig': None, 'chi2': chi2, 'dof': len(r)}
            return res

        # return initial states
        for lc in curves_o:
            lc.tshift = dt_init[lc.Band.Name]
            lc.mshift = dm_init[lc.Band.Name]

        if result.status <= 0:
            print('status = ', result.status)
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
            print("The final params are: tshift=%f+-%f mshift=%f+-%f  chi2: %e" % (
                tshift, tsigma, dmshift, dmsigma, result.fnorm))
        res = {'dt': tshift, 'dtsig': tsigma, 'dm': dmshift, 'dmsig': dmsigma, 'chi2': result.fnorm, 'dof': result.dof}
        return res
