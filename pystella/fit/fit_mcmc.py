import numpy as np

from pystella.fit.fit_lc import FitLc, FitLcResult

import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)
import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)
class FitMCMC(FitLc):
    def __init__(self):
        super().__init__("MCMCFit")

    def fit_lc(self, lc_o, lc_m):
        tshift, tsigma = self.fit_lc_bayesian_1d(lc_o, lc_m, is_debug=self.is_debug)
        fit_result = FitLcResult()
        fit_result.tshift = tshift
        fit_result.tsigma = tsigma
        return fit_result

    def fit_tss(self, tss_m, tss_o):
        # dt0 = first(tss_o).get(tss_o.BandNames[0]).tshift

        sampler = self.fit_curves_bayesian(tss_o, tss_m)

        af = sampler.acceptance_fraction.mean()
        chain = sampler.flatchain
        tshift = chain.mean()
        tsigma = np.std(chain)
        walkers, steps, dim = sampler.chain.shape

        if self.is_info:
            txt = [
                "Walkers: %i" % walkers,
                "Steps:   %i" % steps,
                "Acceptance fraction: %0.2f" % af,
                "-------------------------",
                "tshift      = %0.3f +/- %0.4f" % (tshift, tsigma),
            ]
            print('\n'.join(txt))
        #     sample = sampler.chain  # shape = (nwalkers, nsteps, ndim)
        # sample_dt = sampler.chain[:, nburn:, 0].ravel()  # discard burn-in points
        # sample_mu = sampler.chain[:, nburn:, 1].ravel()  # discard burn-in points
        #     sample = sampler.chain[:, nburn:, :]  # discard burn-in points

        # tshift, tsigma = np.mean(sample_dt), np.std(sample_dt)
        # dm, dmsigma = np.mean(sample_mu), np.std(sample_mu)

        # res = {'dt': (tshift, tsigma), 'measure': np.mean(sampler.acceptance_fraction)}
        # # plot a histogram of the sample
        # if is_info:
        #     # the fraction of steps accepted for each walker.
        #     print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
        # #     plt.hist(sample_dt, bins=50, histtype="stepfilled", alpha=0.3, normed=True)
        #     samples = sampler.chain[:, nburn:, :].reshape((-1, ndim))
        #     # fig = corner.corner(samples, labels=["$delay$", "$dm$"],
        #                     truths=[tshift, dm])

        if self.is_info:
            print("The final params are: tshift=%f tsigma=%f" % (tshift, tsigma))

        fit_result = FitLcResult()
        fit_result.tshift = tshift
        fit_result.tsigma = tsigma
        fit_result.measure = 1 / af
        fit_result.comm = 'result MCMC: np.mean(sampler.acceptance_fraction)).'
        return fit_result

    @staticmethod
    def fit_lc_bayesian_1d(lc_m, lc_o,  err_mdl=0.1, is_debug=True, is_plot=True):
        import emcee
        from sklearn.metrics import mean_squared_error

        def log_prior(theta):
            if theta[0] < -400:
                return 0
            if theta[0] > 400:
                return 0
            return 1  # flat prior

        def log_likelihood_t(theta, lc_o, lc_m):
            t = lc_o.Time
            mo = lc_o.Mag
            e = np.sqrt(lc_o.Err ** 2 + err_mdl ** 2)
            lc_m.tshift = -theta[0]
            # tck = interpolate.splrep(lc_m.Time, lc_m.Mag, s=0)
            # m_mdl = interpolate.splev(t, tck, der=0)
            m_mdl = np.interp(t, lc_m.Time, lc_m.Mag, 0, 0)  # One-dimensional linear interpolation.
            # diff = mo-m_mdl
            rms = np.sqrt(mean_squared_error(mo, m_mdl))
            return -rms
            # return -0.5 * np.sum(np.log(2 * np.pi * (e ** 2))
            #                      + (mo - m_mdl) ** 2 / e ** 2)

        def log_posterior(theta, lc_o, lc_m):
            return log_prior(theta) + log_likelihood_t(theta, lc_o, lc_m)

        ndim = 1  # number of parameters in the model
        nwalkers = 50  # number of MCMC walkers
        nburn = 1000  # "burn-in" period to let chains stabilize
        nsteps = 2000  # number of MCMC steps to take
        if is_debug:
            nburn = 100  # "burn-in" period to let chains stabilize
            nsteps = 200  # number of MCMC steps to take

        # we'll start at random locations between 0 and 2000
        starting_guesses = 100. * np.random.rand(nwalkers, ndim)

        # fit
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(lc_o, lc_m))
        sampler.run_mcmc(starting_guesses, nsteps)

        sample = sampler.chain  # shape = (nwalkers, nsteps, ndim)
        sample = sampler.chain[:, nburn:, :].ravel()  # discard burn-in points

        # plot a histogram of the sample
        # if is_plot:
        #     plt.hist(sample, bins=50, histtype="stepfilled", alpha=0.3, normed=True)

        tshift = np.mean(sample)
        tsigma = np.std(sample)

        return tshift, tsigma

    def fit_curves_bayesian(self, tss_m, tss_o):
        import emcee

        def log_prior(theta):
            return 1  # flat prior

        #
        # def log_likelihood_lc(theta, lc_o, lc_m):
        #     to, mo = lc_o.Time, lc_o.Mag
        #     if lc_o.IsErr:
        #         eo = lc_o.MagErr
        #     else:
        #         eo = np.ones(lc_o.Length)
        #     lc_o.tshift = theta[0]
        #     # tck = interpolate.splrep(lc_m.Time, lc_m.Mag, s=0)
        #     # m_mdl = interpolate.splev(to, tck, der=0)
        #     m_mdl = np.interp(lc_o.Time, lc_m.Time, lc_m.Mag)  # One-dimensional linear interpolation.
        #
        #     return -0.5 * np.sum(np.log(2 * np.pi * (eo ** 2))
        #                          + (mo - m_mdl) ** 2 / eo ** 2)

        def lnprob(theta, ts_m, ts_o):
            # mean, amplitude, stddev = param
            ts_o.tshift = theta[0]
            yp = np.interp(ts_o.Time, ts_m.Time, ts_m.V)
            y = ts_o.V
            diff = (y - yp)
            if ts_o.IsErr:
                eo = ts_o.Err
                diff /= eo
            return -np.dot(diff, diff)

        def log_likelihood_curves(theta, mdl, obs):
            res = 0.
            for ts_o in obs:
                ts_m = mdl.get(ts_o.Name)
                if ts_m is not None:
                    p = lnprob(theta, ts_m, ts_o)
                    res += p
            return res

        def log_posterior(theta, mdl, obs ):
            return log_prior(theta) + log_likelihood_curves(theta, mdl, obs)

        ndim = 1  # number of parameters in the model
        nwalkers = 50  # number of MCMC walkers
        if self.is_debug:
            nburn = 100  # "burn-in" period to let chains stabilize
            nsteps = 200  # number of MCMC steps to take
        else:  # run
            nburn = 1000  # "burn-in" period to let chains stabilize
            nsteps = 2000  # number of MCMC steps to take
        # we'll start at random locations between 0 and 2000
        starting_guesses = np.random.rand(nwalkers, ndim)
        starting_guesses[0] *= 100  # for time delay in [0,100]

        # fit
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(tss_m, tss_o))
        # sampler.run_mcmc(starting_guesses, nsteps)

        # burnin
        pos, prob, state = sampler.run_mcmc(starting_guesses, nburn)
        sampler.reset()
        # run
        sampler.run_mcmc(pos, nsteps)

        return sampler

    def best_curves(self, curves_m, curves_o, dt0=None, dm0=None, is_spline=False,
                    is_plot=False, threads=1, dt_lim=100., dm_lim=5.):
        from scipy.interpolate import InterpolatedUnivariateSpline
        import emcee

        dt_init = {lc.Band.Name: lc.tshift for lc in curves_o}
        dm_init = {lc.Band.Name: lc.mshift for lc in curves_o}

        # заменить модели их сплайном
        curves_m_spline = {}
        if is_spline:
            for lc_m in curves_m:
                curves_m_spline[lc_m.Band.Name] = InterpolatedUnivariateSpline(lc_m.Time, lc_m.Mag, k=1)

        # def log_prior(theta):
        #     #         if (theta[0] < -200) or (theta[0] > 200):
        #     #             return -np.inf
        #     #         if (theta[1] < -3) or (theta[1] > 3):
        #     #             return -np.inf
        #     return 0  # flat prior

        def log_prior(theta, tlim=dt_lim, maglim=dm_lim):
            dt, dm, lnf = theta
            if (dt < -tlim) or (dt > tlim):
                return -np.inf
            if (dm < -maglim) or (dm > maglim):
                return -np.inf
            if (lnf < -10.) or (lnf > 1.):
                return -np.inf
            return 0  # flat prior

        def log_likelihood_curves(theta, cs_m, cs_o):
            dt, dm, lnf = theta
            res = 0.
            for lc_o in cs_o:
                to, mo, eo = lc_o.Time, lc_o.Mag, lc_o.MagErr
                to = to + dt
                mo = mo + dm

                # lc_o.tshift = dt_init[lcm.BName] + dt
                # lc_o.mshift = dm_init[lcm.BName] + dm

                # model interpolation
                lcm = cs_m[lc_o.BName]
                if is_spline:
                    s = curves_m_spline[lcm.BName]
                    model = s(to)
                else:
                    model = np.interp(to, lcm.Time, lcm.Mag)  # One-dimensional linear interpolation

                # to, mo, eo = lc_o.Time, lc_o.Mag, lc_o.MagErr

                inv_sigma2 = 1.0 / (eo ** 2 + model ** 2 * np.exp(2 * lnf))
                chi = -0.5 * (np.sum((mo - model) ** 2 * inv_sigma2 - np.log(inv_sigma2)))
                res += chi
            return res

        def chi2(theta, cs_m, cs_o):
            dt, dm, lnf = theta
            res = 0.
            for lc_o in cs_o:
                to, mo, eo = lc_o.Time, lc_o.Mag, lc_o.MagErr
                to = to + dt
                mo = mo + dm

                # model interpolation
                lcm = cs_m[lc_o.BName]
                if is_spline:
                    s = curves_m_spline[lcm.BName]
                    model = s(to)
                else:
                    model = np.interp(to, lcm.Time, lcm.Mag)  # One-dimensional linear interpolation

                # chi = np.sum((mo - m) ** 2 / (eo ** 2))
                inv_sigma2 = 1.0 / (eo ** 2 + model ** 2 * np.exp(2 * lnf))
                chi = np.sum((mo - model) ** 2 * inv_sigma2)
                res += chi
            return res

        def log_posterior(theta, cs_m, cs_o):
            lp = log_prior(theta)
            if not np.isfinite(lp):
                return -np.inf
            return lp + log_likelihood_curves(theta, cs_m, cs_o)

        if dt0 is not None or dm0 is not None:
            ndim = 3  # number of parameters in the model
            if self.is_debug:
                nwalkers = 100  # number of MCMC walkers
                nburn = 50  # "burn-in" period to let chains stabilize
                nsteps = 100  # number of MCMC steps to take
            else:  # run
                nwalkers = 500  # number of MCMC walkers
                nburn = 100  # "burn-in" period to let chains stabilize
                nsteps = 500  # number of MCMC steps to take
            # we'll start at random locations between 0 and 2000
            # starting_guesses = 1. - 2. * np.random.rand(nwalkers, ndim)  #
            # starting_guesses[:, 0] *= dt_lim  # for time delay in [0,100]
            starting_guesses = np.random.rand(nwalkers, ndim)  #
            # starting_guesses = 1. - 2. * np.random.rand(nwalkers, ndim)  #
            starting_guesses[:, 0] *= dt_lim  # for time delay in [-100,100]
            starting_guesses[:, 1] *= dm_lim  # for  magnification
            starting_guesses[:, 2] = 1. - 2. * np.random.rand(nwalkers)  # for lnf

            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(curves_m, curves_o),
                                            threads=threads, a=3)
            sampler.run_mcmc(starting_guesses, nsteps)

            #     sample = sampler.chain  # shape = (nwalkers, nsteps, ndim)
            sample_dt = sampler.flatchain[nburn:, 0].ravel()  # discard burn-in points
            sample_mu = sampler.flatchain[nburn:, 1].ravel()  # discard burn-in points
            sample_lnf = sampler.flatchain[nburn:, 2].ravel()  # discard burn-in points

            # dt, dtsigma = np.mean(sample_dt), np.std(sample_dt)
            # dm, dmsigma = np.mean(sample_mu), np.std(sample_mu)
            samples = np.array([sample_dt, sample_mu, sample_lnf]).T
            quartiles = list(
                map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]), zip(*np.percentile(samples, (16, 50, 84), axis=0))))
            dt, dtsig1, dtsig2 = quartiles[0][0], quartiles[0][1], quartiles[0][2]
            dm, dmsig1, dmsig2 = quartiles[1][0], quartiles[1][1], quartiles[1][2]
            lnf, lnfsig1, lnfsig2 = quartiles[2][0], quartiles[2][1], quartiles[2][2]

            # res = {'dt': tshift, 'dtsig': tsigma, 'dm': dm, 'dmsig': dmsigma,
            #        'mean': np.mean(sampler.acceptance_fraction)}
            # chi2 = chi2([dt, dm], curves_m, curves_o)
            # dof = np.sum([lc.Length for lc in curves_o]) - 2
            # res = {'dt': dt, 'dtsig': dtsigma, 'dm': dm, 'dmsig': dmsigma,
            #        'chi2': chi2, 'dof': dof}
            chi2 = chi2([dt, dm, lnf], curves_m, curves_o)
            # chi2 = abs(log_likelihood_chi2((dt, dm), curves_m, curves_o))
            dof = np.sum([lc.Length for lc in curves_o]) - 2
            res = {'dt': dt, 'dtsig1': dtsig1, 'dtsig2': dtsig2, 'dtsig': (dtsig1+dtsig2)/2.,
                   'dm': dm, 'dmsig1': dmsig1, 'dmsig2': dmsig2, 'dmsig': (dmsig1+dmsig2)/2.,
                   'lnf': lnf, 'lnfsig1': lnfsig1, 'lnfsig2': lnfsig2, 'lnfsig': (lnfsig1+lnfsig2)/2.,
                   'chi2': chi2, 'dof': dof}
        else:
            chi2 = chi2([0., 0.], curves_m, curves_o)  # just return the weight diff
            dof = np.sum([lc.Length for lc in curves_o]) - 2
            if self.is_info:
                print("No fit: only run least_sq:   chi2={} dof={}".format(chi2, dof))

            res = {'dt': None, 'dtsig': None, 'dm': None, 'dmsig': None, 'chi2': chi2, 'dof': dof}
            return res

        # return initial states
        for lc in curves_o:
            lc.tshift = dt_init[lc.Band.Name]
            lc.mshift = dm_init[lc.Band.Name]

        # plot a histogram of the sample
        if self.is_info:
            # the fraction of steps accepted for each walker.
            print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))

        if self.is_info:
            print("The final params are: dt={:6.2f}+-{:6.4f} dm={:6.2f}+-{:6.4f}  chi2:{:.0f}/{:d}".format(
                res['dt'], res['dtsig'], res['dm'], res['dmsig'], res['chi2'], res['dof']))

        if is_plot:
            from matplotlib import pyplot as plt
            import corner
            plt.hist(sample_dt, bins=50, histtype="stepfilled", alpha=0.3)
            # samples = sampler.chain[:, nburn:, :].reshape((-1, ndim))
            fig = corner.corner(samples, labels=["$delay$", "$dm$", "$lnf$"],
                                truths=[dt, dm, lnf])
            return res, fig
        return res
