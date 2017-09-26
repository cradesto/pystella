import numpy as np

from pystella.fit.fit_lc import FitLc, FitLcResult


class FitLcMcmc(FitLc):
    def __init__(self):
        super().__init__("MCMCFit")

    def fit_lc(self, lc_o, lc_m):
        tshift, tsigma = self.fit_lc_bayesian_1d(lc_o, lc_m, is_debug=super().is_debug)
        fit_result = FitLcResult()
        fit_result.tshift = tshift
        fit_result.tsigma = tsigma
        return fit_result

    def fit_tss(self, tss_o, tss_m):
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
    def fit_lc_bayesian_1d(lc_o, lc_m, err_mdl=0.1, is_debug=True, is_plot=True):
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
            e = np.sqrt(lc_o.Err**2 + err_mdl**2)
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

    def fit_curves_bayesian(self, tss_o, tss_m):
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

        def lnprob(theta, ts_o, ts_m):
            # mean, amplitude, stddev = param
            ts_o.tshift = theta[0]
            yp = np.interp(ts_o.Time, ts_m.Time, ts_m.V)
            y = ts_o.V
            diff = (y - yp)
            if ts_o.IsErr:
                eo = ts_o.Err
                diff /= eo
            return -np.dot(diff, diff)

        def log_likelihood_curves(theta, obs, mdl):
            res = 0.
            for ts_o in obs:
                ts_m = mdl.get(ts_o.Name)
                if ts_m is not None:
                    p = lnprob(theta, ts_o, ts_m)
                    res += p
            return res

        def log_posterior(theta, obs, mdl):
            return log_prior(theta) + log_likelihood_curves(theta, obs, mdl)

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
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(tss_o, tss_m))
        # sampler.run_mcmc(starting_guesses, nsteps)

        # burnin
        pos, prob, state = sampler.run_mcmc(starting_guesses, nburn)
        sampler.reset()
        # run
        sampler.run_mcmc(pos, nsteps)

        return sampler
