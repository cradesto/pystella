import logging
from collections import namedtuple

import numpy as np

from pystella.fit.fit_lc import FitLc, FitLcResult

logger = logging.getLogger(__name__)
# logger.setLevel(logging.INFO)
# logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

try:
    import emcee
except ImportError:
    logging.debug('emcee failed to import', exc_info=True)
    pass

mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)


class FitMCMC(FitLc):
    def __init__(self, nwalkers=200, nburn=200, nsteps=1000):
        super().__init__("MCMCFit")
        self._par = {
            "nwalkers": nwalkers,  # number of MCMC walkers
            "nburn": nburn,  # "burn-in" period to let chains stabilize
            "nsteps": nsteps  # number of MCMC steps to take
        }

    @property
    def nwalkers(self):
        return self.get('nwalkers')

    @nwalkers.setter
    def nwalkers(self, v):
        self._par['nwalkers'] = v

    @property
    def nburn(self):
        return self.get('nburn')

    @nburn.setter
    def nburn(self, v):
        self._par['nburn'] = v

    @property
    def nsteps(self):
        return self.get('nsteps')

    @nsteps.setter
    def nsteps(self, v):
        self._par['nsteps'] = v

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
            print("The final params are: tshift=%f +- %f" % (tshift, tsigma))

        fit_result = FitLcResult()
        fit_result.tshift = tshift
        fit_result.tsigma = tsigma
        fit_result.measure = 1 / af
        fit_result.comm = 'result MCMC: np.mean(sampler.acceptance_fraction)).'
        return fit_result

    @staticmethod
    def fit_lc_bayesian_1d(lc_m, lc_o, err_mdl=0.1, is_debug=True, is_plot=True):
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
            return 0.  # flat prior

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
            # ts_o.tshift = theta[0]
            to, vo = ts_o.Time, ts_o.V
            to += theta[0]

            yp = np.interp(to, ts_m.Time, ts_m.V)
            diff = (vo - yp)
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

        def log_posterior(theta, mdl, obs):
            return log_prior(theta) + log_likelihood_curves(theta, mdl, obs)

        ndim = 1  # number of parameters in the model
        nwalkers = 50  # number of MCMC walkers
        # if self.is_debug:
        #     nburn = 100  # "burn-in" period to let chains stabilize
        #     nsteps = 200  # number of MCMC steps to take
        # else:  # run
        #     nburn = 1000  # "burn-in" period to let chains stabilize
        #     nsteps = 2000  # number of MCMC steps to take
        # we'll start at random locations between 0 and 2000
        starting_guesses = np.random.rand(nwalkers, ndim)
        starting_guesses[0] *= 100  # for time delay in [0,100]

        # fit
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(tss_m, tss_o))
        # sampler.run_mcmc(starting_guesses, nsteps)

        # burnin
        pos, prob, state = sampler.run_mcmc(starting_guesses, self.nburn)
        sampler.reset()
        # run
        sampler.run_mcmc(pos, self.nsteps)

        return sampler

    @staticmethod
    def log_likelihood_curves_Dt(theta, cs_m, cs_o, bnames):
        texp, sigs = FitMCMC.thetaDt2arr(theta, len(bnames))
        # texp, dm0, dts, dms, sigs = theta2arr(theta, len(images), len(bnames))

        res = 0.
        for i, bname in enumerate(bnames):
            lc_o = cs_o[bname]
            to, mo, eo = lc_o.Time, lc_o.Mag, lc_o.MagErr
            # to = to + dt
            # model interpolation
            lcm = cs_m[bname]
            tm, mm = lcm.Time, lcm.Mag
            tm += texp
            em = sigs[i]
            m_fit = np.interp(to, tm, mm)  # One-dimensional linear interpolation

            inv_sigma2 = 1. / (eo ** 2 + em ** 2)
            chi = np.sum((mo - m_fit) ** 2 * inv_sigma2 - np.log(inv_sigma2) + np.log(2*np.pi))
            res += chi
        res = -0.5 * res
        return res

    @staticmethod
    def log_likelihood_curves_DtDm(theta, cs_m, cs_o, bnames):
        texp, dm0, sigs = FitMCMC.thetaDtDm2arr(theta, len(bnames))

        res = 0.
        for i, bname in enumerate(bnames):
            lc_o = cs_o[bname]
            to, mo, eo = lc_o.Time, lc_o.Mag, lc_o.MagErr
            # model interpolation
            lcm = cs_m[bname]
            tm, mm = lcm.Time, lcm.Mag
            mm += dm0
            tm += texp
            em = sigs[i]
            m_fit = np.interp(to, tm, mm)  # One-dimensional linear interpolation

            inv_sigma2 = 1. / (eo ** 2 + em ** 2)
            chi = np.sum((mo - m_fit) ** 2 * inv_sigma2 - np.log(inv_sigma2) + np.log(2*np.pi))
            res += chi
        res = -0.5 * res
        return res

    @staticmethod
    def log_posteriorDtDm(theta, cs_m, cs_o, bnames, kwargs):
        dt_lim = kwargs.get('dt_lim', (-100, 100.))
        dm_lim = kwargs.get('dm_lim', (-3, 3.))
        siglim = kwargs.get('siglim', 3.)
        lp = FitMCMC.log_priorDtDm(theta, tlim=dt_lim, maglim=dm_lim, siglim=siglim)
        if not np.isfinite(lp):
            return -np.inf
        return lp + FitMCMC.log_likelihood_curves_DtDm(theta, cs_m, cs_o, bnames)

    @staticmethod
    def log_posteriorDt(theta, cs_m, cs_o, bnames, kwargs):
        dt_lim = kwargs.get('dt_lim', (-100, 100.))
        siglim = kwargs.get('siglim', 3.)
        lp = FitMCMC.log_priorDt(theta, nb=len(bnames), tlim=dt_lim, siglim=siglim)
        if not np.isfinite(lp):
            return -np.inf
        return lp + FitMCMC.log_likelihood_curves_Dt(theta, cs_m, cs_o, bnames)

    @staticmethod
    def chi2Dt(theta, cs_m, cs_o, bnames):
        # dt, dm, *lnf_bands = theta
        texp, sigs = FitMCMC.thetaDt2arr(theta, len(bnames))
        res = 0.
        N = 0
        for i, bname in enumerate(bnames):
            lc_o = cs_o[bname]
            to, mo, eo = lc_o.Time, lc_o.Mag, lc_o.MagErr
            # model interpolation
            lcm = cs_m[bname]
            tm, mm = lcm.Time, lcm.Mag
            tm += texp
            em = sigs[i]
            m_fit = np.interp(to, tm, mm)  # One-dimensional linear interpolation
            inv_sigma2 = 1. / (eo ** 2 + em ** 2)
            chi = np.sum((mo - m_fit) ** 2 * inv_sigma2)
            res += chi
            N += len(to)
        return res, N

    @staticmethod
    def chi2DtDm(theta, cs_m, cs_o, bnames):
        # dt, dm, *lnf_bands = theta
        texp, dm0, sigs = FitMCMC.thetaDtDm2arr(theta, len(bnames))
        res = 0.
        N = 0
        for i, bname in enumerate(bnames):
            lc_o = cs_o[bname]
            to, mo, eo = lc_o.Time, lc_o.Mag, lc_o.MagErr
            # model interpolation
            lcm = cs_m[bname]
            tm, mm = lcm.Time, lcm.Mag
            mm += dm0
            tm += texp
            em = sigs[i]
            m_fit = np.interp(to, tm, mm)  # One-dimensional linear interpolation
            inv_sigma2 = 1. / (eo ** 2 + em ** 2)
            chi = np.sum((mo - m_fit) ** 2 * inv_sigma2)
            res += chi
            N += len(to)
        return res, N

    def best_curvesDtDm(self, curves_m, curves_o, dt0, dm0,
                        threads=1, dt_lim=(-100., 100.), dm_lim=(-3., 3.), is_sampler=False):
        """
        Find the values of time shift and magnitude shift minimizing the distance between the observational
        and modal light curves
        :param curves_m: modal LCs
        :param curves_o: observational LCs
        :param dt0: initial time shift
        :param dm0: initial magnitude shift
        :param threads: threads for MCMC
        :param dt_lim: limits for  time shift for prior probability
        :param dm_lim: limits for  magnitude shift for prior probability
        :param is_sampler: if True also return samples
        :return: the dictionary with fitting results
        """
        dt_init = {lc.Band.Name: lc.tshift for lc in curves_o}
        dm_init = {lc.Band.Name: lc.mshift for lc in curves_o}

        bnames = curves_o.BandNames
        ndim = len(FitMCMC.thetaDtDmNames(bnames))  # number of parameters in the model

        if dt0 is not None or dm0 is not None:
            # starting_guesses = np.random.rand(self.nwalkers, ndim)  #
            # starting_guesses[:, 0] *= dt_lim  # for time delay in [-100,100]
            # starting_guesses[:, 1] *= dm_lim  # for  magnification

            # starting_guesses = 1. - 2. * np.random.rand(self.nwalkers, ndim)  #
            starting_guesses = np.random.rand(self.nwalkers, ndim)  #
            starting_guesses[:, 0] = dt_lim[0] + starting_guesses[:, 0] * (dt_lim[1] - dt_lim[0])  # delta t_exp
            starting_guesses[:, 1] = dm_lim[0] + starting_guesses[:, 1] * (dm_lim[1] - dm_lim[0])  # delta t_exp
            # starting_guesses[:, 1] = -dm_lim * (1. - 2.*starting_guesses[:, 1])  # dm
            # sig_lambda, model uncertainties for the each band
            # starting_guesses[:, 2:] = np.random.rand(self.nwalkers, len(bnames))

            sampler = emcee.EnsembleSampler(self.nwalkers, ndim, FitMCMC.log_posteriorDtDm,
                                            args=(curves_m, curves_o, bnames, {'dt_lim': dt_lim, 'dm_lim': dm_lim}),
                                            threads=threads, a=3)
            sampler.run_mcmc(starting_guesses, self.nsteps)

            #     sample = sampler.chain  # shape = (nwalkers, nsteps, ndim)
            # sample_dt = sampler.flatchain[self.nburn:, 0].ravel()  # discard burn-in points
            # sample_mu = sampler.flatchain[self.nburn:, 1].ravel()  # discard burn-in points
            # sample_lnf = sampler.flatchain[self.nburn:, 2].ravel()  # discard burn-in points

            # dt, dtsigma = np.mean(sample_dt), np.std(sample_dt)
            # dm, dmsigma = np.mean(sample_mu), np.std(sample_mu)
            # e
            # samples = np.array([sample_dt, sample_mu, sample_lnf]).T
            samples = sampler.chain[:, self.nburn:, :].reshape((-1, ndim))

            quartiles = list(
                map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]), zip(*np.percentile(samples, (16, 50, 84), axis=0))))

            theta_fit = [x[0] for x in quartiles]
            theta_fit_e1 = [x[1] for x in quartiles]
            theta_fit_e2 = [x[2] for x in quartiles]

            if self.is_info:
                print("Theta and it errors: rows[val, +, -]")
                FitMCMC.print_thetaDtDm(theta_fit, bnames)
                FitMCMC.print_thetaDtDm(theta_fit_e1, bnames)
                FitMCMC.print_thetaDtDm(theta_fit_e2, bnames)

            # Results
            # # # we'll start at random locations between 0 and 2000
            # # starting_guesses = np.zeros(self.nwalkers, ndim)  #
            # # starting_guesses[:, 0] = dt0 + np.random.rand(self.nwalkers) * dt_lim  # for time delay in [-100,100]
            # # starting_guesses[:, 1] = np.random.rand(self.nwalkers)
            #
            # sampler = emcee.EnsembleSampler(self.nwalkers, ndim, FitMCMC.log_posteriorDt,
            #                                 args=(curves_m, curves_o, bnames, {'dt_lim': dt_lim}),
            #                                 threads=threads, a=3)
            # sampler.run_mcmc(starting_guesses, self.nsteps)
            #
            # #     sample = sampler.chain  # shape = (nwalkers, nsteps, ndim)
            # sample_dt = sampler.flatchain[self.nburn:, 0].ravel()  # discard burn-in points
            # sample_lnf = sampler.flatchain[self.nburn:, 1].ravel()  # discard burn-in points
            #
            # # dt, dtsigma = np.mean(sample_dt), np.std(sample_dt)
            # # dm, dmsigma = np.mean(sample_mu), np.std(sample_mu)
            # samples = np.array([sample_dt, sample_lnf]).T
            # quartiles = list(
            #     map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]), zip(*np.percentile(samples, (16, 50, 84), axis=0))))
            # dt, dtsig1, dtsig2 = quartiles[0][0], quartiles[0][1], quartiles[0][2]
            # lnf, lnfsig1, lnfsig2 = quartiles[1][0], quartiles[1][1], quartiles[1][2]

            THETA = namedtuple('THETA', FitMCMC.thetaDtDmNames(bnames))
            th, e1, e2 = THETA(*theta_fit), THETA(*theta_fit_e1), THETA(*theta_fit_e2)

            texp, dm0, sigs = FitMCMC.thetaDtDm2arr(th, len(bnames))

            # Statistics
            chi2, N = FitMCMC.chi2DtDm(th, curves_m, curves_o, bnames)
            dof = N - ndim
            l = FitMCMC.log_likelihood_curves_DtDm(theta_fit, curves_m, curves_o, bnames)
            bic = ndim * np.log(N) - 2. * l
            aic = 2 * ndim - 2. * l
            measure = np.sqrt(np.sum(np.array(sigs) ** 2))

            res = th._asdict()
            res['chi2'] = chi2
            res['bic'] = bic
            res['aic'] = aic
            res['dof'] = dof
            res['measure'] = measure
            res['acceptance_fraction'] = np.mean(sampler.acceptance_fraction)

            # res = {'dt': dt, 'dtsig1': dtsig1, 'dtsig2': dtsig2, 'dtsig': (dtsig1 + dtsig2) / 2.,
            #        'lnf': lnf, 'lnfsig1': lnfsig1, 'lnfsig2': lnfsig2, 'lnfsig': (lnfsig1 + lnfsig2) / 2.,
            #        'chi2': chi2, 'bic': bic, 'aic': aic, 'dof': dof, 'measure': measure,
            #        'acceptance_fraction': np.mean(sampler.acceptance_fraction)}

            # return initial states
            for lc in curves_o:
                lc.tshift = dt_init[lc.Band.Name]
                lc.mshift = dm_init[lc.Band.Name]

            # plot a histogram of the sample
            if self.is_info:
                # the fraction of steps accepted for each walker.
                print("Mean acceptance fraction: {0:.3f}".format(res['acceptance_fraction']))

            if self.is_info:
                print(
                    "The final params are: dt={:6.2f}+{:.2f}-{:0.2f}  dm={:6.2f}+{:.2f}-{:0.2f} chi2:{:.0f}/{:d}".format(
                        th.dt, e1.dt, e2.dt, th.dm0, e1.dm0, e2.dm0, res['chi2'], res['dof']))

            if is_sampler:
                return res, (th, e1, e2), sampler  # sampler

            return res, (th, e1, e2)

        #     return res, (th, e1, e2)
        #     dt, dtsig1, dtsig2 = quartiles[0][0], quartiles[0][1], quartiles[0][2]
        #     dm, dmsig1, dmsig2 = quartiles[1][0], quartiles[1][1], quartiles[1][2]
        #     lnf, lnfsig1, lnfsig2 = quartiles[2:][0], quartiles[2:][1], quartiles[2:][2]
        #
        #     # Statistics
        #     chi2, N = FitMCMC.chi2DtDm([dt, dm, lnf], curves_m, curves_o, bnames)
        #     # chi2 = abs(log_likelihood_chi2((dt, dm), curves_m, curves_o))
        #     dof = N - ndim
        #     l = FitMCMC.log_likelihood_curves_DtDm([dt, dm, lnf], curves_m, curves_o, bnames)
        #     bic = ndim * np.log(N) + l
        #     aic = 2 * ndim + l
        #     measure = lnf
        #
        #     res = {'dt': dt, 'dtsig1': dtsig1, 'dtsig2': dtsig2, 'dtsig': (dtsig1 + dtsig2) / 2.,
        #            'dm': dm, 'dmsig1': dmsig1, 'dmsig2': dmsig2, 'dmsig': (dmsig1 + dmsig2) / 2.,
        #            'lnf': lnf, 'lnfsig1': lnfsig1, 'lnfsig2': lnfsig2, 'lnfsig': (lnfsig1 + lnfsig2) / 2.,
        #            'chi2': chi2, 'bic': bic, 'aic': aic, 'dof': dof, 'measure': measure,
        #            'acceptance_fraction': np.mean(sampler.acceptance_fraction)}
        # else:
        #     chi2, N = FitMCMC.chi2DtDm([0.]*ndim, curves_m, curves_o, bnames)  # just return the weight diff
        #     dof = N - ndim
        #     if self.is_info:
        #         print("No fit: only run least_sq:   chi2={} dof={}".format(chi2, dof))
        #
        #     res = {'dt': None, 'dtsig': None, 'dm': None, 'dmsig': None, 'chi2': chi2, 'dof': dof}
        #     return res

        # # restore the initial states
        # for lc in curves_o:
        #     lc.tshift = dt_init[lc.Band.Name]
        #     lc.mshift = dm_init[lc.Band.Name]
        #
        # # plot a histogram of the sample
        # if self.is_info:
        #     # the fraction of steps accepted for each walker.
        #     print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
        #
        # if self.is_info:
        #     print("The final params are: dt={:6.2f}+-{:6.4f} dm={:6.2f}+-{:6.4f}  chi2:{:.0f}/{:d}".format(
        #         res['dt'], res['dtsig'], res['dm'], res['dmsig'], res['chi2'], res['dof']))
        #
        # if is_samples:
        #     return res, samples
        # return res

    def best_curves(self, curves_m, curves_o, dt0, dm0=None, **kwargs):
                    # threads=1, dt_lim=(-100., 100.), dm_lim=(-5., 5.), is_sampler=False):
        """
        Find the values of time shift and magnitude shift minimizing the distance between the observational
        and modal light curves
        :param curves_m: modal LCs
        :param curves_o: observational LCs
        :param dt0: initial time shift
        :param dm0: magnitude shift, default: 0.
        :param kwargs: dictionary with
               threads: threads for MCMC
                dt_lim: limits for  time shift for prior probability
                dm_lim: limits for  magnitude shift for prior probability
                is_sampler: if True also return sampler
                is_fit_sigmas: always True
        :return: the dictionary with fitting results
        """

        threads = kwargs.get("threads", 1)
        dt_lim = kwargs.get("dt_lim", (-300., 300.))
        dm_lim = kwargs.get("dm_lim", (-5., 5.))
        is_sampler = kwargs.get("is_sampler", False)
        is_fit_sigmas = kwargs.get("is_fit_sigmas", True)

        samples = None
        if dm0 is not None:
            result = self.best_curvesDtDm(curves_m, curves_o, dt0, dm0,
                                          dt_lim=dt_lim, dm_lim=dm_lim, threads=threads, is_sampler=is_sampler)
        else:
            result = self.best_curvesDt(curves_m, curves_o, dt0,
                                        dt_lim=dt_lim, threads=threads, is_sampler=is_sampler)

        res = result[0]
        th, e1, e2 = result[1]
        if is_sampler:
            samples = result[2]

        fit_result = FitLcResult()
        fit_result.tshift = res['dt']
        fit_result.tsigma = (e1.dt + e2.dt) / 2.  # res['dtsig']
        fit_result.measure = res['measure']
        fit_result.comm = 'result MCMC: CHI2|BIC|AIC: {:.0f}|{:.0f}|{:.0f} acceptance_fraction: {:.3f}'. \
            format(res['chi2'], res['bic'], res['aic'], res['acceptance_fraction'])

        if is_sampler:
            return fit_result, res, (th, e1, e2), samples
        else:
            return fit_result, res, (th, e1, e2)

    def best_curvesDt(self, curves_m, curves_o, dt0,
                      threads=1, dt_lim=(-100., 100.), is_sampler=False):
        dt_init = {lc.Band.Name: lc.tshift for lc in curves_o}
        bnames = curves_o.BandNames

        ndim = len(FitMCMC.thetaDtNames(bnames))  # number of parameters in the model

        # starting_guesses = np.random.rand(self.nwalkers, ndim)  #
        # starting_guesses[:, 0] *= dt_lim  # for time delay in [-100,100]
        # starting_guesses[:, 1] *= dm_lim  # for  magnification

        starting_guesses = np.random.rand(self.nwalkers, ndim)
        starting_guesses[:, 0] = dt_lim[0] + starting_guesses[:, 0] * (dt_lim[1] - dt_lim[0])  # delta t_exp
        # starting_guesses[:, 0] = starting_guesses[:, 0]*dt_lim + dt0
        # sig_lambda, model uncertainties for the each band
        starting_guesses[:, 1:] = np.random.rand(self.nwalkers, len(bnames))

        # add dt0
        starting_guesses[:, 0] = starting_guesses[:, 0] + dt0

        # initial data
        curves_m.set_tshift(0.)
        curves_m.set_mshift(0.)
        # with Pool(threads) as pool:
        #     sampler = emcee.EnsembleSampler(self.nwalkers, ndim, FitMCMC.log_posteriorDt,
        #                                     args=(curves_m, curves_o, bnames, {'dt_lim': dt_lim}),
        #                                     pool=pool)
        sampler = emcee.EnsembleSampler(self.nwalkers, ndim, FitMCMC.log_posteriorDt,
                                        args=(curves_m, curves_o, bnames, {'dt_lim': dt_lim}),
                                        threads=threads, a=2)
        sampler.run_mcmc(starting_guesses, self.nsteps)
        samples = sampler.chain[:, self.nburn:, :].reshape((-1, ndim))

        quartiles = list(
            map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]), zip(*np.percentile(samples, (16, 50, 84), axis=0))))
        theta_fit = [x[0] for x in quartiles]
        theta_fit_e1 = [x[1] for x in quartiles]
        theta_fit_e2 = [x[2] for x in quartiles]

        if self.is_info:
            print("Theta and it errors: rows[val, +, -]")
            FitMCMC.print_thetaDt(theta_fit, bnames)
            FitMCMC.print_thetaDt(theta_fit_e1, bnames)
            FitMCMC.print_thetaDt(theta_fit_e2, bnames)

        # Results
        # # # we'll start at random locations between 0 and 2000
        # # starting_guesses = np.zeros(self.nwalkers, ndim)  #
        # # starting_guesses[:, 0] = dt0 + np.random.rand(self.nwalkers) * dt_lim  # for time delay in [-100,100]
        # # starting_guesses[:, 1] = np.random.rand(self.nwalkers)
        #
        # sampler = emcee.EnsembleSampler(self.nwalkers, ndim, FitMCMC.log_posteriorDt,
        #                                 args=(curves_m, curves_o, bnames, {'dt_lim': dt_lim}),
        #                                 threads=threads, a=3)
        # sampler.run_mcmc(starting_guesses, self.nsteps)
        #
        # #     sample = sampler.chain  # shape = (nwalkers, nsteps, ndim)
        # sample_dt = sampler.flatchain[self.nburn:, 0].ravel()  # discard burn-in points
        # sample_lnf = sampler.flatchain[self.nburn:, 1].ravel()  # discard burn-in points
        #
        # # dt, dtsigma = np.mean(sample_dt), np.std(sample_dt)
        # # dm, dmsigma = np.mean(sample_mu), np.std(sample_mu)
        # samples = np.array([sample_dt, sample_lnf]).T
        # quartiles = list(
        #     map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]), zip(*np.percentile(samples, (16, 50, 84), axis=0))))
        # dt, dtsig1, dtsig2 = quartiles[0][0], quartiles[0][1], quartiles[0][2]
        # lnf, lnfsig1, lnfsig2 = quartiles[1][0], quartiles[1][1], quartiles[1][2]

        THETA = namedtuple('THETA', FitMCMC.thetaDtNames(bnames))
        th, e1, e2 = THETA(*theta_fit), THETA(*theta_fit_e1), THETA(*theta_fit_e2)

        texp, sigs = FitMCMC.thetaDt2arr(th, len(bnames))

        # Statistics
        chi2, N = FitMCMC.chi2Dt(th, curves_m, curves_o, bnames)
        dof = N - ndim
        l = FitMCMC.log_likelihood_curves_Dt(theta_fit, curves_m, curves_o, bnames)
        bic = ndim * np.log(N) - 2. * l
        aic = 2 * ndim - 2. * l
        measure = np.sum(np.array(sigs) ** 2)

        res = th._asdict()
        res['chi2'] = chi2
        res['bic'] = bic
        res['aic'] = aic
        res['dof'] = dof
        res['measure'] = measure
        res['acceptance_fraction'] = np.mean(sampler.acceptance_fraction)

        # res = {'dt': dt, 'dtsig1': dtsig1, 'dtsig2': dtsig2, 'dtsig': (dtsig1 + dtsig2) / 2.,
        #        'lnf': lnf, 'lnfsig1': lnfsig1, 'lnfsig2': lnfsig2, 'lnfsig': (lnfsig1 + lnfsig2) / 2.,
        #        'chi2': chi2, 'bic': bic, 'aic': aic, 'dof': dof, 'measure': measure,
        #        'acceptance_fraction': np.mean(sampler.acceptance_fraction)}

        # return initial states
        for lc in curves_o:
            lc.tshift = dt_init[lc.Band.Name]

        # plot a histogram of the sample
        if self.is_info:
            # the fraction of steps accepted for each walker.
            print("Mean acceptance fraction: {0:.3f}".format(res['acceptance_fraction']))

        if self.is_info:
            print("The final params are: dt={:6.2f}+{:.2f}-{:0.2f} chi2:{:.0f}/{:d}".format(
                th.dt, e1.dt, e2.dt, res['chi2'], res['dof']))

        if is_sampler:
            return res, (th, e1, e2), sampler  # sampler

        return res, (th, e1, e2)

        # if is_samples:
        #     return res, samples
        # return res
    #
    # def best_curves_sigmas(self, curves_m, curves_o, bnames, dt0, dm0=None,
    #                        threads=1, dt_lim=100., dm_lim=5., is_sampler=False):
    #     """
    #     Find the values of time shift and magnitude shift minimizing the distance between the observational
    #     and modal light curves
    #     :param curves_m: modal LCs
    #     :param curves_o: observational LCs
    #     :param bnames: bands to fit
    #     :param dt0: initial time shift
    #     :param dm0: initial magnitude shift
    #     :param threads: threads for MCMC
    #     :param dt_lim: limits for  time shift for prior probability
    #     :param dm_lim: limits for  magnitude shift for prior probability
    #     :param is_sampler: if True also return sampler
    #     :return: the dictionary with fitting results
    #     """
    #     if dm0 is not None:
    #         return self.best_curves_sigmasDtDm(curves_m, curves_o, bnames, dt0, dm0,
    #                                            threads=threads, dt_lim=dt_lim, dm_lim=dm_lim,
    #                                            is_samples=is_sampler)
    #     else:
    #         return self.best_curves_sigmasDt(curves_m, curves_o, bnames, dt0,
    #                                          threads=threads, dt_lim=dt_lim, is_sampler=is_sampler)
    #
    # def best_curves_sigmasDt(self, curves_m, curves_o, bnames, dt0,
    #                          threads=1, dt_lim=100., is_sampler=False):
    #     dt_init = {lc.Band.Name: lc.tshift for lc in curves_o}
    #
    #     ndim = 1 + len(bnames)  # number of parameters in the model
    #
    #     # we'll start at random locations between 0 and 2000
    #     starting_guesses = np.random.rand(self.nwalkers, ndim)
    #     starting_guesses[:, 0] *= dt_lim  # for time delay in [-100,100]
    #     starting_guesses[:, 0] += dt0  # for time delay in [-100,100]
    #     starting_guesses[:, 1:] = 1. - 2. * starting_guesses[:, 1:]  # for lnf
    #
    #     sampler = emcee.EnsembleSampler(self.nwalkers, ndim, FitMCMC.log_posteriorDt,
    #                                     args=(curves_m, curves_o, bnames, {'dt_lim': dt_lim}),
    #                                     threads=threads, a=3)
    #     sampler.run_mcmc(starting_guesses, self.nsteps)
    #
    #     #     sample = sampler.chain  # shape = (nwalkers, nsteps, ndim)
    #     sample_dt = sampler.flatchain[self.nburn:, 0].ravel()  # discard burn-in points
    #     sample_lnf = sampler.flatchain[self.nburn:, 1:].T  # discard burn-in points
    #
    #     # dt, dtsigma = np.mean(sample_dt), np.std(sample_dt)
    #     # dm, dmsigma = np.mean(sample_mu), np.std(sample_mu)
    #     samples = np.vstack((sample_dt, sample_lnf)).T
    #     # samples = samples.T
    #     # samples = np.array([sample_dt, sample_mu, sample_lnf]).T
    #     quartiles = list(
    #         map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]), zip(*np.percentile(samples, (16, 50, 84), axis=0))))
    #     dt, dtsig1, dtsig2 = quartiles[0][0], quartiles[0][1], quartiles[0][2]
    #     arr_lnf = np.array(quartiles[1:])
    #     lnf, lnfsig1, lnfsig2 = map(np.ravel, (arr_lnf[:, 0], arr_lnf[:, 1], arr_lnf[:, 2]))
    #
    #     # Statistics
    #     chi2, N = FitMCMC.chi2(np.concatenate(([dt, 0.], lnf)), curves_m, curves_o, bnames)
    #     # chi2 = abs(log_likelihood_chi2((dt, dm), curves_m, curves_o))
    #     dof = N - ndim
    #     l = FitMCMC.log_likelihood_curves(np.concatenate(([dt, 0.], lnf)), curves_m, curves_o, bnames)
    #     bic = ndim * np.log(N) - 2. * l
    #     aic = 2 * ndim - 2. * l
    #     measure = lnf
    #
    #     res = {'dt': dt, 'dtsig1': dtsig1, 'dtsig2': dtsig2, 'dtsig': (dtsig1 + dtsig2) / 2.,
    #            'lnf': lnf, 'lnfsig1': lnfsig1, 'lnfsig2': lnfsig2, 'lnfsig': (lnfsig1 + lnfsig2) / 2.,
    #            'signames': bnames,
    #            'chi2': chi2, 'bic': bic, 'aic': aic, 'dof': dof, 'measure': measure,
    #            'acceptance_fraction': np.mean(sampler.acceptance_fraction)
    #            }
    #
    #     # return initial states
    #     for lc in curves_o:
    #         lc.tshift = dt_init[lc.Band.Name]
    #
    #     # plot a histogram of the sample
    #     if self.is_info:
    #         # the fraction of steps accepted for each walker.
    #         print("Mean acceptance fraction: {0:.3f}".format(res['acceptance_fraction']))
    #
    #     if self.is_info:
    #         print("The final params are: dt={:6.2f}+-{:6.4f} chi2:{:.0f}/{:d}".format(
    #             res['dt'], res['dtsig'], res['chi2'], res['dof']))
    #
    #     del sample_dt
    #     del sample_lnf
    #
    #     if is_sampler:
    #         return res, sampler
    #     else:
    #         del samples
    #         return res
    #
    # def best_curves_sigmasDtDm(self, curves_m, curves_o, bnames, dt0, dm0,
    #                            threads=1, dt_lim=100., dm_lim=5., is_samples=False):
    #     """
    #     Find the values of time shift and magnitude shift minimizing the distance between the observational
    #     and modal light curves
    #     :param curves_m: modal LCs
    #     :param curves_o: observational LCs
    #     :param dt0: initial time shift
    #     :param dm0: initial magnitude shift
    #     :param threads: threads for MCMC
    #     :param dt_lim: limits for  time shift for prior probability
    #     :param dm_lim: limits for  magnitude shift for prior probability
    #     :param is_samples: if True also return samples
    #     :return: the dictionary with fitting results
    #     """
    #     dt_init = {lc.Band.Name: lc.tshift for lc in curves_o}
    #     dm_init = {lc.Band.Name: lc.mshift for lc in curves_o}
    #
    #     ndim = 2 + len(bnames)  # number of parameters in the model
    #
    #     starting_guesses = np.random.rand(self.nwalkers, ndim)
    #     starting_guesses[:, 0] *= dt_lim  # for time delay in [-100,100]
    #     starting_guesses[:, 1] *= dm_lim  # for  magnification
    #     starting_guesses[:, 2:] = 1. - 2. * starting_guesses[:, 2:]  # for lnf
    #
    #     sampler = emcee.EnsembleSampler(self.nwalkers, ndim, FitMCMC.log_posteriorDtDm,
    #                                     args=(curves_m, curves_o, bnames, {'dt_lim': dt_lim, 'dm_lim': dm_lim}),
    #                                     threads=threads, a=3)
    #     sampler.run_mcmc(starting_guesses, self.nsteps)
    #
    #     #     sample = sampler.chain  # shape = (nwalkers, nsteps, ndim)
    #     sample_dt = sampler.flatchain[self.nburn:, 0].ravel()  # discard burn-in points
    #     sample_mu = sampler.flatchain[self.nburn:, 1].ravel()  # discard burn-in points
    #     # sample_lnf = sampler.flatchain[self.nburn:, 2:].ravel()  # discard burn-in points
    #     sample_lnf = sampler.flatchain[self.nburn:, 2:].T  # discard burn-in points
    #
    #     # dt, dtsigma = np.mean(sample_dt), np.std(sample_dt)
    #     # dm, dmsigma = np.mean(sample_mu), np.std(sample_mu)
    #     samples = np.vstack((sample_dt, sample_mu, sample_lnf)).T
    #     # samples = samples.T
    #     # samples = np.array([sample_dt, sample_mu, sample_lnf]).T
    #     quartiles = list(
    #         map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]), zip(*np.percentile(samples, (16, 50, 84), axis=0))))
    #     dt, dtsig1, dtsig2 = quartiles[0][0], quartiles[0][1], quartiles[0][2]
    #     dm, dmsig1, dmsig2 = quartiles[1][0], quartiles[1][1], quartiles[1][2]
    #     arr_lnf = np.array(quartiles[2:])
    #     lnf, lnfsig1, lnfsig2 = map(np.ravel, (arr_lnf[:, 0], arr_lnf[:, 1], arr_lnf[:, 2]))
    #
    #     # Statistics
    #     chi2, N = FitMCMC.chi2(np.concatenate(([dt, dm], lnf)), curves_m, curves_o, bnames)
    #     # chi2 = abs(log_likelihood_chi2((dt, dm), curves_m, curves_o))
    #     dof = N - ndim
    #     l = FitMCMC.log_likelihood_curves(np.concatenate(([dt, dm], lnf)), curves_m, curves_o, bnames)
    #     bic = ndim * np.log(N) + l
    #     aic = 2 * ndim + l
    #     measure = np.sum([lnf[i] * curves_o[bname].Length for i, bname in enumerate(bnames)]) / N
    #
    #     res = {'dt': dt, 'dtsig1': dtsig1, 'dtsig2': dtsig2, 'dtsig': (dtsig1 + dtsig2) / 2.,
    #            'dm': dm, 'dmsig1': dmsig1, 'dmsig2': dmsig2, 'dmsig': (dmsig1 + dmsig2) / 2.,
    #            'lnf': lnf, 'lnfsig1': lnfsig1, 'lnfsig2': lnfsig2, 'lnfsig': (lnfsig1 + lnfsig2) / 2.,
    #            'signames': bnames,
    #            'chi2': chi2, 'bic': bic, 'aic': aic, 'dof': dof, 'measure': measure,
    #            'acceptance_fraction': np.mean(sampler.acceptance_fraction)
    #            }
    #     # else:
    #     #     chi2, N = FitMCMC.chi2([0., 0.], curves_m, curves_o)  # just return the weight diff
    #     #     dof = N - ndim
    #     #     if self.is_info:
    #     #         print("No fit: only run least_sq:   chi2={} dof={}".format(chi2, dof))
    #     #
    #     #     res = {'dt': None, 'dtsig': None, 'dm': None, 'dmsig': None, 'chi2': chi2, 'dof': dof}
    #     #     return res
    #
    #     # return initial states
    #     for lc in curves_o:
    #         lc.tshift = dt_init[lc.Band.Name]
    #         lc.mshift = dm_init[lc.Band.Name]
    #
    #     # plot a histogram of the sample
    #     if self.is_info:
    #         # the fraction of steps accepted for each walker.
    #         print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
    #
    #     if self.is_info:
    #         print("The final params are: dt={:6.2f}+-{:6.4f} dm={:6.2f}+-{:6.4f}  chi2:{:.0f}/{:d}".format(
    #             res['dt'], res['dtsig'], res['dm'], res['dmsig'], res['chi2'], res['dof']))
    #
    #     del sample_dt
    #     del sample_mu
    #     del sample_lnf
    #
    #     if is_samples:
    #         return res, samples
    #     else:
    #         del samples
    #         return res

    # @staticmethod
    @staticmethod
    def plot_corner(samples, labels=("dt", 'dm'), bnames=('',), axhist=None, **kwargs):
        """
        To plot  a data set in a multi-dimensional space.

        Parameters
        ----------
        :param samples: The samples. This should be a 1- or 2-dimensional array. For a 1-D
        array this results in a simple histogram. For a 2-D array, the zeroth
        axis is the list of samples and the next axis are the dimensions of
        the space (see corner.corner).
        :param labels: the labels  Default: ("td", 'dm')
        :param bnames: the band names. Used for error labels
        :param bins
        :param axhist: axes to plot hist2d
        :return: fig
        """
        # from matplotlib import pyplot as plt
        import corner
        # plt.hist(sample_dt, bins=50, histtype="stepfilled", alpha=0.3)
        # # samples = sampler.chain[:, nburn:, :].reshape((-1, ndim))
        # fig = corner.corner(samples, labels=["$delay$", "$dm$", "$lnf$"],
        #                     truths=[dt, dm, lnf])
        lw = kwargs.get('lw', 2)
        bins = kwargs.get('bins', 30)
        alpha = kwargs.get('alpha', 1.)
        istex = kwargs.get('istex', False)
        title_fmt = kwargs.get('title_fmt', '.2f')
        quantiles = kwargs.get('quantiles', [0.16, 0.5, 0.84])
        figsize = kwargs.get('figsize', (16, 12))
        verbose = kwargs.get('verbose', False)
        
        if isinstance(labels, str):
            labels = [labels]
            
        if istex:
            labels += tuple(r'$\sigma_m(\mathrm{{{}}})$'.format(bn) for bn in bnames)
        else:            
            labels += tuple('sig+{}'.format(bn) for bn in bnames)
        hist_kwargs = dict(lw=lw, alpha=0.5)
        fig = corner.corner(samples, labels=labels, bins=bins, hist_kwargs=hist_kwargs,
                            quantiles=quantiles,
                            show_titles=True, title_fmt=title_fmt, verbose=verbose,
                            title_kwargs={"fontsize": 12}, label_kwargs={'labelpad': 1, 'fontsize': 14},
                            fontsize=8)

        fig.set_size_inches(figsize)

        if axhist is not None:
            corner.hist2d(samples[:, 0], samples[:, 1], ax=axhist)

        return fig
