import emcee
# import corner
import numpy as np
from scipy import interpolate

# from matplotlib import pyplot as plt


def fit_lc_bayesian_1d(lc_o, lc_m, is_debug=True, is_plot=True):
    def log_prior(theta):
        return 1  # flat prior

    def log_likelihood_t(theta, lc_o, lc_m):
        t = lc_o.Time
        mo = lc_o.Mag
        e = lc_o.MagErr
        lc_m.tshift = theta[0]
        tck = interpolate.splrep(lc_m.Time, lc_m.Mag, s=0)
        m_mdl = interpolate.splev(t, tck, der=0)
        return -0.5 * np.sum(np.log(2 * np.pi * (e ** 2))
                             + (mo - m_mdl) ** 2 / e ** 2)

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


def fit_curves_bayesian_2d(curves_o, curves_m, is_debug=True, is_info=False):
    def log_prior(theta):
        return 1  # flat prior

    def log_likelihood_lc(theta, lc_o, lc_m):
        to, mo, eo = lc_o.Time, lc_o.Mag, lc_o.MagErr
        lc_m.tshift = theta[0]
        lc_m.mshift = theta[1]
        tck = interpolate.splrep(lc_m.Time, lc_m.Mag, s=0)
        m_mdl = interpolate.splev(to, tck, der=0)
        return -0.5 * np.sum(np.log(2 * np.pi * (eo ** 2))
                             + (mo - m_mdl) ** 2 / eo ** 2)

    def log_likelihood_curves(theta, curves_o, curves_m):
        bnames = curves_o.BandNames
        res = 0.
        for bname in bnames:
            lc_o = curves_o.get(bname)
            lc_m = curves_m.get(bname)
            if lc_m is not None:
                p = log_likelihood_lc(theta, lc_o, lc_m)
                res += p
        return res

    def log_posterior(theta, curves_o, curves_m):
        return log_prior(theta) + log_likelihood_curves(theta, curves_o, curves_m)

    ndim = 2  # number of parameters in the model
    nwalkers = 50  # number of MCMC walkers
    if is_debug:
        nburn = 100  # "burn-in" period to let chains stabilize
        nsteps = 200  # number of MCMC steps to take
    else:  # run
        nburn = 1000  # "burn-in" period to let chains stabilize
        nsteps = 2000  # number of MCMC steps to take
    # we'll start at random locations between 0 and 2000
    starting_guesses = np.random.rand(nwalkers, ndim)
    starting_guesses[0] *= 100  # for time delay in [0,100]

    # fit
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(curves_o, curves_m))
    sampler.run_mcmc(starting_guesses, nsteps)

    #     sample = sampler.chain  # shape = (nwalkers, nsteps, ndim)
    sample_dt = sampler.chain[:, nburn:, 0].ravel()  # discard burn-in points
    sample_mu = sampler.chain[:, nburn:, 1].ravel()  # discard burn-in points
    #     sample = sampler.chain[:, nburn:, :]  # discard burn-in points

    tshift, tsigma = np.mean(sample_dt), np.std(sample_dt)
    dm, dmsigma = np.mean(sample_mu), np.std(sample_mu)

    res = {'dt': (tshift, tsigma), 'dm': (dm, dmsigma)}
    # plot a histogram of the sample
    # if is_info:
    #     # the fraction of steps accepted for each walker.
    #     print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
    #     plt.hist(sample_dt, bins=50, histtype="stepfilled", alpha=0.3, normed=True)
    #     samples = sampler.chain[:, nburn:, :].reshape((-1, ndim))
    #     # fig = corner.corner(samples, labels=["$delay$", "$dm$"],
        #                     truths=[tshift, dm])

    return res

