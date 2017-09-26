import unittest
from os.path import dirname, abspath, join

import matplotlib.pyplot as plt

import scipy.optimize as op
import numpy as np
import emcee

from pystella.rf import light_curve_func as lcf
from plugin import sn1999em


class EmceeTests(unittest.TestCase):
    @unittest.skip("just for plot")
    def plot_simple(self):
        def lnprob(x, ivar):
            return -0.5 * np.sum(ivar * x ** 2)

        ndim, nwalkers = 10, 100
        ivar = 1. / np.random.rand(ndim)
        p0 = [np.random.rand(ndim) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[ivar])
        sampler.run_mcmc(p0, 1000)
        for i in range(ndim):
            plt.figure()
            plt.hist(sampler.flatchain[:, i], 100, color="k", histtype="step")
            plt.title("Dimension {0:d}".format(i))
        plt.show()

    def test_fit_linear_model(self):
        """Example: Fitting a Linear Model to Data
        See http://dan.iel.fm/emcee/current/user/line/"""
        # Choose the "true" parameters.
        m_true = -0.9594
        b_true = 4.294
        f_true = 0.534

        # Generate some synthetic data from the model.
        N = 50
        x = np.sort(10 * np.random.rand(N))
        yerr = 0.1 + 0.5 * np.random.rand(N)
        y = m_true * x + b_true
        y += np.abs(f_true * y) * np.random.randn(N)
        y += yerr * np.random.randn(N)

        def lnlike(theta, x, y, yerr):
            m, b, lnf = theta
            model = m * x + b
            inv_sigma2 = 1.0 / (yerr ** 2 + model ** 2 * np.exp(2 * lnf))
            return -0.5 * (np.sum((y - model) ** 2 * inv_sigma2 - np.log(inv_sigma2)))

        nll = lambda *args: -lnlike(*args)
        result = op.minimize(nll, [m_true, b_true, np.log(f_true)], args=(x, y, yerr))
        m_ml, b_ml, lnf_ml = result["x"]

        def lnprior(theta):
            m, b, lnf = theta
            if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < lnf < 1.0:
                return 0.0
            return -np.inf

        def lnprob(theta, x, y, yerr):
            lp = lnprior(theta)
            if not np.isfinite(lp):
                return -np.inf
            return lp + lnlike(theta, x, y, yerr)

        ndim, nwalkers = 3, 100
        pos = [result["x"] + 1e-4 * np.random.randn(ndim) for i in range(nwalkers)]

        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
        sampler.run_mcmc(pos, 500)

        samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

        xl = np.array([0, 10])
        for m, b, lnf in samples[np.random.randint(len(samples), size=100)]:
            plt.plot(xl, m * xl + b, color="k", alpha=0.1)
        plt.plot(xl, m_true * xl + b_true, color="r", lw=2, alpha=0.8)
        plt.errorbar(x, y, yerr=yerr, fmt=".k")
        plt.show()

        samples[:, 2] = np.exp(samples[:, 2])
        m_mcmc, b_mcmc, f_mcmc = map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                     zip(*np.percentile(samples, [16, 50, 84], axis=0)))

        def print_v3(v3):
            print("v = %f + %f - %f" % v3)

        map(print_v3, (m_mcmc, b_mcmc, f_mcmc))

    def test_fit_stella_model(self):
        """Example: Fitting a Stella model to Sn1999em"""

        # obs data
        dm = -29.38
        curves = sn1999em.read_curves()
        lc_obs = curves.get('V')
        lc_obs.mshift = -dm
        lc_obs.tshift = -lc_obs.tmin

        # model light curves
        name = 'cat_R1000_M15_Ni007_E15'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        curves_model = lcf.curves_compute(name, path, bands='V')
        lc_mdl = curves_model.get('V')


        todo
        # Choose the "true" parameters.
        dt_init = 0.  # initial time shift
        dm_init = 0.  # initial time shift
        f_true = 0.534

        # Generate some synthetic data from the model.
        # N = 50
        t = lc_obs.Time
        merr = lc_obs.Err
        m = lc_obs.Mag

        def lnlike(theta, t, m, yerr):
            dt, dm, lnf = theta
            model = lc_mdl.Mag
            inv_sigma2 = 1.0 / (yerr ** 2 + model ** 2 * np.exp(2 * lnf))
            return -0.5 * (np.sum((lc_obs - model) ** 2 * inv_sigma2 - np.log(inv_sigma2)))

        nll = lambda *args: -lnlike(*args)
        result = op.minimize(nll, [dt_init, dm_init, np.log(f_true)], args=(t, m, merr))
        m_ml, b_ml, lnf_ml = result["x"]

        def lnprior(theta):
            m, b, lnf = theta
            if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < lnf < 1.0:
                return 0.0
            return -np.inf

        def lnprob(theta, x, y, yerr):
            lp = lnprior(theta)
            if not np.isfinite(lp):
                return -np.inf
            return lp + lnlike(theta, x, y, yerr)

        ndim, nwalkers = 3, 100
        pos = [result["x"] + 1e-4 * np.random.randn(ndim) for i in range(nwalkers)]

        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(t, m, merr))
        sampler.run_mcmc(pos, 500)

        samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

        xl = np.array([0, 10])
        for m, b, lnf in samples[np.random.randint(len(samples), size=100)]:
            plt.plot(xl, m * xl + b, color="k", alpha=0.1)
        plt.plot(xl, dt_init * xl + dm_init, color="r", lw=2, alpha=0.8)
        plt.errorbar(t, m, yerr=merr, fmt=".k")
        plt.show()

        samples[:, 2] = np.exp(samples[:, 2])
        m_mcmc, b_mcmc, f_mcmc = map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                     zip(*np.percentile(samples, [16, 50, 84], axis=0)))

        def print_v3(v3):
            print("v = %f + %f - %f" % v3)

        map(print_v3, (m_mcmc, b_mcmc, f_mcmc))


def main():
    pass


if __name__ == '__main__':
    main()
