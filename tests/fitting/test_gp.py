import numpy as np
import matplotlib.pyplot as plt
import unittest

import gptools

from plugin import sn1999em, snrefsdal
from pystella.rf import band

__author__ = 'bakl'


class TestFitGaussianProcess(unittest.TestCase):
    def test_GP_example(self):
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

        np.random.seed(1)

        def f(x):
            """The function to predict."""
            return x * np.sin(x)

        # ----------------------------------------------------------------------
        #  First the noiseless case
        X = np.atleast_2d([1., 3., 5., 6., 7., 8.]).T

        # Observations
        y = f(X).ravel()

        # Mesh the input space for evaluations of the real function, the prediction and
        # its MSE
        x = np.atleast_2d(np.linspace(0, 10, 1000)).T

        # Instanciate a Gaussian Process model
        kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
        gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)

        # Fit to data using Maximum Likelihood Estimation of the parameters
        gp.fit(X, y)

        # Make the prediction on the meshed x-axis (ask for MSE as well)
        y_pred, sigma = gp.predict(x, return_std=True)

        # Plot the function, the prediction and the 95% confidence interval based on
        # the MSE
        fig = plt.figure()
        plt.plot(x, f(x), 'r:', label=u'$f(x) = x\,\sin(x)$')
        plt.plot(X, y, 'r.', markersize=10, label=u'Observations')
        plt.plot(x, y_pred, 'b-', label=u'Prediction')
        plt.fill(np.concatenate([x, x[::-1]]),
                 np.concatenate([y_pred - 1.9600 * sigma,
                                 (y_pred + 1.9600 * sigma)[::-1]]),
                 alpha=.5, fc='b', ec='None', label='95% confidence interval')
        plt.xlabel('$x$')
        plt.ylabel('$f(x)$')
        plt.ylim(-10, 20)
        plt.legend(loc='upper left')

        # ----------------------------------------------------------------------
        # now the noisy case
        X = np.linspace(0.1, 9.9, 20)
        X = np.atleast_2d(X).T

        # Observations and noise
        y = f(X).ravel()
        dy = 0.5 + 1.0 * np.random.random(y.shape)
        noise = np.random.normal(0, dy)
        y += noise

        # Instanciate a Gaussian Process model
        gp = GaussianProcessRegressor(kernel=kernel, alpha=(dy / y) ** 2,
                                      n_restarts_optimizer=10)

        # Fit to data using Maximum Likelihood Estimation of the parameters
        gp.fit(X, y)

        # Make the prediction on the meshed x-axis (ask for MSE as well)
        y_pred, sigma = gp.predict(x, return_std=True)

        # Plot the function, the prediction and the 95% confidence interval based on
        # the MSE
        plt.plot(x, f(x), 'r:', label=u'$f(x) = x\,\sin(x)$')
        plt.errorbar(X.ravel(), y, dy, fmt='r.', markersize=10, label=u'Observations')
        plt.plot(x, y_pred, 'b-', label=u'Prediction')
        plt.fill(np.concatenate([x, x[::-1]]),
                 np.concatenate([y_pred - 1.9600 * sigma,
                                 (y_pred + 1.9600 * sigma)[::-1]]),
                 alpha=.5, fc='b', ec='None', label='95% confidence interval')
        plt.xlabel('$x$')
        plt.ylabel('$f(x)$')
        plt.ylim(-10, 20)
        plt.legend(loc='upper left')

        plt.show()

    def test_fit_GP_SN1999em(self):
        # add data
        dm = -29.38  # D = 7.5e6 pc
        # dm = -30.4  # D = 12.e6 pc
        curves = sn1999em.read_curves()
        lc = curves.get('V')
        lc.mshift = dm
        t = lc.Time
        y = lc.Mag
        yerr = lc.MagErr

        # k = gptools.SquaredExponentialKernel()
        # gp = gptools.GaussianProcess(k)
        # k = gptools.SquaredExponentialKernel(param_bounds=[(0, 1e3), (0, 100)])
        k = gptools.SquaredExponentialKernel(param_bounds=[(0., max(np.abs(y))),
                                                           (0, np.std(t))])
        gp = gptools.GaussianProcess(k, mu=gptools.LinearMeanFunction())

        gp.add_data(t, y, err_y=yerr)

        is_mcmc = True
        is_mcmc = False
        if is_mcmc:
            out = gp.predict(t, use_MCMC=True, full_MCMC=True,
                             return_std=False,
                             num_proc=0,
                             nsamp=200,
                             plot_posterior=True,
                             plot_chains=False,
                             burn=100,
                             thin=1)

        else:
            gp.optimize_hyperparameters(verbose=True)
            out = gp.predict(t, use_MCMC=False)

        y_star, err_y_star = out

        # gp.optimize_hyperparameters()
        # y_star, err_y_star = gp.predict(t)

        fig = plt.figure()
        ax = fig.add_axes((0.1, 0.3, 0.8, 0.65))
        ax.invert_yaxis()

        ax.plot(t, y, color='blue', label='L bol', lw=2.5)
        ax.errorbar(t, y, yerr=yerr, fmt='o', color='blue', label='%s obs.')

        #
        # ax.plot(t, y_star, color='red', ls='--', lw=1.5, label='GP')
        # third: plot a constrained function with errors
        ax.plot(t, y_star, '-', color='gray')
        ax.fill_between(t, y_star - 2 * err_y_star, y_star + 2 * err_y_star, color='gray', alpha=0.3)
        # ax.errorbar(t, y_star, err_y_star, fmt='.k', ms=6)

        plt.legend()
        plt.show()

    def test_fit_GP_SNRefsdal(self):
        # k = gptools.SquaredExponentialKernel(param_bounds=[(0, 1e3), (0, 100)])
        # gp = gptools.GaussianProcess(k, mu=gptools.LinearMeanFunction())

        # add data
        dm = -29.38  # D = 7.5e6 pc
        # dm = -30.4  # D = 12.e6 pc
        image = "S1"
        bname = 'F160W'
        curves = snrefsdal.read_curves(snrefsdal.path_data, image)
        lc = curves.get(bname)
        # lc.mshift = dm
        t = lc.Time
        y = lc.Mag
        yerr = lc.MagErr

        #  Gaussian process
        k = gptools.SquaredExponentialKernel(param_bounds=[(0, max(np.abs(y))),
                                                           (0, np.std(t))])
        # k = gptools.SquaredExponentialKernel(param_bounds=[(min(np.abs(y)), max(np.abs(y))),
        #                                                    (0, np.std(t))])
        gp = gptools.GaussianProcess(k)
        # gp = gptools.GaussianProcess(k, mu=gptools.LinearMeanFunction())
        gp.add_data(t, y, err_y=yerr)

        gp.optimize_hyperparameters()
        y_star, err_y_star = gp.predict(t)

        fig = plt.figure()
        ax = fig.add_axes((0.1, 0.3, 0.8, 0.65))
        ax.invert_yaxis()

        ax.plot(t, y, color='blue', label='L bol', lw=2.5)
        ax.errorbar(t, y, yerr=yerr, fmt='o', color='blue', label='%s obs.')

        #
        # ax.plot(t, y_star, color='red', ls='--', lw=1.5, label='GP')
        ax.plot(t, y_star, '-', color='gray')
        ax.fill_between(t, y_star - 2 * err_y_star, y_star + 2 * err_y_star, color='gray', alpha=0.3)

        plt.show()

    def test_fit_GP_SNRefsdal_all_lc(self):
        # k = gptools.SquaredExponentialKernel(param_bounds=[(0, 1e3), (0, 100)])
        # gp = gptools.GaussianProcess(k, mu=gptools.LinearMeanFunction())

        # add data
        dm = -29.38  # D = 7.5e6 pc
        # dm = -30.4  # D = 12.e6 pc
        image = "S1"
        bands = ['F160W', 'F105W', 'F125W']
        # bands = ('F160W','F140W','F105W', 'F125W')
        curves = snrefsdal.read_curves(snrefsdal.path_data, image)
        for bname in bands:
            lc = curves.get(bname)
            # lc.mshift = dm
            t = lc.Time
            y = lc.Mag
            yerr = lc.MagErr

            #  Gaussian process
            k = gptools.SquaredExponentialKernel(param_bounds=[(min(np.abs(y)), max(np.abs(y))),
                                                               (0, np.std(t))])
            # k = gptools.SquaredExponentialKernel(param_bounds=[(min(np.abs(y)), max(np.abs(y))),
            #                                                    (0, np.std(t))])
            gp = gptools.GaussianProcess(k)
            # gp = gptools.GaussianProcess(k, mu=gptools.LinearMeanFunction())
            gp.add_data(t, y, err_y=yerr)

            is_mcmc = True
            if is_mcmc:
                out = gp.predict(t, use_MCMC=True, full_MCMC=True, return_std=True,
                                 num_proc=0, nsamp=100,
                                 plot_posterior=True, plot_chains=False,
                                 burn=10, thin=1)

            else:
                gp.optimize_hyperparameters()
                out = gp.predict(t, use_MCMC=False)
            y_star, err_y_star = out
            # y_star, err_y_star = gp.predict(t)

            fig = plt.figure()
            ax = fig.add_axes((0.1, 0.3, 0.8, 0.65))
            ax.invert_yaxis()

            bcolor = band.bands_colors()[bname]
            ax.plot(t, y, color=bcolor, label='L bol', lw=2.5)
            ax.errorbar(t, y, yerr=yerr, fmt='o', color=bcolor, label='%s obs.')

            #
            # ax.plot(t, y_star, color='red', ls='--', lw=1.5, label='GP')
            ax.plot(t, y_star, '-', color=bcolor)
            ax.fill_between(t, y_star - 2 * err_y_star, y_star + 2 * err_y_star, color=bcolor, alpha=0.3)

        plt.show()

    def test_scikit_GP_SNRefsdal(self):
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

        # add data
        dm = -29.38  # D = 7.5e6 pc
        # dm = -30.4  # D = 12.e6 pc
        image = "S1"
        bname = 'F160W'
        curves = snrefsdal.read_curves(snrefsdal.path_data, image)
        lc = curves.get(bname)
        # lc.mshift = dm
        t = lc.Time
        y = lc.Mag
        yerr = lc.MagErr
        #
        # Instanciate a Gaussian Process model
        # kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
        # Instanciate a Gaussian Process model
        # kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
        # gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)
        gp = GaussianProcessRegressor()

        # Fit to data using Maximum Likelihood Estimation of the parameters
        X = np.atleast_2d(t).T
        gp.fit(X, y)

        # gp = GaussianProcessRegressor()

        # Fit to data using Maximum Likelihood Estimation of the parameters
        # gp.fit(t, y)

        # Make the prediction on the meshed x-axis (ask for MSE as well)
        # y_star, err_y_star = gp.predict(t, return_std=True)
        # Make the prediction on the meshed x-axis (ask for MSE as well)
        y_pred, sigma = gp.predict(t, return_std=True)

        # k = gptools.SquaredExponentialKernel(param_bounds=[(min(np.abs(y)), max(np.abs(y))),
        #                                                    (0, np.std(t))])
        # # k = gptools.SquaredExponentialKernel(param_bounds=[(min(np.abs(y)), max(np.abs(y))),
        # #                                                    (0, np.std(t))])
        # gp = gptools.GaussianProcess(k)
        # # gp = gptools.GaussianProcess(k, mu=gptools.LinearMeanFunction())
        # gp.add_data(t, y, err_y=yerr)
        #
        # gp.optimize_hyperparameters()
        # y_star, err_y_star = gp.predict(t)

        fig = plt.figure()
        ax = fig.add_axes((0.1, 0.3, 0.8, 0.65))
        ax.invert_yaxis()

        ax.plot(t, y, color='blue', label='L bol', lw=2.5)
        ax.errorbar(t, y, yerr=yerr, fmt='o', color='blue', label='%s obs.')

        #
        # ax.plot(t, y_star, color='red', ls='--', lw=1.5, label='GP')
        ax.plot(t, y_pred, '-', color='gray')
        # ax.fill_between(t, y_star - 2 * err_y_star, y_star + 2 * err_y_star, color='gray', alpha=0.3)
        ax.fill(np.concatenate([t, t[::-1]]),
                np.concatenate([y_pred - 1.9600 * sigma,
                                (y_pred + 1.9600 * sigma)[::-1]]),
                alpha=.5, fc='b', ec='None', label='95% confidence interval')

        plt.show()

    def test_GP_brownian_motion(self):
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

        # add data
        t = np.linspace(0, 10, 100)
        #
        # Instanciate a Gaussian Process model
        # kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
        # Instanciate a Gaussian Process model
        kernel = lambda x, y: 1. * min(x, y)
        # kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
        gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)
        # gp = GaussianProcessRegressor()

        # Fit to data using Maximum Likelihood Estimation of the parameters
        X = np.atleast_2d(t).T
        gp.fit(X, y)

        # gp = GaussianProcessRegressor()

        # Fit to data using Maximum Likelihood Estimation of the parameters
        # gp.fit(t, y)

        # Make the prediction on the meshed x-axis (ask for MSE as well)
        # y_star, err_y_star = gp.predict(t, return_std=True)
        # Make the prediction on the meshed x-axis (ask for MSE as well)
        y_pred, sigma = gp.predict(t, return_std=True)

        fig = plt.figure()
        ax = fig.add_axes((0.1, 0.3, 0.8, 0.65))
        ax.invert_yaxis()

        ax.plot(t, y, color='blue', label='L bol', lw=2.5)
        ax.errorbar(t, y, yerr=yerr, fmt='o', color='blue', label='%s obs.')

        #
        # ax.plot(t, y_star, color='red', ls='--', lw=1.5, label='GP')
        ax.plot(t, y_pred, '-', color='gray')
        # ax.fill_between(t, y_star - 2 * err_y_star, y_star + 2 * err_y_star, color='gray', alpha=0.3)
        ax.fill(np.concatenate([t, t[::-1]]),
                np.concatenate([y_pred - 1.9600 * sigma,
                                (y_pred + 1.9600 * sigma)[::-1]]),
                alpha=.5, fc='b', ec='None', label='95% confidence interval')

        plt.show()

    def test_error_func(self):
        def rel_errors(mu, sig, func, num=100):
            x_norm = []
            for x, s in zip(mu, sig):
                x_norm.append(np.random.normal(x, s, num))
                #     x_norm = np.random.normal(mu,sig, num)
            f_dist = func(x_norm)
            return np.mean(f_dist), np.std(f_dist)

        x1, sigX1 = (20.05, 0.20)
        x2, sigX2 = (23.12, 0.11)

        # dS1 = np.random.normal(x1, sigX1, 1000)
        # dS2 = np.random.normal(x2, sigX2, 1000)
        fmu, fsig = rel_errors([x2, x1], [sigX2, sigX1], lambda arg: arg[0] / arg[1], num=10000)
        print 'fmu= %f fsig = %f' % (fmu, fsig)