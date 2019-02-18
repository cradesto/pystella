import numpy as np

from sklearn import gaussian_process
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels \
    import Matern, RBF, WhiteKernel,ConstantKernel, ExpSineSquared, RationalQuadratic

from pystella.rf.lc import LightCurve


class FitGP:
    """
    #  Gaussian process
    """

    @staticmethod
    def fit_lc(lc, Ntime=None, is_RBF=False):
        t = lc.Time
        if is_RBF:
            gp = FitGP.lc2gpRBF(lc)
        else:
            gp = FitGP.lc2gp(lc)

        if Ntime is not None:  # new time points
            t_new = np.linspace(min(t), max(t), Ntime)
            X_samples = t_new.reshape(-1, 1)  # np.array(t_new, ndmin = 2).T
        else:
            t_new = t
            X_samples = t.reshape(-1, 1)

        # Make the prediction on the meshed x-axis (ask for MSE as well)
        y_pred, sigma = gp.predict(X_samples, return_std=True)
        return LightCurve(lc.Band, t_new, y_pred, sigma), gp

    # def fit_lc(lc, Ntime=None):
    #     t = lc.Time
    #     y = lc.Mag
    #     y = np.asarray(y, dtype=np.float64)
    #     yerr = lc.MagErr
    #
    #     kernel = ConstantKernel() + Matern(length_scale=2, nu=3 / 2) + WhiteKernel(noise_level=1)
    #     alpha = yerr ** 2
    #     gp = gaussian_process.GaussianProcessRegressor(kernel=kernel, alpha=alpha)
    #     X_obs = t.reshape(-1, 1)
    #     gp.fit(X_obs, y)
    #
    #     if Ntime is not None:  # new time points
    #         t_new = np.linspace(min(t), max(t), Ntime)
    #         X_samples = t_new.reshape(-1, 1)  # np.array(t_new, ndmin = 2).T
    #     else:
    #         t_new = t
    #         X_samples = X_obs
    #
    #     # Make the prediction on the meshed x-axis (ask for MSE as well)
    #     y_pred, sigma = gp.predict(X_samples, return_std=True)
    #     return LightCurve(lc.Band, t_new, y_pred, sigma), gp

    @staticmethod
    def lc2gp(lc):
        t = lc.Time
        y = lc.Mag
        y = np.asarray(y, dtype=np.float64)
        yerr = lc.MagErr

        time_scale = t[-1] - t[0]
        data_scale = np.max(y) - np.min(y)
        noise_std = np.median(yerr)

        length_scale = 0.01 * time_scale

        kernel = ConstantKernel(0.1) \
                 + Matern(length_scale=length_scale, nu=3 / 2) \
                 + WhiteKernel(noise_level=noise_std**2)
        alpha = (yerr / y) ** 2  # yerr ** 2
        gp = gaussian_process.GaussianProcessRegressor(kernel=kernel, alpha=alpha)
        X_obs = t.reshape(-1, 1)
        gp.fit(X_obs, y)
        return gp

    @staticmethod
    def lc2gpRBF(lc, long_term_length_scale=None, short_term_length_scale=None,
                 noise_level=None):
        """
        See https://github.com/ipashchenko/ogle/blob/master/lc.py
        """

        t = lc.Time
        y = lc.Mag
        y = np.asarray(y, dtype=np.float64)
        yerr = lc.MagErr

        # data = self.data[['mjd', 'mag', 'err']]
        # data = np.atleast_2d(data)
        # time = data[:, 0] - data[0, 0]
        # time = np.atleast_2d(time).T
        #
        time_scale = t[-1] - t[0]
        data_scale = np.max(y) - np.min(y)
        noise_std = np.median(yerr)

        if long_term_length_scale is None:
            long_term_length_scale = 0.5 * time_scale

        if noise_level is None:
            noise_level = noise_std

        # k1 = data_scale ** 2 * RBF(length_scale=long_term_length_scale)
        # k2 = 0.1 * data_scale * \
        #      RBF(length_scale=pre_periodic_term_length_scale) * \
        #      ExpSineSquared(length_scale=periodic_term_length_scale,
        #                     periodicity=periodicity)
        # k3 = WhiteKernel(noise_level=noise_level ** 2,
        #                  noise_level_bounds=(1e-3, 1.))
        # kernel = k1 + k2 + k3
        # gp = GaussianProcessRegressor(kernel=kernel,
        #                               alpha=(yerr / y) ** 2,
        #                               normalize_y=True,
        #                               n_restarts_optimizer=10)

        if long_term_length_scale is None:
            long_term_length_scale = 0.5 * time_scale

        if short_term_length_scale is None:
            short_term_length_scale = 0.05 * time_scale

        if noise_level is None:
            noise_level = noise_std

        k1 = data_scale ** 2 * \
             RationalQuadratic(length_scale=long_term_length_scale)
        k2 = 0.1 * data_scale * RBF(length_scale=short_term_length_scale)
        k3 = WhiteKernel(noise_level=noise_level ** 2,
                         noise_level_bounds=(1e-3, np.inf))
        kernel = k1 + k2 + k3
        gp = GaussianProcessRegressor(kernel=kernel,
                                      alpha=(yerr / y) ** 2,
                                      normalize_y=True)

        X_obs = t.reshape(-1, 1)
        gp.fit(X_obs, y)
        return gp
