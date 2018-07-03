import numpy as np

from sklearn import gaussian_process
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel

from pystella.rf.lc import LightCurve


class FitGP:
    """
    #  Gaussian process
    """
    @staticmethod
    def fit_lc(lc, Ntime=None):
        t = lc.Time
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

        kernel = ConstantKernel() + Matern(length_scale=2, nu=3 / 2) + WhiteKernel(noise_level=1)
        alpha = yerr ** 2
        gp = gaussian_process.GaussianProcessRegressor(kernel=kernel, alpha=alpha)
        X_obs = t.reshape(-1, 1)
        gp.fit(X_obs, y)
        return gp
