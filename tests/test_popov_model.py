import numpy as np
import unittest
import matplotlib.pyplot as plt

from pystella.fit.fit_mcmc import FitLcMcmc
from pystella.model.popov import Popov
import pystella.rf.light_curve_plot as lcp
from plugin import sn87a
from plugin import sn1999em
from plugin import rednova

__author__ = 'bakl'


class TestPopovModel(unittest.TestCase):
    # @unittest.skip("just for plot")
    def test_plot_lc(self):
        n = 100
        start, end = 0.1, 200.
        time = np.exp(np.linspace(np.log(start), np.log(end), n))

        popov = Popov('test', R=400., M=15., Mni=0.04, E=1.)
        popov.plot_Lbol(time)
        plt.show()

    # @unittest.skip("just for plot")
    def test_popov_SN87A(self):
        n = 100
        start, end = 0.1, 1200.
        jd_shift = -2446850.  # moment of explosion SN 1987A, Hamuy 1988, doi:10.1086/114613
        dm = -18.6  # D = 7.5e6 pc
        time = np.exp(np.linspace(np.log(start), np.log(end), n))

        popov = Popov('test', R=50., M=15., Mni=0.07, E=1.)
        ax = popov.plot_Lbol(time)
        sn87a.plot_ubv(ax, path=sn87a.sn_path, jd_shift=jd_shift, mshift=dm)
        plt.show()

    # @unittest.skip("just for plot")
    def test_popov_SN1999em(self):
        n = 100
        start, end = 0.1, 200.
        jd_shift = 20.
        dm = -29.38  # D = 7.5e6 pc
        # dm = -30.4  # D = 12.e6 pc
        time = np.exp(np.linspace(np.log(start), np.log(end), n))
        popov = Popov('test', R=450., M=15., Mni=0.04, E=0.7)
        ax = popov.plot_Lbol(time)
        sn1999em.plot_ubv(ax, path=sn1999em.sn_path, jd_shift=jd_shift, mshift=dm)
        plt.show()

    # @unittest.skip("just for plot")
    def test_popov_rednova(self):
        n = 100
        start, end = 0.1, 90.
        jd_shift = -2500000+42963.
        dm = -24.38  # D = 7.5e5 pc
        # dm = -30.4  # D = 12.e6 pc
        time = np.exp(np.linspace(np.log(start), np.log(end), n))
        popov = Popov('test', R=35., M=1., Mni=0., E=0.002)
        ax = popov.plot_Lbol(time)
        rednova.plot_ubv(ax, path=rednova.sn_path, jd_shift=jd_shift, mshift=dm)
        plt.show()

    # @unittest.skip("just for plot")
    def test_popov_SN1999em_emcee(self):
        n = 100
        start, end = 0.1, 200.
        jd_shift = 20.
        dm = -29.38  # D = 7.5e6 pc
        # dm = -30.4  # D = 12.e6 pc
        time = np.exp(np.linspace(np.log(start), np.log(end), n))
        popov = Popov('test', R=450., M=15., Mni=0.04, E=0.7)
        lc_m = popov.LCBol(time)

        # ax = popov.plot_Lbol(time)
        # sn1999em.plot_ubv(ax, path=sn1999em.sn_path, jd_shift=jd_shift, mshift=dm)
        # plt.show()
        curves_o = sn1999em.read_curves()
        lc_o = curves_o.get('V')
        lc_o.mshift = dm
        print("Run: find tshift with bayesian: obs band %s with %s ..." % (lc_o.Band.Name, popov))
        fitter = FitLcMcmc()
        fitter.is_debug = True
        res = fitter.fit_lc(lc_o, lc_m)
        tshift, tsigma = res.tshift, res.tsigma
        print("Result: tshift= %s tsigma= %s ..." % (tshift, tsigma))

        ax = popov.plot_Lbol(time)
        lcp.lc_plot(lc_m, ax)
        lcp.lc_plot(lc_o, ax, is_line=False)
        # sn1999em.plot_ubv(ax, path=sn1999em.sn_path, jd_shift=-tshift, mshift=dm)
        plt.show()
