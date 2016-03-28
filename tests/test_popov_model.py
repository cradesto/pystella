import numpy as np
import unittest
import matplotlib.pyplot as plt


from pystella.model.popov import Popov
from plugin import sn87a
from plugin import sn1999em

__author__ = 'bakl'


class TestPopovModel(unittest.TestCase):
    def test_plot_lc(self):
        n = 100
        start, end = 0.1, 200.
        time = np.exp(np.linspace(np.log(start), np.log(end), n))

        popov = Popov('test', R=400., M=15., Mni=0.04, E=1.e51)
        popov.plot_Lbol(time)
        plt.show()

    def test_popov_SN87A(self):
        n = 100
        start, end = 0.1, 200.
        jd_shift = -2446850.  # moment of explosion SN 1987A, Hamuy 1988, doi:10.1086/114613
        dm = -18.6  # D = 7.5e6 pc
        time = np.exp(np.linspace(np.log(start), np.log(end), n))

        popov = Popov('test', R=50., M=15., Mni=0.07, E=1.e51)
        ax = popov.plot_Lbol(time)
        sn87a.plot_ubv(ax, path=sn87a.sn_path, jd_shift=jd_shift, mshift=dm)
        plt.show()

    def test_popov_SN1999em(self):
        n = 100
        start, end = 0.1, 200.
        jd_shift = 20.
        dm = -29.38  # D = 7.5e6 pc
        # dm = -30.4  # D = 12.e6 pc
        time = np.exp(np.linspace(np.log(start), np.log(end), n))
        popov = Popov('test', R=450., M=15., Mni=0.07, E=0.7e51)
        ax = popov.plot_Lbol(time)
        sn1999em.plot_ubv(ax, path=sn1999em.sn_path, jd_shift=jd_shift, mshift=dm)
        plt.show()
