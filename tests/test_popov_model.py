import numpy as np
import unittest
import matplotlib.pyplot as plt


from pystella.model.popov import Popov
from pystella.rf.lc import SetLightCurve, LightCurve
from pystella.rf.rad_func import Lum2MagBol
from pystella.util.phys_var import phys

__author__ = 'bakl'


def lc_create(band, m=-19, dt=0.):
    n = 10
    time = np.linspace(0.+dt, 200.+dt, n)
    mags = m * np.ones(n)
    return LightCurve(band, time, mags)


class TestPopovModel(unittest.TestCase):
    def test_plot_lc(self):
        is_plot_Lum = False
        n = 100
        start, end = 0.1, 200.
        # time = np.linspace(start, end, n)
        time = np.exp(np.linspace(np.log(start), np.log(end), n))

        popov = Popov('test', R=900., M=15., Mni=0.04, E=1.e51)

        # lum
        if is_plot_Lum:
            mags = popov.Lbol(time*86400)
            mags_ni = popov.e_rate_ni(time*86400)
        else:  # mag
            mags = popov.MagBol(time*86400)
            mags_ni = Lum2MagBol(popov.e_rate_ni(time*86400))

        plt.plot(time, mags, color='blue', label='L bol')
        plt.plot(time, mags_ni, color='red', ls='-.', label='Eni rate')
        if not is_plot_Lum:
            plt.gca().invert_yaxis()
        plt.xlabel('Time [day]')
        plt.ylabel('Absolute Magnitude')
        plt.legend()
        plt.show()


