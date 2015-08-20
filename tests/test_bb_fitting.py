import unittest
import numpy as np
import pylab as plt
from scipy.optimize import curve_fit

import pystella.rf.spectrum as spectrum
import pystella.util.rf as rf

__author__ = 'bakl'


class TestSpectrumFitting(unittest.TestCase):
    def setUp(self):
        nf = 100
        start, end = 10, 1e5
        wl = np.exp(np.linspace(np.log(start), np.log(end), nf))
        nu = rf.val_to_hz(wl, inp="A")
        T = 5500
        self.sp = spectrum.SpectrumPlanck(nu, T)

    def test_fit_t_color(self):
        Tcol = self.sp.fit_t_color()
        self.assertAlmostEqual(Tcol, self.sp.Tbb,
                               msg="For planck  Tcolor [%f] should be equal Tbb [%f]." % (Tcol, self.sp.Tbb))

    def test_t_wien(self):
        Twien = self.sp.temperature_wien
        self.assertAlmostEqual(Twien, self.sp.Tbb,
                               msg="For planck  Twien [%f] should be equal Tbb [%f]." % (Twien, self.sp.Tbb)
                               , delta=0.1)

    def plot(self):
        # plot the input model and synthetic data
        nu = self.sp.freq
        flux = self.sp.flux_q
        Tcol = self.sp.fit_t_color()
        ybest = rf.planck(nu, Tcol)
        # plot the solution
        plt.plot(nu, flux, 'b*', nu, ybest, 'r-', label='Tinit=%f, Tcol=%f'%(self.sp.Tbb, Tcol))
        plt.xscale('log')
        plt.yscale('log')
        plt.show()

    def test_fit_bb(self):
        nu = self.sp.freq
        flux = self.sp.flux_q
        Tinit = 1.e4

        def func(nu, T):
            return rf.planck(nu, T, inp="Hz", out="freq")

        popt, pcov = curve_fit(func, nu, flux, p0=Tinit)
        bestT = popt
        # bestT, pcov = curve_fit(rf.fit_planck(nu, inp='Hz'), nu, flux, p0=Tinit, sigma=sigma)
        sigmaT = np.sqrt(np.diag(pcov))

        print 'True model values'
        print '  Tbb = %.2f' % self.sp.Tbb

        print 'Parameters of best-fitting model:'
        print '  T = %.2f +/- %.2f' % (bestT, sigmaT)

