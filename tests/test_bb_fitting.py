import numpy as np
import unittest
from scipy.optimize import curve_fit

import pylab as plt

import pystella.rf.rad_func as rf
import pystella.rf.spectrum as spectrum

__author__ = 'bakl'


class TestSpectrumFitting(unittest.TestCase):
    def setUp(self):
        nf = 100
        start, end = 10., 1e5
        wl = np.exp(np.linspace(np.log(start), np.log(end), nf))
        nu = rf.val_to_hz(wl, inp="A")
        T = 5500
        s = spectrum.SpectrumPlanck(nu, T)
        self.sp = s

    def test_fit_t_color(self):
        Tcol = self.sp.T_color
        self.assertAlmostEqual(Tcol, self.sp.T,
                               msg="For planck  Tcolor [%f] should be equal Tbb [%f]." % (Tcol, self.sp.T))

    def test_t_wien(self):
        Twien = self.sp.T_wien
        self.assertAlmostEqual(Twien, self.sp.T,
                               msg="For planck  Twien [%f] should be equal Tbb [%f]." % (Twien, self.sp.T)
                               , delta=Twien*0.05)  # ACCURACY 5 %

    def plot(self):
        # plot the input model and synthetic data
        nu = self.sp._freq
        flux = self.sp.flux_q
        Tcol = self.sp.T_color
        ybest = rf.planck(nu, Tcol)
        # plot the solution
        plt.plot(nu, flux, 'b*', nu, ybest, 'r-', label='Tinit=%f, Tcol=%f'%(self.sp.Tbb, Tcol))
        plt.xscale('log')
        plt.yscale('log')
        plt.show()

    def test_fit_bb(self):
        def func(nu, T):
            return np.pi * rf.planck(nu, T, inp="Hz", out="freq")

        self.sp.cut_flux(max(self.sp.Flux) * 1e-5)

        freq = self.sp.Freq
        flux = self.sp.Flux
        Tinit = 1.e4

        popt, pcov = curve_fit(func, freq, flux, p0=Tinit)
        Tbest = popt
        # bestT, pcov = curve_fit(rf.fit_planck(nu, inp='Hz'), nu, flux, p0=Tinit, sigma=sigma)
        sigmaT = np.sqrt(np.diag(pcov))

        print 'True model values'
        print '  Tbb = %.2f K' % self.sp.T

        print 'Parameters of best-fitting model:'
        print '  T = %.2f +/- %.2f K' % (Tbest, sigmaT)

        # Tcol = self.sp.temp_color
        ybest = func(freq, Tbest) #np.pi * rf.planck(freq, Tbest, inp="Hz", out="freq")
        # plot the solution
        plt.plot(freq, flux, 'b*', label='Spectral T: %f' % self.sp.T)
        plt.plot(freq, ybest, 'r-', label='Best Tcol: %f' % Tbest)
        plt.xscale('log')
        plt.yscale('log')
        plt.legend(loc=3)
        plt.show()

        self.assertAlmostEqual(Tbest, self.sp.T,
                               msg="For planck  Tcolor [%f] should be equal sp.T [%f]." % (Tbest, self.sp.T),
                               delta=Tbest*0.01)
