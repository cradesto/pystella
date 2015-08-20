from pystella.rf.star import Star

__author__ = 'bakl'
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.optimize as opt
import pystella.util.rf as rf
from pystella.util.phys_var import phys


class Spectrum:
    def __init__(self, name, freq=None, flux=None, is_sort_wl=True):
        """Creates a Spectrum instance.  Required parameters:  name."""
        self.name = name
        self.freq = np.array(freq, copy=True)  # frequencies of flux [cm]
        self.flux = np.array(flux, copy=True)  # flux [erg / (cm^2*Hz) ]
        # if is_sort_wl:
        #    freq, flux = [list(x) for x in zip(*sorted(zip(freq, flux), key=itemgetter(0)))]
        if is_sort_wl:
            sorti = np.argsort(self.Wl)
            self.freq = self.freq[sorti]
            self.flux = self.flux[sorti]

    @property
    def Flux_wl(self):  # wavelength of flux [cm]
        return self.flux * self.freq ** 2 / phys.c  # flux [erg/cm^2/cm) ]

    @property
    def Wl(self):
        return phys.c / self.freq

    def plot_spec(self, title=''):
        plt.title('Spectrum: '+title)
        plt.plot(self.Wl*phys.cm_to_angs,  self.Flux_wl)
        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel('Flux')
        plt.xlabel('Wave [A]')
        plt.grid()
        plt.show()

    def fit_t_color(self):
        """
        Fitting Spectrum by planck function and find the color temperature
        :return:   color temperature
        """
        nu = self.freq
        # Tinit = 1.e4
        Tinit = self.temperature_wien

        def func(freq, T):
            return rf.planck(freq, T, inp="Hz", out="freq")
        popt, pcov = opt.curve_fit(func, nu, self.flux, p0=Tinit)
        return popt

    @property
    def temperature_wien(self):
        """
        Find temperature of Spectrum as Wien's law
        :return: temperature
        """
        b = 0.28977721  # [cm]
        wl = self.Wl
        tck = interpolate.splrep(wl, self.Flux_wl)

        def func(wl):
            return -interpolate.splev(wl, tck, der=0)
        wl_max = opt.fmin(func, x0=wl[int(len(wl)/2)])
        Twien = b / wl_max
        return Twien

    @staticmethod
    def flux_to_distance(flux, dl):
        """
        Compute flux at distance dl
        :param flux:  flux
        :param dl: luminosity distance [parsec]
        :return: :raise ValueError:
        """
        return flux / (4*np.pi * dl**2)


class SpectrumPlanck(Spectrum):
    """
        Planck Spectrum
    """

    def __init__(self, freq, temperature, name='bb'):
        flux = math.pi * rf.planck(freq, temperature=temperature)
        Spectrum.__init__(self, name, freq=freq, flux=flux)
        self.T = temperature
        self.zeta = 0.

    @property
    def Tbb(self):
        return self.T

    def correct_zeta(self, zeta):
        self.zeta = zeta
        self.flux *= zeta**2


class SeriesSpectrum:
    def __init__(self, name):
        """Creates a Series of Spectrum instance."""
        self.name = name
        self.times = None
        self.nfreq = None
        self.wl = None  # waves
        self.freq = None
        self.data = None  # array where index -> Spectrum at the times[index]

    def set_times(self, times):
        if times is None or len(times) == 0:
            raise ValueError("times must be array with len > 0.")
        self.times = times

    def is_time(self):
        return self.times is not None and len(self.times) > 0

    def set_freq(self, freqs):
        if freqs is None or len(freqs) == 0:
            raise ValueError("freqs must be array with len > 0.")
        if np.all(freqs == 0):
            raise ValueError("No freqs item equals 0.")
        self.freq = freqs
        self.wl = phys.c / self.freq
        self.nfreq = len(self.freq)

    def set_data(self, sdata):
        if sdata is None or len(sdata) == 0:
            raise ValueError("data must be array with len > 0.")
        self.data = sdata

    def flux_to_mags(self, band, z=0, dl=0):
        if band is None:
            return None
        if not self.is_time():
            return None

        mags = np.zeros(len(self.data))
        for k in range(len(self.data)):
            spec = self.data[k]
            star = Star(k, spec)
            star.set_distance(dl)
            star.set_redshift(z)
            mag = star.flux_to_mag(band)
            mags[k] = mag
            # if self.times[k] > 50:
            #     spec.plot_spec(title="t=%f, mag=%f" % (self.times[k], mag))

        return mags
