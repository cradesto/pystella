from math import log10
import scipy
from pystella.util.phys_var import phys
from operator import itemgetter
from scipy.optimize import curve_fit

__author__ = 'bakl'
from scipy import interpolate
from scipy.integrate import simps as integralfunc
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
import pystella.util.rf as rf

class Spectrum:
    def __init__(self, name, freq=None, flux=None, is_sort_wl=True):
        """Creates a Spectrum instance.  Required parameters:  name."""
        self.name = name
        freq = np.array(freq, copy=True)
        flux = np.array(flux, copy=True)
        # if is_sort_wl:
        #    freq, flux = [list(x) for x in zip(*sorted(zip(freq, flux), key=itemgetter(0)))]
        wl = phys.c / freq
        if is_sort_wl:
            sorti = np.argsort(wl)
            wl = wl[sorti]
            freq = freq[sorti]
            flux = flux[sorti]

        self.wl = wl  # wavelength of flux [cm]
        self.freq = freq  # frequencies of flux [cm]
        self.flux_q = flux  # flux [erg / (cm^2*Hz) ]

    @property
    def flux_wl(self):
        return self.flux_q * phys.c / self.wl ** 2  # flux [erg/cm^2/cm) ]

    def convolution_band(self, band):
        tck = interpolate.splrep(self.wl, self.flux_wl, s=0)
        flux_intrp = interpolate.splev(band.wl, tck, der=0)

        # for i in range(1, len(band.wl)):
        #     print "band.wl=%8g  flux_intrp=%g " %(band.wl[i], flux_intrp[i])
        # for i in range(1, len(self.wl)):
        #     if min(band.wl) < self.wl[i]*phys.cm_to_angs < max(band.wl):
        #         print "wl=%8g    flux_wl=%g" %(self.wl[i],  self.flux_wl[i])
        # a = np.trapz(band.resp * flux_intrp , freq_b)
        a = integralfunc(band.resp * flux_intrp, band.wl)
        # b = np.trapz(band.resp , freq_b)
        b = integralfunc(band.resp, band.wl)
        return a / b

    def response(self, band, z=0, dl=0, is_b_spline=True):
        """
        Compute response flux using provided spectral band
        :param band:  photometric band
        :param z:  redshift
        :param dl: luminosity distance [parsec], default
        :param is_b_spline:  the method of interpolation
        :return: :raise ValueError:
        """
        if min(self.wl) > band.wl[0] or max(self.wl) < band.wl[-1]:
            raise ValueError("Spectrum must be wider then band: " + str(band))

        flux = self.flux_wl
        if dl > 0.:
            flux = Spectrum.flux_to_distance_lum(flux, dl)  # move to lum distance
        flux = flux / phys.cm_to_angs  # to flux [erg/cm^2/A) ]
        wl_s = self.wl * phys.cm_to_angs
        wl_b = band.wl * phys.cm_to_angs

        if z > 0:
            wl_s /= 1. + z  # redshift the flux
            flux *= 1. + z

        if is_b_spline:
            tck = interpolate.splrep(wl_s, flux, s=0)
            flux_intrp = interpolate.splev(wl_b, tck, der=0)
        else:
            flux_intrp = np.interp(wl_b, wl_s, flux, 0, 0)  # One-dimensional linear interpolation.

        a = integralfunc(flux_intrp * band.resp * wl_b, wl_b) / (phys.c * phys.cm_to_angs) / phys.h
        return a

    def flux_to_mag(self, band, z=0, dl=0):
        conv = self.response(band, z, dl=dl)
        if conv <= 0:
            return None
        else:
            mag = -2.5 * np.log10(conv) + band.zp
            return mag

    def flux_to_magAB(self, band, z=0, dl=0):
        conv = self.response(band, z=z, dl=dl)
        if conv <= 0:
            return None
        else:
            mag = -2.5 * np.log10(conv) + phys.ZP_AB
            return mag

    def flux_to_AB(self):
        magAB = -2.5 * np.log10(self.flux_q) + phys.ZP_AB
        return magAB

    def k_cor(self, band_r, band_o, z):
        """
        Compute K-correction for observed and rest-frame bands.

       Args:
          band_r: Rest-frame band.
          band_o: Observed band.
          z:     redshift

       Returns:
          * K: K-correction
          * If failed return None
       """
        # todo make k-correction with b-splinesec
        resp_0 = self.response(band_r, is_b_spline=False)
        resp_z = self.response(band_o, z=z, is_b_spline=False)

        if resp_0 < 0 or resp_z <= 0:
            return None
        else:
            kcor = -2.5*np.log10(resp_z / resp_0 / (1 + z)) + band_r.zp - band_o.zp
            return kcor

    def plot_spec(self, title=''):
        plt.title('Spectrum: '+title)
        plt.plot(self.wl*phys.cm_to_angs,  self.flux_wl)
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
        b = .28977721  # [cm]
        # Tinit = 1.e4
        Tinit = self.t_wien()

        def func(freq, T):
            return rf.planck(freq, T, inp="Hz", out="freq")
        popt, pcov = curve_fit(func, nu, self.flux_q, p0=Tinit)
        return popt

    def t_wien(self):
        """
        Find temperature of Spectrum as Wien's law
        :return: temperature
        """
        b = .28977721  # [cm]
        tck = interpolate.splrep(self.wl, self.flux_wl)
        def func(wl):
            return -interpolate.splev(wl, tck, der=0)
        wl_max = opt.fmin(func, x0=self.wl[len(self.wl)/2])
        Twien = b / wl_max
        return Twien

    @staticmethod
    def flux_to_distance_lum(flux, dl):
        """
        Compute flux at distance dl
        :param flux:  flux
        :param dl: luminosity distance [parsec]
        :return: :raise ValueError:
        """
        return flux / (4*np.pi * rf.pc_to_cm(dl)**2)


class SpectrumPlanck(Spectrum):
    """
        Planck Spectrum
    """

    def __init__(self, freq, temperature, name='bb'):
        flux = rf.planck(freq, temperature=temperature)
        Spectrum.__init__(self, name, freq=freq, flux=flux)
        self.T = temperature

    @property
    def Tbb(self):
        return self.T


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
            mag = spec.flux_to_mag(band, z=z, dl=dl)
            mags[k] = mag
            # if self.times[k] > 50:
            #     spec.plot_spec(title="t=%f, mag=%f" % (self.times[k], mag))

        return mags
