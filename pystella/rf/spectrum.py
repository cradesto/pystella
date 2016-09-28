import math
import numpy as np
try:
    from scipy import interpolate, integrate
    from scipy import ndimage
    from scipy.optimize import fmin, curve_fit
except:
    print("WARNING: pystella: failed to import scipy modules.")


import matplotlib.pyplot as plt

import pystella.rf.rad_func as rf
from pystella.rf import band
from pystella.rf.lc import LightCurve, SetLightCurve
from pystella.rf.star import Star
from pystella.util.phys_var import phys

__author__ = 'bakl'


class Spectrum(object):
    def __init__(self, name, freq, flux, is_sort_wl=True):
        """Creates a Spectrum instance.  Required parameters:  name."""
        self._name = name
        self._freq = np.array(freq, copy=True)  # frequencies of flux [cm]
        self._flux = np.array(flux, copy=True)  # flux [erg / (cm^2*Hz) ]
        # if is_sort_wl:
        #    freq, flux = [list(x) for x in zip(*sorted(zip(freq, flux), key=itemgetter(0)))]
        if is_sort_wl:
            sorti = np.argsort(self.Wl)
            self._freq = self._freq[sorti]
            self._flux = self._flux[sorti]

    @property
    def Name(self):
        return self._name

    @property
    def Freq(self):
        return self._freq

    @property
    def Flux(self):
        return self._flux

    @property
    def Flux_wl(self):  # wavelength of flux [cm]
        return self.Flux * self.Freq**2 / phys.c  # flux [erg/cm^2/cm) ]

    @property
    def Wl(self):
        return phys.c / self.Freq

    def plot_spec(self):
        plt.title('Spectrum: ' % self.Name)
        plt.plot(self.Wl * phys.cm_to_angs, self.Flux_wl)
        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel('Flux')
        plt.xlabel('Wave [A]')
        plt.grid()
        plt.show()

    @property
    def T_color(self):
        """
        Fitting Spectrum by planck function and find the color temperature
        :return:   color temperature
        """
        def func(nu, T):
            return np.pi * rf.planck(nu, T, inp="Hz", out="freq")

        Tinit = self.T_wien
        popt, pcov = curve_fit(func, self.Freq, self.Flux, p0=Tinit)
        return popt[0]

    @property
    def temp_wien_interp(self):
        return rf.temp_wien(self.wl_flux_max_interp)


    @property
    def wl_flux_max_interp(self):
        """
        Find temperature of Spectrum as Wien's law by interpolation
        :return: Wien's temperature
        """
        wl = self.Wl
        tck = interpolate.splrep(wl, self.Flux_wl)

        def func(w):
            return -interpolate.splev(w, tck, der=0)

        wl_max = fmin(func, x0=wl[int(len(wl) / 2)])
        return wl_max

    @property
    def wl_flux_max(self):
        """
        Find max wave length of Spectrum
        :return: the wave length of the spectrum maximum
        """
        idx = np.argmax(self.Flux_wl)
        wl_max = self.Wl[idx]
        return wl_max

    @property
    def T_wien(self):
        """
        Find temperature of Spectrum as Wien's law
        :return: Wien's temperature
        """
        return rf.temp_wien(self.wl_flux_max)

    @property
    def freq_mean(self):
        """
        Find mean frequencies of Spectrum
        :return: mean frequencies
        """
        mean = self.compute_flux_nu_bol() / self.compute_flux_bol()
        return mean

    def compute_flux_nu_bol(self):
        Hnu = integrate.simps(self._flux * self._freq, self._freq)
        return abs(Hnu)  # due to order freq

    def compute_flux_bol(self):
        H = integrate.simps(self.Flux, self.Freq)
        return abs(H)  # due to order freq

    def cut_flux(self, bottom):
        cut = self._flux >= bottom
        self._freq = self._freq[cut]
        self._flux = self._flux[cut]

    def integrate_wl(self, wavelengthMin='min', wavelengthMax='max'):
        """Calculates flux in SED within given wavelength range.

        @type wavelengthMin: float or 'min'
        @param wavelengthMin: minimum of the wavelength range
        @type wavelengthMax: float or 'max'
        @param wavelengthMax: maximum of the wavelength range
        @rtype: float
        @return: relative flux

        """

        if wavelengthMin == 'min':
            wavelengthMin = self.Wl.min()
        if wavelengthMax == 'max':
            wavelengthMax = self.Wl.max()

        mask = np.logical_and(np.greater(self.Wl, wavelengthMin), np.less(self.Wl, wavelengthMax))
        flux = np.trapz(self.Flux_wl[mask], self.Wl[mask])
        return flux

    def smooth(self, smoothPix):
        """Smooths SED.flux with a uniform (boxcar) filter of width smoothPix. Cannot be undone.

        @type smoothPix: int
        @param smoothPix: size of uniform filter applied to SED, in pixels

        """
        smoothed = ndimage.uniform_filter1d(self.flux, smoothPix)
        self.flux = smoothed

    @staticmethod
    def flux_to_distance(flux, dl):
        """
        Compute flux at distance dl
        :param flux:  flux
        :param dl: luminosity distance [parsec]
        :return: :raise ValueError:
        """
        return flux / (4 * np.pi * dl ** 2)


class SpectrumPlanck(Spectrum):
    """
        Planck Spectrum
    """

    def __init__(self, freq, temperature, name='bb'):
        flux = math.pi * rf.planck(freq, temperature=temperature)
        Spectrum.__init__(self, name, freq=freq, flux=flux)
        self._T = temperature
        # self.zeta = None

    @property
    def T(self):
        return self._T

    # def correct_zeta(self, zeta):
    #     self.zeta = zeta
    #     self._flux *= zeta ** 2


class SpectrumDilutePlanck(Spectrum):
    """
        Diluted Planck Spectrum
    """

    def __init__(self, freq, temperature, W, name='wbb'):
        flux = math.pi * W * rf.planck(freq, temperature=temperature)
        Spectrum.__init__(self, name, freq=freq, flux=flux)
        self._W = W
        self._T = temperature

    @property
    def T(self):
        return self._T

    @property
    def W(self):
        return self._W


class SeriesSpectrum(object):
    def __init__(self, name):
        """Creates a Series of Spectrum instance."""
        self.name = name
        self._nfreq = None
        # self._wl = None  # waves
        self._freq = None
        self._times = []  # [day]
        self._data = []  # array where index -> Spectrum at the times[index]

    def get_spec(self, idx):
        return self._data[idx]

    @property
    def is_time(self):
        return len(self._times) > 0

    @property
    def nfreq(self):
        return self._nfreq

    @property
    def Time(self):
        return np.array(self._times)

    @property
    def IdxMax(self):
        return len(self.Time)

    @property
    def Freq(self):
        return self._freq

    @property
    def Freq_min(self):
        res = float("inf")
        for k, t in enumerate(self.Time):
            s = self.get_spec(k)
            res = min(res, np.min(s.Freq))
        return res

    @property
    def Freq_max(self):
        res = float("-inf")
        for k, t in enumerate(self.Time):
            s = self.get_spec(k)
            res = max(res, np.max(s.Freq))
        return res

    @property
    def Wl(self):
        return phys.c / self.Freq

    @property
    def Wl_min(self):
        return phys.c / self.Freq_max

    @property
    def Wl_max(self):
        return phys.c / self.Freq_min

    @property
    def Data(self):
        return self._data

    # def set_times(self, times):
    #     if times is None or len(times) == 0:
    #         raise ValueError("times must be array with len > 0.")
    #     self._times = times

    def set_freq(self, freqs):
        if freqs is None or len(freqs) == 0:
            raise ValueError("freqs must be array with len > 0.")
        if np.all(freqs == 0.):
            raise ValueError("No freqs item equals 0.")
        self._freq = freqs
        # self._wl = phys.c / self._freq
        self._nfreq = len(self._freq)

    def add(self, t, spec):
        self._times.append(t)
        self._data.append(spec)
    # def set_data(self, data):
    #     if data is None or len(data) == 0:
    #         raise ValueError("data must be array with len > 0.")
    #     self._data = data

    def copy(self, nm=None, t_ab=None, wl_ab=None):
        """Copy SeriesSpectrum with times and wave lengths lying inside time and wave intervals.
        t_ab  - time interval [day],
        wl_ab - wave length interval [A]"""
        name = self.name if nm is None else nm
        if t_ab is None:
            t_ab = (float("-inf"), float("inf"))
        if wl_ab is None:
            wl_ab = (0., float("inf"))

        times = self.Time
        res = SeriesSpectrum(name)
        for k, t in enumerate(times):
            if t_ab[0] <= t <= t_ab[1]:
                s = self.get_spec(k)
                freq = s.Freq
                flux = s.Flux
                f_l = phys.c / (wl_ab[1]*phys.angs_to_cm)
                f_h = float("inf") if wl_ab[0] <= 0. else phys.c / (wl_ab[0]*phys.angs_to_cm)
                is_freq = (freq > f_l) & (freq < f_h)
                ss = Spectrum(s.Name, freq[is_freq], flux=flux[is_freq])
                res.add(t, ss)

        return res

    def get_tspec(self, idx):
        return self.Time[idx], self.get_spec(idx)

    def get_spec_by_time(self, time):
        idx = (np.abs(self.Time - time)).argmin()
        return self.get_spec(idx)

    def flux_to_mags(self, b, z=0., d=0., magnification=1.):
        if b is None:
            return None
        if not self.is_time:
            return None

        mags = np.zeros(len(self.Time))
        for k in range(len(self.Time)):
            star = Star(k, self.get_spec(k))
            star.set_distance(d)
            star.set_redshift(z)
            star.set_magnification(magnification)
            # mag = star.flux_to_mag(b)
            mag = star.flux_to_magAB(b)
            mags[k] = mag
            # if self.times[k] > 50:
            #     spec.plot_spec(title="t=%f, mag=%f" % (self.times[k], mag))

        return mags

    def flux_to_curve(self, b, z=0., d=0., magnification=1.):
        if b is None:
            raise ValueError("Band must be defined.")
        if not self.is_time:
            return ValueError("No spectral time points.")

        mags = self.flux_to_mags(b, z, d, magnification)
        # for k in range(len(self.Time)):
        #     star = Star(k, self.get_spec(k))
        #     star.set_distance(d)
        #     star.set_redshift(z)
        #     star.set_magnification(magnification)
        #     # mag = star.flux_to_mag(b)
        #     mag = star.flux_to_magAB(b)
        #     mags[k] = mag

        time = self.Time * (1. + z)
        lc = LightCurve(b, time, mags)
        return lc

    def flux_to_curves(self, bands, z=0., d=rf.pc_to_cm(10.), magnification=1.):
        """

        :param bands:
        :param z:
        :param d:
        :param magnification:
        :return:
        """
        curves = SetLightCurve(self.name)
        for n in bands:

            b = band.band_by_name(n)
            lc = self.flux_to_curve(b, z=z, d=d, magnification=magnification)
            curves.add(lc)
        return curves

    def mags_bands(self, bands, z=0., d=rf.pc_to_cm(10.), magnification=1.):
        mags = dict((k, None) for k in bands)
        for n in bands:
            b = band.band_by_name(n)
            mags[n] = self.flux_to_mags(b, z=z, d=d, magnification=magnification)

        mags['time'] = self.Time * (1. + z)
        return mags

