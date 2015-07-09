from math import log10
from pystella.util.phys_var import phys
from operator import itemgetter

__author__ = 'bakl'
from scipy import interpolate
from scipy.integrate import simps as integralfunc
import numpy as np


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
        self.flux_wl = self.flux_q * phys.c / self.wl**2  # flux [erg/cm^2/cm) ]

    def convolution_band(self, band):
        flux = Spectrum.flux_to_10pc(self.flux_wl)  # move to 10 pc
        # tck = interpolate.splrep(self.freq, self.flux_q, s=0)
        # flux_intrp = interpolate.splev(band.freq, tck, der=0)

        tck = interpolate.splrep(self.wl, flux, s=0)
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

    def response(self, band):
        if min(self.wl) > band.wl[0] or max(self.wl) < band.wl[-1]:
            raise ValueError("Spectrum must be wider then band: "+str(band))

        flux = Spectrum.flux_to_10pc(self.flux_wl)  # move to 10 pc
        flux = flux / phys.cm_to_angs  #  to flux [erg/cm^2/A) ]
        wl_s = self.wl * phys.cm_to_angs
        wl_b = band.wl*phys.cm_to_angs

        tck = interpolate.splrep(wl_s, flux, s=0)
        flux_intrp = interpolate.splev(wl_b, tck, der=0)

        # TODO what about h here?
        a = integralfunc(flux_intrp * band.resp*wl_b, wl_b)/(phys.c*phys.cm_to_angs)/phys.h
        return a

    def flux_to_mag(self, band):
        # conv = self.convolution_band(band)
        conv = self.response(band)
        if conv <= 0:
            return None
        else:
            mag = -2.5 * np.log10(conv) + band.zp
            return mag

    def flux_to_AB(self):
        magAB = -2.5 * np.log10(Spectrum.flux_to_10pc(self.flux_q)) + phys.ZP_AB
        return magAB

    @staticmethod
    def flux_to_10pc(flux):
        return flux / (4 * np.pi * (10 * phys.pc) ** 2)


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

    def flux_to_mags(self, band):
        if band is None:
            return None
        if not self.is_time():
            return None

        mags = np.zeros(len(self.data))
        for k in range(len(self.data)):
            spec = self.data[k]
            mag = spec.flux_to_mag(band)
            mags[k] = mag

        return mags
