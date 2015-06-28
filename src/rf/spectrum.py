from math import log10
from util.phys_var import phys

__author__ = 'bakl'
from scipy import interpolate
import numpy as np


class Spectrum:
    def __init__(self, name, wl=None, flux=None):
        """Creates a Spectrum instance.  Required parameters:  name."""
        self.name = name
        self.wl = wl  # wavelength of flux
        self.flux = flux  # flux

    def convolution_band(self, band):
        x = self.wl
        y = self.flux

        tck = interpolate.splrep(x, y, s=0)
        flux_band = interpolate.splev(band.wl, tck, der=0)

        a = np.trapz(band.wl * band.resp * flux_band, band.wl)
        b = np.trapz(band.wl * band.resp, band.wl)
        return a / b

    def flux_to_mag(self, band):
        conv = self.convolution_band(band)
        if conv <= 0:
            return None
        else:
            mag = -2.5 * np.log10(conv) + band.zp
            return mag


class SeriesSpectrum():
    def __init__(self, name):
        """Creates a Series of Spectrum instance."""
        self.name = name
        self.times = None
        self.nfreq = None
        self.wl = None  # waves
        self.freqs = None
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
        self.freqs = freqs
        self.wl = phys.c / self.freqs
        self.nfreq = len(self.freqs)

    def set_data(self, sdata):
        if sdata is None or len(sdata) == 0:
            raise ValueError("data must be array with len > 0.")
        self.data = sdata

    def flux_to_mag(self, band):
        if band is None:
            return None
        if not self.is_time():
            return None

        mags = []
        for k in range(len(self.data)):
            spec = self.data[k]
            mag = spec.flux_to_mag(band)
            mags[k] = mag

        return mags






