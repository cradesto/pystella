from math import log10

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
        b = np.trapz(band.wl * band.resp,  band.wl)
        return a / b

    def flux_to_mag(self, band):
        conv = self.convolution_band(band)
        if conv <= 0:
            return np.nan
        else:
            mag = (-2.5 * np.log10(conv) + band.zp)
            return mag
