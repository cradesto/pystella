__author__ = 'bakl'


class Spectrum:
    def __init__(self, name, wl=None, flux=None):
        """Creates a Spectrum instance.  Required parameters:  name."""
        self.name = name
        self.wl = wl  # wavelength of flux
        self.flux = flux  # flux

    def convolution_band(self, band):
