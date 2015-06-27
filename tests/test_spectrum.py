import unittest
import src.rf.band as band
import src.rf.spectrum as spectrum
import numpy as np

__author__ = 'bakl'


class TestSpectrum(unittest.TestCase):
    def setUp(self):
        nf = 100
        start, end = 10, 1e5
        wl = np.linspace(start, end, nf)
        flux = np.ones(1, nf)
        self.sp = spectrum.Spectrum('uniform', wl=wl, flux=flux)

    def test_convolution_band(self):
        b = band.band_by_name('U')
        mag = self.sp.convolution_band(b)

    def test_convolution_bands(self):
        bands = ['U', 'B', 'V', 'R', "I"]
        sp = spectrum.Spectrum('test')
        for n in bands:
            b = band.band_by_name(n)
            self.assertTrue(b is not None, "Band %s does not exist." % n)


def main():
    unittest.main()


if __name__ == '__main__':
    main()