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
        flux = np.ones(nf)
        self.sp = spectrum.Spectrum('uniform', wl=wl, flux=flux)

    def test_convolution_band(self):
        b = band.band_by_name('U')
        conv = self.sp.convolution_band(b)
        self.assertAlmostEqual(conv, 1.,
                               msg="For uniform flux=1 it should be mag==1. Now mag is %f for band %s." % (conv, b))

    def test_convolution_bands(self):
        bands = ['U', 'B', 'V', 'R', "I"]
        for n in bands:
            b = band.band_by_name(n)
            conv = self.sp.convolution_band(b)
            self.assertAlmostEqual(conv, 1.,
                                   msg="For uniform flux=1 it should be mag==1. Now mag is %f for band %s." % (conv, b))

    def test_flux_to_mag(self):
        bands = ['U', 'B', 'V', 'R', "I"]
        for n in bands:
            b = band.band_by_name(n)
            mag = self.sp.flux_to_mag(b)
            self.assertTrue(mag > 0.,
                                   msg="For uniform flux=1 it should be mag==1. Now mag is %f for band %s." % (mag, b))


def main():
    unittest.main()


if __name__ == '__main__':
    main()