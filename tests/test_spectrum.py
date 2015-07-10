import unittest
import pystella.rf.band as band
import pystella.rf.spectrum as spectrum
import numpy as np
from pystella.util.phys_var import phys

__author__ = 'bakl'


class TestSpectrum(unittest.TestCase):
    def setUp(self):
        nf = 100
        start, end = 10, 1e5
        wl = np.linspace(start, end, nf)
        freq = phys.c / (wl * 1.e-8)
        # freq = freq[np.argsort(freq)]  # ascending order
        flux = np.ones(nf)
        self.sp = spectrum.Spectrum('uniform', freq=freq, flux=flux)

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

    def test_k_cor_uniform(self):
        b_r = band.band_by_name('U')
        b_o = band.band_by_name('U')
        z = 0
        k_cor, err = self.sp.k_cor(b_r, b_o, z=z)
        self.assertTrue(err > 0., "Return error for k_cor. \
                            Now k_kor is %f for band-rest %s and band-obs %s." % (k_cor, b_r, b_o))
        self.assertAlmostEqual(k_cor, 0.,
                               msg="For uniform flux=1 it should be k_cor==0. \
                            Now k_kor is %f for band-rest %s and band-obs %s." % (k_cor, b_r, b_o))

    def test_k_cor_z(self):
        b_r = band.band_by_name('V')
        b_o = band.band_by_name('V')
        z = 0.1
        k_cor, err = self.sp.k_cor(b_r, b_o, z=z)
        self.assertTrue(err > 0., "Return error for k_cor. \
                            Now k_kor is %f for band-rest %s and band-obs %s." % (k_cor, b_r, b_o))



def main():
    unittest.main()


if __name__ == '__main__':
    main()
