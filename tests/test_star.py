import unittest
import pystella.rf.band as band
import pystella.rf.spectrum as spectrum
import numpy as np
from pystella.rf.star import Star
from pystella.util.phys_var import phys
import pystella.util.rf as rf

__author__ = 'bakl'


class TestStar(unittest.TestCase):
    def setUp(self):
        nf = 100
        start, end = 10, 1e5
        wl = np.exp(np.linspace(np.log(start), np.log(end), nf))
        freq = rf.val_to_hz(wl, inp="A")
        # freq = freq[np.argsort(freq)]  # ascending order
        flux = np.ones(nf)
        self.sp = spectrum.Spectrum('uniform', freq=freq, flux=flux)
        self.distance = 10  # pc

    def test_check_band_zp_UBVRI(self):
        star = Star('test', self.sp)
        # see https://www.astro.umd.edu/~ssm/ASTR620/mags.html#conversions
        delta = {'B': 0.163, 'V': 0.044, 'R': -0.055, "I": -0.309}
        bands = delta.keys()
        for n in bands:
            b = band.band_by_name(n)
            # mag = star.flux_to_mag(b) - phys.ZP_AB_lmb

            mag = star.flux_to_magAB(b) - phys.ZP_AB
            # mag = star.flux_to_magAB(b) - delta[n] - phys.ZP_AB

            self.assertAlmostEqual(mag, 0., delta=0.5,
                                   msg="For uniform flux=1 it should be mag==AB zero point.\n \
                                        Now mag is %f for band %s. ZP is %f" % (mag, b, phys.ZP_AB))

    def test_check_band_zp_ugzri(self):
        bands = ['u', 'g', 'z', 'r', "i"]
        star = Star('test', self.sp)
        # star.set_radius_ph(self.distance)
        # star.set_distance(self.distance)
        for n in bands:
            b = band.band_by_name(n)
            # mag = star.flux_to_mag(b)
            mag = star.flux_to_magAB(b)
            self.assertAlmostEqual(mag, phys.ZP_AB, delta=1.,
                                   msg="For uniform flux=1 it should be mag==AB zero point.\n \
                                        Now mag is %f for band %s. ZP is %f" % (mag, b, phys.ZP_AB))

    def test_check_band_zp_ps1(self):
        bands = ['gps1', 'zps1', 'rps1', "ips1", "y", "w"]
        star = Star('test', self.sp)
        # star.set_radius_ph(self.distance)
        # star.set_distance(self.distance)
        for n in bands:
            b = band.band_by_name(n)
            mag = star.flux_to_magAB(b)
            self.assertAlmostEqual(mag, phys.ZP_AB, delta=1.,
                                   msg="For uniform flux=1 it should be mag==AB zero point.\n \
                                        Now mag is %f for band %s. ZP is %f" % (mag, b, phys.ZP_AB))

    def test_check_band_zp_sdss(self):
        bands = ['UVM2', 'UVW1', 'UVW2']
        star = Star('test', self.sp)
        # star.set_radius_ph(self.distance)
        # star.set_distance(self.distance)
        for n in bands:
            b = band.band_by_name(n)
            mag = star.flux_to_magAB(b)
            self.assertAlmostEqual(mag, phys.ZP_AB, delta=1.,
                                   msg="For uniform flux=1 it should be mag==AB zero point.\n \
                                        Now mag is %f for band %s. ZP is %f" % (mag, b, phys.ZP_AB))

    def test_k_cor_uniform(self):
        b_r = band.band_by_name('U')
        b_o = band.band_by_name('U')
        z = 0.2
        star = Star('test', self.sp)
        k_cor = star.k_cor(b_r, b_o, z=z)
        self.assertIsNotNone(k_cor, "Return error for k_cor. \
                            Band-rest %s and band-obs %s." % (b_r, b_o))
        self.assertAlmostEqual(k_cor, 0.,
                               msg="For uniform flux=1 it should be k_cor==0. \
                            Now k_kor is %f for band-rest %s and band-obs %s." % (k_cor, b_r, b_o), delta=0.0001)

    def test_k_cor_z(self):
        b_r = band.band_by_name('V')
        b_o = band.band_by_name('V')
        z = 0.1
        star = Star('test1', self.sp)
        star.set_redshift(z)
        k_cor = star.k_cor(b_r, b_o)
        self.assertIsNotNone(k_cor, "Return error for k_cor. \
                            Band-rest %s and band-obs %s." % (b_r, b_o))
        self.assertNotAlmostEqual(k_cor, 0.,
                               msg="For z=%f it should be k_cor != 0. \
                            Now k_kor is %f for band-rest %s and band-obs %s." % (z, k_cor, b_r, b_o))


def main():
    unittest.main()


if __name__ == '__main__':
    main()