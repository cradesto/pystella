import numpy as np
import unittest

import pystella as ps

__author__ = 'bakl'


class TestStar(unittest.TestCase):
    def setUp(self):
        nf = 100
        start, end = 10, 1e5
        wl = np.exp(np.linspace(np.log(start), np.log(end), nf))
        freq = ps.rf.val_to_hz(wl, inp="A")
        # freq = freq[np.argsort(freq)]  # ascending order
        flux = np.ones(nf)
        self.sp = ps.Spectrum('uniform', freq=freq, flux=flux)
        self.distance = 10  # pc

    def test_check_band_zp_UBVRI(self):
        star = ps.Star('test', self.sp)
        # see https://www.astro.umd.edu/~ssm/ASTR620/mags.html#conversions
        delta = {'B': 0.163, 'V': 0.044, 'R': -0.055, "I": -0.309}
        bands = delta.keys()
        for n in bands:
            b = ps.band.band_by_name(n)
            # mag = star.flux_to_mag(b) - phys.ZP_AB_lmb

            mag = star.magAB(b) - ps.phys.ZP_AB
            # mag = star.flux_to_magAB(b) - delta[n] - phys.ZP_AB

            self.assertAlmostEqual(mag, 0., delta=0.5,
                                   msg="For uniform flux=1 it should be mag==AB zero point.\n \
                                        Now mag is %f for band %s. ZP is %f" % (mag, b, ps.phys.ZP_AB))

    def test_check_band_zp_ugzri(self):
        bands = ['u', 'g', 'z', 'r', "i"]
        star = ps.Star('test', self.sp)
        # star.set_radius_ph(self.distance)
        # star.set_distance(self.distance)
        for n in bands:
            b = ps.band.band_by_name(n)
            # mag = star.flux_to_mag(b)
            mag = star.magAB(b)
            self.assertAlmostEqual(mag, ps.phys.ZP_AB, delta=1.,
                                   msg="For uniform flux=1 it should be mag==AB zero point.\n \
                                        Now mag is %f for band %s. ZP is %f" % (mag, b, ps.phys.ZP_AB))

    def test_check_band_zp_ps1(self):
        bands = ['PS1g', 'PS1z', 'PS1r', "PS1i", "PS1y", "PS1w"]
        star = ps.Star('test', self.sp)
        for n in bands:
            b = ps.band.band_by_name(n)
            mag = star.magAB(b)
            self.assertAlmostEqual(mag, ps.phys.ZP_AB, delta=1.,
                                   msg="For uniform flux=1 it should be mag==AB zero point.\n \
                                        Now mag is %f for band %s. ZP is %f" % (mag, b, ps.phys.ZP_AB))

    def test_check_band_zp_sdss(self):
        bands = ['UVM2', 'UVW1', 'UVW2']
        star = ps.Star('test', self.sp)
        # star.set_radius_ph(self.distance)
        # star.set_distance(self.distance)
        for n in bands:
            b = ps.band.band_by_name(n)
            mag = star.magAB(b)
            self.assertAlmostEqual(mag, ps.phys.ZP_AB, delta=1.,
                                   msg="For uniform flux=1 it should be mag==AB zero point.\n \
                                        Now mag is %f for band %s. ZP is %f" % (mag, b, ps.phys.ZP_AB))

    def test_k_cor_uniform(self):
        b_r = ps.band.band_by_name('U')
        b_o = ps.band.band_by_name('U')
        z = 0.2
        star = ps.Star('test', self.sp)
        k_cor = star.k_cor(b_r, b_o, z=z)
        self.assertIsNotNone(k_cor, "Return error for k_cor. \
                            Band-rest %s and band-obs %s." % (b_r, b_o))
        self.assertAlmostEqual(k_cor, 0.,
                               msg="For uniform flux=1 it should be k_cor==0. \
                            Now k_kor is %f for band-rest %s and band-obs %s." % (k_cor, b_r, b_o), delta=0.0001)

    def test_k_cor_z(self):
        b_r = ps.band.band_by_name('V')
        b_o = ps.band.band_by_name('V')
        z = 0.1
        star = ps.Star('test1', self.sp)
        star.set_redshift(z)
        k_cor = star.k_cor(b_r, b_o)
        self.assertIsNotNone(k_cor, "Return error for k_cor. \
                            Band-rest %s and band-obs %s." % (b_r, b_o))
        self.assertNotAlmostEqual(k_cor, 0.,
                                  msg="For z=%f it should be k_cor != 0. \
                            Now k_kor is %f for band-rest %s and band-obs %s." % (z, k_cor, b_r, b_o))

    def test_Lum2MagBol(self):
        L_sun = 3.827e33
        m = ps.rf.Lum2MagBol(L_sun)
        self.assertAlmostEqual(m, ps.phys.Mag_sun,
                               msg="Absolute magnitude of Sun is %f. \
                                    You have  m = %f." % (ps.phys.Mag_sun, m), delta=0.05)
        # Vega, see http://iopscience.iop.org/article/10.1088/0004-637X/708/1/71/meta
        L = 40.34 * ps.phys.L_sun
        m = ps.rf.Lum2MagBol(L)
        # absolute visual magnitude, see http://iopscience.iop.org/article/10.1088/0004-6256/136/1/452/meta
        m_vega = 0.582
        self.assertAlmostEqual(m, m_vega, msg="Absolute magnitude of Vega is %4.3f. You have  m = %4.3f." % (m_vega, m),
                               delta=0.05)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
