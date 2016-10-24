import unittest
import numpy as np
import pystella.rf.band as band
from pystella.rf.rad_func import Flux2MagAB, MagAB2Flux

__author__ = 'bakl'


class BandTests(unittest.TestCase):
    def test_available_bands(self):
        bands = ['U', 'B', 'V', 'R', "I"]
        for n in bands:
            b = band.band_by_name(n)
            self.assertTrue(b is not None, "Band %s does not exist." % b)

        bands = ['g', 'i', 'r', 'u', "z"]
        for n in bands:
            b = band.band_by_name(n)
            self.assertTrue(b is not None, "Band %s does not exist." % b)

        bands4 = dict(UVM2="photonUVM2.dat", UVW1="photonUVW1.dat", UVW2="photonUVW2.dat",
                      UVOTU="photonU_UVOT.dat", UVOTB="photonB_UVOT.dat", UVOTV="photonV_UVOT.dat")
        bands = ['UVOTU', 'UVOTB', 'UVOTV', 'UVM2', "UVW1", "UVW2"]
        for n in bands:
            b = band.band_by_name(n)
            self.assertTrue(b is not None, "Band %s does not exist." % b)

    def test_zero_point(self):
        zp = 0.790  # 13.9168580779  for KAIT U
        b = band.band_by_name('U')
        self.assertAlmostEqual(b.zp, zp, msg="Zero points of band %s equals %f. Should be %f" % (b.name, b.zp, zp))

    def test_band_uniform(self):
        b = band.BandUni()
        self.assertTrue(np.any(b.resp == 1), "Response values is equal 1. band: %s" % b.name)

    def test_zp_AB(self):
        # see https://www.gemini.edu/sciops/instruments/magnitudes-and-fluxes
        # see http://ssc.spitzer.caltech.edu/warmmission/propkit/pet/magtojy/
        qq = 1.
        qq1 = MagAB2Flux(Flux2MagAB(qq))
        self.assertAlmostEqual(qq, qq1, msg="MagAB2Flux(Flux2MagAB(x)) %f. Should be %f" % (qq, qq1), delta=0.05)

        # U
        f, ab = 1823., 0.748
        m_ab = Flux2MagAB(f * 1e-23)
        f_ab = MagAB2Flux(ab) * 1e23
        self.assertAlmostEqual(f_ab, f, msg="Zero points of band U equals %f. Should be %f" % (f_ab, f), delta=10.05)
        self.assertAlmostEqual(m_ab, ab, msg="Zero points of band U equals %f. Should be %f" % (m_ab, ab), delta=0.003)

        # B
        f, ab = 4260., -0.174
        m_ab = Flux2MagAB(f * 1e-23)
        f_ab = MagAB2Flux(ab) * 1e23
        self.assertAlmostEqual(f_ab, f, msg="Zero points of band B equals %f. Should be %f" % (f_ab, f), delta=10.05)
        self.assertAlmostEqual(m_ab, ab, msg="Zero points of band B equals %f. Should be %f" % (m_ab, ab), delta=0.003)

        # V
        f, ab = 3640., -0.0028  # https://www.astro.umd.edu/~ssm/ASTR620/mags.html
        # f, ab = 3781., -0.044
        m_ab = Flux2MagAB(f * 1e-23)
        f_ab = MagAB2Flux(ab) * 1e23
        self.assertAlmostEqual(f_ab, f, msg="Zero points of band V equals %f. Should be %f" % (f_ab, f), delta=10.05)
        self.assertAlmostEqual(m_ab, ab, msg="Zero points of band V equals %f. Should be %f" % (m_ab, ab), delta=0.005)

        # R
        f, ab = 3080., 0.18  # https://www.astro.umd.edu/~ssm/ASTR620/mags.html
        m_ab = Flux2MagAB(f * 1e-23)
        f_ab = MagAB2Flux(ab) * 1e23
        self.assertAlmostEqual(f_ab, f, msg="Zero points of band V equals %f. Should be %f" % (f_ab, f), delta=10.05)
        self.assertAlmostEqual(m_ab, ab, msg="Zero points of band V equals %f. Should be %f" % (m_ab, ab), delta=0.005)

        # I
        f, ab = 2550., 0.38  # https://www.astro.umd.edu/~ssm/ASTR620/mags.html
        # f, ab = 3781., -0.044
        m_ab = Flux2MagAB(f * 1e-23)
        f_ab = MagAB2Flux(ab) * 1e23
        self.assertAlmostEqual(f_ab, f, msg="Zero points of band V equals %f. Should be %f" % (f_ab, f), delta=10.05)
        self.assertAlmostEqual(m_ab, ab, msg="Zero points of band V equals %f. Should be %f" % (m_ab, ab), delta=0.005)

        # u
        f, ab = 1545., 0.927
        m_ab = Flux2MagAB(f * 1e-23)
        f_ab = MagAB2Flux(ab) * 1e23
        self.assertAlmostEqual(f_ab, f, msg="Zero points of band u equals %f. Should be %f" % (f_ab, f), delta=10.05)
        self.assertAlmostEqual(m_ab, ab, msg="Zero points of band u equals %f. Should be %f" % (m_ab, ab), delta=0.003)

        # g
        f, ab = 3991., -0.103
        m_ab = Flux2MagAB(f * 1e-23)
        f_ab = MagAB2Flux(ab) * 1e23
        self.assertAlmostEqual(f_ab, f, msg="Zero points of band g equals %f. Should be %f" % (f_ab, f), delta=10.05)
        self.assertAlmostEqual(m_ab, ab, msg="Zero points of band g equals %f. Should be %f" % (m_ab, ab), delta=0.005)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
