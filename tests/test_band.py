import unittest
import numpy as np
import pystella.rf.band as band
from pystella.rf.rad_func import Flux2MagAB, MagAB2Flux
from pystella.util.phys_var import phys

__author__ = 'bakl'


class BandTests(unittest.TestCase):
    def test_load_names(self):
        bands = band.band_load_names()
        self.assertTrue(len(bands) > 0, "You have to see more bands. Not %d" % len(bands))

    def test_band_colors_name(self):
        bands = band.band_load_names()
        for bname in bands:
            self.assertTrue(bname in band.colors(), "You have not color for band: %s" % bname)

    def test_band_by_name(self):
        b = band.band_by_name("BesU")
        self.assertTrue(b.is_load, "The band should be loaded and with data")

    def test_aliases_load(self):
        band.Band.load_settings()
        aliases = band.band_get_aliases()
        self.assertTrue(len(aliases), "Should be more aliases.")

    def test_aliases(self):
        bo = band.band_by_name("BesU")
        ba = band.band_by_name("U")
        self.assertTrue(ba.is_load, "The band should be loaded and with data")
        self.assertCountEqual(bo.wl, ba.wl, msg="The alias wl should be the same as original")
        self.assertCountEqual(bo.resp_wl, ba.resp_wl, msg="The alias wl should be the same as original")

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
                      SwiftU="photonU_Swift.dat", SwiftB="photonB_Swift.dat", SwiftV="photonV_Swift.dat")
        bands = ['SwiftU', 'SwiftB', 'SwiftV', 'UVM2', "UVW1", "UVW2"]
        for n in bands:
            b = band.band_by_name(n)
            self.assertTrue(b is not None, "Band %s does not exist." % n)

    def test_zero_point(self):
        zp = 0.748  # See filters.ini
        b = band.band_by_name('U')
        self.assertAlmostEqual(b.zp, zp, msg="Zero points of band %s equals %f. Should be %f" % (b.Name, b.zp, zp))

    def test_band_ubvri(self):
        import pylab as plt

        b = band.band_by_name(band.Band.NameUBVRI)
        plt.plot(b.wl * phys.cm_to_angs, b.resp_wl, band.colors(b.Name), label=b.Name, linewidth=2)
        plt.legend(loc=4)
        plt.ylabel('Amplitude Response')
        plt.xlabel('Wave [A]')
        plt.grid(linestyle=':')
        plt.show()

    def test_wl_eff(self):
        # wl_eff = {'U': 3650, 'B': 4450, 'V': 5510, 'R': 6580, 'I': 8060}
        wl_eff = {'u': 3560, 'g': 4830, 'r': 6260, 'i': 7670, 'z': 8890,
                  'U': 3600, 'B': 4380, 'V': 5450, 'R': 6410, 'I': 7980, 'J': 12200, 'H': 16300,
                  'K': 21900}
        for bname, wl in wl_eff.items():
            b = band.band_by_name(bname)
            res = b.wl_eff_angs
            print('{} {:.0f} VS {:.0f}'.format(bname, res, wl))
            self.assertAlmostEqual(res, wl, delta=wl * 0.03, msg="The effective wavelength of band %s equals %f. "
                                                                 "Should be %f" % (b.Name, res, wl))

    def test_fwhm(self):
        # wl_eff = {'U': 3650, 'B': 4450, 'V': 5510, 'R': 6580, 'I': 8060}
        wl_fwhm = {'U': 660., 'B': 940., 'V': 880, 'R': 1380., 'I': 1490.,
                   'J': 2130., 'H': 3070., 'K': 3900.}  # AA
        # convert to cm
        wl_fwhm = {bn: wl * 1e-8 for bn, wl in wl_fwhm.items()}
        for bname, wl in wl_fwhm.items():
            b = band.band_by_name(bname)
            res = b.fwhm
            print('{} {:.3e} VS {:.3e}'.format(bname, res, wl))
            self.assertAlmostEqual(res, wl, delta=wl * 0.1,
                                   msg="The fwhm of band {} equals {:.3e}.  Should be {:.3e}".format(b.Name, res, wl))

    def test_band_uniform(self):
        b = band.BandUni()
        self.assertTrue(np.any(b.resp_wl == 1), "Response values is equal 1. band: %s" % b.name)
        self.assertTrue(np.any(b.resp_fr == 1), "Response values is equal 1. band: %s" % b.name)

    def test_band_zp_vs_Jy(self):
        bands = band.band_load_names()
        for bname in bands:
            b = band.band_by_name(bname)
            if b.is_zp and b.is_Jy:
                m_ab = Flux2MagAB(b.Jy * phys.jy_to_erg)
                self.assertAlmostEqual(m_ab, b.zp, msg="Band [%s] zp and Jy should be coincide each other. "
                                                       "zp=%f,  m_zp(Jy) = %f, Jy = %f"
                                                       % (b.Name, b.zp, m_ab, b.Jy),
                                       delta=0.01)

    def test_zp_AB(self):
        # see https://www.gemini.edu/sciops/instruments/magnitudes-and-fluxes
        # see http://ssc.spitzer.caltech.edu/warmmission/propkit/pet/magtojy/
        qq = 1.
        qq1 = MagAB2Flux(Flux2MagAB(qq))
        self.assertAlmostEqual(qq, qq1, msg="MagAB2Flux(Flux2MagAB(x)) %f. Should be %f" % (qq, qq1), delta=0.05)

        # U
        f, ab = 1823., 0.748
        m_ab = Flux2MagAB(f * phys.jy_to_erg)
        f_ab = MagAB2Flux(ab) / phys.jy_to_erg
        self.assertAlmostEqual(f_ab, f, msg="Zero points of band U equals %f. Should be %f" % (f_ab, f), delta=10.05)
        self.assertAlmostEqual(m_ab, ab, msg="Zero points of band U equals %f. Should be %f" % (m_ab, ab), delta=0.003)

        # B
        f, ab = 4260., -0.174
        m_ab = Flux2MagAB(f * phys.jy_to_erg)
        f_ab = MagAB2Flux(ab) / phys.jy_to_erg
        self.assertAlmostEqual(f_ab, f, msg="Zero points of band B equals %f. Should be %f" % (f_ab, f), delta=10.05)
        self.assertAlmostEqual(m_ab, ab, msg="Zero points of band B equals %f. Should be %f" % (m_ab, ab), delta=0.003)

        # V
        f, ab = 3640., -0.0028  # https://www.astro.umd.edu/~ssm/ASTR620/mags.html
        # f, ab = 3781., -0.044
        m_ab = Flux2MagAB(f * phys.jy_to_erg)
        f_ab = MagAB2Flux(ab) / phys.jy_to_erg
        self.assertAlmostEqual(f_ab, f, msg="Zero points of band V equals %f. Should be %f" % (f_ab, f), delta=10.05)
        self.assertAlmostEqual(m_ab, ab, msg="Zero points of band V equals %f. Should be %f" % (m_ab, ab), delta=0.005)

        # R
        f, ab = 3080., 0.18  # https://www.astro.umd.edu/~ssm/ASTR620/mags.html
        m_ab = Flux2MagAB(f * phys.jy_to_erg)
        f_ab = MagAB2Flux(ab) / phys.jy_to_erg
        self.assertAlmostEqual(f_ab, f, msg="Zero points of band V equals %f. Should be %f" % (f_ab, f), delta=10.05)
        self.assertAlmostEqual(m_ab, ab, msg="Zero points of band V equals %f. Should be %f" % (m_ab, ab), delta=0.005)

        # I
        f, ab = 2550., 0.38  # https://www.astro.umd.edu/~ssm/ASTR620/mags.html
        # f, ab = 3781., -0.044
        m_ab = Flux2MagAB(f * phys.jy_to_erg)
        f_ab = MagAB2Flux(ab) / phys.jy_to_erg
        self.assertAlmostEqual(f_ab, f, msg="Zero points of band V equals %f. Should be %f" % (f_ab, f), delta=10.05)
        self.assertAlmostEqual(m_ab, ab, msg="Zero points of band V equals %f. Should be %f" % (m_ab, ab), delta=0.005)

        # J
        f, ab = 1600., 0.88970004336
        m_ab = Flux2MagAB(f * phys.jy_to_erg)
        f_ab = MagAB2Flux(ab) / phys.jy_to_erg
        print("Flux of Zero points of band u equals %f. m_zp = %f" % (f, m_ab))
        self.assertAlmostEqual(f_ab, f, msg="Zero points of band u equals %f. Should be %f" % (f_ab, f), delta=10.05)
        self.assertAlmostEqual(m_ab, ab, msg="Zero points of band u equals %f. Should be %f" % (m_ab, ab), delta=0.003)

        # g
        f, ab = 3991., -0.103
        m_ab = Flux2MagAB(f * phys.jy_to_erg)
        f_ab = MagAB2Flux(ab) / phys.jy_to_erg
        self.assertAlmostEqual(f_ab, f, msg="Zero points of band g equals %f. Should be %f" % (f_ab, f), delta=10.05)
        self.assertAlmostEqual(m_ab, ab, msg="Zero points of band g equals %f. Should be %f" % (m_ab, ab), delta=0.005)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
