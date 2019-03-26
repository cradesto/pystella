import numpy as np
import unittest

import pystella.rf.rad_func as rf
import pystella.rf.spectrum as spectrum
from pystella.util.phys_var import phys
from   pystella.rf.rad_func import kcorrection, kcorrection_spec
from   pystella.rf.band import band_is_exist, band_by_name

__author__ = 'bakl'


class TestSpec(unittest.TestCase):
    def setUp(self):
        nf, start, end = 100, 10., 1e5
        wl = np.exp(np.linspace(np.log(end), np.log(start), nf))
        freq = rf.val_to_hz(wl, inp="A")
        flux = np.ones(nf)
        self.sp = spectrum.Spectrum('test_spectrum', freq=freq, flux=flux, is_sort_wl=True)

    def test_wl_sort(self):
        freq = self.sp.Freq
        for i in range(len(freq) - 1):
            self.assertTrue(freq[i] > freq[i + 1],
                            "Freq should be ordered, but freq[%d] > freq[%d]." % (i, i + 1))

    def test_spec_kcorr_z0(self):
        nf, start, end = 100, 10., 1e5
        wl = np.exp(np.linspace(np.log(end), np.log(start), nf))
        freq = rf.val_to_hz(wl, inp="A")
        W = 0.5
        Tcolor = 7e3
        sp = spectrum.SpectrumDilutePlanck(freq, Tcolor, W, 'test')

        z = 0.
        bn_rest = 'V'
        # bn_obs = 'V'
        for bn_rest in ['U', 'B', 'V', 'R', 'g']:
            bn_obs = bn_rest
            br = band_by_name(bn_rest)
            bo = band_by_name(bn_obs)

            kcorr = kcorrection_spec(sp, z, br, bo)

            self.assertAlmostEqual(kcorr, 0., 8, "{}: For z={} kcorr should be 0. But kcorr={:.5f}.".
                                   format(bn_rest, z, kcorr))

    def test_spec_kcorr_z0_bad(self):  # todo WHY???
        nf, start, end = 100, 10., 1e5
        wl = np.exp(np.linspace(np.log(end), np.log(start), nf))
        freq = rf.val_to_hz(wl, inp="A")
        W = 0.5
        Tcolor = 7e3
        sp = spectrum.SpectrumDilutePlanck(freq, Tcolor, W, 'test')

        z = 0.
        bn_rest = 'I'
        # bn_obs = 'V'
        # for bn_rest in ['I', 'r']:
        bn_obs = bn_rest
        br = band_by_name(bn_rest)
        bo = band_by_name(bn_obs)

        kcorr = kcorrection_spec(sp, z, br, bo)

        self.assertAlmostEqual(kcorr, 0., 8, "{}: For z={} kcorr should be 0. But kcorr={:.5f}.".
                               format(bn_rest, z, kcorr))


def main():
    unittest.main()


if __name__ == '__main__':
    main()
