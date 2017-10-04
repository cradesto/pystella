import numpy as np
import unittest

import pystella.rf.rad_func as rf
import pystella.rf.spectrum as spectrum
from pystella.util.phys_var import phys

__author__ = 'bakl'


class TestSpectrum(unittest.TestCase):
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

    def test_spectrum_T_color_zeta(self):
        nf, start, end = 100, 10., 1e5
        wl = np.exp(np.linspace(np.log(end), np.log(start), nf))
        freq = rf.val_to_hz(wl, inp="A")
        W = 0.5
        Tcolor = 7e3
        sp = spectrum.SpectrumDilutePlanck(freq, Tcolor, W, 'test')

        Tbb, zeta = sp.T_color_zeta()
        self.assertAlmostEqual(Tcolor, Tbb, 5, "Init Tcolor={:.2f} should be as Tbb={:.2f}.".format(Tcolor, Tbb))
        self.assertAlmostEqual(W, zeta, 5, "Init W={:.2f} should be as zeta={:.2f}.".format(W, zeta))


def main():
    unittest.main()


if __name__ == '__main__':
    main()
