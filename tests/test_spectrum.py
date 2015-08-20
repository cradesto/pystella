import unittest
import numpy as np

import pystella.rf.spectrum as spectrum
import pystella.util.rf as rf

__author__ = 'bakl'


class TestSpectrum(unittest.TestCase):
    def test_wl_sort(self):
        nf, start, end = 100, 10., 1e5
        wl = np.exp(np.linspace(np.log(end), np.log(start), nf))
        freq = rf.val_to_hz(wl, inp="A")
        flux = np.ones(nf)
        sp = spectrum.Spectrum('uniform', freq=freq, flux=flux, is_sort_wl=True)
        freq = sp.freq
        for i in range(len(freq) - 1):
            self.assertTrue(freq[i] > freq[i + 1],
                            "Freq should be ordered, but freq[%d] > freq[%d]." % (i, i + 1))


def main():
    unittest.main()


if __name__ == '__main__':
    main()
