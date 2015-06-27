__author__ = 'bakl'

import unittest

import src.rf.band as band


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

    def test_zero_point(self):
        b = band.band_by_name('U')
        self.assertAlmostEqual(b.zp, 13.9168580779, "Band %s does not exist." % b)

def main():
    unittest.main()


if __name__ == '__main__':
    main()