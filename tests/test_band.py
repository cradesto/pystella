__author__ = 'bakl'

import unittest

import src.rf.band as band


class BandTests(unittest.TestCase):
    def test_available_bands(self):
        bands = ['U', 'B', 'V', 'R', "I"]
        for n in bands:
            b = band.band_by_name(n)
            self.assertTrue(b is not None, "Band %s does not exist." % n)

        bands = ['g', 'i', 'r', 'u', "z"]
        for n in bands:
            b = band.band_by_name(n)
            self.assertTrue(b is not None, "Band %s does not exist." % n)


def main():
    unittest.main()


if __name__ == '__main__':
    main()