__author__ = 'bakl'

import unittest

import src.band.band as band


class BandTests(unittest.TestCase):
    def testAvailableBands(self):
        bands = ['U', 'B', 'V', 'R', "I"]
        for n in bands:
            b = band.band_by_name(n)
            self.assertTrue(b is not None, "Band %s does not exist." % n)


def main():
    unittest.main()


if __name__ == '__main__':
    main()