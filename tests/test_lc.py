import numpy as np
import unittest

from pystella.rf.lc import SetLightCurve, LightCurve

__author__ = 'bakl'


def lc_create(band):
    n = 10
    time = np.linspace(0., 200., n)
    mags = -19. * np.ones(n)
    return LightCurve(band, time, mags)


class TestLightCurve(unittest.TestCase):
    def test_BandName(self):
        band = 'U'
        lc = lc_create(band)
        self.assertEqual(band, lc.Band.Name,
                         msg="It should be equal band names.\n \
                         Now band is %s but  lc.Band.Name is  %s." % (band, lc.Band.Name))

    def test_SetLightCurve_BandNames(self):
        bands = ['U', 'B', 'V']
        curves = SetLightCurve()
        for b in bands:
            curves.add(lc_create(b))

        self.assertItemsEqual(bands, curves.BandNames,
                              msg="Error for band names.\n \
                Now band is %s but  lc.Band.Name is  %s." % (' '.join(bands), ' '.join(curves.BandNames)))


def main():
    unittest.main()


if __name__ == '__main__':
    main()
