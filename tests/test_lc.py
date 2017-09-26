import numpy as np
import unittest

from pystella.rf.band import Band
from pystella.rf.lc import SetLightCurve, LightCurve
from pystella.rf.light_curve_func import curves_save

__author__ = 'bakl'


def lc_create(bname, m=-19, dt=0., n=10, is_err=False):
    time = np.linspace(0.+dt, 200.+dt, n)
    mags = m * np.linspace(0.1, 1., n)
    band = Band(bname)
    if is_err:
        errs = m * np.linspace(0.01, 0.3, n)
        return LightCurve(band, time, mags, errs)
    else:
        return LightCurve(band, time, mags)


class TestLightCurve(unittest.TestCase):
    def test_BandName(self):
        band = 'U'
        lc = lc_create(band)
        self.assertEqual(band, lc.Band.Name,
                         msg="It should be equal band names.\n \
                         Now band is %s but  lc.Band.Name is  %s." % (band, lc.Band.Name))

    def test_lc_leastsq(self):
        dt_init = 10.
        lc1 = lc_create('U', dt=0.)
        lc2 = lc_create('U', dt=0.)

        # self.assertItemsEqual(bands, curves.BandNames,
        #                       msg="Error for band names.\n \
        #         Now band is %s but  lc.Band.Name is  %s." % (' '.join(bands), ' '.join(curves.BandNames)))


class TestSetLightCurve(unittest.TestCase):
    def test_SetLightCurve_BandNames(self):
        bands = ['U', 'B', 'V']
        curves = SetLightCurve()
        for b in bands:
            curves.add(lc_create(b))

        self.assertCountEqual(bands, curves.BandNames,
                              msg="Error for band names.\n \
                Now band is %s but  lc.Band.Name is  %s." % (' '.join(bands), ' '.join(curves.BandNames)))

    def test_SetLightCurve_save_true(self):
            bands = ['U', 'B', 'V']
            curves = SetLightCurve()
            for b in bands:
                curves.add(lc_create(b))
            res = curves_save(curves, 'tmp_curves')
            self.assertTrue(res, msg="Error:  curves_save should return True")

    def test_SetLightCurve_save_true_with_errors(self):
            bands = ['U', 'B', 'V']
            curves = SetLightCurve()
            for b in bands:
                curves.add(lc_create(b, is_err=True))
            curves.add(lc_create('I'))
            res = curves_save(curves, 'tmp_curves')
            self.assertTrue(res, msg="Error:  curves_save should return True")

    def test_SetLightCurve_save_NoIsCommonTime(self):
            bands = ['U', 'B', 'V']
            curves = SetLightCurve()
            for b in bands:
                curves.add(lc_create(b))
            curves.add(lc_create('TimeDif', dt=1.))

            res = curves_save(curves, 'tmp_curves_2')
            self.assertFalse(res, msg="Error:  curves_save should return False due to IsCommonTime=False")


def main():
    unittest.main()


if __name__ == '__main__':
    main()
