# coding=utf-8
import unittest

import numpy as np

from pystella.rf import extinction
from pystella.rf.lc import LightCurve
from pystella.rf.light_curve_func import curves_reddening

__author__ = 'bakl'


def lc_create(band, m=-19, dt=0.):
    n = 10
    time = np.linspace(0.+dt, 200.+dt, n)
    mags = m * np.ones(n)
    return LightCurve(band, time, mags)


class TestReddening(unittest.TestCase):
    def test_names(self):
        law = 'Rv3.1'
        data, laws = extinction.read_data()
        self.assertTrue(law in laws,
                        msg="There is no  %s  in laws: %s." % (law, ' '.join(laws)))

    def test_reddening_shift(self):
        """
        see extinction_schlafly.csv
        Bandpass       λeff      2.1     3.1     4.1     5.1
        Landolt_U     3508.2    5.614   4.334   3.773   3.460
        Landolt_B     4329.0    4.355   3.626   3.290   3.096
        Landolt_V     5421.7    2.953   2.742   2.645   2.589
        Landolt_R     6427.8    2.124   2.169   2.189   2.201
        Landolt_I     8048.4    1.410   1.505   1.548   1.573
        :return:
        """
        # band, mshift = 'V', 2.742
        band, mshift = 'U', 4.334
        ebv = 1.
        lc = lc_create(band)
        magExpect = lc.Mag + mshift  # Av for e=1.
        curves_reddening(lc, ebv=ebv)
        magRed = lc.Mag
        res = list(filter(lambda x: abs(x) > 1e-4, magExpect - magRed))
        self.assertTrue(len(res) == 0,
                        msg="Some mags were shifted more then %f." % mshift)

    def test_reddening_shift_with_z(self):
        """
        see extinction_schlafly.csv
        Bandpass       λeff      2.1     3.1     4.1     5.1
        Landolt_U     3508.2    5.614   4.334   3.773   3.460
        Landolt_B     4329.0    4.355   3.626   3.290   3.096
        Landolt_V     5421.7    2.953   2.742   2.645   2.589
        Landolt_R     6427.8    2.124   2.169   2.189   2.201
        Landolt_I     8048.4    1.410   1.505   1.548   1.573
        :return:
        """
        # band, mshift = 'V', 2.742
        band, mshift = 'V', 4.334  # shifted to U
        z = 5421.7/3508.2 - 1.
        ebv = 1.
        lc = lc_create(band)
        magExpect = lc.Mag + mshift  # Av for e=1.
        curves_reddening(lc, ebv=ebv, z=z)
        magRed = lc.Mag
        res = list(filter(lambda x: abs(x) > 1e-4, magExpect - magRed))
        self.assertTrue(len(res) == 0,
                        msg="Some mags were shifted more then %f." % mshift)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
