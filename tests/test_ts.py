import numpy as np
import unittest

from pystella.rf.ts import TimeSeries, SetTimeSeries

__author__ = 'bakl'


def ts_create(bname, m=-19, dt=0., n=10):
    time = np.linspace(0.+dt, 200.+dt, n)
    mags = m * np.ones(n)
    return TimeSeries(bname, time, mags)


class TestSetTimeSeries(unittest.TestCase):
    def test_SetTimeSeries_TimeCommon(self):
        bands = ['U', 'B', 'V']
        curves = SetTimeSeries()
        for b in bands:
            curves.add(ts_create(b))

        self.assertTrue(curves.IsCommonTime, msg="Error: IsCommonTime should be True")

        curves.add(ts_create('TT', dt=1.))
        self.assertFalse(curves.IsCommonTime, msg="Error: IsCommonTime should be False")


def main():
    unittest.main()


if __name__ == '__main__':
    main()
