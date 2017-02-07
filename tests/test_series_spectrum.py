import numpy as np
import unittest

import pystella.rf.rad_func as rf
import pystella.rf.spectrum as spectrum
from pystella.util.phys_var import phys

__author__ = 'bakl'


class TestSeriesSpectrum(unittest.TestCase):
    def setUp(self):
        T = 5500
        nf, start, end = 100, 10., 1e5
        wl = np.exp(np.linspace(np.log(end), np.log(start), nf))
        freq = rf.val_to_hz(wl, inp="A")

        times = np.linspace(0., 200., 20)
        ss = spectrum.SeriesSpectrum("test_SeriesSpectrum")

        for k, t in enumerate(times):
            s = spectrum.SpectrumPlanck(freq=freq, temperature=T, name='test_spectrum')
            ss.add(t, s)

        self.series = ss

    def test_series_spectrum_copy_time(self):
        t_l, t_h = 0., 10.
        ss = self.series.copy(t_ab=(t_l, t_h))
        self.assertTrue(np.min(ss.Time) >= 0.,
                        "Low time boundary should be equal 0, but np.min(cp.Time) = %f" % np.min(ss.Time))
        self.assertTrue(np.max(ss.Time) <= t_h,
                        "Max time boundary should be less %f, but np.max(cp.Time) = %f" % (t_h, np.max(ss.Time)))

    def test_series_spectrum_copy_wl(self):
        wl_l, wl_h = 110., 10e3
        ss = self.series.copy(wl_ab=(wl_l, wl_h))
        for i, t in enumerate(ss.Time):
            s = ss.get_spec(i)
            self.assertTrue(np.min(s.Wl) * phys.cm_to_angs >= wl_l,
                            "Low wave length boundary for t=%f (idx=%d) should be equal 0, "
                            "but np.min(cp.Wl) = %f" % (t, i, np.min(s.Wl)))
            self.assertTrue(np.max(s.Wl) * phys.cm_to_angs <= wl_h,
                            "Max wave length boundary should be less %f, but np.max(cp.Wl) = %f" % (wl_h, np.max(s.Wl)))

    def test_series_get_T_color(self):
        aT_color = self.series.get_T_color()

        self.assertEqual(len(self.series.Time), len(aT_color), 'You should get array of T_color the same length as Time')

        for k, t in enumerate(self.series.Time):
            s = self.series.get_spec(k)
            self.assertAlmostEqual(s.T_color, aT_color[k],
                                   "The temperatures shoild be the same "
                                   "but k,time= [{},{}] s.T_color = {:f} ".format(k, t, s.T_color, aT_color[k]))

    def test_series_get_T_wien(self):
        aT_wien = self.series.get_T_wien()

        self.assertEqual(len(self.series.Time), len(aT_wien), 'You should get array of T_color the same length as Time')

        for k, t in enumerate(self.series.Time):
            s = self.series.get_spec(k)
            self.assertAlmostEqual(s.T_wien, aT_wien[k],
                                   "The temperatures shoild be the same "
                                   "but k,time= [{},{}] s.T_color = {:f} ".format(k, t, s.T_wien, aT_wien[k]))


def main():
    unittest.main()


if __name__ == '__main__':
    main()
