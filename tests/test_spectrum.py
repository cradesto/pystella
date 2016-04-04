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


class TestSeriesSpectrum(unittest.TestCase):
    def setUp(self):
        nf, start, end = 100, 10., 1e5
        wl = np.exp(np.linspace(np.log(end), np.log(start), nf))
        freq = rf.val_to_hz(wl, inp="A")

        times = np.linspace(0., 200., 20)
        ss = spectrum.SeriesSpectrum("test_SeriesSpectrum")

        for k, t in enumerate(times):
            flux = np.ones(nf)  # * np.sin(t)
            s = spectrum.Spectrum('test_spectrum', freq=freq, flux=flux)
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
            self.assertTrue(np.min(s.Wl)*phys.cm_to_angs >= wl_l,
                            "Low wave length boundary for t=%f (idx=%d) should be equal 0, "
                            "but np.min(cp.Wl) = %f" % (t, i, np.min(s.Wl)))
            self.assertTrue(np.max(s.Wl)*phys.cm_to_angs <= wl_h,
                            "Max wave length boundary should be less %f, but np.max(cp.Wl) = %f" % (wl_h, np.max(s.Wl)))


def main():
    unittest.main()


if __name__ == '__main__':
    main()
