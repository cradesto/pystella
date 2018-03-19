import numpy as np
import unittest

import pystella as ps
# from pystella.rf.band import Band
# from pystella.rf.lc import SetLightCurve, LightCurve, LC_interp
# from pystella.rf.light_curve_func import curves_save, curves_read

__author__ = 'bakl'


def lc_create(bname, m=-19, dt=0., n=10, tend=200., is_err=False):
    time = np.linspace(0. + dt, tend + dt, n)
    mags = m * np.linspace(0.1, 1., n)
    band = ps.Band(bname)
    if is_err:
        errs = m * np.linspace(0.01, 0.3, n)
        return ps.LightCurve(band, time, mags, errs)
    else:
        return ps.LightCurve(band, time, mags)


class TestLightCurve(unittest.TestCase):
    def test_BandName(self):
        band = 'U'
        lc = lc_create(band)
        self.assertEqual(band, lc.Band.Name,
                         msg="It should be equal band names.\n \
                         Now band is %s but  lc.Band.Name is  %s." % (band, lc.Band.Name))

    def test_LC_interp(self):
        lc = lc_create('U', dt=0.)

        time = np.linspace(10, 50, 5)
        lc_interp = ps.rf.lc.LC_interp(lc, time)

        self.assertEqual(len(time), lc_interp.Length, msg='The lenght of Interp LC should be equal len(time)')

        # plot
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots(1, 1)
        ax.plot(lc.Time, lc.Mag, label='Origin')
        ax.plot(lc_interp.Time, lc_interp.Mag, marker='o', label='Interpolated')
        ax.invert_yaxis()
        ax.legend()
        plt.show()

    def test_lc_leastsq(self):
        dt_init = 10.
        lc1 = lc_create('U', dt=0.)
        lc2 = lc_create('U', dt=0.)

    def test_lc_bol(self):
        import matplotlib.pyplot as plt
        from scipy.integrate import simps

        m1 = ps.Stella('cat_R500_M15_Ni006_E12', path='data/stella')
        tt1 = m1.get_tt().read()
        curves = m1.curves(bands=['bol'])
        ax = ps.light_curve_plot.curves_plot(curves, xlim=(-10, 155), ylim=(-9, -20), is_line=False)
        t = tt1['time']
        ax.plot(t, tt1['Mbol'], label='tt-bolometric LC ', color='red', lw=2, ls=':')
        # ph
        if False:
            ph = m1.get_ph()
            m_bol = []
            for spec in ph:
                lum = simps(spec.Flux[::-1], spec.Freq[::-1])
                bol = 4.75 - 2.5 * np.log10(np.abs(lum)/3.86e33)
                m_bol.append(bol)
            ax.plot(ph.Time, m_bol, label='ph-bolometric LC ', color='green', lw=2, ls='-.')
        ax.legend()
        plt.show()


class TestSetLightCurve(unittest.TestCase):
    def test_SetLightCurve_BandNames(self):
        bands = ['U', 'B', 'V']
        curves = ps.SetLightCurve()
        for b in bands:
            curves.add(lc_create(b))

        self.assertCountEqual(bands, curves.BandNames,
                              msg="Error for band names.\n \
                Now band is %s but  lc.Band.Name is  %s." % (' '.join(bands), ' '.join(curves.BandNames)))

    def test_SetLightCurve_save_true(self):
        bands = ['U', 'B', 'V']
        curves = ps.SetLightCurve()
        for b in bands:
            curves.add(lc_create(b))
        res = ps.curves_save(curves, 'tmp_curves')
        self.assertTrue(res, msg="Error:  curves_save should return True")

    def test_SetLightCurve_save_read(self):
        bands = ['U', 'B', 'V']
        curves = ps.SetLightCurve()
        for b in bands:
            curves.add(lc_create(b))
        ps.curves_save(curves, 'tmp_curves')
        read = ps.curves_read('tmp_curves')
        self.assertTrue((np.array(curves.BandNames == read.BandNames)).all(),
                        msg="Error for the initial band names [%s] "
                            "VS secondary BandNames are %s."
                            % (' '.join(curves.BandNames), ' '.join(read.BandNames)))
        self.assertTrue(np.allclose(curves.TimeCommon, read.TimeCommon),
                        msg="Error for the initial TimeCommon of Bands.\n \
                                  Now band were %s but BandNames are %s."
                            % (' '.join(curves.BandNames), ' '.join(read.BandNames)))
        # todo correct testing
        # self.assertSequenceEqual(curves.TimeCommon, read.TimeCommon, msg="The time columns are not equal")

    def test_SetLightCurve_save_true_with_errors(self):
        bands = ['U', 'B', 'V']
        curves = ps.SetLightCurve()
        for b in bands:
            curves.add(lc_create(b, is_err=True))
        curves.add(lc_create('I'))
        res = ps.curves_save(curves, 'tmp_curves')
        self.assertTrue(res, msg="Error:  curves_save should return True")

    def test_SetLightCurve_save_NoIsCommonTime(self):
        bands = ['U', 'B', 'V']
        curves = ps.SetLightCurve()
        for b in bands:
            curves.add(lc_create(b))
        curves.add(lc_create('TimeDif', dt=1.))

        res = ps.curves_save(curves, 'tmp_curves_2')
        self.assertTrue(res, msg="Error:  curves_save should return True")


def main():
    unittest.main()


if __name__ == '__main__':
    main()
