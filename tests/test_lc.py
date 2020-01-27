import numpy as np
import unittest

import pystella as ps

__author__ = 'bakl'


def lc_create(bname, m=-19, tbeg=0., tend=200., n=10, is_err=False):
    time = np.linspace(0. + tbeg, tend + tbeg, n)
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
        lc = lc_create('U', tbeg=0.)

        time = np.linspace(10, 50, 5)
        lc_interp = ps.rf.lc.LC_interp(lc, time)

        self.assertEqual(len(time), lc_interp.Length, msg='The length of Interp LC should be equal len(time)')

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
        lc1 = lc_create('U', tbeg=0.)
        lc2 = lc_create('U', tbeg=0.)

    def test_lc_copy(self):
        ps.band.Band.load_settings()
        lc1 = lc_create('U', tbeg=0.)
        lc2 = lc1.copy()
        self.assertEqual(lc1.Length, lc2.Length, msg='The length of copy  should be equal the original length')
        self.assertEqual(lc1.Band.Name, lc2.Band.Name, msg='The lc of copy  should be equal the original lc')
        np.testing.assert_array_equal(lc1.Time, lc2.Time)
        np.testing.assert_array_equal(lc1.Mag, lc2.Mag)

    def test_lc_copy_filter(self):
        tlim = (10., 99.)
        ps.band.Band.load_settings()
        # Time
        lc1 = lc_create('V', m=-19, tbeg=1., tend=200., n=10,  is_err=False)
        lc2 = lc1.copy(f=lambda x: (tlim[0] <= x.Time) & (x.Time <= tlim[1]))
        self.assertGreater(lc1.Length, lc2.Length, msg='The length of copy  should be equal the original length')
        self.assertEqual(lc1.Band.Name, lc2.Band.Name, msg='The lc of copy  should be equal the original lc')
        self.assertTrue(np.any(lc2.Time >= tlim[0]), msg='The lc.Time should be greater the lower limit')
        self.assertTrue(np.any(lc2.Time <= tlim[1]), msg='The lc.Time should be smaller the lower limit')
        # Magnitude
        maglim = (-18., -10.)
        lc3 = lc1.copy(f=lambda x: (maglim[0] <= x.Mag) & (x.Mag <= maglim[1]))
        self.assertGreater(lc1.Length, lc3.Length, msg='The length of copy  should be equal the original length')
        self.assertEqual(lc1.Band.Name, lc3.Band.Name, msg='The lc of copy  should be equal the original lc')
        self.assertTrue(np.any(lc3.Mag >= maglim[0]), msg='The lc.Mag should be greater the lower limit')
        self.assertTrue(np.any(lc3.Mag <= maglim[1]), msg='The lc.Mag should be smaller the lower limit')

    def test_lc_clone(self):
        lc1 = lc_create('U', tbeg=0.)
        lc2, tshift, mshift = lc1.clone()
        self.assertEqual(lc1.Length, lc2.Length, msg='The length of clone  should be equal the original length')
        self.assertEqual(lc1.Band.Name, lc2.Band.Name, msg='The band of clone  should be equal the original band')
        np.testing.assert_array_equal(lc1.Time, lc2.Time)
        np.testing.assert_array_equal(lc1.Mag, lc2.Mag)

    def test_lc_clone_add_err(self):
        lc1 = lc_create('U', tbeg=0.)
        err = [1] * lc1.Length
        lc2, tshift, mshift = lc1.clone(err=err)
        self.assertEqual(lc1.Length, lc2.Length, msg='The length of clone  should be equal the original length')
        np.testing.assert_array_equal(err, lc2.MagErr)
        np.testing.assert_array_equal(lc1.Mag, lc2.Mag)

    def test_lc_bol(self):
        import matplotlib.pyplot as plt
        from scipy.integrate import simps

        m1 = ps.Stella('cat_R500_M15_Ni006_E12', path='data/stella')
        curves = m1.curves(bands=['bol'], t_diff=1.0000001)
        # ax = ps.light_curve_plot.curves_plot(curves, xlim=(0.7, 1), ylim=(-14, -24), is_line=False)
        ax = ps.lcp.curves_plot(curves, xlim=(-10, 155), ylim=(-14, -24), is_line=False)
        # tt
        tt1 = m1.get_tt().load()
        t = tt1['time']
        ax.plot(t, tt1['Mbol'], label='tt-bolometric LC ', color='red', lw=2, ls=':')
        # ph
        if True:
            ph = m1.get_ph()
            m_bol = []
            for t, spec in ph:
                lum = simps(spec.Flux[::-1], spec.Freq[::-1])
                bol = 4.75 - 2.5 * np.log10(np.abs(lum) / 3.86e33)
                m_bol.append(bol)
            ax.plot(ph.Time, m_bol, label='ph-bolometric LC ', color='green', lw=2, ls='-.')
        ax.legend()
        plt.show()
        import warnings
        warnings.warn("Should be check for shorck breakout")

    def test_bol_Uni(self):
        import matplotlib.pyplot as plt
        m1 = ps.Stella('cat_R500_M15_Ni006_E12', path='data/stella')
        fig, ax = plt.subplots()

        # Bol
        curves1 = m1.curves(bands=['bol'], wlrange=(1e0, 42.), is_nfrus=False)
        for lc in curves1:
            color = 'blue'
            ax.plot(lc.Time, lc.Mag, label=lc.Band.Name, color=color, linewidth=2, ls='--')

        band03kEv = ps.BandUni(name='bol', wlrange=(1e0, 42.), length=300)
        wl_ab = np.min(band03kEv.wl2args), np.max(band03kEv.wl2args)
        curves2 = m1.curves(bands=[band03kEv], is_nfrus=False, wl_ab=wl_ab)
        for lc in curves2:
            color = 'red'
            ax.plot(lc.Time, lc.Mag, label=lc.Band.Name, color=color, linewidth=2, ls=':')

        ax.invert_yaxis()
        #
        ax.legend()
        # ax.set_ylim(-14, -24)
        plt.show()
        import warnings
        warnings.warn("Should be check for shorck breakout")


class TestSetLightCurve(unittest.TestCase):
    def test_SetLightCurve_BandNames(self):
        bands = ['U', 'B', 'V']
        curves = ps.SetLightCurve()
        for b in bands:
            curves.add(lc_create(b))

        self.assertCountEqual(bands, curves.BandNames,
                              msg="Error for band names.\n \ Now band is %s but  lc.Band.Name is  %s."
                                  % (' '.join(bands), ' '.join(curves.BandNames)))

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
        curves.add(lc_create('TimeDif', tbeg=1.))

        res = ps.curves_save(curves, 'tmp_curves_2')
        self.assertTrue(res, msg="Error:  curves_save should return True")

    def test_SetLightCurve_copy_tmlim(self):
        ps.band.Band.load_settings()
        bands = ['U', 'B', 'V']
        curves = ps.SetLightCurve()
        for b in bands:
            curves.add(lc_create(b, m=-19, tbeg=0., tend=200., n=10, is_err=False))
        curves.add(lc_create('R', tbeg=1.))

        tlim = (10, 99)
        mlim = (10, -18)
        curves_cut = curves.copy_tmlim(tlim=tlim, mlim=mlim)
        self.assertTrue(curves_cut.TimeMin >= tlim[0])
        self.assertTrue(curves_cut.TimeMax <= tlim[1])

    def test_SetLightCurve_clone_add_err(self):
        bands = ['U', 'B', 'V']
        bname = bands[1]
        curves = ps.SetLightCurve()
        for b in bands:
            curves.add(lc_create(b))
        curves.add(lc_create('TimeDif', tbeg=1.))
        lc = curves[bname]
        # Time
        t = np.ones(lc.Length)
        curves_clone = curves.clone(t=t)
        self.assertEqual(curves_clone.Length, curves.Length,
                         msg=f'The length of clone{curves_clone.Length}  should be equal the original length {curves.Length}')
        lc_clone = curves_clone[bname]
        np.testing.assert_array_equal(t, lc_clone.Time)

        # Mag
        mag = np.ones(lc.Length)
        curves_clone = curves.clone(m=mag)
        self.assertEqual(curves_clone.Length, curves.Length,
                         msg=f'The length of clone{curves_clone.Length}  should be equal the original length {curves.Length}')
        lc_clone = curves_clone[bname]
        np.testing.assert_array_equal(mag, lc_clone.Mag)
        # Err
        err = np.ones(lc.Length)
        curves_clone = curves.clone(err=err)
        self.assertEqual(curves_clone.Length, curves.Length,
                         msg=f'The length of clone{curves_clone.Length}  should be equal the original length {curves.Length}')
        lc_clone = curves_clone[bname]
        np.testing.assert_array_equal(err, lc_clone.MagErr)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
