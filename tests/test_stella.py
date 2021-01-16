import unittest
from os.path import join, dirname, abspath

import matplotlib.pyplot as plt
import numpy as np

import pystella as ps

__author__ = 'bakl'


class TestStellaLightCurves(unittest.TestCase):
    def test_stella_curves(self):
        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        bands = ('U', 'B', 'V')

        mdl = ps.Stella(name, path=path)
        curves = mdl.curves(bands)

        print(ps.first(curves).Time[:300])
        self.assertTrue((np.array(sorted(curves.BandNames) == sorted(bands))).all(),
                        msg="Error for the initial band names [%s] "
                            "VS secondary BandNames are %s."
                            % (' '.join(bands), ' '.join(curves.BandNames)))

    def test_stella_curves_tbeg(self):
        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        bands = ('U', 'B', 'V')

        mdl = ps.Stella(name, path=path)
        t_beg = 1.
        curves = mdl.curves(bands, t_beg=t_beg)

        print(ps.first(curves).Time[:3])
        self.assertTrue(np.any(ps.first(curves).Time >= t_beg),
                        msg="There ara values Time less then t_beg = {}".format(t_beg))

    def test_stella_curves_VS_tt_plot(self):
        from pystella.rf.band import colors
        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        bands = ('U', 'B', 'V', 'R', 'I')

        mdl = ps.Stella(name, path=path)
        curves = mdl.curves(bands)

        tt = mdl.get_tt().load()

        ax = ps.lcp.curves_plot(curves)
        for bname in bands:
            ax.plot(tt['time'], tt['M'+bname], label="tt "+bname, color=colors(bname),
                    marker='*', markersize=3, ls='')
        ax.legend()

        plt.grid(linestyle=':', linewidth=1)
        plt.show()

    def test_stella_curves_reddening_plot(self):
        from matplotlib import gridspec

        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        bands = ('UVW1', 'UVW2', 'UVM2')
        # bands = ('UVW1', 'UVW2', 'UVM2', 'U', 'B', 'R', 'I')
        ebv = 1

        # mags reddening
        cs = ps.lcf.curves_compute(name, path, bands, t_diff=1.05)

        mdl = ps.Stella(name, path=path)
        is_SMC = False
        if is_SMC:
            curves_mags = ps.lcf.curves_reddening(cs, ebv=ebv, law='Rv2.1')
            curves = mdl.curves(bands, ebv=ebv, t_diff=1.05, mode=ps.ReddeningLaw.SMC)  # best SMC MW
        else:
            curves_mags = ps.lcf.curves_reddening(cs, ebv=ebv, law=ps.extinction.law_default)
            curves = mdl.curves(bands, ebv=ebv, t_diff=1.05,  mode=ps.ReddeningLaw.MW)
        # curves = mdl.curves(bands, ebv=ebv, law=LawFitz, mode=ReddeningLaw.SMC)  # best SMC

        self.assertTrue((np.array(sorted(curves.BandNames) == sorted(curves_mags.BandNames))).all(),
                        msg="Error for the initial band names [%s] "
                            "VS secondary BandNames are %s."
                            % (' '.join(curves_mags.BandNames), ' '.join(curves.BandNames)))

        # plot reddening with mags
        fig = plt.figure(figsize=(12, 12))
        gs1 = gridspec.GridSpec(4, 1)
        axUbv = fig.add_subplot(gs1[:-1, 0])
        axDM = fig.add_subplot(gs1[3, 0])
        lt = {lc.Band.Name: 'o' for lc in curves_mags}
        ax = ps.lcp.curves_plot(curves_mags, ax=axUbv, lt=lt, markersize=2, is_legend=False)
        xlim = ax.get_xlim()

        ps.lcp.curves_plot(curves, ax=axUbv)

        x = curves.TimeCommon
        for b in bands:
            y = curves.get(b).Mag - curves_mags.get(b).Mag
            axDM.plot(x, y, label="Delta {}".format(b))
        axDM.set_xlim(xlim)
        axDM.legend()

        plt.grid(linestyle=':', linewidth=1)
        plt.show()
