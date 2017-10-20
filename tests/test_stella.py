import unittest
from os.path import join, dirname, abspath
import numpy as np

from pystella.model.stella import Stella
import pystella.rf.light_curve_func as lcf
import pystella.rf.light_curve_plot as lcp
from pystella.rf.reddening import ReddeningLaw

__author__ = 'bakl'


class TestStellaLightCurves(unittest.TestCase):
    def test_stella_curves(self):
        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        bands = ('U', 'B', 'V')

        mdl = Stella(name, path=path)
        curves = mdl.curves(bands)

        self.assertTrue((np.array(sorted(curves.BandNames) == sorted(bands))).all(),
                        msg="Error for the initial band names [%s] "
                            "VS secondary BandNames are %s."
                            % (' '.join(bands), ' '.join(curves.BandNames)))

    def test_stella_curves_reddening(self):
        import matplotlib.pyplot as plt
        from matplotlib import gridspec

        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        bands = ('UVW1', 'UVW2', 'UVM2', 'U', 'B', 'R', 'I')
        ebv = 1

        # mags reddening
        cs = lcf.curves_compute(name, path, bands)
        curves_mags = lcf.curves_reddening(cs, ebv=ebv)

        mdl = Stella(name, path=path)
        curves = mdl.curves(bands, ebv=ebv)  # best SML
        # curves = mdl.curves(bands, ebv=ebv, law=ReddeningLaw.SMC)  # best SML

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
        ax = lcp.curves_plot(curves_mags, ax=axUbv, lt=lt, markersize=2, is_legend=False)
        xlim = ax.get_xlim()

        lcp.curves_plot(curves, ax=axUbv)

        x = curves.TimeCommon
        for b in bands:
            y = curves.get(b).Mag - curves_mags.get(b).Mag
            axDM.plot(x, y, label="Delta {}".format(b))
        axDM.set_xlim(xlim)
        axDM.legend()

        plt.grid(linestyle=':', linewidth=1)
        plt.show()
