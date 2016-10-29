import matplotlib.pyplot as plt
import re
import unittest
from os.path import dirname, abspath, join

import numpy as np

from pystella.model.stella import Stella

from pystella.rf import light_curve_func as lcf

__author__ = 'bakl'


class TestStellaTt(unittest.TestCase):
    def setUp(self):
        # name = 'ccsn2007bi1dNi6smE23bRlDC'
        name = 'cat_R500_M15_Ni006_E12'
        # name = 'cat_R1000_M15_Ni007_E15'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        self.stella = Stella(name, path=path)
        self.tt = self.stella.get_tt()

    def test_info_parse(self):
        tt_header = """
                              <===== HYDRODYNAMIC RUN OF MODEL cat_R1000_M15_Ni007_E15.prf                                                     =====>

                              MASS(SOLAR)=15.000       RADIUS(SOLAR)= 1000.000
                              EXPLOSION ENERGY(10**50 ERG)= 1.50000E+01

                              <=====                                 =====>


  INPUT PARAMETERS
 EPS   =         0.00300          Rce   =     1.00000E-02 SOL.Rad.
 HMIN  =     1.00000E-11 S        AMht  =     1.00000E-02 SOL.MASS
 HMAX  =     5.00000E+04 S        Tburst=     1.00000E-01 S
 THRMAT=     1.00000E-30          Ebstht=     1.50000E+01 1e50 ergs
 METH  =               3          CONV  =               F
 JSTART=               0          EDTM  =               T
 MAXORD=               4          CHNCND=               T
 NSTA  =           -4500          FitTau=     3.00000E+02
 NSTB  =              -1          TauTol=     1.30000E+00
 NOUT  =             100          IOUT  =              -1
 TcurA =         0.00000          Rvis   =        1.00000
 TcurB =       200.00000          BQ    =     1.00000E+00
 NSTMAX=          360000          DRT   =     1.00000E-01
 XMNI  =     7.00000E-02 SOL.MASS NRT   =               1
 XNifor=     1.16561E-01
 MNicor=     1.16999E-01 SOL.MASS SCAT  =               T

               """
        pattern = filter(lambda x: len(x) > 0, tt_header.splitlines())
        pattern = map(str.strip, pattern)
        p = r"(.*?)\s*=\s+([-+]?\d*\.\d+|\d+)"
        for line in pattern:
            res = re.findall(p, line)
            if len(res) > 0:
                for k, v in res:
                    print "key: %s  v: %f " % (k, float(v))

    def test_info_parse(self):
        info = self.tt.Info.parse()
        info.show()

        tmp = 1000.
        self.assertEquals(info.R, tmp, "Radius [%f] should be %f" % (info.R, tmp))

        tmp = 15.
        self.assertEquals(info.M, tmp, "Mass [%f] should be %f" % (info.M, tmp))

        tmp = 15.
        self.assertEquals(info.E, tmp, "Ebstht [%f] should be %f" % (info.E, tmp))

    def test_tt_vs_ph(self):
        curves_tt = self.tt.read_curves_tt()
        bands = curves_tt.BandNames

        serial_spec = self.stella.read_series_spectrum(t_diff=1.00001)
        curves_ph = serial_spec.flux_to_curves(bands)

        models_dic = {'tt': curves_tt, 'ph': curves_ph}
        lc_types = {'tt': '--', 'ph': '-'}

        plt.matplotlib.rcParams.update({'font.size': 14})
        fig, ax = plt.subplots(1, 1)

        lcf.plot_models_curves_fixed_bands(ax, models_dic, bands, lc_types=lc_types, ylim=(-10, -24))
        # lcf.plot_models_curves(ax, models_dic, bands, lc_types=lc_types, ylim=(-10, -24), xlim=(0, 20))
        plt.legend()
        plt.show()

    def test_gri_vs_ph(self):
        curves_gri = self.tt.read_curves_gri()
        bands = curves_gri.BandNames
        # bands = ('J','H','K')

        serial_spec = self.stella.read_series_spectrum(t_diff=1.)
        curves_ph = serial_spec.flux_to_curves(bands)

        models_dic = {'gri': curves_gri, 'ph': curves_ph}
        lc_types = {'gri': '--', 'ph': ':'}

        plt.matplotlib.rcParams.update({'font.size': 14})
        fig, ax = plt.subplots(1, 1)

        # lcf.plot_models_curves_fixed_bands(ax, models_dic, bands=('r','B','V'), lc_types=lc_types, ylim=(-13, -23), lw=3)
        lcf.plot_models_curves(ax, models_dic, lc_types=lc_types, ylim=(-10, -23), lw=3)
        plt.legend()
        plt.show()

    def test_tt_vs_gri_vs_ph(self):
        curves_tt = self.tt.read_curves_tt()
        bands_tt = curves_tt.BandNames

        curves_gri = self.tt.read_curves_gri()
        bands_gri = curves_gri.BandNames

        bands = np.unique(np.array(bands_tt + bands_gri))

        serial_spec = self.stella.read_series_spectrum(t_diff=1.00001)
        curves_ph = serial_spec.flux_to_curves(bands)

        models_dic = {'tt': curves_tt, 'gri': curves_gri, 'ph': curves_ph}
        lc_types = {'tt': ':', 'gri': '--', 'ph': '-'}

        plt.matplotlib.rcParams.update({'font.size': 14})
        fig, ax = plt.subplots(1, 1)

        lcf.plot_models_curves(ax, models_dic, lc_types=lc_types, ylim=(-10, -19), lw=3)
        # lcf.plot_models_curves_fixed_bands(ax, models_dic, bands=('B', 'V'), lc_types=lc_types, ylim=(-13, -23), lw=3)
        plt.legend()
        plt.show()


def main():
    unittest.main()


if __name__ == '__main__':
    main()
