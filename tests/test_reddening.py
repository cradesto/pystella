# coding=utf-8
import unittest

import numpy as np

from pystella.rf import extinction, band
from pystella.rf.lc import LightCurve
from pystella.rf.light_curve_func import curves_reddening, flux_reddening_wl, flux_reddening

__author__ = 'bakl'


def lc_create(bname, m=-19, dt=0.):
    n = 10
    time = np.linspace(0.+dt, 200.+dt, n)
    mags = m * np.ones(n)
    return LightCurve(bname, time, mags)


class TestReddening(unittest.TestCase):
    def setUp(self):
        band.Band.load_settings()

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
        bname, mshift = 'U', 4.334
        ebv = 1.
        lc = lc_create(bname)
        magExpect = lc.Mag + mshift  # Av for e=1.
        lc = curves_reddening(lc, ebv=ebv)
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
        bname, mshift = 'V', 4.334  # shifted to U
        z = 5421.7/3508.2 - 1.
        ebv = 1.
        lc = lc_create(bname)
        magExpect = lc.Mag + mshift  # Av for e=1.
        lc = curves_reddening(lc, ebv=ebv, z=z)
        magRed = lc.Mag
        res = list(filter(lambda x: abs(x) > 1e-4, magExpect - magRed))
        self.assertTrue(len(res) == 0,
                        msg="Some mags were shifted more then %f." % mshift)

    def test_ReddeningLaw_xi(self):
        from pystella.rf.reddening import ReddeningLaw
        import matplotlib.pyplot as plt

        wl = np.logspace(3.0, np.log10(5e4), num=100)
        for law in ReddeningLaw.Modes:
            xi = ReddeningLaw.xi(wl, law=law)
            x = 1e4 / wl
            plt.plot(x, xi, label=ReddeningLaw.Names[law])
        plt.legend()
        plt.xlabel(r'$1/\lambda\; [\mu m^{-1}]$')
        plt.ylabel(r'$\xi(\lambda) = A_\lambda / A_V$')
        plt.grid(linestyle=':', linewidth=1)
        plt.show()

    def test_LawPei_Almd(self):
        from pystella.rf.reddening import LawPei
        import matplotlib.pyplot as plt

        wl = np.logspace(3.0, np.log10(5e4), num=100)
        ebv = 0.074
        for mode in LawPei.Modes:
            Almd = LawPei.Almd(wl, ebv, Rv=LawPei.Rv[mode])
            x = 1e4 / wl
            plt.plot(x, Almd, label=LawPei.Names[mode])
        plt.legend()
        plt.xlabel(r'LawPei: $1/\lambda\; [\mu m^{-1}]$')
        plt.ylabel(r'$A_\lambda$')
        plt.grid(linestyle=':', linewidth=1)
        plt.show()

    def test_LawFitz_Almd(self):
        from pystella.rf.reddening import LawFitz
        import matplotlib.pyplot as plt

        wl = np.logspace(3.0, np.log10(5e4), num=100)
        ebv = 0.074
        for mode in LawFitz.Modes:
            Almd = LawFitz.Almd(wl, ebv, Rv=LawFitz.Rv[mode])
            x = 1e4 / wl
            plt.plot(x, Almd, label=LawFitz.Names[mode])
        plt.legend()
        plt.xlabel(r'LawFitz: $1/\lambda\; [\mu m^{-1}]$')
        plt.ylabel(r'$A_\lambda$')
        plt.grid(linestyle=':', linewidth=1)
        plt.show()

    def test_flux_freq_with_ReddeningLaw(self):
        import matplotlib.pyplot as plt
        from pystella.rf.reddening import ReddeningLaw
        import pystella.rf.rad_func as rf
        from pystella.util.phys_var import phys

        wl = np.logspace(3.0, np.log10(5e4), num=100)
        freq = phys.c / (wl * phys.angs_to_cm)
        T = 15500.
        flux = rf.planck(freq, temperature=T)

        ebv = 0.074
        plt.plot(freq, flux, color='black', label='Origin planck T={}'.format(T))
        for mode in ReddeningLaw.Modes:
            x = freq
            y = flux_reddening(freq, flux, ebv, mode=mode)
            plt.plot(x, y, label=ReddeningLaw.Names[mode])
        plt.legend()
        plt.xscale("log", nonposx='clip')
        plt.yscale("log", nonposy='clip')
        plt.xlabel(r'$Frequency\; [Hz]$')
        plt.ylabel(r'$F_\nu\; [ergs s^{-1} cm^{-2} Hz^{-1}]$')
        plt.grid(linestyle=':', linewidth=1)
        plt.show()

    def test_flux_wl_with_ReddeningLaw(self):
        import matplotlib.pyplot as plt
        from pystella.rf.reddening import ReddeningLaw
        import pystella.rf.rad_func as rf
        from pystella.util.phys_var import phys

        wl = np.logspace(3.0, np.log10(5e4), num=100)
        freq = phys.c / (wl * phys.angs_to_cm)
        T = 15500.
        flux = rf.planck(freq, temperature=T)
        flux_wl = rf.Fnu2Fwl(freq, flux)
        # flux_wl = flux * freq**2 / phys.c

        ebv = 0.074
        plt.plot(wl, flux_wl, color='black', label='Origin planck T={}'.format(T))
        for mode in ReddeningLaw.Modes:
            x = wl
            y = flux_reddening_wl(wl, flux_wl, ebv, mode=mode)
            plt.plot(x, y, label=ReddeningLaw.Names[mode])
        plt.legend()
        plt.xscale("log", nonposx='clip')
        plt.yscale("log", nonposy='clip')
        plt.xlabel(r'$\lambda\; [A]$')
        plt.ylabel(r'$F_\lambda\; [ergs s^{-1} cm^{-2} A^{-1}]$')
        plt.grid(linestyle=':', linewidth=1)
        plt.show()

    def test_flux_ones_with_LawFitz(self):
        import matplotlib.pyplot as plt
        from pystella.rf.reddening import LawFitz
        from pystella.util.phys_var import phys

        wl = np.logspace(3.0, np.log10(5e4), num=100)
        freq = phys.c / (wl * phys.angs_to_cm)
        flux = np.ones_like(freq)
        # flux = np.ones(len(freq))
        ebv = 1.
        x = 1e4 / wl
        plt.plot(x, flux, color='black', label='Origin flux')
        law = LawFitz
        for mode in law.Modes:
            y = flux_reddening_wl(wl, flux, ebv, law=law, mode=mode)
            plt.plot(x, y, label='Flux {}'.format(law.Names[mode]))
            Almd = law.Almd(wl, ebv, Rv=law.Rv[mode])
            yy = 10.**(-0.4*Almd)
            plt.plot(x, yy, label=r'$A_\lambda$ {}'.format(law.Names[mode]),
                     marker='o', ls='', markersize='1')

        plt.legend()
        plt.xscale("log", nonposx='clip')
        plt.xlabel(r'$1/\lambda\; [\mu m^{-1}]$')
        plt.ylabel(r'$F_\nu\; [ergs s^{-1} cm^{-2} Hz^{-1}]$')
        plt.grid(linestyle=':', linewidth=1)
        plt.show()

    def test_band_EXT_LAWS(self):
        import matplotlib.pyplot as plt
        from pystella.rf.reddening import LawPei, LawFitz

        # ebv = .074
        ebv = 1.
        data, laws = extinction.read_cache_data()
        x = data['lambdaeff']
        sorti = np.argsort(x)
        x = x[sorti]

        lawBands = ['Rv2.1', 'Rv3.1', 'Rv4.1', 'Rv5.1']
        for lawBand in lawBands:
            # lawBand = 'Rv2.1'
            # lawBand = 'Rv3.1'
            # lawBand = 'Rv4.1'
            # lawBand = 'Rv5.1'

            y = data[lawBand]
            # sort
            y = y[sorti] * ebv

            plt.plot(x, y, label='Bands {}'.format(lawBand), ls='--', marker='o', markersize=3)

        markers = ['<', '>', 'd']
        for i, mode in enumerate(LawPei.Modes):
            Rv = LawPei.Rv[mode]
            # nm = LawPei.Names[law]
            yy = LawPei.Almd(x, ebv, Rv=Rv)
            # yy = LawPei.Almd(x, ebv, Rv=Rv, law=LawPei.MW)
            plt.plot(x, yy, label=r'$A_\lambda$: Pei Rv={}'.format(Rv),
                     marker=markers[i], ls='', markersize=2)

        for i, mode in enumerate(LawFitz.Modes):
            Rv = LawFitz.Rv[mode]
            # nm = LawFitz.Names[law]
            yy = LawFitz.Almd(x, ebv, Rv=Rv)
            plt.plot(x, yy, label=r'$A_\lambda$: Fitz Rv={}'.format(Rv),
                     marker=markers[i], ls='', markersize=2)
        # Rv = LawFitz.Rv[LawFitz.MW]
        # yy = LawFitz.Almd(x, ebv, Rv=Rv, law=LawFitz.MW)
        # plt.plot(x, yy, label=r'$A_\lambda$ {}'.format(LawFitz.Names[LawFitz.MW]),
        #          marker='d', ls=':', markersize=2)

        # yy = ReddeningLaw.Almd(x, ebv, law=ReddeningLaw.LMC)
        # plt.plot(x, yy, label=r'$A_\lambda$ {}'.format(ReddeningLaw.Names[ReddeningLaw.LMC]),
        #          marker='<', ls=':', markersize=2)
        #
        # yy = ReddeningLaw.Almd(x, ebv, law=ReddeningLaw.SMC)
        # plt.plot(x, yy, label=r'$A_\lambda$ {}'.format(ReddeningLaw.Names[ReddeningLaw.SMC]),
        #          marker='>', ls=':', markersize=2)

        plt.legend()
        plt.xscale("log", nonposx='clip')
        # plt.xlabel(r'$1/\lambda\; [\mu m^{-1}]$')
        # plt.ylabel(r'$F_\nu\; [ergs s^{-1} cm^{-2} Hz^{-1}]$')
        plt.grid(linestyle=':', linewidth=1)
        plt.show()


def main():
    unittest.main()


if __name__ == '__main__':
    main()
