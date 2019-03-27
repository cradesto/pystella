import numpy as np
import unittest

import pystella.rf.rad_func as rf
import pystella.rf.spectrum as spectrum
from pystella.util.phys_var import phys
from pystella.rf.rad_func import kcorrection, kcorrection_spec
from pystella.rf.band import band_is_exist, band_by_name

__author__ = 'bakl'


class TestSpec(unittest.TestCase):
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

    def test_spec_kcorr_z0(self):
        nf, start, end = 100, 10., 1e5
        wl = np.exp(np.linspace(np.log(end), np.log(start), nf))
        freq = rf.val_to_hz(wl, inp="A")
        Tcolor = 7e3
        sp = spectrum.SpectrumPlanck(freq, Tcolor, 'test')

        z = 0.
        bn_rest = 'V'
        # bn_obs = 'V'
        for bn_rest in ['U', 'B', 'V', 'R', 'g']:
            bn_obs = bn_rest
            br = band_by_name(bn_rest)
            bo = band_by_name(bn_obs)

            kcorr = kcorrection_spec(sp, z, br, bo)

            self.assertAlmostEqual(kcorr, 0., 8, "{}: For z={} kcorr should be 0. But kcorr={:.5f}.".
                                   format(bn_rest, z, kcorr))

    def test_spec_kcorr_with_curve(self):
        nf, start, end = 100, 10., 1e5
        wl = np.exp(np.linspace(np.log(end), np.log(start), nf))
        freq = rf.val_to_hz(wl, inp="A")
        Tcolor = 7e3
        sp = spectrum.SpectrumPlanck(freq, Tcolor, 'test')

        z = 0.1
        distance = 1e6  # phys.cosmology_D_by_z(z)
        bn_rest = 'g'
        # bn_obs = 'V'
        bn_obs = bn_rest
        br = band_by_name(bn_rest)
        bo = band_by_name(bn_obs)

        M_Q = sp.to_mag(br, z=0., d=phys.pc2cm(10.))
        m_r = sp.to_mag(bo, z=z, d=phys.pc2cm(distance))
        DM = phys.dist2MD(distance)
        K_QR = m_r - M_Q - DM

        kcorr = kcorrection_spec(sp, z, br, bo)

        print("{}: z= {} K_QR= {:.5f} kcorr= {:.5f}".format(bn_rest, z, K_QR, kcorr))

        self.assertAlmostEqual(K_QR, kcorr, 8, "{}: For z={} kcorr should be like K_QR. But kcorr={:.5f} K_QR= {:.5f}.".
                               format(bn_rest, z, kcorr, K_QR))

    def test_spec_kcorr_positive(self):
        nf = 300
        start, end = 100., 1e5
        wl = np.exp(np.linspace(np.log(end), np.log(start), nf))
        freq = rf.val_to_hz(wl, inp="A")
        Tcolor = 7e3
        sp = spectrum.SpectrumPlanck(freq, Tcolor, 'test')

        z = 0.1
        for bn_rest in ['U', 'B', 'V', 'R', 'u', 'g', 'r']:
            bn_obs = bn_rest
            br = band_by_name(bn_rest)
            bo = band_by_name(bn_obs)

            kcorr = kcorrection_spec(sp, z, br, bo)
            print('{}: z= {}  kcorr= {}'.format(bn_rest, z, kcorr))
            self.assertGreater(kcorr, 0., "{}: For z={} kcorr should be > 0. But kcorr={:.5f}.".
                               format(bn_rest, z, kcorr))

    def test_spec_kcorr_zplot(self):
        import matplotlib.pyplot as plt

        nf, start, end = 100, 10., 1e5
        wl = np.exp(np.linspace(np.log(end), np.log(start), nf))
        freq = rf.val_to_hz(wl, inp="A")
        Tcolor = 7.e3
        sp = spectrum.SpectrumPlanck(freq, Tcolor, 'test')

        bn_rest = 'B'
        bn_obs = bn_rest
        br = band_by_name(bn_rest)
        bo = band_by_name(bn_obs)

        z_arr = np.linspace(0, 0.3, 33)
        kcoors = []
        for z in z_arr:
            k = kcorrection_spec(sp, z, br, bo)
            kcoors.append(k)

        fig, ax = plt.subplots(figsize=(4, 16))
        ax.plot(z_arr, kcoors, marker='o', ls='', markersize=3)
        # ax.set_ylim(-0., .3)
        ax.grid()
        ax.set_title("{} band, Tcol= {:.0f} K".format(bn_rest, Tcolor))
        ax.set_xlabel('Redshift')
        ax.set_ylabel('k-correction')
        plt.show()
        fig.savefig('k-correction.pdf')


def main():
    unittest.main()


if __name__ == '__main__':
    main()
