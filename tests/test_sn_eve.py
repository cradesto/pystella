import os
from os.path import dirname, abspath, join
import unittest
import matplotlib.pyplot as plt

import numpy as np

from pystella.model import snec
import pystella as ps
from pystella.model import sn_eve
from pystella.model.sn_eve import PreSN

__author__ = 'bakl'


class SnEveTests(unittest.TestCase):
    def test_rho_load(self):
        name = 'cat_R1000_M15_Ni007'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        rho_file = os.path.join(path, name + '.rho')

        presn = sn_eve.load_rho(rho_file)
        self.assertTrue(presn.nzon > 0, "Rho-file have not been loaded: %s" % rho_file)

    def test_rho_el(self):
        name = 'cat_R1000_M15_Ni007'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        rho_file = os.path.join(path, name + '.rho')

        presn = sn_eve.load_rho(rho_file)
        h = presn.lg_el('H')
        self.assertTrue(len(h) == len(presn.m), "No data for el: %s" % 'H')

    # @staticmethod
    # @unittest.skip("just for plot")
    def test_rho_plot(self):
        name = 'cat_R1000_M15_Ni007'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        rho_file = os.path.join(path, name + '.rho')
        presn = sn_eve.load_rho(rho_file)
        ax = presn.plot_chem(ylim=(1e-12, 1.))
        plt.show()

    def test_hyd_load(self):
        name = 'm030307mhh'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        # eve = StellaHydAbn(name, path=path).load_hyd()
        eve = sn_eve.load_hyd_abn(name, path=path)
        self.assertTrue(eve.is_set('Rho'), "hyd-file have not been loaded: %s" % name)

    def test_hyd_load_hyd_abn_no_dum(self):
        name = 'u50fenv'
        path = '/home/bakl/Sn/my/papers/18/fischer/model'
        # eve = StellaHydAbn(name, path=path).load_hyd()
        eve = sn_eve.load_hyd_abn(name, path=path, is_dum=False)
        self.assertTrue(eve.is_set('Rho'), "hyd-file have not been loaded: %s" % name)

    def test_presn_reshape(self):
        name = 'u50fenv'
        path = '/home/bakl/Sn/my/papers/18/fischer/model'
        # eve = StellaHydAbn(name, path=path).load_hyd()
        eve = sn_eve.load_hyd_abn(name, path=path, is_dum=False)

        # for same size must be the same
        nstart = 1
        nzon = eve.nzon
        evenew = eve.reshape(nz=eve.nzon)
        self.assertEqual(evenew.nzon, eve.nzon,
                         "Reshape PreSn: you have {} zone, but it should be {}".format(eve.nzon, eve.nzon))
        # todo check range for Rho, V, T
        k = 0
        self.assertEqual(eve.m[k], evenew.m[k],
                         "Mass PreSn: you have the first zone where old mass {} = new {}".format(eve.m[k], evenew.m[k]))
        k = -1
        self.assertEqual(eve.m[k], evenew.m[k],
                         "Mass PreSn: you have the last zone where old mass {} = new {}".format(eve.m[k], evenew.m[k]))

        nzon = 300
        nstart = 309
        evenew = eve.reshape(nz=nzon, nstart=nstart, nend=None)
        evenew.plot_chem()
        plt.show()
        self.assertEqual(evenew.nzon, nstart + nzon,
                         "Reshape PreSn: you have {} zone, but it should be {}".format(eve.nzon, nstart + nzon))

    def test_presn_zone_reduce(self):
        name = 'u50fenv'
        path = '/home/bakl/Sn/my/papers/18/fischer/model'
        # eve = StellaHydAbn(name, path=path).load_hyd()
        eve = sn_eve.load_hyd_abn(name, path=path, is_dum=False)

        # evenew = eve.zone_reduce(diff=1.005, start=110, end=-20)
        evenew = eve.zone_reduce(by='M', diff=1.01, start=309, end=-100)

        ii = np.arange(evenew.nzon)
        if not np.all(np.diff(evenew.m) > 0.):
            print("Not np.all(np.diff(evenew.m) > 0.): {}".format(ii[np.diff(evenew.m) < 0.]))
        self.assertTrue(np.all(np.diff(evenew.m) > 0.), " Any M[i] < M[i+1]")
        self.assertTrue(np.all(np.diff(evenew.r) > 0.), " Any R[i] < R[i+1]")

        # plot
        plt.plot(evenew.m / ps.phys.M_sun, np.log10(evenew.rho), marker='d', ls='', label='eveinterp',
                 linewidth=4)
        plt.plot(evenew.m / ps.phys.M_sun, np.log10(evenew.rho), ls='', marker='<', label='evenew', markersize='3')
        plt.plot(eve.m / ps.phys.M_sun, np.log10(eve.rho), ls='--', label='u50fenv', linewidth=2)
        plt.xlabel('M, [Msun]')
        plt.ylabel(r'$log_{10}(\rho)$')

        plt.legend(loc='best')
        plt.show()

    def test_presn_chem_norm(self):
        nzon = 10
        presn = PreSN('sol', nzon)
        data = np.ones(presn.nzon)
        for ename in PreSN.presn_elements:
            presn.set_chem(ename, data)

        for k in np.arange(1, presn.nzon):
            self.assertTrue(sum(presn.abun(k)) > 1., msg=" Norm should be large 1 ")
        # normalize
        for k in np.arange(1, presn.nzon):
            presn.chem_norm(k)
        # check
        for k in np.arange(1, presn.nzon):
            self.assertAlmostEqual(sum(presn.abun(k)), 1., msg="  Norm should be 1 ")

    def test_presn_set_composition(self):
        nzon = 10
        presn = PreSN('sol', nzon)
        presn.set_hyd(PreSN.sM, np.arange(1, presn.nzon+1)*ps.phys.M_sun)
        # set chemical
        data = np.ones(presn.nzon)
        for ename in PreSN.presn_elements:
            presn.set_chem(ename, data)

        for k in np.arange(1, presn.nzon+1):
            presn.chem_norm(k)

        idxs = np.arange(presn.nzon)
        zones = idxs[presn.el('S') > 0.01]+1
        # Remove S
        sample = {'S': 1e-4}
        presn.set_composition(zones, sample=sample, is_add=False, is_normalize=True)
        ax = presn.plot_chem(ylim=(1e-6, 1.1))

        # set solar composition
        zones = range(3, presn.nzon-3)
        presn.set_composition(zones, is_add=False)

        # plot
        presn.plot_chem()
        plt.show()

    def test_hyd_write(self):
        name = 'm030307mhh'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        presn = sn_eve.load_hyd_abn(name, path=path)
        fname = 'tmp.hyd'
        res = presn.write_hyd(fname)
        self.assertTrue(res, "Fail to write in %s" % fname)

    def test_abn_write(self):
        name = 'm030307mhh'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        eve = sn_eve.load_hyd_abn(name, path=path)
        fname = 'tmp.abn'
        res = eve.write_abn(fname)
        self.assertTrue(res, "Fail to write in %s" % fname)

    def test_rho_write_hyd_abn(self):
        name = 'cat_R1000_M15_Ni007'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        rho_file = os.path.join(path, name + '.rho')
        presn = sn_eve.load_rho(rho_file)

        fname = 'tmprho.hyd'
        res = presn.write_hyd(fname)
        self.assertTrue(res, "Fail to write in %s" % fname)
        fname = 'tmprho.abn'
        res = presn.write_abn(fname)
        self.assertTrue(res, "Fail to write in %s" % fname)

    def test_hyd_abn_load(self):
        name = 'm030307mhh'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        eve = sn_eve.load_hyd_abn(name, path=path)

        self.assertTrue(eve.nzon > 0,
                        "Zones numbers should be more 0 [%d]." % eve.nzon)

    # def test_hyd_abn_load_fischer(self):
    #     name = 'u50fenv'
    #     path = join(dirname(abspath(__file__)), 'data', 'stella')
    #     eve = sn_eve.load_hyd_abn(name, path='/home/bakl/Sn/my/papers/18/fischer/model', is_dum=False)
    #
    #     self.assertTrue(eve.nzon > 0,
    #                     "Zones numbers should be more 0 [%d]." % eve.nzon)

    def test_hyd_abn_load_H(self):
        name = 'm030307mhh'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        eve = sn_eve.load_hyd_abn(name, path=path)
        self.assertTrue(len(eve.el('H')) > 0, "abn-file have not been loaded: %s" % name)

    def test_snec_to_presn(self):
        path = join(dirname(abspath(__file__)), 'data', 'snec')
        prb_snec = snec.Problem('s15s7b2').load_chem(os.path.join(path, 's15s7b2_15isotope.dat'))
        prb_snec.load_profile(os.path.join(path, 's15s7b2.short'))

        presn = snec.to_presn(prb_snec)
        self.assertTrue(presn.nzon > 0, "Fail to convert snec to presn in %s" % prb_snec.Name)

    def test_snec_to_presn_plot(self):
        path = join(dirname(abspath(__file__)), 'data', 'snec')
        prb_snec = snec.Problem('s15s7b2').load_chem(os.path.join(path, 's15s7b2_15isotope.dat'))
        prb_snec.load_profile(os.path.join(path, 's15s7b2.short'))

        presn = snec.to_presn(prb_snec)
        self.assertTrue(presn.nzon > 0, "Fail to convert snec to presn in %s" % prb_snec.Name)

        presn.plot_chem(ylim=(1e-10, 1.))
        plt.show()

    def test_snec_write_hyd_abn(self):
        path = join(dirname(abspath(__file__)), 'data', 'snec')
        prb_snec = snec.Problem('s15s7b2').load_chem(os.path.join(path, 's15s7b2_15isotope.dat'))
        prb_snec.load_profile(os.path.join(path, 's15s7b2.short'))

        presn = snec.to_presn(prb_snec)

        fname = presn.Name + '.hyd'
        res = presn.write_hyd(fname)
        self.assertTrue(res, "Fail to write in %s" % fname)
        fname = presn.Name + '.abn'
        res = presn.write_abn(fname)
        self.assertTrue(res, "Fail to write in %s" % fname)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
