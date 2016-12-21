import os
from os.path import dirname, abspath, join
import unittest
from pystella.model import sn_eve
import matplotlib.pyplot as plt

__author__ = 'bakl'


class SnEveTests(unittest.TestCase):
    def test_eve_load_rho(self):
        name = 'cat_R1000_M15_Ni007'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        rho_file = os.path.join(path, name + '.rho')

        presn = sn_eve.load_rho(rho_file)
        self.assertTrue(presn.nzon > 0, "Rho-file have not been loaded: %s" % rho_file)

    def test_eve_el(self):
        name = 'cat_R1000_M15_Ni007'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        rho_file = os.path.join(path, name + '.rho')

        presn = sn_eve.load_rho(rho_file)
        h = presn.lg_el('H')
        self.assertTrue(len(h) == len(presn.m), "No data for el: %s" % 'H')

    # @staticmethod
    # @unittest.skip("just for plot")
    def test_eve_plot(self):
        name = 'cat_R1000_M15_Ni007'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        rho_file = os.path.join(path, name + '.rho')
        presn = sn_eve.load_rho(rho_file)
        ax = presn.plot_chem(ylim=(-8,0))
        plt.show()

    def test_eve_load_hyd(self):
        name = 'm030307mhh'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        # eve = StellaHydAbn(name, path=path).load_hyd()
        eve = sn_eve.load_hyd_abn(name, path=path)
        self.assertTrue(eve.is_set('Rho'), "hyd-file have not been loaded: %s" % name)

    def test_eve_write_hyd(self):
        name = 'm030307mhh'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        presn = sn_eve.load_hyd_abn(name, path=path)
        fname = 'tmp.hyd'
        res = presn.write_hyd(fname)
        self.assertTrue(res, "Fail to write in %s" % fname)

    def test_eve_write_abn(self):
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

    def test_eve_load_zones(self):
        name = 'm030307mhh'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        eve = sn_eve.load_hyd_abn(name, path=path)

        self.assertTrue(eve.nzon > 0,
                        "Zones numbers should be more 0 [%d]." % eve.nzon)

    def test_eve_load_abn(self):
        name = 'm030307mhh'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        eve = sn_eve.load_hyd_abn(name, path=path)
        self.assertTrue(len(eve.el('H')) > 0, "abn-file have not been loaded: %s" % name)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
