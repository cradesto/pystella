from os.path import dirname, abspath, join
import unittest
from pystella.model.sn_eve import StellaEve, StellaHydAbn

__author__ = 'bakl'


class SnEveTests(unittest.TestCase):
    def test_eve_load_rho(self):
        name = 'cat_R1000_M15_Ni007'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        eve = StellaEve(name, path=path).load_rho()
        self.assertTrue(eve.is_load_rho, "Rho-file have not been loaded: %s" % eve.rho_file)

    def test_eve_el(self):
        name = 'cat_R1000_M15_Ni007'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        eve = StellaEve(name, path=path).load()
        h = eve.lg_el('H')
        self.assertTrue(len(h) == len(eve.mass), "No data for el: %s" % 'H')

    @staticmethod
    @unittest.skip("just for plot")
    def test_eve_plot():
        name = 'cat_R1000_M15_Ni007'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        eve = StellaEve(name, path=path).load()
        eve.plot_chem()

    def test_eve_load_hyd(self):
        name = 'm030307mhh'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        eve = StellaHydAbn(name, path=path).load_hyd()
        self.assertTrue(eve.is_load_hyd, "hyd-file have not been loaded: %s" % eve.hyd_file)

    def test_eve_write_hyd(self):
        name = 'm030307mhh'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        eve = StellaHydAbn(name, path=path).load_hyd()
        fname = 'tmp.hyd'
        res = eve.write_hyd(fname)
        self.assertTrue(res, "Fail to write in %s" % fname)

    def test_eve_load_abn(self):
        name = 'm030307mhh'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        eve = StellaHydAbn(name, path=path).load_abn()
        self.assertTrue(eve.is_load_chem, "abn-file have not been loaded: %s" % eve.abn_file)

    def test_eve_write_abn(self):
        name = 'm030307mhh'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        eve = StellaHydAbn(name, path=path).load_abn()
        fname = 'tmp.abn'
        res = eve.write_abn(fname)
        self.assertTrue(res, "Fail to write in %s" % fname)

    def test_eve_load_zones(self):
        name = 'm030307mhh'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        eve = StellaHydAbn(name, path=path).load_hyd().load_abn()

        self.assertTrue(eve.nzon_abn == eve.nzon_hyd,
                        "Zones numbers should be the same in hyd-[%d] and abn-[%d] files." \
                        % (eve.nzon_hyd, eve.nzon_abn))

        # self.assertAlmostEqual(b.zp, zp, "Zero points of band %s equals %f. Should be %f" % (b, b.zp, zp))


def main():
    unittest.main()


if __name__ == '__main__':
    main()
