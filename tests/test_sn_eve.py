from os.path import dirname, abspath, join
import unittest
from pystella.model.sn_eve import StellaEve

__author__ = 'bakl'


class SnEveTests(unittest.TestCase):
    def test_eve_load(self):
        name = 'cat_R1000_M15_Ni007'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        eve = StellaEve(name, path=path).load()
        self.assertTrue(eve.is_chem_load, "Rho-file have not been loaded: %s" % eve.rho_file)

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

        # self.assertAlmostEqual(b.zp, zp, "Zero points of band %s equals %f. Should be %f" % (b, b.zp, zp))


def main():
    unittest.main()


if __name__ == '__main__':
    main()
