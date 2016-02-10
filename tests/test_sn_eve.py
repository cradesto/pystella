from os.path import dirname, abspath, join
import unittest
from pystella.model.sn_eve import StellaEve

__author__ = 'bakl'


class SnEveTests(unittest.TestCase):
    def test_load_eve(self):
        name = 'cat_R1000_M15_Ni007'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        eve = StellaEve(name, path=path)
        data = eve.rho_load()
        self.assertIsNotNone(data, "Rho-file have not been loaded: %s." % eve.rho_file)

        # self.assertAlmostEqual(b.zp, zp, "Zero points of band %s equals %f. Should be %f" % (b, b.zp, zp))


def main():
    unittest.main()


if __name__ == '__main__':
    main()
