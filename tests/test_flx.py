import unittest
from os.path import dirname, abspath, join

import os


import matplotlib.pyplot as plt

import pystella.model.sn_flx as flx
from pystella.model.stella import Stella

__author__ = 'bakl'


class TestStellaFlx(unittest.TestCase):
    def setUp(self):
        pass
        # name = 'cat_R1000_M15_Ni007_E15'
        # path = join(dirname(abspath(__file__)), 'data', 'stella')
        # stella = Stella(name, path=path)
        # self.flx = stella.get_flx()

    def test_flx_reader(self):
        name = 'cat_R1000_M15_Ni007_E15'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        fname = os.path.join(path, name + '.flx')
        flx_data = flx.flx_reader(fname)
        flx_data.show_emergent_Fl(logy=False)
        plt.show()

    def test_stella_get_flx(self):
        name = 'cat_R1000_M15_Ni007_E15'
        path = join(dirname(abspath(__file__)), 'data', 'stella')

        flx_data = Stella(name, path=path).get_flx()
        flx_data.show_emergent_Fl(logy=False)
        plt.show()


def main():
    unittest.main()


if __name__ == '__main__':
    main()
