import os
import unittest

import pystella as ps

__author__ = 'bakl'


class TestTau(unittest.TestCase):
    def setUp(self):
        name = 'levJ_R450_M15_Ni004_E10'
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'stella')
        stella = ps.Stella(name, path=path)
        self.tau = stella.get_tau()

    def test_find_start_end_blocks(self):
        se, times = self.tau.parse_start_end_blocks(self.tau.FName)
        s, e = se[0]
        for i in (0, 1, 2, -2, -1):
            print(i, se[i], times[i])

        self.assertEqual(len(se), len(times),
                         "Error in find_start_end_blocks: len(se)=len(times), "
                         "but should be len(se)= {}, but len(times)= {}".format(len(se), len(times)))

    def test_load(self):
        nzon = 100
        ntimes = 202

        self.tau.load()

        self.assertEqual(self.tau.Nzon, nzon, "Error in zon setting: tau.Nzon=%d, but should be %d"
                         % (self.tau.Nzon, nzon))

        self.assertEqual(self.tau.Ntimes, ntimes, "Error in time reading: tau.Ntimes=%d, but should be %d"
                         % (self.tau.Ntimes, ntimes))

    def test_block(self):
        nfreqs = 50
        self.tau.load()
        b = self.tau[0]

        self.assertEqual(b.Nzon, self.tau.Nzon, "Error in block zones: b.Nzon=%d, but should be tau.Nzon= %d"
                         % (b.Nzon, self.tau.Nzon))
        self.assertEqual(b.NFreq, nfreqs, "Error in block reading: b.NFreq=%d, but should be %d"
                         % (b.NFreq, nfreqs))

    def test_vars(self):
        import matplotlib.pyplot as plt

        self.tau.load()
        b = self.tau.block_nearest(4.)

        plt.semilogx(b.R, b.V8)
        plt.xlabel('R [cm]')
        plt.ylabel('V8')
        plt.show()


def main():
    unittest.main()


if __name__ == '__main__':
    main()
