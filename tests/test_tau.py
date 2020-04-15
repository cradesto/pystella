import os
import sys
import unittest

import pystella as ps

__author__ = 'bakl'


class TestTau(unittest.TestCase):
    def setUp(self):
        name = 'levJ_R450_M15_Ni004_E10'
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'stella')
        stella = ps.Stella(name, path=path)
        self.tau = stella.get_tau()

    def test_find_startend_blocks(self):
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
        nfreqs = 50

        self.tau.load()

        self.assertEqual(self.tau.Nzon, nzon, "Error in zon setting: tau.Nzon=%d, but should be %d"
                         % (self.tau.Nzon, nzon))

        self.assertEqual(self.tau.Ntimes, ntimes, "Error in time reading: tau.Ntimes=%d, but should be %d"
                         % (self.tau.Ntimes, ntimes))

        self.assertEqual(self.tau.NFreq, nfreqs, "Error in time reading: tau.NFreq=%d, but should be %d"
                         % (self.tau.NFreq, nfreqs))

    def test_block(self):
        self.tau.load()
        b = self.tau[0]

        nzon = b.Size
        self.assertEqual(self.tau.Nzon, nzon, "Error in block zones: nzon=%d, but should be tau.Nzon= %d"
                         % (nzon, self.tau.Nzon))


def main():
    unittest.main()


if __name__ == '__main__':
    main()
