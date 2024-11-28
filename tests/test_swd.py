from os.path import dirname, abspath, join
import unittest

import matplotlib.pyplot as plt
from matplotlib import gridspec

from pystella.model import sn_swd
from pystella.model.stella import Stella

__author__ = 'bakl'


class TestStellaShockWaveDetail(unittest.TestCase):
    def setUp(self):
        name = 'rednova_R3.2_M6_Ni0_E0.25'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        stella = Stella(name, path=path)
        self.swd = stella.get_swd()
        self.swd.load()

    def test_reading(self):
        nzon = 100
        ntimes = 20
        self.assertEquals(self.swd.NzonMax, nzon, "Error in zon setting: swd.Nzon=%d, but should be %d"
                          % (self.swd.NzonMax, nzon))
        self.assertEquals(self.swd.Ntimes, ntimes, "Error in time reading: swd.Ntimes=%d, but should be %d"
                          % (self.swd.Ntimes, ntimes))

    def test_reading_block(self):
        time = 2.
        idx, t = self.swd.time_nearest(time)
        b = self.swd.block_nearest(time)

        self.assertEquals(self.swd.NzonMax, len(b.R), "Error in block size: len(block)[%d] != Nzon [%d]"
                          % (len(b.R), self.swd.NzonMax))

        # self.assertTrue(np.all([b.Tau[i] < b.Tau[i+1] for i in range(0, b.Nzon-1)]),
        #                 "The Tau should be a monotonically increasing function")

    # @unittest.skip("just for plot")
    def test_block_plot(self):
        time = 2.
        idx, t = self.swd.time_nearest(time)
        b = self.swd.block_nearest(time)

        fig = plt.figure(num=None, figsize=(12, 8), dpi=100, facecolor='w', edgecolor='k')
        gs1 = gridspec.GridSpec(4, 1)
        plt.matplotlib.rcParams.update({'font.size': 14})
        # if is_vel:
        #     axUbv = fig.add_subplot(gs1[:-1, 0])
        #     axVel = fig.add_subplot(gs1[3, 0])
        # else:
        #     axUbv = fig.add_subplot(gs1[:, 0])

        ax = fig.add_subplot(gs1[:, 0])
        # ax = fig.add_axes((0.1, 0.3, 0.8, 0.65))
        sn_swd.plot_swd(ax, b)
        # ax.legend()
        plt.show()
        # self.assertTrue(np.all([b.Tau[i] < b.Tau[i+1] for i in range(0, b.Nzon-1)]),
        #                 "The Tau should be a monotonically increasing function")

    # @unittest.skip("just for plot")
    def test_block_plot_many(self):
        times = [1., 2., 3.]

        fig = plt.figure(num=None, figsize=(8, 12), dpi=100, facecolor='w', edgecolor='k')
        gs1 = gridspec.GridSpec(len(times), 1)
        plt.matplotlib.rcParams.update({'font.size': 14})

        i = 0
        for t in times:
            ax = fig.add_subplot(gs1[i, 0])
            b = self.swd.block_nearest(t)
            sn_swd.plot_swd(ax, b, is_xlabel=i==len(times)-1, vnorm=1e7, lumnorm=1e36, is_legend=False)
            i += 1
        plt.show()

    def test_swd_evolution_exeption_nz(self):
        def func():
            self.swd.load()
            times = []
            rhos = []
            for t, y in self.swd.evolution('Rho', self.swd.Nzon+1):
                times.append(t)
                rhos.append(y)
            # # t, y = self.swd.evolution('Rho', self.swd.Nzon+1)
            # print('Len(times)= {}'.format(len(times)))
        self.assertRaises(ValueError, func)

    def test_swd_evolution_exeption_val(self):
        def func():
            t, y = self.swd.evolution('bad_value', self.swd.Nzon)

        self.assertRaises(ValueError, func)

    def test_swd_evolution(self):
        times = []
        rhos = []
        for t, y in self.swd.evolution('Rho', 20):
            times.append(t)
            rhos.append(y)

        # t, y = self.swd.evolution('Rho', 20)
        self.assertCountEqual(times, self.swd.Times, "The times should be the same size.")

    def test_swd_params_ph(self):
        params_ph = self.swd.params_ph()
        self.assertCountEqual(params_ph['time'], self.swd.Times, "The times should be the same size.")


def main():
    unittest.main()


if __name__ == '__main__':
    main()
