import os
import unittest
from os.path import dirname, abspath, join

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import interactive

from pystella import velocity as vel

# interactive(True)


__author__ = 'bakl'


class TestVelocity(unittest.TestCase):
    def setUp(self):
        self.mname = 'cat_R1000_M15_Ni007_E15'
        self.path = join(dirname(abspath(__file__)), 'data', 'stella')

    # @unittest.skip("just for plot")
    def test_velocity_ttres(self):
        fig = plt.figure(num=None, figsize=(12, 8), dpi=100, facecolor='w', edgecolor='k')
        gs1 = gridspec.GridSpec(4, 1)
        plt.matplotlib.rcParams.update({'font.size': 14})

        nm, path = self.mname, self.path
        # nm = 'rn300_R10_M3Mht3t30_Ni0_E003wAq2e3'
        # path = '~/Sn/Release/svn_kepler/stella/branches/lucy/run/res/sncurve/rednovaM31/tt/'
        print(nm, os.path.expanduser(path))
        vels = vel.compute_vel_res_tt(nm, os.path.expanduser(path), is_new_std=False)

        ax = fig.add_subplot(gs1[:, 0])
        vel.plot_vel(ax, vels)
        plt.show()    # @unittest.skip("just for plot")

    def test_velocity_ttres_is_new_std(self):
        fig = plt.figure(num=None, figsize=(12, 8), dpi=100, facecolor='w', edgecolor='k')
        gs1 = gridspec.GridSpec(4, 1)
        plt.matplotlib.rcParams.update({'font.size': 14})

        nm, path = self.mname, self.path
        # nm = 'rn300_R10_M3Mht3t30_Ni0_E003wAq2e3'
        # path = '~/Sn/Release/svn_kepler/stella/branches/lucy/run/res/sncurve/rednovaM31/tt/'
        nm, path = 'nirefE5R50M26Ni3m2b2m4Z01', '/home/bakl/Sn/Release/seb_git/run/87a/2fit/tmp'
        print(nm, os.path.expanduser(path))
        vels = vel.compute_vel_res_tt(nm, os.path.expanduser(path), is_new_std=True)

        ax = fig.add_subplot(gs1[:, 0])
        vel.plot_vel(ax, vels)
        plt.show()

    def test_velocity_swd(self):
        fig = plt.figure(num=None, figsize=(12, 8), dpi=100, facecolor='w', edgecolor='k')
        gs1 = gridspec.GridSpec(4, 1)
        plt.matplotlib.rcParams.update({'font.size': 14})

        nm, path = self.mname, self.path
        # nm = 'rn300_R10_M3Mht3t30_Ni0_E003wAq2e3'
        # path = '~/Sn/Release/svn_kepler/stella/branches/lucy/run/res/sncurve/rednovaM31/tt/'
        vels = vel.compute_vel_swd(nm, os.path.expanduser(path))

        ax = fig.add_subplot(gs1[:, 0])
        vel.plot_vel(ax, vels)
        plt.show()

    def test_velocity_compare_ttresVSswd(self):
        nm, path = self.mname, self.path
        # nm = 'rn300_R10_M3Mht3t30_Ni0_E003wAq2e3'
        # path = '~/Sn/Release/svn_kepler/stella/branches/lucy/run/res/sncurve/rednovaM31/tt/'
        vels_tt = vel.compute_vel_res_tt(nm, os.path.expanduser(path))
        vels_swd = vel.compute_vel_swd(nm, os.path.expanduser(path))

        plt.matplotlib.rcParams.update({'font.size': 14})
        fig = plt.figure(num=None, figsize=(12, 8), dpi=100, facecolor='w', edgecolor='k')
        gs1 = gridspec.GridSpec(4, 1)
        ax = fig.add_subplot(gs1[:, 0])
        # fig = plt.figure(figsize=(20, 10))
        # ax = fig.add_axes((0.1, 0.3, 0.8, 0.65))

        vel.plot_vel(ax, vels_tt, label='tt')
        vel.plot_vel(ax, vels_swd, color='green', label='swd')

        ax.legend()
        plt.grid(linestyle=':', linewidth=1)
        plt.show()


def main():
    unittest.main()


if __name__ == '__main__':
    main()
