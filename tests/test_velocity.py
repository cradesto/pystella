import os
import unittest

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import interactive

from pystella import velocity as vel

interactive(True)


__author__ = 'bakl'


class TestVelocity(unittest.TestCase):
    # @unittest.skip("just for plot")
    def test_velocity(self):
        fig = plt.figure(num=None, figsize=(12, 8), dpi=100, facecolor='w', edgecolor='k')
        gs1 = gridspec.GridSpec(4, 1)
        plt.matplotlib.rcParams.update({'font.size': 14})

        nm = 'rn300_R10_M3Mht3t30_Ni0_E003wAq2e3'
        path = '~/Sn/Release/svn_kepler/stella/branches/lucy/run/res/sncurve/rednovaM31/tt/'
        vels = vel.compute_vel_res_tt(nm, os.path.expanduser(path))

        ax = fig.add_subplot(gs1[:, 0])
        vel.plot_vel(ax, vels)
        plt.show()
        # plt.savefig("temp.png")
        # self.assertEquals(self.swd.Nzon, len(b.R), "Error in block size: len(block)[%d] != Nzon [%d]"
        #                   % (len(b.R), self.swd.Nzon))


def main():
    unittest.main()


if __name__ == '__main__':
    main()
