##
#  Callbacks
##
import numpy as np

import pystella.rf.light_curve_func as lc
from pystella.rf import band


def plot(ax, arg):
    lw = 2.
    band_shift = arg
    colors = band.bands_colors()
    fname = "~/Desktop/Downloads/2/100z0E60Ni_6.ph.hsc.2"
    data = np.loadtxt(fname, comments='#')
    fs = list('grizy')
    x = data[:, 0]
    for i in range(len(fs)):
        y = data[:, i + 1]
        bcolor = colors[fs[i]]
        ax.plot(x, y, label='%s Tolstov' % lc.lbl(fs[i], band_shift),
                color=bcolor, ls="-.", linewidth=lw)
