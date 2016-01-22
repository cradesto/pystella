import os

import numpy as np

from pystella.rf import band


def plot(ax, arg):
    #     plot_snrefsdal(ax=ax,  arg=arg)
    #
    #
    # def plot_snrefsdal(ax, arg=None):
    print "Plot Sn Refsdal, "
    d = os.path.expanduser('~/Sn/my/papers/2016/snrefsdal/data')
    if arg is None:
        jd_shift = -57000
    else:
        jd_shift = float(arg[0])

    colors = band.bands_colors()

    # kelly's data from plot
    if False:
        fs = {'F125W': os.path.join(d, 'snrefsdal_F125W_S2.csv'), 'F160W': d + 'snrefsdal_F160W_S2.csv'}

        for b, fname in fs.items():
            data = np.loadtxt(fname, comments='#')
            x = data[:, 0] + jd_shift
            y = data[:, 1]
            bcolor = colors[b]
            ax.plot(x, y, label='%s Sn, Kelly' % b, ls=".", color=bcolor, markersize=8, marker="o")

    # from Rodney_tbl4
    rodney = np.loadtxt(os.path.join(d, 'rodney_all.csv'), comments='#', skiprows=3)
    bands = np.unique(rodney[:, 0])
    colS = 2  # S1 - 2 col  S4 - 8 col
    for b in bands:
        bn = 'F%sW' % int(b)
        data = rodney[rodney[:, 0] == b, ]
        x = data[:, 1] + jd_shift
        y = data[:, colS]
        yerr = data[:, colS + 1]
        bcolor = colors[bn]
        # ax.plot(x, y, label='%s Sn, Rodney' % lbl(bn, band_shift), ls="-.", color=bcolor, markersize=8, marker="*")
        ax.errorbar(x, y, yerr=yerr, fmt='o', color=bcolor, label='%s Sn, Rodney' % bn)
        # print max
        t_min = data[np.argmin(y), 1]
        print "t_max( %s) = %8.1f, shifted=%8.1f" % (bn, t_min, t_min + jd_shift)

    #
    if len(arg) > 1 is not None:
        band_max = 'F160W'
        data = rodney[rodney[:, 0] == int(band_max[1:4]), ]  # extract maximum for F160W
        t_min = data[np.argmin(data[:, colS]), 1]
        print "Plot Sn Refsdal velocities: t_max( %s) = %8.1f" % (band_max, t_min)
        axVel = arg[-1]
        data = np.loadtxt(os.path.join(d, 'kelly_vel_halpha.txt'), comments='#', skiprows=2)
        x = data[:, 0] + t_min + jd_shift
        y = data[:, 1]
        yerr = data[:, 2]
        bcolor = 'red'
        axVel.errorbar(x, y, yerr=yerr, fmt='o', color=bcolor, label=r'$H_{\alpha}$, Kelly')
