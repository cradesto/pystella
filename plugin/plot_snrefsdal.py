import os

import numpy as np

from pystella.rf import band


def plot(ax, arg):
    d = os.path.expanduser('~/Sn/my/papers/2016/snrefsdal/data')
    band_max = 'F160W'
    print "Plot Sn Refsdal, path: %s, band_max=%s " % (d, band_max)

    z = 0.

    if arg is None:
        jd_shift = -56950
    else:
        jd_shift = float(arg.pop(0))

    t_min = plot_ubv(ax=ax, path=d, jd_shift=jd_shift, band_max=band_max)

    if len(arg) > 0 is not None:
        axVel = arg[-1]
        if len(arg) > 0:
            z = float(arg.pop(0))
        print "Plot Sn Refsdal velocities: jd_shift=%8.1f t_max( %s) = %8.1f" % (jd_shift, band_max, t_min)
        plot_vel(ax=axVel, path=d, jd_shift=jd_shift+t_min, z=z)


def plot_ubv(ax, path, jd_shift, band_max):
    colors = band.bands_colors()

    # kelly's data from plot
    if False:
        fs = {'F125W': os.path.join(path, 'snrefsdal_F125W_S2.csv'), 'F160W': path + 'snrefsdal_F160W_S2.csv'}

        for b, fname in fs.items():
            data = np.loadtxt(fname, comments='#')
            x = data[:, 0] + jd_shift
            y = data[:, 1]
            bcolor = colors[b]
            ax.plot(x, y, label='%s Sn, Kelly' % b, ls=".", color=bcolor, markersize=8, marker="o")

    # from Rodney_tbl4
    rodney = np.loadtxt(os.path.join(path, 'rodney_all.csv'), comments='#', skiprows=3)
    bands = np.unique(rodney[:, 0])
    colS = 4  # S1 - 2, col  S2 - 4, S3 - 6, S4 - 8 col
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

    data = rodney[rodney[:, 0] == int(band_max[1:4]), ]  # extract maximum for F160W
    t_min = data[np.argmin(data[:, colS]), 1]
    return t_min


def plot_vel(ax, path, jd_shift, z=0.):
        data = np.loadtxt(os.path.join(path, 'kelly_vel_halpha.txt'), comments='#', skiprows=2)
        x = data[:, 0]*(1.+z) + jd_shift
        y = data[:, 1]
        yerr = data[:, 2]
        bcolor = 'red'
        ax.errorbar(x, y, yerr=yerr, fmt='o', color=bcolor, label=r'$H_{\alpha}$, Kelly')
