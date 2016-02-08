import os

import numpy as np

from pystella.rf import band


def coef_time_delay(model):
    a = {'oguri':
             {'S1': 0, 'S2': 9.4, 'S3': 5.6, 'S4': 20.9, 'SX': 335.6}
         }
    return a[model]


def coef_magnification(model):
    a = {'oguri':  # Masamune Oguri
              {'S1': 1., 'S2': 17.7 / 15.4, 'S3': 18.3 / 15.4, 'S4': 9.8 / 15.4, 'SX': 4.2 / 15.4}
         }
    return a[model]


def plot(ax, dic=None):
    d = os.path.expanduser('~/Sn/my/papers/2016/snrefsdal/data')
    band_max = 'F160W'
    print "Plot Sn Refsdal, path: %s, band_max=%s " % (d, band_max)

    z = 0.
    arg = None
    im = 'S1'

    if 'args' in dic:
        arg = dic['args']
    if 'image' in dic:
        im = dic['image']

    if arg is None:
        jd_shift = -56950
    else:
        jd_shift = float(arg.pop(0))

    t_min = plot_ubv(ax=ax, path=d, jd_shift=jd_shift, band_max=band_max, image=im)

    if 'ax2' in dic and dic['ax2'] is not None:
        axVel = dic['ax2']
        if len(arg) > 0:
            z = float(arg.pop(0))
        print "Plot Sn Refsdal velocities: jd_shift=%8.1f t_max( %s) = %8.1f dt_exp = %8.1f"\
              % (jd_shift, band_max, t_min, t_min+jd_shift)
        plot_vel(ax=axVel, path=d, jd_shift=jd_shift + t_min, z=z)

        # show spectral obs.
        dt = 1
        ts_obs = np.array([57033, 57158])
        for t in ts_obs+jd_shift:
             axVel.axvspan(t-dt, t+dt, facecolor='r', alpha=0.5)


def plot_ubv(ax, path, jd_shift, band_max, image='S1'):
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
    sn_images = {'S1': 2, 'S2': 4, 'S3': 6, 'S4': 8}

    time_delay = coef_time_delay('oguri')

    colS = sn_images[image]  # S1 - 2, col  S2 - 4, S3 - 6, S4 - 8 col

    for b in bands:
        bn = 'F%sW' % int(b)
        data = rodney[rodney[:, 0] == b,]
        x = data[:, 1] + jd_shift + time_delay[image]
        y = data[:, colS]  # compensate the magnification
        yerr = data[:, colS + 1]
        bcolor = colors[bn]
        lowerlimits = np.zeros(x.shape, dtype=bool)
        lowerlimits[yerr == -1] = True

        # plot values with errors
        witherrors = np.invert(lowerlimits)
        xx = x[witherrors]
        if len(xx) > 0:
            yy = y[witherrors]
            yyerr = yerr[witherrors]
            ax.errorbar(xx, yy, yerr=yyerr, fmt='o', color=bcolor, label='%s Sn, Rodney' % bn)

        # plot upper limits
        xxx = x[lowerlimits]
        if len(xx) > 0:
            yyy = y[lowerlimits]
            yerr = [-0.5*np.ones(xxx.shape), np.zeros(xxx.shape)]
            # print '%s: xxx=%s & yyy=%s lowerlimits=%s' % (bn, xxx.shape, yyy.shape, lowerlimits.shape)
            ax.errorbar(xxx, yyy, yerr=yerr, lolims=True, fmt='kv', color=bcolor, lw=1.2)

        # ax.errorbar(x, y, yerr=yerr, lolims=lowerlimits, fmt='o', color=bcolor, label='%s Sn, Rodney' % bn)

        # print max
        t_min = data[np.argmin(y), 1]
        print "t_max( %s) = %8.1f, shifted=%8.1f" % (bn, t_min, t_min + jd_shift)

    data = rodney[rodney[:, 0] == int(band_max[1:4]),]  # extract maximum for F160W
    t_min = data[np.argmin(data[:, colS]), 1]
    return t_min


def plot_vel(ax, path, jd_shift, z=0.):
    data = np.loadtxt(os.path.join(path, 'kelly_vel_halpha.txt'), comments='#', skiprows=2)
    x = data[:, 0] * (1. + z) + jd_shift
    xerr = data[:, 1] * (1. + z)
    y = data[:, 2]
    yerr = data[:, 3]
    bcolor = 'red'
    ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', color=bcolor, label=r'$H_{\alpha}$, Kelly')

    # for xlim in zip(x-xerr, x+xerr):
    #     ax.axvspan(xlim[0], xlim[1], facecolor='g', alpha=0.5)
    # # ylim = ax.get_ylim()
    # fb = ax.fill_between([x-errx, x+errx], ylim[0], ylim[1], facecolor='blue', alpha=0.2, label='Visible')
