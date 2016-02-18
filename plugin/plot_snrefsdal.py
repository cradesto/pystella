import os

import numpy as np
from scipy import interpolate

from pystella.rf import band

grav_lens_def = 'ogu-a'

marker_glens = {'ogu-g': 'o', 'ogu-a': 's', 'gri-g': '+', 'sha-a': '*', 'die-a': 'd'}


def coef_glens():  # see http://arxiv.org/abs/1510.05750
    a = {'ogu-g': {'time': {'S1': 0,  'S2': 8.7, 'S3': 5.1, 'S4': 18.8, 'SX': 311.},
                   'mag':  {'S1': 1., 'S2': 1.14, 'S3': 1.22, 'S4': 0.67, 'SX': 0.27}},
         'ogu-a': {'time': {'S1': 0,  'S2': 9.4, 'S3': 5.6, 'S4': 20.9, 'SX': 335.6},
                   'mag':  {'S1': 1., 'S2': 17.7 / 15.4, 'S3': 18.3 / 15.4, 'S4': 9.8 / 15.4, 'SX': 4.2 / 15.4}},
         'gri-g': {'time': {'S1': 0,  'S2': 10.6, 'S3': 4.8, 'S4': 25.9, 'SX': 361.},
                   'mag':  {'S1': 1., 'S2': 0.92, 'S3': 0.99, 'S4': 0.42, 'SX': 0.36}},
         'sha-a': {'time': {'S1': 0,  'S2': 8.,   'S3': 5., 'S4': 17., 'SX': 233.},
                   'mag':  {'S1': 1., 'S2': 0.84,  'S3': 1.46, 'S4': 0.44, 'SX': 0.19}},
         'die-a': {'time': {'S1': 0,  'S2': -17., 'S3': -4., 'S4': 74., 'SX': 262.},
                   'mag':  {'S1': 1., 'S2': 1.89, 'S3': 0.64, 'S4': 0.35, 'SX': 0.31}}
         }
    return a


def coef_time_delay(model):  # see http://arxiv.org/abs/1510.05750
    a = coef_glens()
    return a[model]['time']


def coef_magnification(model):
    a = coef_glens()
    return a[model]['mag']


def plot(ax, dic=None):
    d = os.path.expanduser('~/Sn/my/papers/2016/snrefsdal/data')
    band_max = 'F160W'

    z = 0.
    arg = None
    im = 'S2'
    glens = grav_lens_def

    if 'args' in dic:
        arg = dic['args']
    if 'glens' in dic:
        glens = dic['glens']
    if 'image' in dic:
        im = dic['image']

    if arg is None:
        jd_shift = -56950
    else:
        jd_shift = float(arg.pop(0))

    print "Plot Sn Refsdal %s, glens: %s, path: %s, band_max=%s " % (im, glens, d, band_max)

    # plot B-V
    if 'bv' in dic:
        plot_BV(ax=ax, path=d, jd_shift=jd_shift, glens=glens, image=im)
    elif ax is not None:  # plot ubv
        t_lc_max = plot_ubv(ax=ax, path=d, jd_shift=jd_shift, band_max=band_max, glens=glens, image=im)

    # plot velocities
    if 'ax2' in dic and dic['ax2'] is not None:
        axVel = dic['ax2']
        if len(arg) > 0:
            z = float(arg.pop(0))
        print "Plot Sn Refsdal %s velocities: jd_shift=%8.1f t_max( %s) = %8.1f dt_exp = %8.1f" \
              % (im, jd_shift, band_max, t_lc_max, t_lc_max + jd_shift)
        plot_vel(ax=axVel, path=d, jd_shift=jd_shift + t_lc_max, z=z)

        # show spectral obs.
        dt = 1
        # ts_obs = np.array([57021, 57158])  # see Kelly ?
        ts_obs = np.array([57021, 57178])  # see Kelly ?
        for t in ts_obs + jd_shift:
            axVel.axvspan(t - dt, t + dt, facecolor='r', alpha=0.5)


def get_xy(d, b, colS):
    d2 = d[d[:, 0] == b, ]
    x = d2[:, 1]
    y = d2[:, colS]
    is_good = y != 0
    return x[is_good], y[is_good]


def plot_BV(ax, path, jd_shift, glens, image):
    colors = band.bands_colors()

    # from Rodney_tbl4
    sn_images = {'S1': 2, 'S2': 4, 'S3': 6, 'S4': 8, 'SX': 10}
    # if image in sn_images.keys():
    lc_data = np.loadtxt(os.path.join(path, 'rodney_all.csv'), comments='#', skiprows=3)

    # Kelly for SX
    z = np.zeros((np.shape(lc_data)[0], 2))
    lc_data = np.append(lc_data, z, axis=1)
    fs = {'F125W': 'kelly_F125W_SX.csv', 'F160W': 'kelly_F160W_SX.csv'}
    fs = dict((k, os.path.join(path, v)) for k, v in fs.items())
    for bn, fname in fs.items():
        # bn = 'F%sW' % int(b)
        b = int(bn[1:4])
        data = np.loadtxt(os.path.join(path, fname), comments='#')
        for row in data:
            r = np.zeros(np.shape(lc_data)[1])
            r[0] = int(b)
            r[1] = row[0]
            r[sn_images['SX']] = row[1]
            r[sn_images['SX'] + 1] = row[2]

            lc_data = np.vstack([lc_data, r])

    bands = np.unique(lc_data[:, 0])
    time_delay = coef_time_delay(glens)

    colS = sn_images[image]  # S1 - 2, col  S2 - 4, S3 - 6, S4 - 8 col

    for b1, b2 in zip(bands[:-1], bands[1:]):
        bn1 = 'F%sW' % int(b1)
        bn2 = 'F%sW' % int(b2)
        x1, y1, = get_xy(lc_data, b1, colS)
        x1 = x1 + jd_shift - time_delay[image]
        if np.all(y1 <= 0):
            continue
        x2, y2, = get_xy(lc_data, b2, colS)
        x2 = x2 + jd_shift - time_delay[image]
        if np.all(y2 <= 0):
            continue

        if x1[-1] > x2[-1]:
            y_spline = interpolate.splrep(x1, y1, s=0)
            xx = x2
            yy1 = interpolate.splev(xx, y_spline)
            bv = yy1 - y2
        else:
            y_spline = interpolate.splrep(x2, y2, s=0)
            xx = x1
            yy2 = interpolate.splev(xx, y_spline)
            bv = y1 - yy2

        ax.plot(xx, bv,  color=colors[bn1], label='%s-%s .' % (bn1, bn2))


def plot_ubv(ax, path, jd_shift, band_max, glens, image):
    colors = band.bands_colors()

    # from Rodney_tbl4
    sn_images = {'S1': 2, 'S2': 4, 'S3': 6, 'S4': 8, 'SX': 10}
    # if image in sn_images.keys():
    lc_data = np.loadtxt(os.path.join(path, 'rodney_all.csv'), comments='#', skiprows=3)

    # Kelly for SX
    z = np.zeros((np.shape(lc_data)[0], 2))
    lc_data = np.append(lc_data, z, axis=1)
    fs = {'F125W': 'kelly_F125W_SX.csv', 'F160W': 'kelly_F160W_SX.csv'}
    fs = dict((k, os.path.join(path, v)) for k, v in fs.items())
    for bn, fname in fs.items():
        # bn = 'F%sW' % int(b)
        b = int(bn[1:4])
        data = np.loadtxt(os.path.join(path, fname), comments='#')
        for row in data:
            r = np.zeros(np.shape(lc_data)[1])
            r[0] = int(b)
            r[1] = row[0]
            r[sn_images['SX']] = row[1]
            r[sn_images['SX'] + 1] = row[2]

            lc_data = np.vstack([lc_data, r])

    bands = np.unique(lc_data[:, 0])
    time_delay = coef_time_delay(glens)

    colS = sn_images[image]  # S1 - 2, col  S2 - 4, S3 - 6, S4 - 8 col

    for b in bands:
        bn = 'F%sW' % int(b)
        data = lc_data[lc_data[:, 0] == b, ]
        x = data[:, 1] + jd_shift - time_delay[image]
        y = data[:, colS]
        if np.all(y <= 0):
            continue

        yerr = data[:, colS + 1]
        bcolor = colors[bn]
        lower_limits = np.zeros(x.shape, dtype=bool)
        lower_limits[yerr == -1] = True

        # plot values with errors
        with_errors = np.invert(lower_limits)
        xx = x[with_errors]
        if len(xx) > 0:
            yy = y[with_errors]
            yyerr = yerr[with_errors]
            # ax.errorbar(xx, yy, yerr=yyerr, fmt=marker_glens[glens], color=bcolor, label='%s obs.' % bn)
            ax.errorbar(xx, yy, yerr=yyerr, fmt='o', color=bcolor, label='%s obs.' % bn)

        # plot upper limits
        xxx = x[lower_limits]
        if len(xx) > 0:
            yyy = y[lower_limits]
            yerr = [-0.5 * np.ones(xxx.shape), np.zeros(xxx.shape)]
            ax.errorbar(xxx, yyy, yerr=yerr, lolims=True, fmt='kv', color=bcolor, lw=1.2)

    data = lc_data[lc_data[:, 0] == int(band_max[1:4]), ]  # extract maximum for F160W
    data = data[data[:, colS] > 0., :]  # remove 0
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
