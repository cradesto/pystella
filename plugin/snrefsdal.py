import numpy as np
import os
from scipy import interpolate

from pystella.rf.lc import SetLightCurve, LightCurve

path_data = os.path.expanduser('~/Sn/my/papers/2016/snrefsdal/data')
grav_lens_def = 'ogu-g'
marker_glens = {'ogu-g': 'o', 'ogu-a': 's', 'gri-g': '+', 'sha-g': '*', 'sha-a': '*', 'die-a': 'd',
                'bakl': 'D', 'obs-pol': '>', 'obs-tmp': '<'}

colors_gl = {'ogu-g': 'blue', 'gri-g': 'magenta', 'sha-g': 'orange',
             'ogu-a': 'skyblue', 'sha-a': 'red', 'die-a': 'olive',
             'obs-tmp': 'chocolate', 'obs-sn87a': 'brown', 'obs-pol': 'black'}
colors_band = dict(F105W="magenta", F435W="skyblue", F606W="cyan", F125W="g",
                   F140W="orange", F160W="r", F814W="blue")


def coef_glens():  # see http://arxiv.org/abs/1510.05750
    a = {'ogu-g': {'time': {'S1': 0, 'S2': 8.7, 'S3': 5.1, 'S4': 18.8, 'SX': 311.},
                   'etime': {'S1': (0., 0.), 'S2': (.7, .7), 'S3': (.5, .5), 'S4': (1.7, 1.7), 'SX': (24., 24.)},
                   'mag': {'S1': 1., 'S2': 1.14, 'S3': 1.22, 'S4': 0.67, 'SX': 0.27},
                   'emag': {'S1': (0, 0), 'S2': (.24, .24), 'S3': (.24, .24), 'S4': (.17, .17), 'SX': (.05, .05)}
                   },
         'ogu-a': {'time': {'S1': 0, 'S2': 9.4, 'S3': 5.6, 'S4': 20.9, 'SX': 335.6},
                   'etime': {'S1': (0., 0.), 'S2': (1.1, 1.1), 'S3': (.5, .5), 'S4': (2.0, 2.0), 'SX': (21, 21)},
                   'mag': {'S1': 1., 'S2': 17.7 / 15.4, 'S3': 18.3 / 15.4, 'S4': 9.8 / 15.4, 'SX': 4.2 / 15.4},
                   'emag': {'S1': (0, 0), 'S2': (.17, .17), 'S3': (.17, .17), 'S4': (.11, .11), 'SX': (0.03, 0.03)}
                   },
         'gri-g': {'time': {'S1': 0, 'S2': 10.6, 'S3': 4.8, 'S4': 25.9, 'SX': 361.},
                   'etime': {'S1': (0., 0.), 'S2': (3.0, 6.2), 'S3': (1.8, 3.2), 'S4': (4.3, 8.1), 'SX': (27, 19)},
                   'mag': {'S1': 1., 'S2': 0.92, 'S3': 0.99, 'S4': 0.42, 'SX': 0.36},
                   'emag': {'S1': (0, 0), 'S2': (.52, .43), 'S3': (.33, .52), 'S4': (.2, .19), 'SX': (0.09, 0.11)}
                   },
         'sha-g': {'time': {'S1': 0, 'S2': 6., 'S3': -1., 'S4': 12., 'SX': 277.},
                   'etime': {'S1': (0., 0.), 'S2': (5, 6), 'S3': (5, 7), 'S4': (3, 3), 'SX': (21, 11)},
                   'mag': {'S1': 1., 'S2': 0.84, 'S3': 1.68, 'S4': 0.57, 'SX': 0.25},
                   'emag': {'S1': (0, 0), 'S2': (.06, .18), 'S3': (.21, .55), 'S4': (.04, .11), 'SX': (0.02, 0.05)}
                   },
         'sha-a': {'time': {'S1': 0, 'S2': 8., 'S3': 5., 'S4': 17., 'SX': 233.},
                   'etime': {'S1': (0., 0.), 'S2': (5, 7), 'S3': (7, 10), 'S4': (5, 6), 'SX': (13, 46)},
                   'mag': {'S1': 1., 'S2': 0.84, 'S3': 1.46, 'S4': 0.44, 'SX': 0.19},
                   'emag': {'S1': (0, 0), 'S2': (.19, .2), 'S3': (.49, .07), 'S4': (.1, .05), 'SX': (0.04, 0.01)}
                   },
         'die-a': {'time': {'S1': 0, 'S2': -17., 'S3': -4., 'S4': 74., 'SX': 262.},
                   'etime': {'S1': (0., 0.), 'S2': (19, 19), 'S3': (27, 27), 'S4': (43, 43), 'SX': (55, 55)},
                   'emag': {'S1': (0, 0), 'S2': (.79, .79), 'S3': (.19, .19), 'S4': (.11, .11), 'SX': (0.1, 0.1)},
                   'mag': {'S1': 1., 'S2': 1.89, 'S3': 0.64, 'S4': 0.35, 'SX': 0.31}
                   },
         'obs-tmp': {'time': {'S1': 0, 'S2': 4., 'S3': 2., 'S4': 24., 'SX': 0},  # Rodney
                     'etime': {'S1': (0., 0.), 'S2': (4, 4), 'S3': (5, 5), 'S4': (7, 7), 'SX': (0, 0)},
                     'mag': {'S1': 1., 'S2': 1.15, 'S3': 1.01, 'S4': 0.34, 'SX': 0},
                     'emag': {'S1': (0, 0), 'S2': (.05, .05), 'S3': (.04, .04), 'S4': (.02, .02), 'SX': (0, 0)}
                     },
         # 'obs-sn87a': {'time': {'S1': 0, 'S2': -1., 'S3': 0.4, 'S4': 14.1, 'SX': 0},  # Rodney
         #               'mag': {'S1': 1., 'S2': 1.14, 'S3': 1.05, 'S4': 0.34, 'SX': 0}},
         'obs-pol': {'time': {'S1': 0, 'S2': 7., 'S3': 0.6, 'S4': 27., 'SX': 0},  # Rodney
                     'etime': {'S1': (0., 0.), 'S2': (2, 2), 'S3': (3, 3), 'S4': (8, 8), 'SX': (0, 0)},
                     'mag': {'S1': 1., 'S2': 1.17, 'S3': 1., 'S4': 0.38, 'SX': 0},
                     'emag': {'S1': (0, 0), 'S2': (.02, .02), 'S3': (.01, .01), 'S4': (.02, .02), 'SX': (0, 0)}
                     },
         'bakl': {'time': {'S1': 0, 'S2': 4.56, 'S3': -1.79, 'S4': 15.07, 'SX': 0},  # sn_obs_refsdal_grav-lens.ipynb
                  'etime': {'S1': (0., 0.), 'S2': (.25, .25), 'S3': (1., 1.), 'S4': (2.7, 2.7), 'SX': (0, 0)},
                  'mag': {'S1': 1., 'S2': 1.15, 'S3': 1.02, 'S4': 0.37, 'SX': 10},
                  'emag': {'S1': (0, 0), 'S2': (.006, .006), 'S3': (.006, .006), 'S4': (.007, .007), 'SX': (0, 0)}
                  },
         'zero': {'time': {'S1': 0, 'S2': 0., 'S3': 0., 'S4': 0., 'SX': 0.},  # sn_obs_refsdal_grav-lens.ipynb
                  'etime': {'S1': (0., 0.), 'S2': (0., 0.), 'S3': (0., 0.), 'S4': (0., 0.), 'SX': (0., 0.)},
                  'mag': {'S1': 1., 'S2': 1., 'S3': 1., 'S4': 1., 'SX': 1.},
                  'emag': {'S1': (0, 0), 'S2': (0, 0), 'S3': (0, 0), 'S4': (0, 0), 'SX': (0, 0)}

                  }

         # 'obs-tmp': {'time': {'S1': 0, 'S2': -2.1, 'S3': 5.6, 'S4': 22., 'SX': 0},  # Treu
         #             'mag': {'S1': 1., 'S2': 1.09, 'S3': 1.04, 'S4': 0.35, 'SX': 0}}, # Treu
         # 'obs-sn87a': {'time': {'S1': 0, 'S2': -0.8, 'S3': -0.9, 'S4': 14.9, 'SX': 0},
         #               'mag': {'S1': 1., 'S2': 1.13, 'S3': 1.03, 'S4': 0.34, 'SX': 0}},
         # 'obs-pol': {'time': {'S1': 0, 'S2': 8., 'S3': -0.4, 'S4': 30.7, 'SX': 0},  # Treu
         #             'mag': {'S1': 1., 'S2': 1.17, 'S3': 1.01, 'S4': 0.38, 'SX': 0}}
         }
    return a


def coef_time_delay(model, img=None):  # see http://arxiv.org/abs/1510.05750
    a = coef_glens()
    if img is None:
        return a[model]['time']
    return a[model]['time'][img]


def coef_etime(model, img=None):  # see http://arxiv.org/abs/1510.05750
    a = coef_glens()
    if img is None:
        return a[model]['etime']
    return a[model]['etime'][img]


def coef_magnification(model, img=None):
    a = coef_glens()
    if img is None:
        return a[model]['mag']
    return a[model]['mag'][img]


def coef_emagnification(model, img=None):
    a = coef_glens()
    if img is None:
        return a[model]['emag']
    return a[model]['emag'][img]


def coef_time_mag(model, img=None):
    return coef_time_delay(model, img), coef_magnification(model, img)


def plot(ax, dic=None):
    d = path_data
    band_max = 'F160W'

    z = 0.
    arg = []
    bnames = None
    is_lens_shift = False
    jd_shift = None
    im = 'S1'
    glens = grav_lens_def
    t_lc_max = 0.

    if dic is not None:
        arg = dic.get('args', [])
        glens = dic.get('glens', grav_lens_def)
        im = dic.get('image', 'S1')
        jd_shift = dic.get('jd_shift', None)
        bnames = dic.get('bands', None)
        is_lens_shift = dic.get('shift', '') == 'lens'

    if jd_shift is None:
        jd_shift = float(arg.pop(0)) if arg is not None else -56950

    print "Plot Sn Refsdal %s, glens: %s, path: %s, band_max=%s " % (im, glens, d, band_max)

    # plot B-V
    if 'bv' in dic:
        plot_BV(ax=ax, path=d, jd_shift=jd_shift, glens=glens, image=im)
    elif ax is not None:  # plot ubv
        t_lc_max = plot_ubv(ax=ax, path=d, jd_shift=jd_shift, band_max=band_max, glens=glens, image=im, bnames=bnames,
                            is_lens_shift=is_lens_shift)

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

    if False:
        for a, l in coef_glens().items():
            print '\hline %s  & %s  & %s  & %s  & %s  \\ ' % (
                a, l['mag']['S1'], l['mag']['S2'], l['mag']['S3'], l['mag']['S4'])


def get_xy(d, b, colS):
    d2 = d[d[:, 0] == b,]
    x = d2[:, 1]
    y = d2[:, colS]
    is_good = y != 0
    return x[is_good], y[is_good]


def get_band_num(name):
    return float(name[1:4])


def plot_BV(ax, path, jd_shift, glens, image):
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

    # bands = np.unique(lc_data[:, 0])
    bands = str('F814W-F105W-F125W-F160W').split('-')
    time_delay = coef_time_delay(glens)

    colS = sn_images[image]  # S1 - 2, col  S2 - 4, S3 - 6, S4 - 8 col

    for bn1, bn2 in zip(bands[:-1], bands[1:]):
        b1 = get_band_num(bn1)
        b2 = get_band_num(bn2)
        x1, y1, = get_xy(lc_data, b1, colS)
        x1 = x1 + jd_shift - time_delay[image]
        if np.all(y1 <= 0):
            continue
        x2, y2, = get_xy(lc_data, b2, colS)
        x2 = x2 + jd_shift - time_delay[image]
        if np.all(y2 <= 0):
            continue

        if x1[-1] > x2[-1]:
            xx = x2
            if len(x1) > 4:
                y_spline = interpolate.splrep(x1, y1, s=0)
                yy1 = interpolate.splev(xx, y_spline)
            else:
                yy1 = np.interp(xx, x1, y1, 0, 0)
            bv = yy1 - y2
        else:
            xx = x1
            if len(x2) > 4:
                y_spline = interpolate.splrep(x2, y2, k=min(len(x2) - 1, 3), s=0.)
                yy2 = interpolate.splev(xx, y_spline)
            else:
                yy2 = np.interp(xx, x2, y2, 0, 0)
            bv = y1 - yy2

        ax.plot(xx, bv, color=colors_band[bn1], label='%s-%s .' % (bn1, bn2), marker='o', markersize=9, lw=1.5, ls='')


def plot_ubv(ax, path, jd_shift, band_max, glens, image, bnames=None, is_lens_shift=False):
    bands_excluded = ['F435W']
    lc_data, sn_images = read_lc(path)

    bands = np.unique(lc_data[:, 0])
    time_delay = coef_time_delay(glens)
    mu_ratio = coef_magnification(glens)

    marker_s = {'S1': 'o', 'S2': 's', 'S3': 'd', 'S4': '*', 'SX': '+'}
    colS = sn_images[image]  # S1 - 2, col  S2 - 4, S3 - 6, S4 - 8 col

    def m_mu(x):
        return 10 ** (-x * 0.4)

    def mu_m(x):
        return -2.5 * np.log10(x)

    for b in bands:
        bn = 'F%sW' % int(b)
        if bnames is not None:
            if bn not in bnames:
                continue
        if bn in bands_excluded:
            continue

        data = lc_data[lc_data[:, 0] == b,]
        x = data[:, 1] + jd_shift - time_delay[image]
        y = data[:, colS]
        if is_lens_shift:
            y += -mu_m(mu_ratio[image])
        if np.all(y <= 0):
            continue

        yerr = data[:, colS + 1]
        bcolor = colors_band[bn]
        # print "bn: %s bcolor: %s" % (bn, bcolor)
        lower_limits = np.zeros(x.shape, dtype=bool)
        lower_limits[yerr == -1] = True

        # plot values with errors
        with_errors = np.invert(lower_limits)
        xx = x[with_errors]
        if len(xx) > 0:
            yy = y[with_errors]
            yyerr = yerr[with_errors]
            # ax.errorbar(xx, yy, yerr=yyerr, fmt=marker_glens[glens], color=bcolor, label='%s obs.' % bn)
            ax.errorbar(xx, yy, yerr=yyerr, fmt=marker_s[image], color=bcolor, label='%s obs. %s' % (bn, image))

        # plot upper limits
        xxx = x[lower_limits]
        if len(xx) > 0:
            yyy = y[lower_limits]
            yerr = [-0.5 * np.ones(xxx.shape), np.zeros(xxx.shape)]
            ax.errorbar(xxx, yyy, yerr=yerr, lolims=True, fmt='kv', color=bcolor, lw=1.2)

    data = lc_data[lc_data[:, 0] == int(band_max[1:4]),]  # extract maximum for F160W
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


def read_lc(path=path_data):
    # from Rodney_tbl4
    lc_data = np.loadtxt(os.path.join(path, 'rodney_all.csv'), comments='#', skiprows=3)

    # Kelly for SX
    z = np.zeros((np.shape(lc_data)[0], 2))
    lc_data = np.append(lc_data, z, axis=1)
    fs = {'F125W': 'kelly_F125W_SX.csv', 'F160W': 'kelly_F160W_SX.csv'}
    fs = dict((k, os.path.join(path, v)) for k, v in fs.items())
    col_pos = {'S1': 2, 'S2': 4, 'S3': 6, 'S4': 8, 'SX': 10}
    for bn, fname in fs.items():
        # bn = 'F%sW' % int(b)
        b = get_band_num(bn)
        data = np.loadtxt(os.path.join(path, fname), comments='#')
        for row in data:
            r = np.zeros(np.shape(lc_data)[1])
            r[0] = int(b)
            r[1] = row[0]
            r[col_pos['SX']] = row[1]
            r[col_pos['SX'] + 1] = row[2]

            lc_data = np.vstack([lc_data, r])
    return lc_data, col_pos


def read_curves(path=path_data, image='S1', bnames=None):
    curves = SetLightCurve("SN Refsdal, image: %s" % image)

    lc_data, sn_images = read_lc(path)
    bands = np.unique(lc_data[:, 0])
    colS = sn_images[image]  # S1 - 2, col  S2 - 4, S3 - 6, S4 - 8 col

    for ib in bands:
        bname = 'F%sW' % int(ib)
        if bnames is not None and bname not in bnames:
            continue
        data = lc_data[lc_data[:, 0] == ib,]
        t = data[:, 1]
        m = data[:, colS]
        e = data[:, colS + 1]
        if np.all(m <= 0):
            continue
        # filter data
        is_good = m > 0
        time = t[is_good]
        mags = m[is_good]
        yerr = e[is_good]
        lc = LightCurve(bname, time, mags, errs=yerr)
        curves.add(lc)
    return curves
