#!/usr/bin/env python3

import getopt
import os
import sys
from os.path import dirname

import numpy as np
from matplotlib import gridspec
from scipy import interpolate
from scipy.optimize import fmin

import pystella as ps

__author__ = 'bakl'

t_fit_zeta_max = 999.

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))

_colors = ["blue", "cyan", "brown", 'darkseagreen', 'tomato', 'olive', 'orange',
           'skyblue', 'darkviolet']
# colors = {"B-V": "blue", 'B-V-I': "cyan", 'V-I': "brown"}
lntypes = {"B-V": "-", 'B-V-I': "-.", 'V-I': "--",
           "u-g": "-", 'g-r-i': "-.", 'r-i': "--"}
markers = {u'D': u'diamond', 6: u'caretup', u's': u'square', u'x': u'x',
           5: u'caretright', u'^': u'triangle_up', u'd': u'thin_diamond', u'h': u'hexagon1',
           u'+': u'plus', u'*': u'star', u'o': u'circle', u'p': u'pentagon', u'3': u'tri_left',
           u'H': u'hexagon2', u'v': u'triangle_down', u'8': u'octagon', u'<': u'triangle_left'}
markers = list(markers.keys())

epm_coef = {
    'dessart': {
        'B-V': [0.47188, -0.25399, 0.32630],
        'B-V-I': [0.63241, -0.38375, 0.28425],
        'V-I': [0.81662, -0.62896, 0.33852],
        'J-H-K': [0.10786, 1.12374, 0.]
    },
    'eastman': {
        'B-V': [0.674, -0.741, 0.456, 0.11],
        'B-V-I': [0.686, -0.577, 0.316, 0.05],
        'V-I': [0.445, 0.0136, 0., 0.08],
        'J-H-K': [1.45, -0.45, 0.]
    },
    'hamuy': {
        'B-V': [0.7557, -0.8997, 0.5199, 0.048],
        'B-V-I': [0.7336, -0.6942, 0.3740, 0.027],
        'V-I': [0.7013, -0.5304, 0.2646, 0.029],
        'J-H-K': [1.4787, -0.4799, 0., 0.046]
    },
    'bakl': {  # dir /home/bakl/Sn/Release/seb_git/res/tt/tcolor/r500/1
        'B-V': [0.64, -0.69, 0.51],
        'B-V-I': [0.65, -0.55, 0.4],
        'V-I': [0.83, -0.7, 0.41],
        'J-H-K': [-0.7, 2.76, -0.99],
        'g-r': [0.64, -0.65, 0.47],
        'g-r-i': [0.68, -0.62, 0.43],
        'r-i': [0.95, -0.86, 0.46]
    }
}


# def plot_zeta_oneframe(models_dic, set_bands, t_cut=4.9, is_fit=False, is_plot_Tcolor=True,
#                        is_plot_Tnu=True, is_time_points=False):
#     import matplotlib.pyplot as plt
#
#     t_points = [1, 2, 3, 4, 5, 10, 20, 40, 80, 150]
#
#     xlim = [0, 25000]
#     ylim = [0, 2.5]
#     # xlim = [2500, 15000]
#     # ylim = [0, 1.5]
#
#     # setup figure
#     plt.matplotlib.rcParams.update({'font.size': 14})
#     # plt.rc('text', usetex=True)
#     fig = plt.figure(num=None, figsize=(7, 7), dpi=100, facecolor='w', edgecolor='k')
#     gs1 = gridspec.GridSpec(1, 1)
#     gs1.update(wspace=0.3, hspace=0.3, left=0.1, right=0.95)
#     ax = fig.add_subplot(gs1[0, 0])
#     lw = 2.5
#
#     mi = 0
#     ib = 0
#     for mname, mdic in models_dic.items():
#         mi += 1
#         for bset in set_bands:
#             ib += 1
#             if is_plot_Tcolor:
#                 x = mdic[bset]['Tcol']
#                 y = mdic[bset]['zeta']
#                 z = mdic[bset]['time']
#                 x = x[z > t_cut]
#                 y = y[z > t_cut]
#                 z = z[z > t_cut]
#                 bcolor = _colors[ib % (len(_colors) - 1)]
#                 ax.plot(x, y, marker=markers[mi % (len(markers) - 1)], label='%s T_mag %s' % (bset, mname),
#                         markersize=4, color=bcolor, ls="", linewidth=lw)
#                 if is_fit:  # dessart & eastman
#                     xx = x[z > 2.]
#                     yd = zeta_fit(xx, bset, "dessart")
#                     if yd is not None:
#                         ax.plot(xx, yd, color=bcolor, ls="-", linewidth=lw, label='%s Dessart 05' % bset)
#                     yd = zeta_fit(xx, bset, "eastman")
#                     if yd is not None:
#                         ax.plot(xx, yd, color=bcolor, ls="--", linewidth=lw, label='%s Eastman 96' % bset)
#                 if is_time_points:
#                     integers = [np.abs(z - t).argmin() for t in t_points]  # set time points
#                     for (X, Y, Z) in zip(x[integers], y[integers], z[integers]):
#                         ax.annotate('{:.0f}'.format(Z), xy=(X, Y), xytext=(-10, 20), ha='right',
#                                     textcoords='offset points', color=bcolor,
#                                     arrowprops=dict(arrowstyle='->', shrinkA=0))
#                 t_min = z[y[x > 5000.].argmin()]
#                 print("t_min( %s) = %f" % (bset, t_min))
#
#             if is_plot_Tnu:
#                 z = mdic[bset]['time']
#                 xTnu = mdic[bset]['Tnu']
#                 zTeff = mdic[bset]['Teff']
#                 yW = mdic[bset]['W']
#                 xTnu = xTnu[z > t_cut]
#                 yW = yW[z > t_cut]
#                 zTeff = zTeff[z > t_cut]
#                 z = z[z > t_cut]
#                 ax.plot(xTnu, yW, label='T_nu ' + mname, markersize=5, color="magenta",
#                         ls="-.", linewidth=lw)
#                 ax.plot(zTeff, yW, label='T_eff ' + mname, markersize=5, color="red",
#                         ls="-.", linewidth=lw)
#                 if is_time_points:
#                     integers = [np.abs(z - t).argmin() for t in t_points]  # set time points
#                     for (X, Y, Z) in zip(xTnu[integers], yW[integers], z[integers]):
#                         ax.annotate('{:.0f}'.format(Z), xy=(X, Y), xytext=(-10, 20), ha='right',
#                                     textcoords='offset points', color='magenta',
#                                     arrowprops=dict(arrowstyle='->', shrinkA=0))
#                     for (X, Y, Z) in zip(zTeff[integers], yW[integers], z[integers]):
#                         ax.annotate('{:.0f}'.format(Z), xy=(X, Y), xytext=(-10, 20), ha='right',
#                                     textcoords='offset points', color='red',
#                                     arrowprops=dict(arrowstyle='->', shrinkA=0))
#
#     ax.legend(prop={'size': 8})
#     ax.set_xlim(xlim)
#     ax.set_ylim(ylim)
#     ax.set_ylabel(r'$\zeta$')
#     ax.set_xlabel(r'$T_{color}$')
#     # ax.set_title(bset)
#
#     #     plt.title('; '.join(set_bands) + ' filter response')
#     plt.grid()
#     plt.show()


def plot_zeta(models_dic, set_bands, theta, t_points=None, is_time_points_only=False,
              is_plot_Tcolor=True, is_plot_Tnu=True,
              is_fit_bakl=False, fits=None,
              xlim=(0, 18000), ylim=(0, 2.5), tcut=None, t_fit_lim=None):
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    # t_points = [1, 2, 3, 4, 5, 10, 30, 80, 150]

    # setup figure
    plt.matplotlib.rcParams.update({'font.size': 14})
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    fig = plt.figure(num=len(set_bands), figsize=(9, 9), dpi=100, facecolor='w', edgecolor='k')
    # fig.set_tight_layout(False)
    gs1 = gridspec.GridSpec(len(set_bands) // 2 + len(set_bands) % 2, 2)
    gs1.update(wspace=0., hspace=0., left=0.1, right=0.9)

    ax_cache = {}

    # create the grid of figures
    for ib, bset in enumerate(set_bands):
        # ib += 1
        icol = int(ib % 2)
        irow = int(ib / 2)
        ax = fig.add_subplot(gs1[irow, icol])
        ax_cache[bset] = ax

        if icol > 0:
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position("right")
        ax.set_ylabel(r'$\zeta$')

        if irow == 1:
            ax.set_xlabel(r'$T_{color}$')

        ax.set_xlim(xlim)
        xstart, xend = 0, 20000.
        ax.xaxis.set_ticks(np.arange(5000, xend, (xend - xstart) / 4.))
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
        # y
        ax.set_ylim(ylim)
        ax.set_yticks(np.arange(0.5, ylim[1], 0.5))
        # ax.text(.5, .9, bset, horizontalalignment='center', transform=ax.transAxes)
        # ax.set_title(bset)
        ax.set_title(bset, x=0.5, y=0.9)

    if len(models_dic) > 0:
        for ib, bset in enumerate(set_bands):
            mi = 0
            for mname, tbl in models_dic[bset].items():
                ax = ax_cache[bset]
                mi += 1
                # filter data
                if tcut is not None:
                    tbl = table_cut_by_col(tbl, tcut, 'time')

                if is_plot_Tcolor:
                    x = tbl['Tcol']
                    y = tbl['zeta']
                    z = tbl['time']
                    bcolor = "grey"
                    # bcolor = _colors[ib % (len(_colors) - 1)]

                    if t_points is not None:
                        integers = [np.abs(z - t).argmin() for t in t_points]  # set time points
                        if mi % 10 == 1:
                            # integers = [np.abs(z - t).argmin() for t in t_points]  # set time points
                            for (X, Y, Z) in zip(x[integers], y[integers], z[integers]):
                                ax.annotate('{:.0f}'.format(Z), xy=(X, Y), xytext=(20, 20), ha='right',
                                            textcoords='offset points', color=bcolor,
                                            arrowprops=dict(arrowstyle='->', shrinkA=0))
                                # t_min = z[y[x > 5000.].argmin()]
                                # print "t_min( %s) = %f" % (bset, t_min)
                        if is_time_points_only:
                            ax.plot(x[integers], y[integers], marker=markers[mi % (len(markers) - 1)], label='',
                                    markersize=5, color=bcolor, ls="", linewidth=1.5)
                        else:
                            ax.plot(x, y, marker=markers[mi % (len(markers) - 1)], label='',
                                    # mname, # label='T_mag ' + mname,
                                    markersize=5, color=bcolor, ls="", linewidth=1.5)
                    else:
                        ax.plot(x, y, marker=markers[mi % (len(markers) - 1)], label='',
                                # mname, # label='T_mag ' + mname,
                                markersize=5, color=bcolor, ls="", linewidth=1.5)

                if is_plot_Tnu:
                    z = tbl['time']
                    xTnu = tbl['Tnu']
                    zTeff = tbl['Teff']
                    yW = tbl['W']
                    ax.plot(xTnu, yW, label='T_nu ' + mname, markersize=5, color="magenta",
                            ls="-", linewidth=2.5)
                    ax.plot(zTeff, yW, label='T_eff ' + mname, markersize=5, color="red",
                            ls="-", linewidth=2.5)
                    if t_points is not None:
                        if is_time_points_only:
                            integers = [np.abs(z - t).argmin() for t in t_points]  # set time points
                            for (X, Y, Z) in zip(xTnu[integers], yW[integers], z[integers]):
                                ax.annotate('{:.0f}'.format(Z), xy=(X, Y), xytext=(-10, 20), ha='right',
                                            textcoords='offset points', color='magenta',
                                            arrowprops=dict(arrowstyle='->', shrinkA=0))
                        else:
                            for (X, Y, Z) in zip(zTeff, yW, z):
                                ax.annotate('{:.0f}'.format(Z), xy=(X, Y), xytext=(-10, 20), ha='right',
                                            textcoords='offset points', color='red',
                                            arrowprops=dict(arrowstyle='->', shrinkA=0))

    # find & plot fit zeta-Tcol for Stella
    if fits is not None:  # dessart, eastman, hamuy
        xx = np.linspace(max(100, xlim[0]), xlim[1], num=50)
        for bset in set_bands:
            ax = ax_cache[bset]
            if "hamuy" in fits:
                yh = zeta_fit(xx, bset, "hamuy")
                if yh is not None:
                    bcolor = "skyblue"
                    ax.plot(xx, yh, color=bcolor, ls="-.", linewidth=2.5, label='Hamuy 01')

            if "eastman" in fits:
                ye = zeta_fit(xx, bset, "eastman")
                if ye is not None:
                    bcolor = "tomato"
                    ax.plot(xx, ye, color=bcolor, ls="-.", linewidth=2.5, label='Eastman 96')

            if "dessart" in fits:
                yd = zeta_fit(xx, bset, "dessart")
                if yd is not None:
                    bcolor = "darkviolet"
                    ax.plot(xx, yd, color=bcolor, ls="--", linewidth=2.5, label='Dessart 05')

            if "bakl" in fits:
                yb = zeta_fit(xx, bset, "bakl")
                if yb is not None:
                    bcolor = "orange"
                    ax.plot(xx, yb, color=bcolor, ls="--", linewidth=2.5, label='Baklanov')

            # PRINT coef
            print(bset + "  a_i coefficients")
            if zeta_fit_coef_exists(bset, 'dessart'):
                for nm in ["eastman", "hamuy", "dessart", 'bakl']:
                    if nm in fits:
                        print("  {:8s}: {}".format(nm, ' '.join(map(str, zeta_fit_coef(bset, nm)))))

    # new fit
    if theta is not None:
        xx = np.linspace(max(100, xlim[0]), xlim[1], num=50)
        for bset in set_bands:
            ax = ax_cache[bset]
            yf = zeta_fit_rev_temp(xx, theta[bset]['v'])
            bcolor = "orange"
            ax.plot(xx, yf, color=bcolor, dashes=[12, 6, 12, 6, 3, 6], linewidth=2.5, label=r'MCMC $\zeta-T$')

            # PRINT coef
            print(bset + "  a_i coefficients")
            print("  {:8s}: {}".format('MCMC', ', '.join(map(str, np.round(theta[bset]['v'], 2)))))
            # print("    MCMC zeta-T %s: %s " % (bset, ' '.join(map(str, np.round(theta_dic[bset]['v'], 4)))))

    if is_fit_bakl:  # bakl fit
        # find a_coef
        a = {}
        err = {}
        total_zt = models_join(models_dic)
        for bset, tbl in total_zt.items():
            # filter data
            if t_fit_lim is not None:
                tbl = table_cut_by_col(tbl, t_fit_lim, 'time')
                if t_fit_lim[1] == t_fit_zeta_max:
                    idx = np.argmax(tbl['zeta'])
                    cut = np.zeros(len(tbl['zeta']), dtype=bool)
                    cut[:idx] = True
                    tbl = tbl[cut]

            a[bset], err[bset] = zeta_fit_coef_my(tbl)
            # print "%s & %s " % (bset, ', '.join([str(round(x, 4)) for x in a[bset]]))
            print(" Baklan zeta-T  {}: {} : err {:.3f}".
                  format(bset, ' '.join([str(round(x, 4)) for x in a[bset]]), err[bset]))
            # print " Baklan errors  %s: %s " % (bset, ' '.join([str(round(x, 4)) for x in ]))
            # print("")

            # show fit
            xx = np.linspace(max(100, xlim[0]), xlim[1], num=50)
            bcolor = "orange"
            ax = ax_cache[bset]
            yb = zeta_fit_rev_temp(xx, a[bset])
            if yb is not None:
                ax.plot(xx, yb, color=bcolor, ls="--", linewidth=2., label=' Baklan')

    # legend
    for bset in set_bands:
        ax_cache[bset].legend(prop={'size': 6})

    # plt.title('; '.join(set_bands) + ' filter response')
    # plt.grid()
    return fig


def plot_fits(set_bands, is_grid=True, used=('dessart', 'eastman', 'hamuy', 'bakl'),
              xlim=(0, 20000), ylim=(0, 2)):
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    plt.matplotlib.rcParams.update({'font.size': 14})

    fig = plt.figure(num=len(set_bands), figsize=(9, 9), dpi=100, facecolor='w', edgecolor='k')
    if is_grid:
        gs1 = gridspec.GridSpec(len(set_bands) // 2 + len(set_bands) % 2, 2)
    else:
        gs1 = gridspec.GridSpec(1, 1)
    gs1.update(wspace=0.3, hspace=0.3, left=0.15, right=0.95)

    ax = fig.add_subplot(gs1[0, 0])
    for ib, bset in enumerate(set_bands):
        if is_grid and ib > 0:
            icol = int(ib % 2)
            irow = int(ib / 2)
            ax = fig.add_subplot(gs1[irow, icol])

        bcolor = _colors[(ib + 1) % (len(_colors) - 1)]

        # figure parameters
        xstart, xend = xlim[0], 20000.
        ax.xaxis.set_ticks(np.arange(2000, xend, (xend - xstart) / 4.))
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_ylabel(r'$\zeta(' + bset + ')$')
        ax.set_xlabel(r'$T_{color}$')
        ax.set_title(bset)

        # lines
        opts = {
            'dessart': {'c': "darkviolet", 'ls': "--", 'lbl': 'Dessart 05'},
            'eastman': {'c': "tomato", 'ls': "-.", 'lbl': 'Eastman 96'},
            'hamuy': {'c': "skyblue", 'ls': "-.", 'lbl': 'Hamuy 01'},
            'bakl': {'c': "orange", 'ls': "-", 'lbl': 'Baklanov 18'},
        }
        xx = np.linspace(max(100, xlim[0]), xlim[1], num=50)
        for mn in used:
            yy = zeta_fit(xx, bset, mn)
            if yy is not None:
                if is_grid:
                    bcolor = opts[mn]['c']
                ax.plot(xx, yy, color=bcolor, ls=opts[mn]['ls'], linewidth=2.5, label=opts[mn]['lbl'])

        # mn = "dessart"
        # if mn in used:
        #     xx = np.linspace(max(100, xlim[0]), xlim[1], num=50)
        #     yd = zeta_fit(xx, bset, mn)
        #     if yd is not None:
        #         if is_grid:
        #             bcolor = "darkviolet"
        #         ax.plot(xx, yd, color=bcolor, ls="--", linewidth=2.5, label='Dessart 05')
        #
        # ye = zeta_fit(xx, bset, "eastman")
        # if ye is not None:
        #     if is_grid:
        #         bcolor = "tomato"
        #     ax.plot(xx, ye, color=bcolor, ls="-.", linewidth=2.5, label='Eastman 96')
        #
        # if False:
        #     yh = zeta_fit(xx, bset, "hamuy")
        #     if yh is not None:
        #         if is_grid:
        #             bcolor = "skyblue"
        #         ax.plot(xx, yh, color=bcolor, ls="-.", linewidth=2.5, label='Hamuy 01')
        #
        # yb = zeta_fit(xx, bset, "bakl")
        # if yb is not None:
        #     if is_grid:
        #         bcolor = "orange"
        #     ax.plot(xx, yb, color=bcolor, ls="-", linewidth=2.5, label='Baklanov 18')

        ax.legend(prop={'size': 6})
        # if is_grid:
        #     icol = (ib+1) % 2
        #     irow = int((ib) / 2)
        #     ax = fig.add_subplot(gs1[irow, icol])

    return fig


def zeta_fit(Tcol, bset, src):
    """
    Zeta fit from Dessart, L., & Hillier, D. J. (2005). doi:10.1051/0004-6361:20053217
    :param Tcol:
    :param bset:
    :param src:
    :return:
    """
    a = zeta_fit_coef(bset, src)
    if a is not None:
        zeta = zeta_fit_rev_temp(Tcol, a[:3])
        return zeta
    return None


def zeta_fit_coef_exists(bset, src=None):
    if (src is not None) and (src not in epm_coef):
        return False

    if src is not None:
        return bset in epm_coef[src].keys()

    for src, v in epm_coef.items():
        if bset in v.keys():
            return True
    return False


def zeta_fit_coef(bset, src):
    """
    Coefficients for zeta fit from Dessart, L., & Hillier, D. J. (2005). doi:10.1051/0004-6361:20053217
    :param src:
    :param bset:
    :return:
    """

    if src not in epm_coef:
        return None
    if bset not in epm_coef[src]:
        return None
    return epm_coef[src][bset]


def models_join(bands_models):
    res = {}
    for bset, models in bands_models.items():
        total_bset = None
        for mname, tbl in models.items():
            if total_bset is None:
                total_bset = tbl
            else:
                total_bset = np.concatenate((total_bset, tbl))
        # sort by time
        total_bset = np.sort(total_bset, order='time')
        res[bset] = total_bset
    return res


def zeta_fit_rev_temp(T, a_coef):
    z = 0
    i = 0
    for ai in a_coef:
        z += ai * (1.e4 / T) ** i
        i += 1
    return z


def zeta_fit_coef_my(tbl_zt):
    """
    Zeta fit for Stella model data
    :param tbl_zt:  table with data
    :return:
    """
    a_init = [0.5, -0.5, 0.]
    # a = fmin(epsilon_fit_zeta, x0=a_init, args=tuple([tbl_zt]), disp=0)
    a = fmin(epsilon_fit_zeta, x0=a_init, args=tuple([tbl_zt]), disp=False)
    err = epsilon_fit_zeta(a, tbl_zt)
    return a, err


def epsilon_fit_zeta(x, tbl):
    Tcol = tbl['Tcol']
    zeta = tbl['zeta']

    z_fit = zeta_fit_rev_temp(Tcol, x)
    e = np.sum((zeta - z_fit) ** 2) / len(zeta)
    # print "epsilon: err=%f" % e
    return e


def epsilon(theta, freq, mag, bands, radius, dist, z):
    temp_color, zeta = theta
    e = 0
    if temp_color < 0 or zeta < 0:
        for b in bands:
            e += mag[b] ** 2
        return e
    # sp = spectrum.SpectrumDilutePlanck(freq, temp_color, W=zeta**2)
    sp = ps.SpectrumDilutePlanck(freq, temp_color, W=zeta ** 2)
    # sp.correct_zeta(zeta)

    star = ps.Star("bb", sp)
    star.set_radius_ph(radius)
    star.set_distance(dist)
    star.set_redshift(z)
    mag_bb = {b: star.magAB(ps.band.band_by_name(b)) for b in bands}
    for b in bands:
        e += (mag[b] - mag_bb[b]) ** 2
    return e


def compute_Tcolor_zeta(mags, tt, bands, freq, d, z):
    from scipy.interpolate import InterpolatedUnivariateSpline

    temp = list()
    zeta_radius = list()
    times = list()
    # Rph_spline = interpolate.splrep(tt['time'], tt['Rph'], s=0)
    Rph_spline = InterpolatedUnivariateSpline(tt['time'], tt['Rph'], k=1)
    lc_time = mags["time"]

    for nt in range(len(lc_time)):
        t = lc_time[nt]
        if t < min(tt['time']):
            continue
        if t > max(tt['time']):
            break
        mag = {b: mags[b][nt] for b in bands}
        # radius = interpolate.splev(t, Rph_spline)
        radius = Rph_spline(t)
        # res = minimize(lambda x: epsilon(x, freq, mag, bands, radius, d, z),
        #                x0=[1.e4, 1], method='Nelder-Mead', tol=1e-4)
        # tcolor, w = res.x
        tcolor, w = fmin(epsilon, x0=np.array([1.e4, 1.]), args=(freq, mag, bands, radius, d, z), disp=False)
        temp.append(tcolor)
        zeta_radius.append(w)
        times.append(t)
    return temp, zeta_radius, times


def compute_Tnu_w(serial_spec, tt):
    temp_nu = list()
    temp_eff = list()
    W = list()
    Rph_spline = interpolate.splrep(tt['time'], tt['Rph'], s=0)
    x_bb = ps.rf.compute_x_bb()
    for nt in range(len(serial_spec.Time)):
        t, spec = serial_spec.get_tspec(nt)
        if t < min(tt['time']):
            continue
        if t > max(tt['time']):
            break
        Hnu = spec.compute_flux_nu_bol()
        H = spec.compute_flux_bol()
        nu_bb = Hnu / H
        radius = interpolate.splev(t, Rph_spline)
        H /= 4. * np.pi * radius ** 2

        Tnu = ps.phys.h / ps.phys.k * nu_bb / x_bb
        Teff = (H / ps.phys.sigma_SB) ** 0.25
        dilution = (Teff / Tnu) ** 4

        temp_nu.append(Tnu)
        temp_eff.append(Teff)
        W.append(dilution)
    return temp_nu, temp_eff, W


def compute_tcolor(name, path, bands, d=ps.phys.pc2cm(10.), z=0., t_cut=(1., np.inf), t_diff=1.1, ):
    model = ps.Stella(name, path=path)

    if not model.is_ph:
        print("No ph-data for: " + str(model))
        return None

    if not model.is_tt:
        print("No tt-data for: " + str(model))
        return None

    # serial_spec = model.read_series_spectrum(t_diff=1.)
    # curves = serial_spec.flux_to_curves(bands, d=distance)
    serial_spec = model.get_ph(t_diff=t_diff, t_beg=t_cut[0], t_end=t_cut[1])
    # # curves = serial_spec.
    mags = serial_spec.mags_bands(bands, z=z, d=d)

    # curves = model.curves(bands, z=z, distance=d)
    # read R_ph
    tt = model.get_tt().load()
    tt = tbl_rm_equal_el(tt, 'time')
    tt = tt[np.logical_and(t_cut[0] <= tt['time'], tt['time'] <= t_cut[1])]  # time cut  days

    # compute Tnu, W
    Tnu, Teff, W = compute_Tnu_w(serial_spec, tt=tt)

    # fit mags by B(T_col) and get \zeta\theta & T_col
    Tcolors, zetaR, times = compute_Tcolor_zeta(mags, tt=tt, bands=bands, freq=serial_spec.Freq, d=d, z=z)

    # show results
    res = np.zeros(len(Tcolors), dtype=np.dtype({'names': ['time', 'Tcol', 'zeta', 'Tnu', 'Teff', 'W'],
                                                 'formats': [np.float64] * 6}))
    res['time'] = times
    res['Tcol'] = Tcolors
    res['zeta'] = zetaR
    res['Tnu'] = Tnu
    res['Teff'] = Teff
    res['W'] = W

    return res


def log_prior(theta):
    if np.any(theta > 5):
        return -np.inf
    if np.any(theta < -5):
        return -np.inf
    return 0  # flat prior


def log_likelihood_t(theta, tbl_zt):
    Tcol = tbl_zt['Tcol']
    zeta = tbl_zt['zeta']
    z_fit = zeta_fit_rev_temp(Tcol, theta)
    e = 0.01 * zeta
    lhood = -0.5 * np.sum(np.log(2 * np.pi * (e ** 2)) + (zeta - z_fit) ** 2 / e ** 2)
    # print "epsilon: err=%f" % l
    return lhood


def log_posterior(theta, zt):
    p = log_prior(theta)
    for k, v in zt.items():
        p += log_likelihood_t(theta, v)
    return p


def fit_bayesian_models(models_zt, is_debug=True, is_info=False, title='', threads=2):
    import emcee

    # if is_fit_zeta_max:
    #     models = {}
    #     for mname, tbl in models_zt.items():
    #         idx = np.argmax(tbl['zeta'])
    #         cut = np.zeros(len(tbl['zeta']), dtype=bool)
    #         cut[:idx] = True
    #         models[mname] = tbl[cut]
    # else:
    #     models = models_zt

    ndim = 3  # number of parameters in the model
    # run
    nwalkers = 50  # number of MCMC walkers
    nburn = 100  # "burn-in" period to let chains stabilize
    nsteps = 1000  # number of MCMC steps to take
    if is_debug:  # debug
        nwalkers = 20
        nburn = 50  # "burn-in" period to let chains stabilize
        nsteps = 300  # number of MCMC steps to take

    # we'll start at random locations between 0 and 2000
    # starting_guesses = np.zeros(nwalkers, ndim)
    starting_guesses = np.random.rand(nwalkers, ndim)
    starting_guesses[0] = 1.

    # fit
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[models_zt], threads=threads)
    sampler.run_mcmc(starting_guesses, nsteps)
    # plot
    if is_info:
        import matplotlib.pyplot as plt
        # Plot
        fig = plt.figure(num=None, figsize=(6, ndim * 2), dpi=100, facecolor='w', edgecolor='k')
        gs1 = gridspec.GridSpec(ndim, 1)
        plt.matplotlib.rcParams.update({'font.size': 10})
        fig.suptitle(title)
        for i in range(ndim):
            ax = fig.add_subplot(gs1[i, 0])
            ax.hist(sampler.flatchain[:, i], 100, color="k", histtype="step",
                    alpha=0.3, normed=True, label='a_%d' % i)
            ax.legend()

            # the fraction of steps accepted for each walker.
            print("Mean acceptance fraction: {0:.3f}"
                  .format(np.mean(sampler.acceptance_fraction)))

    a = []
    sigma = []
    for i in range(ndim):
        sample = sampler.chain[:, nburn:, i].ravel()
        a.append(np.mean(sample))
        sigma.append(np.std(sample))

    res = {'v': a, 'e': sigma}
    return res


def fit_bayesian(total_zt, is_debug=True, is_info=False, title='', threads=2):
    import emcee
    ndim = 3  # number of parameters in the model
    # run
    nwalkers = 50  # number of MCMC walkers
    nburn = 200  # "burn-in" period to let chains stabilize
    nsteps = 1000  # number of MCMC steps to take
    if is_debug:  # debug
        nwalkers = 20
        nburn = 50  # "burn-in" period to let chains stabilize
        nsteps = 300  # number of MCMC steps to take

    # we'll start at random locations between 0 and 2000
    # starting_guesses = np.zeros(nwalkers, ndim)
    starting_guesses = np.random.rand(nwalkers, ndim)
    starting_guesses[0] = 1.

    # fit
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood_t, args=[total_zt], threads=threads)
    sampler.run_mcmc(starting_guesses, nsteps)
    # plot
    if is_info:
        import matplotlib.pyplot as plt
        # Plot
        fig = plt.figure(num=None, figsize=(6, ndim * 2), dpi=100, facecolor='w', edgecolor='k')
        gs1 = gridspec.GridSpec(ndim, 1)
        plt.matplotlib.rcParams.update({'font.size': 10})
        fig.suptitle(title)
        for i in range(ndim):
            ax = fig.add_subplot(gs1[i, 0])
            ax.hist(sampler.flatchain[:, i], 100, color="k", histtype="step",
                    alpha=0.3, normed=True, label='a_%d' % i)
            ax.legend()

            # the fraction of steps accepted for each walker.
            print("Mean acceptance fraction: {0:.3f}"
                  .format(np.mean(sampler.acceptance_fraction)))

    a = []
    sigma = []
    for i in range(ndim):
        sample = sampler.chain[:, nburn:, i].ravel()
        a.append(np.mean(sample))
        sigma.append(np.std(sample))

    res = {'v': a, 'e': sigma}
    return res


def table_cut_by_col(tbl, xlim, cname):
    x = tbl[cname]
    x_beg = xlim[0]
    x_end = xlim[1]
    if x_beg is None:
        cut = x <= x_end
    elif x_end is None:
        cut = x >= x_beg
    else:
        cut = (x >= x_beg) & (x <= x_end)

    res = np.zeros(np.sum(cut), dtype=tbl.dtype)
    for col in res.dtype.names:
        res[col] = tbl[col][cut]
    return res


def usage():
    print("Usage:")
    print("  zeta.py [params]")
    print("  -b <set_bands>: delimiter '_'. Default: B-V-I_B-V_V-I")
    print("  -d <distance> [pc].  Default: 10 pc")
    print("  -i <model name>.  Example: cat_R450_M15_Ni007_E7")
    print("  -p <model path(directory)>, default: ./")
    print("  -e <model extension> is used to define model name, default: tt ")
    print("  -s  silence mode: no info, no plot")
    print("  -f  force mode: rewrite zeta-files even if it exists")
    print(
        "  -o  options: <dessart:eastman:hamuy:bakl:fitb:ubv:Tnu> - fit E&D&H&B: fit bakl: show time points: plot UBV")
    print("  -x  <tbeg:tend> - time range for fitting. Special value tend={0}, used as tend=t(zeta_max)"
          " Default: 1:{0} , used all days.".format(int(t_fit_zeta_max)))
    print("  -w  write magnitudes to file, default 'False'")
    print("  -z <redshift>.  Default: 0")
    print("  -h  print usage")
    print("  ---  ")

    ps.band.print_bands()


def print_coef(theta):
    strs = []
    for v, e in zip(theta['v'], theta['e']):
        strs.append("{0:.2f}+/-{1:.4f}".format(v, e))

    # print(""" Results: """)
    for i, s in enumerate(strs):
        print("{0}: {1}".format(i, s))


def fitcoef2str(theta):
    strs = []
    for v, e in zip(theta['v'], theta['e']):
        strs.append("{0:.2f}+/-{1:.4f}".format(v, e))

    # print(""" Results: """)
    res = ''
    for i, s in enumerate(strs):
        res += "{0}: {1}  ".format(i, s)
    return res


def tcut2str(tcut):
    s = ''
    if tcut is not None:
        s = "_tcut"
        s += str(tcut[0]).replace('.', 'p')
        if tcut[1] == np.inf:
            s += '-inf'
        else:
            s += str(tcut[1]).replace('.', 'p')
    return s


def gen_check_increasing(a):
    prev = -np.inf
    for i, t in enumerate(a):
        if prev >= t:
            # print("{} {}".format(i, prev))
            yield i
        prev = t


def tbl_rm_equal_el(tbl, col):
    idx = [i for i in gen_check_increasing(tbl[col])]
    res = np.delete(tbl, idx, None)
    return res


def cache_name(name, path, bands, z=0., tcut=None):
    fname = os.path.join(path, "zeta_%s.%s" % (name, bands))
    if z > 0.:
        fname += ".z%.2f" % z

    fname += tcut2str(tcut)
    #     print("cshe name: {}".format(fname))
    return fname


def main(name='', path='./', is_force=False, is_save=False, is_plot_Tnu=False):
    import pickle as pickle

    fits = None
    is_info = False
    is_save_plot = False
    is_fit_bakl = False
    is_fit_only = False
    model_ext = '.tt'
    theta = None
    distance = ps.phys.pc2cm(10.)  # pc
    z = 0.
    t_cut = (1.5, np.inf)
    xlim = (8., t_fit_zeta_max)
    is_fit_zeta_max = True
    # if bset == 'J-H-K':
    #     is_fit_zeta_max = False

    ps.Band.load_settings()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "fhswtb:d:e:i:p:o:x:z:")
    except getopt.GetoptError as err:
        print(str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(args) > 0:
        path, name = os.path.split(str(args[0]))
        path = os.path.expanduser(path)
    elif len(opts) == 0:
        usage()
        sys.exit(2)

    # if not name:
    #     for opt, arg in opts:
    #         if opt == '-i':
    #             name = os.path.splitext(os.path.basename(str(arg)))[0]
    #             break

    set_bands = ['B-V', 'B-V-I', 'V-I', 'J-H-K']
    # set_bands = ['B-V', 'B-V-I', 'V-I', 'J-H-K']
    # set_bands = ['U-B-V-I', 'U-B-V-R-I', 'U-B', 'V-R', 'B-V-I', 'B-V', 'V-I']

    for opt, arg in opts:
        if opt == '-e':
            model_ext = '.' + arg
            continue
        if opt == '-b':
            set_bands = str(arg).split('_')
            for bset in set_bands:
                for b in bset.split('-'):
                    if not ps.band.is_exist(b):
                        print('No such band: ' + b)
                        sys.exit(2)
            continue
        if opt == '-z':
            z = float(arg)
            continue
        if opt == '-d':
            distance = ps.phys.pc2cm(float(arg))
            continue
        if opt == '-s':
            is_save_plot = True
            continue
        if opt == '-o':
            ops = str(arg).split(':')
            fits = ops
            is_plot_Tnu = "Tnu" in ops
            is_fit_only = "only" in ops
            is_fit_bakl = "fitb" in ops
            continue
        if opt == '-x':
            xlim = ps.str2interval(arg, llim=0, rlim=float('inf'))
            continue
        if opt == '-w':
            is_save = True
            continue
        if opt == '-f':
            is_force = True
            is_save = True
            continue
        if opt == '-p':
            path = os.path.expanduser(str(arg))
            if not (os.path.isdir(path) and os.path.exists(path)):
                print("No such directory: " + path)
                sys.exit(2)
            continue
        elif opt == '-h':
            usage()
            sys.exit(2)

    # Plot only fits
    # , used=('dessart', 'eastman', 'hamuy', 'bakl')
    if is_fit_only:
        fig = plot_fits(set_bands, is_grid=True, used=('hamuy', 'eastman', 'dessart'))
        if is_save_plot:
            fsave = "epm_fit_{}.pdf".format('_'.join(set_bands))
            print("Save plot in %s" % fsave)
            fig.savefig(fsave, bbox_inches='tight', format='pdf')
        else:
            import matplotlib.pyplot as plt
            plt.show()
        return

    # Set model names
    names = []
    if not name:
        for opt, arg in opts:
            if opt == '-i':
                name = os.path.splitext(os.path.basename(str(arg)))[0]
                names.append(name)
    else:
        names.append(name)

    if len(names) == 0:  # run for all files in the path
        names = ps.path.get_model_names(path, model_ext)

    if len(names) > 0:
        results = {}  # dict((k, None) for k in names)
        ib = 0
        is_err = False
        for bset in set_bands:
            ib += 1
            im = 0
            dic = {}
            print("\nRun: %s [%d/%d], z=%.2f, d=%.2e, %s" % (bset, ib, len(set_bands), z, distance, tcut2str(t_cut)))
            for name in names:
                im += 1
                is_err = False
                fname = cache_name(name, path, bset, z, tcut=t_cut)
                fname_pkl = fname + ".pkl"
                fname_csv = fname + ".csv"
                if not is_force and os.path.exists(fname_pkl):
                    print("Load: %s [%d/%d]" % (name, im, len(names)))
                    # res = ps.util.cache_load(fname)
                    with open(fname_pkl, 'rb') as f:
                        res = pickle.load(f)
                else:
                    print("Run: %s [%d/%d]" % (name, im, len(names)))
                    res = compute_tcolor(name, path, bset.split('-'), d=distance, z=z, t_cut=t_cut)
                    if is_save and res is not None:
                        import pandas as pd
                        print("Save Tcolor & Zeta for %s in %s" % (bset, fname_pkl))
                        with open(fname_pkl, 'wb') as output:
                            pickle.dump(res, output, pickle.HIGHEST_PROTOCOL)
                        # ps.util.cache_save(res, fname=fname)
                        df = pd.DataFrame(res)
                        df.to_csv(fname_csv, float_format='%.4f', index=False)
                        print("Saved Tcolor & Zeta for {} in {}".format(bset, fname_csv))

                dic[name] = res
                # check for errors
                idx = np.argmin(np.abs(res['Tcol'] - 1.e4))
                if abs(res['zeta'][idx] - 1.) <= 1e-3:
                    is_err = True
            if is_err:
                print(" ERROR for %s in %s" % (name, bset))

            results[bset] = dic
            print("Finish: %s" % name)

        if fits is not None:
            # join all data
            total_zt = models_join(results)

            theta = {}
            print("Fitting with MCMC")
            if xlim is not None:
                print("   tbeg={}  tend={}".format(xlim[0], xlim[1]))
            for bset, tbl in total_zt.items():
                # filter data
                if xlim is not None:
                    tbl = table_cut_by_col(tbl, xlim, 'time')

                if is_fit_zeta_max:
                    idx = np.argmax(tbl['zeta'])
                    cut = np.zeros(len(tbl['zeta']), dtype=bool)
                    cut[:idx] = True
                    tbl = tbl[cut]

                # print("\nFit: %s [%d/%d]" % (bset, im, len(set_bands)))
                # theta = fit_bayesian(models, is_debug=True, is_info=is_info, title=bset)
                theta_mcmc = fit_bayesian(tbl, is_debug=True, is_info=is_info, title=bset)
                # results_filter[bset] = total_zt[bset]
                theta[bset] = theta_mcmc
                # print_coef(theta)
                print("{}: {}".format(bset, fitcoef2str(theta_mcmc)))

        fig = plot_zeta(results, set_bands, theta=theta, tcut=t_cut, t_fit_lim=xlim, fits=fits,
                        is_fit_bakl=is_fit_bakl, is_plot_Tnu=is_plot_Tnu)

        if is_save_plot and len(results) > 0:
            n = max([len(results[b]) for b in set_bands])
            fsave = "~/zeta_{}__Len{}.pdf".format('_'.join(set_bands), n)
            fsave = os.path.expanduser(fsave)
            print("Save plot to %s " % fsave)
            fig.savefig(fsave, bbox_inches='tight', format='pdf')
        else:
            import matplotlib.pyplot as plt
            plt.show()

    else:
        print("There are no models in the directory: %s with extension: %s " % (path, model_ext))


if __name__ == '__main__':
    main()

    # -b g-r_g-r-i_r-i
    # set_fit_bands = ['B-V', 'B-V-I', 'V-I']
    # set_fit_bands = ['B-V', 'B-V-I', 'V-I', 'J-H-K']
    # # , used=('dessart', 'eastman', 'hamuy', 'bakl')
    # fig = plot_fits(set_fit_bands, is_grid=True, used=('hamuy', 'eastman', 'dessart'))
    # import matplotlib.pyplot as plt
    # plt.show()

# main(name="cat_R1000_M15_Ni007_E15", path="/home/bakl/Sn/Release/seb_git/res/tt",
#      is_force=False, is_save=True, is_plot_time_points=True)
