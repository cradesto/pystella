#!/usr/bin/python3
# -*- coding: utf-8 -*-
import getopt
import os
import sys
from os.path import dirname

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib import gridspec
from scipy import interpolate
from scipy.optimize import fmin

import pystella.rf.rad_func as rf
from pystella.model.stella import Stella
from pystella.rf import band, spectrum
from pystella.rf.star import Star
from pystella.util.path_misc import get_model_names
from pystella.util.phys_var import phys
from pystella.util.string_misc import cache_load, cache_save

__author__ = 'bakl'

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
markers = markers.keys()

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
        'B-V': [0.5948, -0.5891, 0.4784],
        'B-V-I': [0.6416, -0.3788, 0.2955],
        'V-I': [0.9819, -0.7699, 0.3919],
        'J-H-K': [1.331, -0.4201, 0.0891],
        'g-r': [0.69, -0.58, 0.42],
        'g-r-i': [0.72, -0.49, 0.33],
        'r-i': [0.82, -0.49, 0.26]
    }
}


def plot_zeta_oneframe(models_dic, set_bands, t_cut=4.9, is_fit=False, is_plot_Tcolor=True,
                       is_plot_Tnu=True, is_time_points=False):
    t_points = [1, 2, 3, 4, 5, 10, 20, 40, 80, 150]

    xlim = [0, 25000]
    ylim = [0, 2.5]
    # xlim = [2500, 15000]
    # ylim = [0, 1.5]

    # setup figure
    plt.matplotlib.rcParams.update({'font.size': 14})
    # plt.rc('text', usetex=True)
    fig = plt.figure(num=None, figsize=(7, 7), dpi=100, facecolor='w', edgecolor='k')
    gs1 = gridspec.GridSpec(1, 1)
    gs1.update(wspace=0.3, hspace=0.3, left=0.1, right=0.95)
    ax = fig.add_subplot(gs1[0, 0])
    lw = 2.5

    mi = 0
    ib = 0
    for mname, mdic in models_dic.iteritems():
        mi += 1
        for bset in set_bands:
            ib += 1
            if is_plot_Tcolor:
                x = mdic[bset]['Tcol']
                y = mdic[bset]['zeta']
                z = mdic[bset]['time']
                x = x[z > t_cut]
                y = y[z > t_cut]
                z = z[z > t_cut]
                bcolor = _colors[ib % (len(_colors) - 1)]
                ax.plot(x, y, marker=markers[mi % (len(markers) - 1)], label='%s T_mag %s' % (bset, mname),
                        markersize=4, color=bcolor, ls="", linewidth=lw)
                if is_fit:  # dessart & eastman
                    xx = x[z > 2.]
                    yd = zeta_fit(xx, bset, "dessart")
                    if yd is not None:
                        ax.plot(xx, yd, color=bcolor, ls="-", linewidth=lw, label='%s Dessart 05' % bset)
                    yd = zeta_fit(xx, bset, "eastman")
                    if yd is not None:
                        ax.plot(xx, yd, color=bcolor, ls="--", linewidth=lw, label='%s Eastman 96' % bset)
                if is_time_points:
                    integers = [np.abs(z - t).argmin() for t in t_points]  # set time points
                    for (X, Y, Z) in zip(x[integers], y[integers], z[integers]):
                        ax.annotate('{:.0f}'.format(Z), xy=(X, Y), xytext=(-10, 20), ha='right',
                                    textcoords='offset points', color=bcolor,
                                    arrowprops=dict(arrowstyle='->', shrinkA=0))
                t_min = z[y[x > 5000.].argmin()]
                print("t_min( %s) = %f" % (bset, t_min))
        if is_plot_Tnu:
            z = mdic[bset]['time']
            xTnu = mdic[bset]['Tnu']
            zTeff = mdic[bset]['Teff']
            yW = mdic[bset]['W']
            xTnu = xTnu[z > t_cut]
            yW = yW[z > t_cut]
            zTeff = zTeff[z > t_cut]
            z = z[z > t_cut]
            ax.plot(xTnu, yW, label='T_nu ' + mname, markersize=5, color="magenta",
                    ls="-.", linewidth=lw)
            ax.plot(zTeff, yW, label='T_eff ' + mname, markersize=5, color="red",
                    ls="-.", linewidth=lw)
            if is_time_points:
                integers = [np.abs(z - t).argmin() for t in t_points]  # set time points
                for (X, Y, Z) in zip(xTnu[integers], yW[integers], z[integers]):
                    ax.annotate('{:.0f}'.format(Z), xy=(X, Y), xytext=(-10, 20), ha='right',
                                textcoords='offset points', color='magenta',
                                arrowprops=dict(arrowstyle='->', shrinkA=0))
                for (X, Y, Z) in zip(zTeff[integers], yW[integers], z[integers]):
                    ax.annotate('{:.0f}'.format(Z), xy=(X, Y), xytext=(-10, 20), ha='right',
                                textcoords='offset points', color='red',
                                arrowprops=dict(arrowstyle='->', shrinkA=0))

    ax.legend(prop={'size': 8})
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_ylabel(r'$\zeta$')
    ax.set_xlabel(r'$T_{color}$')
    # ax.set_title(bset)

    #     plt.title('; '.join(set_bands) + ' filter response')
    plt.grid()
    plt.show()


def plot_zeta(models_dic, set_bands, theta_dic, t_cut=4.9,
              is_plot_Tcolor=True, is_plot_Tnu=True, is_fit=False,
              is_fit_bakl=False, is_time_points=False):
    t_points = [1, 5, 10, 30, 80, 140]
    # t_points = [1, 2, 3, 4, 5, 10, 30, 80, 150]

    xlim = [0, 18000]
    ylim = [0, 2.5]

    # setup figure
    plt.matplotlib.rcParams.update({'font.size': 14})
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    fig = plt.figure(num=len(set_bands), figsize=(9, 9), dpi=100, facecolor='w', edgecolor='k')
    gs1 = gridspec.GridSpec(len(set_bands) // 2 + len(set_bands) % 2, 2)
    # gs1 = gridspec.GridSpec(2, 4, width_ratios=(8, 1, 8, 1))
    gs1.update(wspace=0., hspace=0., left=0.1, right=0.9)

    ax_cache = {}

    # create the grid of figures
    ib = 0
    for bset in set_bands:
        ib += 1
        icol = (ib - 1) % 2
        irow = (ib - 1) / 2
        ax = fig.add_subplot(gs1[irow, icol])
        ax_cache[bset] = ax
        # ax.legend(prop={'size': 6})
        # x
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

    # plot data
    # i = 0
    ib = 0
    for bset in set_bands:
        ib += 1
        mi = 0
        for mname, tbl in models_dic[bset].iteritems():
            ax = ax_cache[bset]
            mi += 1
            # i += 1

            if is_plot_Tcolor:
                x = tbl['Tcol']
                y = tbl['zeta']
                z = tbl['time']
                x = x[z > t_cut]
                y = y[z > t_cut]
                z = z[z > t_cut]
                # bcolor = "black"
                bcolor = _colors[ib % (len(_colors) - 1)]

                if is_time_points:
                    integers = [np.abs(z - t).argmin() for t in t_points]  # set time points
                    ax.plot(x[integers], y[integers], marker=markers[mi % (len(markers) - 1)], label='',
                            markersize=5, color=bcolor, ls="", linewidth=1.5)
                else:
                    ax.plot(x, y, marker=markers[mi % (len(markers) - 1)], label='',  # mname, # label='T_mag ' + mname,
                            markersize=5, color=bcolor, ls="", linewidth=1.5)

                if is_time_points and mi % 10 == 1:
                    # integers = [np.abs(z - t).argmin() for t in t_points]  # set time points
                    for (X, Y, Z) in zip(x[integers], y[integers], z[integers]):
                        ax.annotate('{:.0f}'.format(Z), xy=(X, Y), xytext=(20, 20), ha='right',
                                    textcoords='offset points', color=bcolor,
                                    arrowprops=dict(arrowstyle='->', shrinkA=0))
                        # t_min = z[y[x > 5000.].argmin()]
                        # print "t_min( %s) = %f" % (bset, t_min)
            if is_plot_Tnu:
                z = tbl['time']
                xTnu = tbl['Tnu']
                zTeff = tbl['Teff']
                yW = tbl['W']
                xTnu = xTnu[z > t_cut]
                yW = yW[z > t_cut]
                yW = np.sqrt(yW)
                zTeff = zTeff[z > t_cut]
                z = z[z > t_cut]
                ax.plot(xTnu, yW, label='T_nu ' + mname, markersize=5, color="magenta",
                        ls="-", linewidth=2.5)
                ax.plot(zTeff, yW, label='T_eff ' + mname, markersize=5, color="red",
                        ls="-", linewidth=2.5)
                if is_time_points:
                    integers = [np.abs(z - t).argmin() for t in t_points]  # set time points
                    for (X, Y, Z) in zip(xTnu[integers], yW[integers], z[integers]):
                        ax.annotate('{:.0f}'.format(Z), xy=(X, Y), xytext=(-10, 20), ha='right',
                                    textcoords='offset points', color='magenta',
                                    arrowprops=dict(arrowstyle='->', shrinkA=0))
                    for (X, Y, Z) in zip(zTeff[integers], yW[integers], z[integers]):
                        ax.annotate('{:.0f}'.format(Z), xy=(X, Y), xytext=(-10, 20), ha='right',
                                    textcoords='offset points', color='red',
                                    arrowprops=dict(arrowstyle='->', shrinkA=0))

    # find & plot fit zeta-Tcol for Stella
    if is_fit:  # dessart, eastman, hamuy
        xx = np.linspace(max(100, xlim[0]), xlim[1], num=50)
        for bset in set_bands:
            ax = ax_cache[bset]
            yd = zeta_fit(xx, bset, "dessart")
            if yd is not None:
                bcolor = "darkviolet"
                ax.plot(xx, yd, color=bcolor, ls="--", linewidth=2.5, label='Dessart 05')

            ye = zeta_fit(xx, bset, "eastman")
            if ye is not None:
                bcolor = "tomato"
                ax.plot(xx, ye, color=bcolor, ls="-.", linewidth=2.5, label='Eastman 96')

            if True:
                yh = zeta_fit(xx, bset, "hamuy")
                if yh is not None:
                    bcolor = "skyblue"
                    ax.plot(xx, yh, color=bcolor, ls="-.", linewidth=2.5, label='Hamuy 01')

            # yb = zeta_fit(xx, bset, "bakl")
            # if yb is not None:
            #     bcolor = "orange"
            #     ax.plot(xx, yb, color=bcolor, ls="-", linewidth=2.5, label='Baklanov 16')

            # new fit
            if theta_dic is not None:
                yf = zeta_fit_rev_temp(xx, theta_dic[bset]['v'])
                bcolor = "orange"
                ax.plot(xx, yf, color=bcolor, dashes=[12, 6, 12, 6, 3, 6], linewidth=2.5, label='Baklanov 16')

        # PRINT coef
        for bset in set_bands:
            if is_fit:
                if zeta_fit_coef_exists(bset, 'dessart'):
                    print("Dessart zeta-T  %s: %s " % (bset, ' '.join(map(str, zeta_fit_coef(bset, "dessart")))))
                if zeta_fit_coef_exists(bset, 'eastman'):
                    print("Eastman zeta-T  %s: %s " % (bset, ' '.join(map(str, zeta_fit_coef(bset, "eastman")))))
                if zeta_fit_coef_exists(bset, 'hamuy'):
                    print("Hamuy01 zeta-T  %s: %s " % (bset, ' '.join(map(str, zeta_fit_coef(bset, "hamuy")))))
                if zeta_fit_coef_exists(bset, 'bakl'):
                    print("Baklanov zeta-T %s: %s " % (bset, ' '.join(map(str, zeta_fit_coef(bset, "bakl")))))
                if theta_dic is not None:
                    print("     New zeta-T %s: %s " % (bset, ' '.join(map(str, np.round(theta_dic[bset]['v'], 4)))))

            if is_fit_bakl:  # bakl fit
                # find a_coef
                a = {}
                err = {}
                # t_beg = t_cut
                t_beg, t_end = 5., 110  # None  # 100.
                a[bset], err[bset] = zeta_fit_coef_my(models_dic, bset, t_beg=t_beg, t_end=t_end)  # todo check t_end
                # print "%s & %s " % (bset, ', '.join([str(round(x, 4)) for x in a[bset]]))
                print(" Baklan zeta-T  %s: %s : err %f" % (
                    bset, ' '.join([str(round(x, 4)) for x in a[bset]]), err[bset]))
                # print " Baklan errors  %s: %s " % (bset, ' '.join([str(round(x, 4)) for x in ]))
                print("")

                # show fit
                xx = np.linspace(max(100, xlim[0]), xlim[1], num=50)
                bcolor = "orange"
                ax = ax_cache[bset]
                yb = zeta_fit_rev_temp(xx, a[bset])
                if yb is not None:
                    ax.plot(xx, yb, color=bcolor, ls="--", linewidth=2., label='baklan fit')

    # legend
    for bset in set_bands:
        ax_cache[bset].legend(prop={'size': 6})

    # plt.title('; '.join(set_bands) + ' filter response')
    # plt.grid()
    plt.show()
    return fig


def plot_fits(set_bands, is_grid=True):
    xlim = [4000, 20000]
    ylim = [0, 1.5]

    plt.matplotlib.rcParams.update({'font.size': 14})
    fig = plt.figure(num=len(set_bands), figsize=(9, 9), dpi=100, facecolor='w', edgecolor='k')
    if is_grid:
        gs1 = gridspec.GridSpec(len(set_bands) // 2 + len(set_bands) % 2, 2)
    else:
        gs1 = gridspec.GridSpec(1, 1)
    gs1.update(wspace=0.3, hspace=0.3, left=0.15, right=0.95)

    ax = fig.add_subplot(gs1[0, 0])
    ib = 1
    for bset in set_bands:
        bcolor = _colors[ib % (len(_colors) - 1)]

        # figure parameters
        xstart, xend = 0, 20000.
        ax.xaxis.set_ticks(np.arange(5000, xend, (xend - xstart) / 4.))
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_ylabel(r'$\zeta(' + bset + ')$')
        ax.set_xlabel(r'$T_{color}$')
        ax.set_title(bset)

        # lines
        xx = np.linspace(max(100, xlim[0]), xlim[1], num=50)

        yd = zeta_fit(xx, bset, "dessart")
        if yd is not None:
            if is_grid:
                bcolor = "darkviolet"
            ax.plot(xx, yd, color=bcolor, ls="--", linewidth=2.5, label='Dessart 05')

        ye = zeta_fit(xx, bset, "eastman")
        if ye is not None:
            if is_grid:
                bcolor = "tomato"
            ax.plot(xx, ye, color=bcolor, ls="-.", linewidth=2.5, label='Eastman 96')

        if False:
            yh = zeta_fit(xx, bset, "hamuy")
            if yh is not None:
                if is_grid:
                    bcolor = "skyblue"
                ax.plot(xx, yh, color=bcolor, ls="-.", linewidth=2.5, label='Hamuy 01')

        yb = zeta_fit(xx, bset, "bakl")
        if yb is not None:
            if is_grid:
                bcolor = "orange"
            ax.plot(xx, yb, color=bcolor, ls="-", linewidth=2.5, label='Baklanov 15')

        ax.legend(prop={'size': 6})
        ib += 1
        if is_grid and ib < len(set_bands):
            icol = (ib - 1) % 2
            irow = (ib - 1) / 2
            ax = fig.add_subplot(gs1[irow, icol])

    plt.show()


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


def zeta_fit_rev_temp(T, a_coef):
    z = 0
    i = 0
    for ai in a_coef:
        z += ai * (1.e4 / T) ** i
        i += 1
    return z


def zeta_fit_coef_my(models_dic, bset, t_beg, t_end=None):
    """
    Zeta fit for Stella model data
    :param models_dic:  model dictionary with data
    :param bset: band set, ex. B-V
    :param t_beg:
    :param t_end:
    :return:
    """
    a_init = zeta_fit_coef(bset, src="bakl")
    if a_init is None:
        a_init = [0.5, 0.5, 0.]
    a = fmin(epsilon_fit_zeta, x0=a_init, args=(models_dic, bset, t_beg, t_end), disp=0)
    err = epsilon_fit_zeta(a, models_dic, bset, t_beg, t_end)
    return a, err


def epsilon_fit_zeta(x, models_dic, bset, t_beg, t_end=None):
    e = 0
    for mname, tbl in models_dic[bset].iteritems():
        Tcol = tbl['Tcol']
        zeta = tbl['zeta']
        z = tbl['time']
        if t_end is None:
            cut = z >= t_beg
        else:
            cut = (z >= t_beg) & (z <= t_end)

        Tcol = Tcol[cut]
        zeta = zeta[cut]

        z_fit = zeta_fit_rev_temp(Tcol, x)
        e += np.sum((zeta - z_fit) ** 2) / len(zeta)
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
    sp = spectrum.SpectrumDilutePlanck(freq, temp_color, W=zeta ** 2)
    # sp.correct_zeta(zeta)

    star = Star("bb", sp)
    star.set_radius_ph(radius)
    star.set_distance(dist)
    star.set_redshift(z)
    mag_bb = {b: star.magAB(band.band_by_name(b)) for b in bands}
    for b in bands:
        e += (mag[b] - mag_bb[b]) ** 2
    return e


def compute_Tcolor_zeta(mags, tt, bands, freq, d, z):
    temp = list()
    zeta_radius = list()
    times = list()
    Rph_spline = interpolate.splrep(tt['time'], tt['Rph'], s=0)
    lc_time = mags["time"]

    for nt in range(len(lc_time)):
        t = lc_time[nt]
        if t < min(tt['time']):
            continue
        if t > max(tt['time']):
            break
        mag = {b: mags[b][nt] for b in bands}
        radius = interpolate.splev(t, Rph_spline)
        # res = minimize(lambda x: epsilon(x, freq, mag, bands, radius, d, z),
        #                x0=[1.e4, 1], method='Nelder-Mead', tol=1e-4)
        # tcolor, w = res.x
        tcolor, w = fmin(epsilon, x0=np.array([1.e4, 1]), args=(freq, mag, bands, radius, d, z), disp=False)
        temp.append(tcolor)
        zeta_radius.append(w)
        times.append(t)
    return temp, zeta_radius, times


def compute_Tnu_w(serial_spec, tt):
    temp_nu = list()
    temp_eff = list()
    W = list()
    Rph_spline = interpolate.splrep(tt['time'], tt['Rph'], s=0)
    x_bb = rf.compute_x_bb()
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

        Tnu = phys.h / phys.k * nu_bb / x_bb
        Teff = (H / phys.sigma_SB) ** 0.25
        dilution = (Teff / Tnu) ** 4

        temp_nu.append(Tnu)
        temp_eff.append(Teff)
        W.append(dilution)
    return temp_nu, temp_eff, W


def compute_tcolor(name, path, bands, d=rf.pc_to_cm(10.), z=0., t_cut=1.):
    model = Stella(name, path=path)

    if not model.is_ph_data:
        print("No ph-data for: " + str(model))
        return None

    if not model.is_tt_data:
        print("No tt-data for: " + str(model))
        return None

    # serial_spec = model.read_series_spectrum(t_diff=1.)
    # curves = serial_spec.flux_to_curves(bands, d=distance)
    serial_spec = model.read_series_spectrum(t_diff=1.05, t_beg=t_cut)
    # curves = serial_spec.
    mags = serial_spec.mags_bands(bands, z=z, d=d)

    # read R_ph
    tt = model.read_tt_data()
    tt = tt[tt['time'] > t_cut]  # time cut  days

    # compute Tnu, W
    Tnu, Teff, W = compute_Tnu_w(serial_spec, tt=tt)

    # fit mags by B(T_col) and get \zeta\theta & T_col
    Tcolors, zetaR, times = compute_Tcolor_zeta(mags, tt=tt, bands=bands, freq=serial_spec.Freq, d=d, z=z)

    # show results
    res = np.array(np.zeros(len(Tcolors)),
                   dtype=np.dtype({'names': ['time', 'Tcol', 'zeta', 'Tnu', 'Teff', 'W'],
                                   'formats': [np.float64] * 6}))
    res['time'] = times
    res['Tcol'] = Tcolors
    res['zeta'] = zetaR
    res['Tnu'] = Tnu
    res['Teff'] = Teff
    res['W'] = W

    return res


def fit_bayesian(models_zt, is_debug=True, is_info=False, title=''):
    import emcee

    def log_prior(theta):
        if np.any(theta > 5):
            return -np.inf
        if np.any(theta < -5):
            return -np.inf
        return 1  # flat prior

    def log_likelihood_t(theta, tbl_zt):
        Tcol = tbl_zt['Tcol']
        zeta = tbl_zt['zeta']
        z_fit = zeta_fit_rev_temp(Tcol, theta)
        e = 0.01 * zeta
        l = -0.5 * np.sum(np.log(2 * np.pi * (e ** 2)) + (zeta - z_fit) ** 2 / e ** 2)
        # print "epsilon: err=%f" % l
        return l

    def log_posterior(theta, models):
        p = log_prior(theta)
        for mname, tbl in models.iteritems():
            p += log_likelihood_t(theta, tbl)
        return p

    ndim = 3  # number of parameters in the model
    nwalkers = 50  # number of MCMC walkers
    if is_debug:  # debug
        nwalkers = 6
    # run
    nburn = 1000  # "burn-in" period to let chains stabilize
    nsteps = 2000  # number of MCMC steps to take
    if is_debug:  # debug
        nburn = 100  # "burn-in" period to let chains stabilize
        nsteps = 200  # number of MCMC steps to take

    # we'll start at random locations between 0 and 2000
    # starting_guesses = np.zeros(nwalkers, ndim)
    starting_guesses = np.random.rand(nwalkers, ndim)
    starting_guesses[0] = 1.

    # fit
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[models_zt])
    sampler.run_mcmc(starting_guesses, nsteps)
    # plot
    if is_info:
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
    print("  -o  options: <fit:fitb:time:ubv:Tnu> - fit E&D: fit bakl:  show time points: plot UBV")
    print("  -w  write magnitudes to file, default 'False'")
    print("  -z <redshift>.  Default: 0")
    print("  -h  print usage")
    print("  ---  ")

    band.print_bands()


def print_coef(theta):
    strs = []
    for v, e in zip(theta['v'], theta['e']):
        strs.append("{0:.2f}+/-{1:.4f}".format(v, e))

    # print(""" Results: """)
    for i, s in enumerate(strs):
        print("{0}: {1}".format(i, s))


def cache_name(name, path, bands, z=0.):
    fname = os.path.join(path, "%s.%s" % (name, bands))
    if z > 0.:
        fname += ".z%.2f" % z
    fname += ".zeta"
    return fname


def main(name='', path='./', is_force=False, is_save=False, is_plot_Tnu=False, is_plot_time_points=False):
    is_info = False
    is_fit = False
    is_save_plot = False
    is_fit_bakl = False
    model_ext = '.tt'
    theta_dic = None
    distance = rf.pc_to_cm(10.)  # pc
    z = 0.

    band.Band.load_settings()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "fhswtb:d:e:i:p:o:z:")
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

    if not name:
        for opt, arg in opts:
            if opt == '-i':
                name = os.path.splitext(os.path.basename(str(arg)))[0]
                break

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
                    if not band.band_is_exist(b):
                        print('No such band: ' + b)
                        sys.exit(2)
            continue
        if opt == '-z':
            z = float(arg)
            continue
        if opt == '-d':
            distance = rf.pc_to_cm(float(arg))
            continue
        if opt == '-s':
            is_save_plot = True
            continue
        if opt == '-o':
            ops = str(arg).split(':')
            is_plot_Tnu = "Tnu" in ops
            is_plot_time_points = "time" in ops
            is_fit = "fit" in ops
            is_fit_bakl = "fitb" in ops
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

    names = []
    if name != '':
        names.append(name)
    else:  # run for all files in the path
        names = get_model_names(path, model_ext)

    if len(names) > 0:
        results = {}  # dict((k, None) for k in names)
        ib = 0
        is_err = False
        for bset in set_bands:
            ib += 1
            im = 0
            dic = {}
            print("\nRun: %s [%d/%d], z=%.2f, d=%.2e" % (bset, ib, len(set_bands), z, distance))
            for name in names:
                im += 1
                print("Run: %s [%d/%d]" % (name, im, len(names)))
                is_err = False
                fname = cache_name(name, path, bset)
                if not is_force and os.path.exists(fname):
                    res = cache_load(fname)
                else:
                    res = compute_tcolor(name, path, bset.split('-'), d=distance, z=z, t_cut=0.9)
                    if is_save and res is not None:
                        print("Save Tcolor & Zeta for %s in %s" % (bset, fname))
                        cache_save(res, fname=fname)

                dic[name] = res
                # check for errors
                idx = np.argmin(np.abs(res['Tcol'] - 1.e4))
                if abs(res['zeta'][idx] - 1.) <= 1e-3:
                    is_err = True
            if is_err:
                print(" ERROR for %s in %s" % (name, bset))

            results[bset] = dic
            print("Finish: %s" % name)
        results_filter = {}
        if is_fit:
            theta_dic = {}
            im = 0
            for bset in set_bands:
                im += 1
                models = results[bset]
                tlim = (7., 80)
                if tlim is not None:
                    models_zt_new = {}
                    for mname, tbl in models.iteritems():
                        models_zt_new[mname] = table_cut_by_col(tbl, tlim, 'time')
                    models = models_zt_new

                # templim = (4.e3, 12e3)
                # if templim is not None:
                #     models_zt_new = {}
                #     for mname, tbl in models.iteritems():
                #         models_zt_new[mname] = table_cut_by_col(tbl, templim, 'Tcol')
                #     models = models_zt_new

                print("\nFit: %s [%d/%d]" % (bset, im, len(set_bands)))
                theta = fit_bayesian(models, is_debug=True, is_info=is_info, title=bset)
                results_filter[bset] = models
                theta_dic[bset] = theta
                print_coef(theta)

        # fig = plot_zeta(results, set_bands, theta_dic, t_cut=1.9, is_fit=is_fit, is_fit_bakl=is_fit_bakl,
        fig = plot_zeta(results_filter, set_bands, theta_dic, t_cut=1.9, is_fit=is_fit,
                        is_fit_bakl=is_fit_bakl, is_plot_Tnu=is_plot_Tnu,
                        is_time_points=is_plot_time_points)
        if is_save_plot and len(results) > 0:
            fsave = os.path.join(os.path.expanduser('~/'), 'epm_' + '_'.join(set_bands) + '.pdf')
            print("Save plot in %s" % fsave)
            fig.savefig(fsave, bbox_inches='tight', format='pdf')

    else:
        print("There are no models in the directory: %s with extension: %s " % (path, model_ext))


if __name__ == '__main__':
    main()

# -b g-r_g-r-i_r-i
# set_bands = ['B-V', 'B-V-I', 'V-I']
# set_bands = ['B-V', 'B-V-I', 'V-I', 'J-H-K']
# plot_fits(set_bands, is_grid=False)

# main(name="cat_R1000_M15_Ni007_E15", path="/home/bakl/Sn/Release/seb_git/res/tt",
#      is_force=False, is_save=True, is_plot_time_points=True)
