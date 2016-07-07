#!/usr/bin/python
# -*- coding: utf-8 -*-
import getopt
import numpy as np
import os
import sys
from os.path import isfile, join, dirname
from scipy import interpolate
from scipy.optimize import fmin, minimize

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import gridspec

import pystella.rf.rad_func as rf
from pystella.model.stella import Stella
from pystella.rf import band, spectrum
from pystella.rf.star import Star
from pystella.util.phys_var import phys
from pystella.util.string_misc import cache_load, cache_name, cache_save

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


def plot_zeta_oneframe(models_dic, set_bands, t_cut=4.9, is_fit=False, is_fit_bakl=False,
                       is_plot_Tcolor=True, is_plot_Tnu=True, is_time_points=False):
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
                print "t_min( %s) = %f" % (bset, t_min)
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


def plot_zeta(models_dic, set_bands, t_cut=4.9,
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
    gs1 = gridspec.GridSpec(len(set_bands) / 2 + len(set_bands) % 2, 2)
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
    mi = 0
    i = 0
    for mname, mdic in models_dic.iteritems():
        mi += 1
        ib = 0
        for bset in set_bands:
            ib += 1
            i += 1
            if bset in ax_cache:
                ax = ax_cache[bset]
            else:
                icol = (ib - 1) % 2
                irow = (ib - 1) / 2
                ax = fig.add_subplot(gs1[irow, icol])
                ax_cache[bset] = ax

            if is_plot_Tcolor:
                x = mdic[bset]['Tcol']
                y = mdic[bset]['zeta']
                z = mdic[bset]['time']
                x = x[z > t_cut]
                y = y[z > t_cut]
                z = z[z > t_cut]
                # bcolor = "black"
                bcolor = _colors[ib % (len(_colors) - 1)]

                if is_time_points:
                    integers = [np.abs(z - t).argmin() for t in t_points]  # set time points
                    ax.plot(x[integers], y[integers], marker=markers[mi % (len(markers) - 1)], label='',  # mname,   # label='T_mag ' + mname,
                            markersize=5, color=bcolor, ls="", linewidth=1.5)
                else:
                    ax.plot(x, y, marker=markers[mi % (len(markers) - 1)], label='',  # mname,   # label='T_mag ' + mname,
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
                z = mdic[bset]['time']
                xTnu = mdic[bset]['Tnu']
                zTeff = mdic[bset]['Teff']
                yW = mdic[bset]['W']
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

    if is_fit:  # dessart, eastman, hamuy
        xx = np.linspace(max(100, xlim[0]), xlim[1], num=50)
        for bset in set_bands:
            ax = ax_cache[bset]
            yd = zeta_fit(xx, bset, "dessart")
            bcolor = "darkviolet"
            if yd is not None:
                ax.plot(xx, yd, color=bcolor, ls="--", linewidth=2.5, label='Dessart 05')

            ye = zeta_fit(xx, bset, "eastman")
            bcolor = "tomato"
            if ye is not None:
                ax.plot(xx, ye, color=bcolor, ls="-.", linewidth=2.5, label='Eastman 96')

            if True:
                yh = zeta_fit(xx, bset, "hamuy")
                bcolor = "skyblue"
                if yh is not None:
                    ax.plot(xx, yh, color=bcolor, ls="-.", linewidth=2.5, label='Hamuy 01')

            yb = zeta_fit(xx, bset, "bakl")
            bcolor = "orange"
            if yb is not None:
                ax.plot(xx, yb, color=bcolor, ls="-", linewidth=2.5, label='Baklanov 16')

    # find & plot fit zeta-Tcol for Stella
    if is_fit_bakl:  # bakl fit
        # find a_coef
        a = {}
        err = {}
        # t_beg = t_cut
        t_beg, t_end = 5., 110  # None  # 100.
        for bset in set_bands:
            a[bset], err[bset] = zeta_fit_coef_my(models_dic, bset, t_beg=t_beg, t_end=t_end)  # todo check t_end

        # PRINT coef
        for bset in set_bands:
            # print "%s & %s " % (bset, ', '.join([str(round(x, 4)) for x in a[bset]]))
            print " Baklan zeta-T  %s: %s : err %f" % (bset, ' '.join([str(round(x, 4)) for x in a[bset]]), err[bset])
            # print " Baklan errors  %s: %s " % (bset, ' '.join([str(round(x, 4)) for x in ]))
            if is_fit:
                print "Dessart zeta-T  %s: %s " % (bset, ' '.join(map(str, zeta_fit_coef(bset, "dessart"))))
                print "Eastman zeta-T  %s: %s " % (bset, ' '.join(map(str, zeta_fit_coef(bset, "eastman"))))
                print "Hamuy01 zeta-T  %s: %s " % (bset, ' '.join(map(str, zeta_fit_coef(bset, "hamuy"))))
                print "Baklanov zeta-T %s: %s " % (bset, ' '.join(map(str, zeta_fit_coef(bset, "bakl"))))
                print ""

        # show fit
        xx = np.linspace(max(100, xlim[0]), xlim[1], num=50)
        bcolor = "orange"
        for bset in set_bands:
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


def plot_fits(set_bands, is_grid=True):
    xlim = [4000, 20000]
    ylim = [0, 1.5]

    plt.matplotlib.rcParams.update({'font.size': 14})
    fig = plt.figure(num=len(set_bands), figsize=(9, 9), dpi=100, facecolor='w', edgecolor='k')
    if is_grid:
        gs1 = gridspec.GridSpec(len(set_bands) / 2 + len(set_bands) % 2, 2)
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


def zeta_fit_coef(bset, src):
    """
    Coefficients for zeta fit from Dessart, L., & Hillier, D. J. (2005). doi:10.1051/0004-6361:20053217
    :param src:
    :param bset:
    :return:
    """
    a = {
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
            'J-H-K': [1.331, -0.4201, 0.0891]
        }
    }
    if src not in a:
        return None
    if bset not in a[src]:
        return None
    return a[src][bset]


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
    for mname, mdic in models_dic.iteritems():
        Tcol = mdic[bset]['Tcol']
        zeta = mdic[bset]['zeta']
        z = mdic[bset]['time']
        if t_end is None:
            cut = z >= t_beg
        else:
            cut = (z >= t_beg) & (z <= t_end)

        Tcol = Tcol[cut]
        zeta = zeta[cut]

        z_fit = zeta_fit_rev_temp(Tcol, x)
        e += np.sum((zeta - z_fit)**2) / len(zeta)
    # print "epsilon: err=%f" % e
    return e


def epsilon(x, freq, mag, bands, radius, dist):
    temp_color, zeta = x
    e = 0
    if temp_color < 0 or zeta < 0:
        for b in bands:
            e += mag[b] ** 2
        return e
    sp = spectrum.SpectrumDilutePlanck(freq, temp_color, W=zeta**2)
    # sp.correct_zeta(zeta)

    star = Star("bb", sp)
    star.set_radius_ph(radius)
    star.set_distance(dist)
    mag_bb = {b: star.flux_to_magAB(band.band_by_name(b)) for b in bands}
    for b in bands:
        e += (mag[b] - mag_bb[b])**2
    return e


def compute_Tcolor_zeta(mags, tt, bands, freq, d):
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
        # tcolor, w = minimize(epsilon, x0=np.array([1.e4, 1]), args=(freq, mag, bands, radius, d))
        tcolor, w = fmin(epsilon, x0=np.array([1.e4, 1]), args=(freq, mag, bands, radius, d), disp=False)
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


def compute_tcolor(name, path, bands, t_cut=1.):
    model = Stella(name, path=path)

    if not model.is_ph_data:
        print "No ph-data for: " + str(model)
        return None

    if not model.is_tt_data:
        print "No tt-data for: " + str(model)
        return None

    distance = rf.pc_to_cm(10.)  # pc for Absolute magnitude
    # serial_spec = model.read_series_spectrum(t_diff=1.)
    # curves = serial_spec.flux_to_curves(bands, d=distance)
    serial_spec = model.read_series_spectrum(t_diff=1.05, t_beg=t_cut)
    # curves = serial_spec.
    mags = serial_spec.mags_bands(bands, z=0., d=distance)

    # read R_ph
    tt = model.read_tt_data()
    tt = tt[tt['time'] > t_cut]  # time cut  days

    # compute Tnu, W
    Tnu, Teff, W = compute_Tnu_w(serial_spec, tt=tt)

    # fit mags by B(T_col) and get \zeta\theta & T_col
    Tcolors, zetaR, times = compute_Tcolor_zeta(mags, tt=tt, bands=bands, freq=serial_spec.Freq, d=distance)

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


def usage():
    bands = band.band_get_names().keys()
    print "Usage:"
    print "  zeta.py [params]"
    print "  -b <set_bands>: delimiter '_'. Default: B-V-I_B-V_V-I.\n" \
          "     Available: " + '-'.join(sorted(bands))
    print "  -i <model name>.  Example: cat_R450_M15_Ni007_E7"
    print "  -p <model path(directory)>, default: ./"
    print "  -e <model extension> is used to define model name, default: tt "
    print "  -s  silence mode: no info, no plot"
    print "  -f  force mode: rewrite zeta-files even if it exists"
    # print "  -a  plot the Eastman & Dessart fits"
    print "  -o  options: <fit:fitb:time:ubv:Tnu> - fit E&D: fit bakl:  show time points: plot UBV"
    print "  -w  write magnitudes to file, default 'False'"
    print "  -h  print usage"


def main(name='', path='./', is_force=False, is_save=False, is_plot_Tnu=False, is_plot_time_points=False):
    is_silence = False
    is_fit = False
    is_plot_ubv = False
    is_fit_bakl = False
    model_ext = '.tt'
    ubv_args = ''

    try:
        opts, args = getopt.getopt(sys.argv[1:], "fhswtp:e:i:b:o:")
    except getopt.GetoptError as err:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if not name:
        if len(opts) == 0:
            usage()
            sys.exit(2)
        for opt, arg in opts:
            if opt == '-i':
                path = ROOT_DIRECTORY
                name = str(arg)
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
                        print 'No such band: ' + b
                        sys.exit(2)
            continue
        if opt == '-s':
            is_silence = True
            ubv_args += opt + ' '
            continue
        if opt == '-o':
            ops = str(arg).split(':')
            is_plot_ubv = "ubv" in ops
            is_plot_Tnu = "Tnu" in ops
            is_plot_time_points = "time" in ops
            is_fit = "fit" in ops
            is_fit_bakl = "fitb" in ops
            ubv_args += " %s %s ".format(opt, arg)
            continue
        if opt == '-w':
            is_save = True
            ubv_args += opt + ' '
            continue
        if opt == '-f':
            is_force = True
            is_save = True
            continue
        if opt == '-p':
            path = os.path.expanduser(str(arg))
            if not (os.path.isdir(path) and os.path.exists(path)):
                print "No such directory: " + path
                sys.exit(2)
            continue
        elif opt == '-h':
            usage()
            sys.exit(2)

    names = []
    if name != '':
        names.append(name)
    else:  # run for all files in the path
        files = [f for f in os.listdir(path) if isfile(join(path, f)) and f.endswith(model_ext)]
        for f in files:
            names.append(os.path.splitext(f)[0])

    im = 0
    if len(names) > 0:
        dic_results = {}  # dict((k, None) for k in names)
        for name in names:
            im += 1
            dic = {}  # dict((k, None) for k in set_bands)
            print "\nRun: %s [%d/%d]" % (name, im, len(names))
            for bset in set_bands:
                fname = cache_name(name, path, bset)
                if not is_force and os.path.exists(fname):
                    dic[bset] = cache_load(fname)
                else:
                    res = compute_tcolor(name, path, bset.split('-'), t_cut=0.9)
                    if res is not None:
                        dic[bset] = res
                        if is_save:
                            print "Save Tcolor & Zeta for %s in %s" % (bset, fname)
                            cache_save(dic[bset], fname=fname)

            dic_results[name] = dic
            print "Finish: %s" % name
        if is_plot_ubv:
            os.system("./ubv.py -i %s -p %s %s & " % (name, path, ubv_args))
        if not is_silence and len(dic_results) > 0:
            plot_zeta(dic_results, set_bands, t_cut=1.9, is_fit=is_fit, is_fit_bakl=is_fit_bakl,
                      is_plot_Tnu=is_plot_Tnu, is_time_points=is_plot_time_points)
            # plot_zeta_oneframe(dic_results, set_bands, t_cut=1.9, is_fit=is_fit,
            #                    is_plot_Tnu=is_plot_Tnu, is_time_points=is_plot_time_points)

    else:
        print "There are no models in the directory: %s with extension: %s " % (path, model_ext)


if __name__ == '__main__':
    main()

     # set_bands = ['B-V', 'B-V-I', 'V-I']
     # set_bands = ['B-V', 'B-V-I', 'V-I', 'J-H-K']
     # plot_fits(set_bands, is_grid=False)

    # main(name="cat_R1000_M15_Ni007_E15", path="/home/bakl/Sn/Release/seb_git/res/tt",
    #      is_force=False, is_save=True, is_plot_time_points=True)
