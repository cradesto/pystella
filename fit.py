#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
This program finds the time offset that minimizes chi-square.


"""
import argparse
import logging
import os
import sys
from itertools import cycle
from collections import OrderedDict
from concurrent import futures

import numpy as np
import scipy as sci

import pystella.util.callback as cb
from pystella import velocity
from pystella.fit.fit_mcmc import FitLcMcmc
from pystella.fit.fit_mpfit import FitMPFit
from pystella.rf import band
from pystella.rf import light_curve_func as lcf
from pystella.rf.band import band_is_exist
from pystella.rf.lc import SetLightCurve
from pystella.rf.ts import SetTimeSeries
from pystella.util.arr_dict import first
from pystella.util.path_misc import get_model_names
from pystella.util.string_misc import str2interval
from pystella.velocity import VelocityCurve, SetVelocityCurve

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

markers = {u'x': u'x', u'd': u'thin_diamond',
           u'+': u'plus', u'*': u'star', u'o': u'circle', u'v': u'triangle_down', u'<': u'triangle_left'}
markers_style = list(markers.keys())


def rel_errors(mu, sig, func, num=10000):
    x_norm = []
    for x, s in zip(mu, sig):
        x_norm.append(np.random.normal(x, s, num))
    # x_norm = np.random.normal(mu,sig, num)
    f_dist = func(x_norm)
    return np.mean(f_dist), np.std(f_dist)


def m_mu(x):
    return 10 ** (-x * 0.4)


def get_parser():
    parser = argparse.ArgumentParser(description='Process light curve fitting.')
    print(" Observational data could be loaded with plugin, ex: -c lcobs:filedata:tshift:mshift")

    parser.add_argument('-b', '--band',
                        required=False,
                        dest="bnames",
                        help="-b <bands>: string, default: U-B-V-R-I, for example g-i-r-UVW2")
    parser.add_argument('-c', '--call',
                        required=True,
                        nargs='+',
                        type=str,
                        dest='call',
                        help='Call observational data')
    parser.add_argument('-g', '--engine',
                        required=False,
                        type=str,
                        default='mpfit',
                        dest='engine',
                        help='The fitter algorithm-engine: [{}]. Default: mpfit'.format(', '.join(engines())))
    parser.add_argument('-d',
                        required=False,
                        type=float,
                        default=10,
                        dest="distance",
                        help="Distance to the model [pc].  Default: 10 pc")
    parser.add_argument('-e',
                        required=False,
                        type=float,
                        default=0.,
                        dest="color_excess",
                        help="Color excess E(B-V) is used to compute extinction via  A_nu = R_nu*E(B-V).  Default: 0")
    parser.add_argument('-z',
                        required=False,
                        type=float,
                        default=0,
                        dest="redshift",
                        help="Redshift for the model .  Default: 0")
    parser.add_argument('-i', '--input',
                        required=False,
                        dest="input",
                        help="Model name, example: cat_R450_M15_Ni007")
    parser.add_argument('-p', '--path',
                        required=False,
                        type=str,
                        default='./',
                        dest="path",
                        help="Model directory")
    parser.add_argument('-t', '--time',
                        required=False,
                        type=str,
                        default=None,
                        dest="times",
                        help="The range of fitting in model LC. Default: None (all points). Format: {0}".format('2:50'))
    parser.add_argument('-dt', '--dtshift',
                        required=False,
                        type=str,
                        default=None,
                        dest="dtshift",
                        help="The range of tshift in model LC. Default: None (any time). Format: {0}".format('2:50'))
    parser.add_argument('-nq', '--no-quiet',
                        action='store_const',
                        const=True,
                        dest="is_not_quiet",
                        help="Result with additional information")
    parser.add_argument('-n', '--node',
                        required=False,
                        type=int,
                        default=1,
                        dest="nodes",
                        help="-n <nodes>: number of processes ")
    return parser


def engines(nm=None):
    switcher = {
        'mpfit': FitMPFit(),
        'mcmc': FitLcMcmc(),
    }
    if nm is not None:
        return switcher.get(nm)
    return list(switcher.keys())


def plot_curves(curves_o, res_models, res_sorted, **kwargs):
    from pystella.rf import light_curve_plot as lcp
    from matplotlib import pyplot as plt

    font_size = kwargs.get('font_size', 8)

    xlim = None
    ylim = None
    num = len(res_sorted)
    nrow = int(num / 2.1) + 1
    ncol = 2 if num > 1 else 1
    fig = plt.figure(figsize=(12, nrow * 4))
    plt.matplotlib.rcParams.update({'font.size': font_size})

    tshift0 = first(curves_o).tshift
    i = 0
    for k, v in res_sorted.items():
        i += 1
        ax = fig.add_subplot(nrow, ncol, i)

        tshift_best = v.tshift
        curves = res_models[k]
        lcp.curves_plot(curves, ax=ax, figsize=(12, 8), linewidth=1, is_legend=False)
        xlim = ax.get_xlim()
        lt = {lc.Band.Name: 'o' for lc in curves_o}
        curves_o.set_tshift(tshift0+tshift_best)
        lcp.curves_plot(curves_o, ax, xlim=xlim, lt=lt, markersize=2, is_legend=False)

        ax.text(0.99, 0.94, k, horizontalalignment='right', transform=ax.transAxes)
        ax.text(0.98, 0.85, "dt={:.2f}".format(tshift_best), horizontalalignment='right', transform=ax.transAxes)
        ax.text(0.01, 0.05, "$\chi^2: {:.2f}$".format(v.measure), horizontalalignment='left', transform=ax.transAxes,
                bbox=dict(facecolor='green', alpha=0.3))
        # ax.text(0.9, 0.9, "{:.2f}".format(v.measure), horizontalalignment='right', transform=ax.transAxes, bbox=dict(facecolor='green', alpha=0.3))

        # fix axes
        if i % ncol == 0:
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position("right")
            # ax.set_ylabel('')
            # ax.set_yticklabels([])
            # ax2.set_ylabel('Magnitude')
        else:
            ax.yaxis.tick_left()
            ax.yaxis.set_label_position("left")
            # ax.set_ylabel('Magnitude')
        ax.yaxis.set_ticks_position('both')
        # legend
        ax.legend(curves.BandNames, loc='lower right', frameon=False, ncol=min(5, len(curves.BandNames)),
                  fontsize='small', borderpad=1)
        # lc_colors = band.bands_colors()

    fig.subplots_adjust(wspace=0, hspace=0)
    plt.subplots_adjust(left=0.07, right=0.96, top=0.97, bottom=0.06)
    plt.show()


def plot_curves_vel(curves_o, vels_o, res_models, res_sorted, vels_m, **kwargs):
    from pystella.rf import light_curve_plot as lcp
    from matplotlib import pyplot as plt

    font_size = kwargs.get('font_size', 10)
    linewidth = kwargs.get('linewidth', 2.0)
    markersize = kwargs.get('markersize', 5)

    xlim = None
    ylim = None
    num = len(res_sorted)
    nrow = int(num/2.1) + 1
    ncol = 2
    fig = plt.figure(figsize=(12, nrow * 4))
    plt.matplotlib.rcParams.update({'font.size': font_size})
    tshift0 = 0
    if curves_o is not None:
        tshift0 = first(curves_o).tshift
    elif vels_o is not None:
        tshift0 = first(vels_o).tshift

    i = 0
    for k, v in res_sorted.items():
        i += 1
        #  Plot UBV
        axUbv = fig.add_subplot(nrow, ncol, i)

        tshift_best = v.tshift
        axUbv.text(0.99, 0.94, k, horizontalalignment='right', transform=axUbv.transAxes)
        axUbv.text(0.98, 0.85, "dt={:.2f}".format(tshift_best), horizontalalignment='right', transform=axUbv.transAxes)
        axUbv.text(0.01, 0.05, "$\chi^2: {:.2f}$".format(v.measure), horizontalalignment='left', transform=axUbv.transAxes,
                   bbox=dict(facecolor='green', alpha=0.3))

        curves = res_models[k]
        if curves is not None:
            lcp.curves_plot(curves, ax=axUbv, figsize=(12, 8), linewidth=1, is_legend=False)
            lt = {lc.Band.Name: 'o' for lc in curves_o}
            curves_o.set_tshift(tshift0 + tshift_best)
            lcp.curves_plot(curves_o, axUbv, xlim=axUbv.get_xlim(), lt=lt, markersize=2, is_legend=False)
            # legend
            axUbv.legend(curves.BandNames, loc='lower right', frameon=False, ncol=min(5, len(curves.BandNames)),
                         fontsize='small', borderpad=1)

        if i % ncol == 0:
            axUbv.yaxis.tick_right()
            axUbv.set_ylabel('')
            axUbv.yaxis.set_ticks([])
        else:
            axUbv.yaxis.tick_left()
            axUbv.yaxis.set_label_position("left")
        axUbv.grid(linestyle=':')
        #  Plot Vel
        axVel = axUbv.twinx()
        axVel.set_ylim((0., 29))
        axVel.set_ylabel('Velocity [1000 km/s]')

        ts_vel = vels_m[k]
        x = ts_vel.Time
        y = ts_vel.V
        axVel.plot(x, y, label='Vel  %s' % k, color='blue', ls="-", linewidth=linewidth)
        # Obs. vel

        markers_cycler = cycle(markers_style)
        for vel_o in vels_o:
            vel_o.tshift = tshift0 + tshift_best
            axVel.plot(vel_o.Time, vel_o.V, label=vel_o.Name, color='blue', ls='', marker=next(markers_cycler),
                       markersize=markersize)

    fig.subplots_adjust(wspace=0, hspace=0)
    plt.subplots_adjust(left=0.07, right=0.96, top=0.97, bottom=0.06)
    plt.show()


def plot_squared_grid(res_sorted, path='./', **kwargs):
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    font_size = kwargs.get('font_size', 10)
    num = 4
    nrow = int(num / 2.1) + 1
    ncol = 2 if num > 1 else 1
    fig = plt.figure(figsize=(12, nrow * 6))
    plt.matplotlib.rcParams.update({'font.size': font_size})

    # R-M
    axRM = fig.add_subplot(nrow, ncol, 1)
    plot_squared(axRM, res_sorted, path, p=('R', 'M'), **kwargs)

    # R-E
    axRE = fig.add_subplot(nrow, ncol, 2)
    plot_squared(axRE, res_sorted, path, p=('R', 'E'), **kwargs)

    # M-E
    axME = fig.add_subplot(nrow, ncol, 3)
    plot_squared(axME, res_sorted, path, p=('M', 'E'), **kwargs)

    # R-M
    ax = fig.add_subplot(nrow, ncol, 4, projection='3d')
    plot_squared_3d(ax, res_sorted, path, p=('R', 'M', 'E'), is_surface=False, **kwargs)

    #
    # ax = fig.add_subplot(nrow, ncol, 4, projection='polar')
    # plot_squared_3d(ax, res_sorted, path, p=('R', 'M', 'E'), is_polar=True, **kwargs)
    #
    plt.show()


def plot_squared(ax, res_sorted, path='./', p=('R', 'M'), **kwargs):
    from matplotlib import pyplot as plt
    from pystella.model.stella import Stella

    is_rbf = kwargs.get('is_rbf', True)
    is_surface = kwargs.get('is_surface', True)
    is_scatter = kwargs.get('is_scatter', not is_surface and False)
    is_not_quiet = False
    # is_not_quiet = kwargs.get('is_not_quiet', True)

    # graph
    # ax.set_title('-'.join(p))
    ax.set_xlabel(p[0])
    ax.set_ylabel(p[1])

    # find parameters
    i = 0
    data = []
    chi = []
    models = []
    for name, res_chi in res_sorted.items():
        i += 1
        stella = Stella(name, path=path)
        if stella.is_tt_data:
            try:
                info = stella.get_tt().Info
                v = [getattr(info, pp) for pp in p]

                # print info
                if is_not_quiet:
                    if i == 1:
                        print("| %40s |  %7s |  %6s" % ('Model', p[0], p[1]))
                    print("| {:40s} | ".format(info.Name) + ' '.join("{0:6.2f}".format(vv) for vv in v))

                k = -1
                for vo in data:
                    k += 1
                    if np.array_equal(v, vo):
                        if is_not_quiet:
                            print("|   | " + ' '.join("{0:6.2f}".format(vv) for vv in v)
                                  + " |  {:40s} | chi_saved={:6.2f}  chi_new={:6.2f}".
                                  format('This is not a unique point', chi[k], res_chi.measure))
                        if res_chi.measure < chi[k]:
                            print("| %50s | k = %5d  Chi [%7.2f] < [%6.2f]" %
                                  (info.Name + ' '.join("{0:6.2f}".format(vv) for vv in v),
                                   k, res_chi.measure, chi[k]))
                            chi[k] = res_chi.measure
                        break
                else:
                    models.append(name)
                    data.append(v)
                    chi.append(res_chi.measure)
            except KeyError as ex:
                print("Error for model {}. Message: {} ".format(name, ex))

    if len(models) == 0:
        print('There are no tt-data for any models.')
        return

    # plot
    # x, y = map(np.array, zip(data))
    x = [v[0] for v in data]
    y = [v[1] for v in data]
    x = np.array(x)
    y = np.array(y)
    chi = np.array(chi)

    if is_surface:
        # Set up a regular grid of interpolation points
        xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
        xi, yi = np.meshgrid(xi, yi)

        # Interpolate
        if is_rbf:
            rbf = sci.interpolate.Rbf(x, y, chi, function='linear')
            zi = rbf(xi, yi)
        else:
            zi = sci.interpolate.griddata((x, y), chi, (xi, yi), method="linear")

        # im = ax.imshow(zi, cmap=plt.cm.RdBu, vmin=chi.min(), vmax=chi.max(), origin='lower',
        im = ax.imshow(zi, cmap=plt.cm.bone, vmin=chi.min(), vmax=chi.max(), origin='lower',
                       extent=[x.min(), x.max(), y.min(), y.max()], interpolation='none',
                       aspect='auto', alpha=0.5)
        # try:
        #     from skimage import measure
        #     # Find contours at a constant value
        #     levels = np.linspace(np.min(chi), np.max(chi), 10)
        #     levels = levels[1:len(levels)-2]
        #     for level in levels:
        #         contours = measure.find_contours(zi, level)
        #         for n, contour in enumerate(contours):
        #             lx, ly = contour[:, 1], contour[:, 0]
        #             l, = ax.plot(lx, ly, linewidth=2, label="%.1f"%level)
        #             pos = [(lx[-2] + lx[-1]) / 2., (ly[-2] + ly[-1]) / 2.]
        #             # xscreen = ax.transData.transform(zip(lx[-2::], ly[-2::]))
        #             # rot = np.rad2deg(np.arctan2(*np.abs(np.gradient(xscreen)[0][0][::-1])))
        #             rot = 0
        #             ltex = plt.text(pos[0], pos[1], "%.1f"%level, size=9, rotation=rot,
        #                             color=l.get_color(), ha="center", va="center",
        #                             bbox=dict(ec='1', fc='1'))
        #
        #     # ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, fancybox=True)
        # except ImportError:
        cset = ax.contour(xi, yi, zi, linewidths=1., cmap=plt.cm.bone)
        # cset = ax.contour(xi, yi, zi, linewidths=2, cmap=plt.cm.Set2)
        ax.clabel(cset, inline=True, fmt='%1.1f', fontsize=9)
        cbar = plt.colorbar(im)
        cbar.ax.set_ylabel(r'$\chi^2$')
        # plt.setp(cb.ax.get_yticklabels(), visible=False)

        ax.scatter(x, y, c=chi / np.max(chi), cmap=plt.cm.bone, picker=True)

        # ax.set_picker(True)

        def on_pick(event):
            try:
                ind = event.ind[0]
                print('{} {}: R={} M={} chi^2={:.2f}'.format(ind, models[ind], x[ind], y[ind], chi[ind]))
            except AttributeError:
                pass

        ax.figure.canvas.mpl_connect('pick_event', on_pick)
        # def onclick(event):
        #     print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
        #           (event.button, event.x, event.y, event.xdata, event.ydata))

        # ax.figure.canvas.mpl_connect('button_press_event', onclick)
        # ax.scatter(x, y, c=chi, cmap=plt.cm.RdBu)

    elif is_scatter:
        # Sort the points by density, so that the densest points are plotted last
        idx = chi.argsort()[::-1]
        x, y, chi = x[idx], y[idx], chi[idx]

        area = np.pi * (100 * np.log10(chi))**2  # 0 to 15 point radii
        # plt.scatter(x, y, s=area, c=colors, alpha=0.5)
        cax = plt.scatter(x, y, s=area, c=chi, cmap='gray', edgecolor='', alpha=0.5)
        plt.colorbar(cax)
    else:
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt

        # # Make data.
        xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
        xi, yi = np.meshgrid(xi, yi)
        z = chi
        #
        # # Interpolate
        if is_rbf:
            rbf = sci.interpolate.Rbf(x, y, chi, function='linear')
            zi = rbf(xi, yi)
        else:
            zi = sci.interpolate.griddata((x, y), chi, (xi, yi), method="linear")

        surf = ax.plot_trisurf(x, y, z, linewidth=0.2, antialiased=True, cmap='gray')
        # ax.plot_trisurf(xi, yi, zi, linewidth=0.2, antialiased=True)
        # # Plot the surface.
        # surf = ax.plot_surface(xi, yi, zi, cmap=plt.cm.coolwarm,
        #                        linewidth=0, antialiased=False)
        #
        # # Customize the z axis.
        # ax.set_zlim(-1.01, 1.01)
        # ax.zaxis.set_major_locator(LinearLocator(10))
        # ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        #
        # # Add a color bar which maps values to colors.
        plt.colorbar(surf, shrink=0.5, aspect=5)
        plt.subplots_adjust(left=0.07, right=0.06, top=0.97, bottom=0.06)
        # plt.show()


def plot_squared_3d(ax, res_sorted, path='./', p=('R', 'M', 'E'), is_rbf=True, **kwargs):
    from matplotlib import pyplot as plt
    from pystella.model.stella import Stella

    is_not_quiet = False
    is_polar = kwargs.get('is_polar', False)
    is_show = False

    # find parameters
    i = 0
    # info_models = {}
    data = []
    # data = np.empty(len(p), len(res_sorted))
    chi = []
    for name, res_chi in res_sorted.items():
        i += 1
        stella = Stella(name, path=path)
        if stella.is_tt_data:
            try:
                info = stella.get_tt().Info
                v = [getattr(info, pp) for pp in p]

                # print info
                if is_not_quiet:
                    if i == 1:
                        print("| %40s |  %7s |  %6s" % ('Model', p[0], p[1]))
                    # print("| %40s |  %7.2f |  %6.2f" % (info.Name) + v)
                    print("| {:40s} | ".format(info.Name) + ' '.join("{0:6.2f}".format(vv) for vv in v))

                k = -1
                for vo in data:
                    k += 1
                    if np.array_equal(v, vo):
                        print("|   | " + ' '.join("{0:6.2f}".format(vv) for vv in v)
                              + " |  {:40s} | chi_saved={:6.2f}  chi_new={:6.2f}".
                              format('This is not a unique point', chi[k], res_chi.measure))
                        if res_chi.measure < chi[k]:
                            print("| %40s | k = %5d  Chi [%7.2f] => [%6.2f]" %
                                  (info.Name, k, res_chi.measure, chi[k]))
                            chi[k] = res_chi.measure
                        break
                # if -1 < k < len(x)-1:  # This is not a unique point
                #     # if v1 in x and v2 in y:
                #         todo Условие на то что chi меньше старого
                #     print("| %40s |  %7.2f |  %6.2f |  %s" %
                #           ('   ', v1, v2, 'This is not a unique point'))
                else:
                    data.append(v)
                    # y.append(v2)
                    chi.append(res_chi.measure)
            except KeyError as ex:
                print("Error for model {}. Message: {} ".format(name, ex))

    if len(data) == 0:
        print('There are no tt-data for any models.')
        return

    # plot
    x = [v[0] for v in data]
    y = [v[1] for v in data]
    z = [v[2] for v in data]

    x = np.array(x[::-1])
    y = np.array(y[::-1])
    z = np.array(z[::-1])
    chi = np.array(chi[::-1])

    if ax is None:
        fig = plt.figure(figsize=(12, 8))
        if is_polar:
            ax = fig.add_subplot(1, 1, 1, projection='polar')
        else:
            ax = fig.add_subplot(1, 1, 1)
        is_show = True

    if is_polar:
        C = 1.9
        theta = (x - np.min(x)) / (np.max(x) - np.min(x)) * C * np.pi  # R
        width = y / np.max(y) * np.pi / 8  # M
        # radii = np.log10(chi+1) * 10 # chi
        radii = 10 + (chi - np.min(chi)) / (np.max(chi) - np.min(chi)) * 100  # chi

        bars = ax.bar(theta, radii, width=width, bottom=0.0)

        # Use custom colors and opacity
        labels = z
        for z, bar, l in zip(z / np.max(z), bars, labels):
            # bar.set_facecolor(plt.cm.jet(r))  # E
            bar.set_facecolor(plt.cm.plasma(z))
            bar.set_alpha(0.5)
            bar.set_label(l)

        xlabel_max = np.min(x) + 7. / 4. / C * (np.max(x) - np.min(x))
        xlabel = ["{:6.0f} R".format(xx) for xx in np.linspace(np.min(x), xlabel_max, 8)]
        ax.set_xticklabels(xlabel)
        ax.legend(bbox_to_anchor=(1.3, 1.05))
    else:
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        # graph
        # ax.set_title('-'.join(p))
        ax.set_xlabel(p[0])
        ax.set_ylabel(p[1])
        ax.set_zlabel(p[2])

        N = chi / np.max(chi)
        surf = ax.scatter(x, y, z, c=N, cmap="gray")
        from matplotlib.ticker import LinearLocator
        from matplotlib.ticker import FormatStrFormatter

        ax.yaxis.set_major_locator(LinearLocator(10))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

        plt.colorbar(surf, shrink=0.5, aspect=5)
        ax.grid(True)
    if is_show:
        plt.show()


def fit_mfl(args, curves_o, vels_o, bnames, fitter, name, path, t_diff, times):
    tss_m = SetTimeSeries("Models")
    tss_o = SetTimeSeries("Obs")
    # tss_o = {}

    curves_m = None
    # tshifts = {}
    # light curves
    if curves_o is not None:
        curves_m = lcf.curves_compute(name, path, bnames, z=args.redshift, distance=args.distance,
                                      t_beg=times[0], t_end=times[1], t_diff=t_diff)
        if args.color_excess:
            curves_m = lcf.curves_reddening(curves_m, ebv=args.color_excess, z=args.redshift, is_info=args.is_not_quiet)

        for lc in curves_m:
            tss_m.add(lc)
            # tss_m[lc.Band.Name] = lc

        for lc in curves_o:
            l, tshift, mshift = lc.shifted()
            tss_o.add(l)
            # tss_o[lc.Band.Name], tshift, mshift = lc.shifted()
            # tshifts[lc.Band.Name] = tshift

    vel_m = None
    if vels_o is not None:
        # compute model velocities
        try:
            tbl = velocity.compute_vel_swd(name, path)
            # tbl = velocity.compute_vel_res_tt(name, path)
            Vnorm = 1e8
            vel_m = VelocityCurve('Vel', tbl['time'], tbl['vel']/Vnorm)
        except ValueError as ext:
            print(ext)

    # velocity
    if vel_m is not None:
        if curves_m is not None:
            # To increase the weight of Velocities with fitting
            for i in range(curves_m.Length):
                for vel_o in vels_o:
                    key = 'Vel{:d}{}'.format(i,vel_o.Name)
                    tss_m.add(vel_m.copy(name=key))
                    tss_o.add(vel_o.copy(name=key))
                    # tss_o[key] = vel_o
                    # tss_m[key] = vel_m
        else:
            tss_m.add(vel_m.copy(name='Vel'))
            tss_o.add(vels_o.copy(name='Vel'))
            # tss_m['Vel'] = vel_m
            # tss_o['Vel'] = vels_o

    # fit
    res = fitter.fit_tss(tss_o, tss_m)

    return curves_m, vel_m, res


def main():
    model_ext = '.ph'
    t_diff = 1.01
    name = None
    # z = 0.
    # distance = 10  # pc
    # bnames = ['U', 'B', 'V', 'R', "I"]
    Nbest = 33
    NbestPlot = 6
    band.Band.load_settings()

    parser = get_parser()
    args, unknownargs = parser.parse_known_args()

    # Set model names
    names = []

    if len(unknownargs) > 0:
        path, name = os.path.split(unknownargs[0])
        path = os.path.expanduser(path)
        name = name.replace(model_ext, '')
    else:
        if args.path:
            path = os.path.expanduser(args.path)
        else:
            path = os.getcwd()
        if args.input is not None:
            name = os.path.splitext(os.path.basename(args.input))[0]  # remove extension

    # if len(unknownargs) == 0:
    #     parser.print_help()
    #     sys.exit(2)

    if name is None:
        names = get_model_names(path, model_ext)  # run for all files in the path
    else:
        names.append(name)

    # Set band names
    bnames = []
    if args.bnames:
        for bname in args.bnames.split('-'):
            if not band.band_is_exist(bname):
                print('No such band: ' + bname)
                parser.print_help()
                sys.exit(2)
            bnames.append(bname)

    # if args.distance:
    # distance = args.distance

    if args.call:
        if len(args.call) > 1:
            a = []
            for line in args.call:
                c = cb.lc_wrapper(line, method='load')
                a.append(c)
            callback = cb.CallBackArray(a)
        elif len(args.call) == 1:
            callback = cb.lc_wrapper(args.call[0], method='load')
        else:
            callback = None
    else:
        print('No obs data. Use key: -c: ')
        parser.print_help()
        sys.exit(2)

    # Get observations
    curves_o = None
    vels_o = None
    obs = callback.run({'is_debug': args.is_not_quiet})
    if isinstance(obs, list):
        for o in obs:
            if isinstance(o, SetLightCurve):
                curves_o = SetLightCurve.Merge(curves_o, o)
            if isinstance(o, SetVelocityCurve):
                vels_o = SetVelocityCurve.Merge(vels_o, o)
    else:
        if isinstance(obs, SetLightCurve):
            curves_o = obs
        if isinstance(obs, SetVelocityCurve):
            vels_o = obs

    is_vel = vels_o is not None and vels_o.Length > 0
    is_curves_o = curves_o is not None
    # curves_o = callback.load({'is_debug': not args.is_quiet})
    # todo Information about tshift

    if is_curves_o and len(bnames) == 0:
        bnames = [bn for bn in curves_o.BandNames if band_is_exist(bn)]

    # Time limits for models
    times = (0, None)

    if args.times:
        times = list(map(float, args.times.split(':')))

    print('Time limits for models: {}'.format(':'.join(map(str, times))))

    # The fit engine
    fitter = engines(args.engine)
    fitter.is_info = args.is_not_quiet  # fitter = FitMPFit(is_debug=args.is_not_quiet)
    fitter.is_debug = args.is_not_quiet

    # The filter results by tshift
    if args.dtshift:
        dtshift = str2interval(args.dtshift, llim=float("-inf"), rlim=float('inf'))
    else:
        dtshift = (float("-inf"), float("inf"))

    # tshift = 0.
    res_models = {}
    vels_m = {}
    res_chi = {}
    if len(names) == 1:
        name = names[0]
        if args.is_not_quiet:
            if times is not None:
                print("Fitting for model %s %s for %s moments" % (path, name, times))
            else:
                print("Fitting for model %s %s " % (path, name))
        # curves_m = lcf.curves_compute(name, path, bnames, z=args.redshift, distance=args.distance,
        #                               t_beg=times[0], t_end=times[1], t_diff=t_diff)
        # res = fitter.fit_curves(curves_o, curves_m)
        curves_m, vel_m, res = fit_mfl(args, curves_o, vels_o, bnames, fitter, name, path, t_diff, times)

        print("{}: time shift  = {:.2f}+/-{:.4f} Measure: {:.4f}".format(name, res.tshift, res.tsigma, res.measure))
        # best_tshift = res.tshift
        res_models[name] = curves_m
        vels_m[name] = vel_m
        res_chi[name] = res
        res_sorted = res_chi
    elif len(names) > 1:
        if args.nodes > 1:
            print("Run parallel fitting: nodes={}, models  {}".format(args.nodes, len(names)))

            with futures.ProcessPoolExecutor(max_workers=args.nodes) as executor:
                future_to_name = {
                    executor.submit(fit_mfl, args, curves_o, vels_o, bnames, fitter, n, path, t_diff, times):
                        n for n in names
                }
                i = 0
                for future in futures.as_completed(future_to_name):
                    i += 1
                    name = future_to_name[future]
                    try:
                        data = future.result()
                    except Exception as exc:
                        print('%r generated an exception: %s' % (name, exc))
                    else:
                        res_models[name] = data[0]
                        vels_m[name] = data[1]
                        v = data[2]
                        res_chi[name] = v
                        print("[{}/{}] {:30s} -> {}".format(i, len(names), name, v.comm))
        else:
            i = 0
            for name in names:
                i += 1
                txt = "Fitting [{}] for model {:30s}  [{}/{}]".format(fitter.Name, name, i, len(names))
                if args.is_not_quiet:
                    print(txt)
                else:
                    sys.stdout.write(u"\u001b[1000D" + txt)
                    sys.stdout.flush()
                curves_m, vel_m, res = fit_mfl(args, curves_o, vels_o, bnames, fitter, name, path, t_diff, times)
                res_models[name] = curves_m
                vels_m[name] = vel_m
                res_chi[name] = res

        # select with dtshift
        res_chi_sel = {}
        for k, v in res_chi.items():
            if dtshift[0] < v.tshift < dtshift[1]:
                res_chi_sel[k] = v
        res_chi = res_chi_sel
        # sort with measure
        res_sorted = OrderedDict(sorted(res_chi.items(), key=lambda kv: kv[1].measure))
    else:
        print("No any data about models. Path: {}".format(path))
        parser.print_help()
        sys.exit(2)

    # print results
    print("\n Results (tshift in range:{:.2f} -- {:.2f}".format(dtshift[0], dtshift[1]))
    print("{:40s} ||{:18s}|| {:10}".format('Model', 'dt+-t_err', 'Measure'))
    for k, v in res_sorted.items():
        print("{:40s} || {:7.2f}+/-{:7.4f} || {:.4f}".format(k, v.tshift, v.tsigma, v.measure))

    if len(res_sorted) >= Nbest:
        # plot chi squared
        plot_squared_grid(res_sorted, path, par=('M', 'E'), is_not_quiet=args.is_not_quiet)

        # plot only Nbest modeles
        while len(res_sorted) > Nbest:
            res_sorted.popitem()

        plot_squared_3d(None, res_sorted, path, p=('R', 'M', 'E'), is_polar=True)

    best_mdl, res = first(res_sorted.items())
    # res = first(res_sorted.values())[0]
    print("Best fit model:")
    print("{}: time shift  = {:.2f}+/-{:.4f} Measure: {:.4f}".format(best_mdl, res.tshift, res.tsigma, res.measure))

    # shift observational data
    # curves_o.set_tshift(best_tshift)

    # plot only NbestPlot modeles
    while len(res_sorted) > NbestPlot:
        res_sorted.popitem()

    if is_vel:
        # vel_o.tshift = best_tshift
        plot_curves_vel(curves_o, vels_o, res_models, res_sorted, vels_m)
    else:
        plot_curves(curves_o, res_models, res_sorted)


if __name__ == '__main__':
    main()
