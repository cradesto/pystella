#!/usr/bin/python
# -*- coding: utf-8 -*-
import getopt
import numpy as np
import os
import sys
from os.path import dirname
from scipy import interpolate
from scipy.optimize import fmin

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import gridspec
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.mlab import griddata
from matplotlib.ticker import LinearLocator
from mpl_toolkits.mplot3d import Axes3D

import pystella.rf.rad_func as rf
from pystella.model.stella import Stella
from pystella.rf import band, spectrum
from pystella.rf.star import Star
from pystella.util.phys_var import phys

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))

colors_band = dict(U="blue", B="cyan", V="black", R="red", I="magenta",
                   J="green", H="cyan", K="black",
                   UVM2="skyblue", UVW1="orange", UVW2="blue",
                   g="g", r="red", i="magenta", u="blue", z="chocolate",
                   y='olive', w='tomato')

_colors = ["blue", "cyan", "brown", 'darkseagreen', 'tomato', 'olive', 'orange',
           'skyblue', 'darkviolet']
# colors = {"B-V": "blue", 'B-V-I': "cyan", 'V-I': "brown"}
lntypes = {"B-V": "-", 'B-V-I': "-.", 'V-I': "--"}
markers = {u'D': u'diamond', 6: u'caretup', u's': u'square', u'x': u'x',
           5: u'caretright', u'^': u'triangle_up', u'd': u'thin_diamond', u'h': u'hexagon1',
           u'+': u'plus', u'*': u'star', u'o': u'circle', u'p': u'pentagon', u'3': u'tri_left',
           u'H': u'hexagon2', u'v': u'triangle_down', u'8': u'octagon', u'<': u'triangle_left'}
markers = markers.keys()


def plot_spec(dic_stars, times, set_bands, is_planck=False, is_filter=False):
    xlim = [1000., 20000]
    # xlim = [3000., 6000]
    ylim = [0, 3.]

    # setup figure
    plt.matplotlib.rcParams.update({'font.size': 14})
    fig = plt.figure(num=len(set_bands), figsize=(9, 9), dpi=100, facecolor='w', edgecolor='k')
    gs1 = gridspec.GridSpec(len(set_bands), 1)
    gs1.update(wspace=0.3, hspace=0.3, left=0.1, right=0.9)

    ax_cache = {}
    ax2_cache = {}
    ib = 0
    for bset in set_bands:
        ib += 1
        irow = ib - 1
        if bset in ax_cache:
            ax = ax_cache[bset]
        else:
            ax = fig.add_subplot(gs1[irow, 0])
            ax_cache[bset] = ax

        for it in range(len(times)):
            star = dic_stars[it]
            x = star.Wl * phys.cm_to_angs
            y = star.FluxAB
            # y = star.FluxWlObs
            bcolor = _colors[it % (len(_colors) - 1)]
            ax.plot(x, y, marker=markers[it % (len(markers) - 1)],
                    label='%.0f: %0.2e, %0.2f' % (times[it], star.get_Tcol(bset), star.get_zeta(bset)),
                    markersize=5, color=bcolor, ls="", linewidth=1.5)
            if is_planck:
                star_bb = planck_fit(star, bset)
                xx = star_bb.Wl * phys.cm_to_angs
                yy = star_bb.FluxAB
                # yy = star_bb.FluxWlObs
                if yy is not None:
                    ax.plot(xx, yy, color=bcolor, ls="-", linewidth=2.5)
                    # ax.plot(xx, yy, color=bcolor, ls="--", linewidth=2.5, label='Planck')
        if is_filter:
            if bset in ax2_cache:
                ax2 = ax2_cache[bset]
            else:
                ax2 = fig.add_subplot(gs1[irow, 0], sharex=ax, frameon=False)
                # ax2.invert_yaxis()
                ax2_cache[bset] = ax2
            ax2.yaxis.tick_right()
            ax2.yaxis.set_label_position("right")
            ax2.set_ylabel(r'AB Magnitudes')
            # ax2.set_ylabel("Filter Response")
            bands = bset.split('-')
            for n in bands:
                b = band.band_by_name(n)
                xx = b.wl * phys.cm_to_angs
                yy = b.resp
                ax2.plot(xx, yy, color=colors_band[n], ls="--", linewidth=1.5, label=n)
            ax2.set_xscale('log')
            ax2.set_xlim(xlim)
            ax2.legend(prop={'size': 9}, loc=2, borderaxespad=0.)

    ib = 0
    for bset in set_bands:
        ib += 1
        ax = ax_cache[bset]
        ax.set_xlim(xlim)
        # ax.set_ylim(ylim)
        ax.invert_yaxis()
        ax.set_ylabel(r'Absolute AB Magnitude')
        # ax.set_ylabel(r'$F_\lambda, \, [erg\, s^{-1} cm^2]$')
        ax.set_xscale('log')
        # ax.set_yscale('log')
        if ib == len(set_bands):
            ax.set_xlabel(r'$\lambda, \, [\AA]$')
        ax.set_title(bset)
        # if is_filter:
        #     ax2 = ax2_cache[bset]
        #     ax2.legend(prop={'size': 9}, loc=4)
        # else:
        ax.legend(prop={'size': 9}, loc=4, borderaxespad=0.)
    # plt.title('; '.join(set_bands) + ' filter response')
    plt.grid()
    plt.show()


def color_map_temp():
    cdict = {'blue': ((0.0, 0.0, 0.0),
                      (0.25, 0.0, 0.0),
                      (0.5, 0.8, 1.0),
                      (0.75, 1.0, 1.0),
                      (1.0, 0.4, 1.0)),

             'green': ((0.0, 0.0, 0.0),
                       (0.25, 0.0, 0.0),
                       (0.5, 0.9, 0.9),
                       (0.75, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),

             'red': ((0.0, 0.0, 0.4),
                     (0.25, 1.0, 1.0),
                     (0.5, 1.0, 0.8),
                     (0.75, 0.0, 0.0),
                     (1.0, 0.0, 0.0))
             }
    cdict = {'blue': ((0.0, 0.0, 0.0),
                      (0.25, 0.0, 0.0),
                      (0.25, 0.0, 0.0),
                      (0.5, 0.3, 0.3),
                      (0.75, 0.75, 1.0),
                      (1.0, 0.9, 1.0)),

             'green': ((0.0, 0.0, 0.0),
                       (0.25, 0., 0.0),
                       (0.5, 0.9, 0.9),
                       (0.75, 0., 0.0),
                       (1.0, 0.0, 0.0)),

             'red': ((0.0, 0.0, 0.4),
                     (0.15, 0.9, 0.9),
                     (0.25, 0.8, 0.8),
                     (0.5, 0.5, 0.5),
                     (0.75, 0.0, 0.0),
                     (1.0, 0.0, 0.0))
             }
    # cdict = {'red': ((0.0, 1.0, 1.0),
    #                  (0.1, 1.0, 1.0),  # red
    #                  (0.4, 1.0, 1.0),  # violet
    #                  (1.0, 0.0, 0.0)),  # blue
    #
    #          'green': ((0.0, 0.0, 0.0),
    #                    (1.0, 0.0, 0.0)),
    #
    #          'blue': ((0.0, 0.0, 0.0),
    #                   (0.1, 0.0, 0.0),  # red
    #                   (0.5, 1.0, 1.0),  # violet
    #                   (1.0, 1.0, 0.0))  # blue
    #          }

    return LinearSegmentedColormap('UWR', cdict, 256)
    #
    # cm.register_cmap(name='UWR', cmap=cmap1)
    # UWR = cm.get_cmap('UWR')
    # return UWR


def plot_spec_poly(series, moments=None, fcut=1.e-20):
    # moments = moments or (1., 3., 5, 10., 15., 25., 35., 50., 70., 100., 120., 150., 200.)
    # moments = moments or np.arange(0., 200., 3)
    moments = moments or np.exp(np.linspace(np.log(0.5), np.log(400.), 40))

    # init graph
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    pos = 0
    t_data = []
    for i in range(len(series.Time)):
        t = series.Time[i]
        if t > moments[pos]:
            t_data.append(t)
            pos += 1

    verts = []
    T_wiens = []
    x_lim = [float("inf"), float("-inf")]
    z_lim = [float("inf"), float("-inf")]
    for i in range(len(t_data)):
        spec = series.get_spec_by_time(t_data[i])
        spec.cut_flux(fcut * max(spec.Flux))  # cut flux
        ys = spec.Flux
        # ys = np.log10(ys)
        wl = spec.Wl * phys.cm_to_angs
        verts.append(list(zip(wl, ys)))

        x_lim[0] = min(x_lim[0], np.min(wl))
        x_lim[1] = max(x_lim[1], np.max(wl))

        z_lim[0] = min(z_lim[0], np.min(ys))
        z_lim[1] = max(z_lim[1], np.max(ys))

        T_wiens.append(spec.temp_wien)
        print "time: %f  T_wien=%f" % (t_data[i], spec.temp_wien)

    T_wiens_lg = np.log10(T_wiens)
    # cc = lambda arg: colorConverter.to_rgba(arg, alpha=0.6)
    # facecolors = [cm.jet(np.sqrt(x)/10.) for x in t_data]
    # facecolors = UWR(heatmap)

    T_wiens_c = T_wiens_lg - np.min(T_wiens_lg)

    Tmap = T_wiens_c/np.max(T_wiens_c)
    color_map = color_map_temp()
    # cm.register_cmap(name='UWR', cmap=color_map)
    # UWR = cm.get_cmap('UWR')
    # Create a variable for the colorbar
    m = cm.ScalarMappable(cmap=color_map)
    m.set_array(Tmap)
    # m.set_array([min(T_wiens_lg), max(T_wiens_lg)])
    # m.set_clim(vmin=min(T_wiens), vmax=max(T_wiens))
    # m.set_array(T_wiens)
    facecolors = m.to_rgba(Tmap*1.e-4)
    # facecolors = color_map(Tmap)
    # facecolors = color_map(T_wiens_c)
    # facecolors = color_map(5e3/T_wiens)

    poly = PolyCollection(verts, facecolors=facecolors, linewidths=1.5)  # np.ones(len(t_data)))
    poly.set_alpha(0.5)
    ax.add_collection3d(poly, zs=t_data, zdir='y')

    # Create a color bar with 11 ticks
    cbar = plt.colorbar(m, ticks=LinearLocator(10), shrink=0.85)
    cbar.ax.yaxis.set_label_text("Wien temperature [K]")
    ticks = np.round(np.exp(np.linspace(np.log(1e3), np.log(max(T_wiens)), 10)), -2)
    # ticks = np.round(np.exp(np.linspace(np.log(min(T_wiens)), np.log(max(T_wiens)), 5)), -2)
    cbar.ax.set_yticklabels(ticks)

    # cbar = plt.colorbar(m, ticks=LinearLocator(11), shrink=0.85)
    # Make the tick label go from 0 to 1 in steps of 0.1
    # cbar.ax.set_yticklabels(arange(0, 1.01, 0.1))

    # ax.set_xlabel('Wave')
    ax.set_xlim3d(x_lim)
    # ax.set_xlim3d(min(wl), max(wl))

    # ax.set_ylabel('Time')
    ax.set_ylim3d(min(t_data), max(t_data))

    # ax.set_zlabel('Flux')
    ax.set_zlim3d(z_lim)

    # ax.set_xscale('log')
    ax.set_zscale('log')

    ax.xaxis.set_label_text('Wave [A]')
    ax.yaxis.set_label_text('Time [days]')
    ax.zaxis.set_label_text('Flux')

    plt.show()


def plot_spec_t(series, wl_lim=None, moments=None):
    wl_lim = wl_lim or (1e1, 5e4)
    moments = moments or (1., 3., 5, 10., 15., 25., 35., 50., 70., 100., 120., 150., 200.)
    # moments = np.arange(0.5, 200., 3)
    pos = 0
    t_data = []
    spec_array = []
    for i in range(len(series.Time)):
        t = series.Time[i]
        if t > moments[pos]:
            t_data.append(t)
            spec_array.append(series.get_spec(i).Flux)
            pos += 1

    y_data = series.Wl * phys.cm_to_angs
    spec_array = np.array(spec_array)

    x_data, y_data = np.meshgrid(t_data, y_data)

    x = x_data.flatten()
    y = y_data.flatten()
    z = spec_array.flatten()

    ## filters
    # z
    zlim = np.max(z) * 1.e-20
    is_z = z > zlim
    # y
    is_y = (y > wl_lim[0]) & (y < wl_lim[1])

    is_good = is_y & is_z
    x = x[is_good]
    y = y[is_good]
    z = z[is_good]

    fig = plt.figure()

    # scatter3D
    if True:
        ax = Axes3D(fig)
        ax.scatter3D(x, y, z, c=z, cmap=cm.jet)
        plt.show()
        return None

    # plot_surface
    is_plot_surface = False
    if is_plot_surface:
        xi = np.linspace(min(x), max(x))
        yi = np.linspace(min(y), max(y))
        X, Y = np.meshgrid(xi, yi)
        # interpolation
        Z = griddata(x, y, z, xi, yi)
        ax = Axes3D(fig)
        ax.scatter3D(x, y, z, c=z, cmap=cm.jet)
        ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)
        plt.show()
        return None

        #
        # surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.coolwarm,
        #                        linewidth=0, antialiased=False)
        #
        # ax = fig.add_subplot(111, projection='3d')
        #
        # X, Y = np.meshgrid(x_data, y_data)
        # z_data = griddata(x, y, data_array, x_data, y_data)
        #
        # if is_spline:
        #     spline = sp.interpolate.Rbf(x, y, z, function='thin-plate')
        #     xi = np.linspace(min(x), max(x))
        #     yi = np.linspace(min(y), max(y))
        #     X, Y = np.meshgrid(xi, yi)
        #     # interpolation
        #     Z = spline(X, Y)

        # x_data, y_data = np.meshgrid(np.arange(data_array.shape[1]),
        #                              np.arange(data_array.shape[0]))
        #
        # Flatten out the arrays so that they may be passed to "ax.bar3d".
        # Basically, ax.bar3d expects three one-dimensional arrays:
        # x_data, y_data, z_data. The following call boils down to picking
        # one entry from each array and plotting a bar to from
        # (x_data[i], y_data[i], 0) to (x_data[i], y_data[i], z_data[i]).
        #
        # x_data = x_data.flatten()
        # y_data = y_data.flatten()
        # z_data = data_array.flatten()
        # ax.bar3d(x_data,
        #          y_data,
        #          np.zeros(len(z_data)),
        #          1, 1, z_data)
        # surf = ax.plot_surface(x_data, y_data, z_data, rstride=1, cstride=1, cmap=cm.coolwarm,
        #                        linewidth=0, antialiased=False)
        # plt.show()


def planck_fit(star, bset):
    sp = spectrum.SpectrumPlanck(star.Freq, star.get_Tcol(bset))
    star_bb = Star("bb", sp)
    sp.correct_zeta(star.get_zeta(bset))
    star_bb.set_radius_ph(star.radius_ph)
    star_bb.set_distance(star.distance)
    return star_bb


def epsilon(x, freq, mag, bands, radius, dist):
    temp_color, zeta = x
    sp = spectrum.SpectrumPlanck(freq, temp_color)
    sp.correct_zeta(zeta)

    star = Star("bb", sp)
    star.set_radius_ph(radius)
    star.set_distance(dist)
    mag_bb = {b: star.flux_to_magAB(band.band_by_name(b)) for b in bands}
    e = 0
    for b in bands:
        e += abs(mag[b] - mag_bb[b])
    return e


def compute_tcolor(star, bands):
    mags = {}
    for n in bands:
        b = band.band_by_name(n)
        mags[n] = star.flux_to_magAB(b)

    Tcol, zeta = fmin(epsilon, x0=np.array([1.e4, 1]),
                      args=(star.Freq, mags, bands, star.radius_ph, star.distance), disp=0)
    return Tcol, zeta


def usage():
    bands = band.band_get_names().keys()
    print "Usage: show F(t,nu) from ph-file"
    print "  plot_spec.py [params]"
    print "  -b <set_bands>: delimiter '_'. Default: B-V.\n" \
          "     Available: " + '-'.join(sorted(bands))
    print "  -i <model name>.  Example: cat_R450_M15_Ni007_E7"
    print "  -o  options: [fit] - plot spectral fit in bands"
    print "  -p <model path(directory)>, default: ./"
    print "  -s  silence mode: no info, no plot"
    print "  -f  force mode: rewrite tcolor-files even if it exists"
    print "  -w  write magnitudes to file, default 'False'"
    print "  -h  print usage"


def main():
    is_silence = False
    is_fit = False

    try:
        opts, args = getopt.getopt(sys.argv[1:], "fhsuwp:i:b:o:")
    except getopt.GetoptError as err:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    name = ''
    path = os.path.expanduser('./')
    # name = 'cat_R500_M15_Ni006_E12'

    if not name:
        if len(opts) == 0:
            usage()
            sys.exit(2)
        for opt, arg in opts:
            if opt == '-i':
                path = ROOT_DIRECTORY
                name = str(arg)
                break

    # set_bands = ['B-V']
    # set_bands = ['B-V', 'B-V-I']
    # set_bands = ['U-B', 'U-B-V', 'B-V']
    set_bands = ['B-V', 'B-V-I', 'V-I']
    # set_bands = ['B-V', 'B-V-I', 'V-I', 'J-H-K']
    times = [15., 30., 60., 120.]

    for opt, arg in opts:
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
            continue
        if opt == '-w':
            is_save = True
            continue
        if opt == '-f':
            is_force = True
            is_save = True
            continue
        if opt == '-o':
            ops = str(arg).split(':')
            is_fit = "fit" in ops
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

    if not name:
        print "No model. Use key -i."
        sys.exit(2)

    model = Stella(name, path=path)

    if not model.is_ph_data:
        print "No ph-data for: " + str(model)
        return None

    series_spec = model.read_series_spectrum(t_diff=1.05)

    if not is_fit:
        plot_spec_poly(series_spec)
        # plot_spec_t(series_spec)
        print "Plot spectral F(t,nu): " + str(model)
        sys.exit(2)

    if not model.is_tt_data:
        print "No tt-data for: " + str(model)
        return None

    tt = model.read_tt_data()
    tt = tt[tt['time'] > min(times) - 1.]  # time cut  days

    Rph_spline = interpolate.splrep(tt['time'], tt['Rph'], s=0)

    distance = rf.pc_to_cm(10.)  # pc for Absolute magnitude
    dic_results = {}  # dict((k, None) for k in names)
    it = 0
    for time in times:
        print "\nRun: %s t=%f [%d/%d]" % (name, time, it, len(times))
        spec = series_spec.get_spec_by_time(time)
        spec.cut_flux(max(spec.Flux) * .1e-5)  # cut flux

        star = Star("%s: %f" % (name, time), spec=spec, is_flux_eq_luminosity=True)
        star.set_distance(distance)
        radius = interpolate.splev(time, Rph_spline)
        star.set_radius_ph(radius)
        star.set_redshift(0.)

        for bset in set_bands:
            Tcol, zeta = compute_tcolor(star, bset.split('-'))
            # save results
            star.set_Tcol(Tcol, bset)
            star.set_zeta(zeta, bset)
            print "\nRun: %s Tcol=%f zeta=%f " % (bset, Tcol, zeta)

        dic_results[it] = star
        it += 1
    if not is_silence:
        plot_spec(dic_results, times, set_bands, is_planck=True, is_filter=True)


if __name__ == '__main__':
    main()
