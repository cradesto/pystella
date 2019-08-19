#!/usr/bin/env python3
# #!/usr/bin/python3

import getopt
import os
import sys
from os.path import dirname

import numpy as np
from scipy import interpolate
from scipy.optimize import fmin

import matplotlib.pyplot as plt
# from matplotlib import cm
from matplotlib import gridspec
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D

import pystella as ps

import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(os.path.abspath(__file__))

colors_band = {'U': "blue", 'B': "cyan", 'V': "black", 'R': "red", 'I': "magenta",
               'J': "green", 'H': "cyan", 'K': "black",
               'UVM2': "skyblue", 'UVW1': "orange", 'UVW2': "blue",
               'g': "g", 'r': "red", 'i': "magenta",
               'u': "blue", 'z': "chocolate", 'y': 'olive', 'w': 'tomato'}

_colors = ["blue", "cyan", "brown", 'darkseagreen', 'tomato', 'olive', 'orange',
           'skyblue', 'darkviolet']
# colors = {"B-V": "blue", 'B-V-I': "cyan", 'V-I': "brown"}
lntypes = {"B-V": "-", 'B-V-I': "-.", 'V-I': "--"}
markers = {u'D': u'diamond', 6: u'caretup', u's': u'square', u'x': u'x',
           5: u'caretright', u'^': u'triangle_up', u'd': u'thin_diamond', u'h': u'hexagon1',
           u'+': u'plus', u'*': u'star', u'o': u'circle', u'p': u'pentagon', u'3': u'tri_left',
           u'H': u'hexagon2', u'v': u'triangle_down', u'8': u'octagon', u'<': u'triangle_left'}
markers = list(markers.keys())


def plot_spec(dic_stars, times, set_bands, is_planck=False, is_filter=False):
    xlim = [1000., 20000]
    # xlim = [3000., 6000]

    # setup figure
    plt.matplotlib.rcParams.update({'font.size': 14})
    fig = plt.figure(num=len(set_bands), figsize=(8, len(set_bands) * 4), dpi=100, facecolor='w', edgecolor='k')
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
            x = star.Wl * ps.phys.cm_to_angs
            y = star.FluxAB
            # y = star.FluxWlObs
            bcolor = _colors[it % (len(_colors) - 1)]
            ax.plot(x, y, marker=markers[it % (len(markers) - 1)],
                    label=r'$%4.0f^d: T=%7.1e\,K, \zeta=%0.2f$' % (times[it], star.get_Tcol(bset), star.get_zeta(bset)),
                    markersize=5, color=bcolor, ls="", linewidth=1.5)
            if is_planck:
                star_bb = planck_fit(star, bset)
                xx = star_bb.Wl * ps.phys.cm_to_angs
                yy = star_bb.FluxAB
                # yy = star_bb.FluxWlObs
                if yy is not None:
                    ax.plot(xx, yy, color=bcolor, ls=":", linewidth=2.5)
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
            ax2.set_ylabel(r'Filter transmission')
            # ax2.set_ylabel("Filter Response")
            bands = bset.split('-')
            for n in bands:
                b = ps.band.band_by_name(n)
                xx = b.wl * ps.phys.cm_to_angs
                yy = b.resp_wl
                # ax2.plot(xx, yy, color=colors_band[n], ls="--", linewidth=1.5, label=n)
                ax2.fill_between(xx, 0, yy, color=colors_band[n], alpha=0.2)
            ax2.set_xscale('log')
            ax2.set_xlim(xlim)
            ax2.legend(prop={'size': 9}, loc=2, borderaxespad=0., fontsize='large')

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
        ax.legend(prop={'size': 11}, loc=4, borderaxespad=0.)
    # plt.title('; '.join(set_bands) + ' filter response')
    # plt.grid()
    return fig


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


def plot_spec_poly(series, moments=None, fcut=1.e-20, is_info=False):
    # moments = moments or (1., 3., 5, 10., 15., 25., 35., 50., 70., 100., 120., 150., 200.)
    # moments = moments or np.arange(0., 200., 3)
    moments = moments or np.exp(np.linspace(np.log(0.5), np.log(400.), 40))

    # init graph
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    pos = 0
    t_data = []
    for i, t in enumerate(series.Time):
        if t > moments[pos]:
            t_data.append(t)
            pos += 1

    verts = []
    T_cols = []
    x_lim = [float("inf"), float("-inf")]
    z_lim = [float("inf"), float("-inf")]
    for i, t in enumerate(t_data):
        spec = series.get_spec_by_time(t_data[i])
        spec.cut_flux(fcut * max(spec.Flux))  # cut flux
        ys = spec.Flux
        # ys = np.log10(ys)
        wl = spec.Wl * ps.phys.cm_to_angs
        verts.append(list(zip(wl, ys)))

        x_lim[0] = min(x_lim[0], np.min(wl))
        x_lim[1] = max(x_lim[1], np.max(wl))

        z_lim[0] = min(z_lim[0], np.min(ys))
        z_lim[1] = max(z_lim[1], np.max(ys))

        T_cols.append(spec.T_color)
        if is_info:
            print("time: %f  T_color=%f" % (t, spec.T_color))
            # print "time: %f  T_wien=%f wl_max=%f" % (t_data[i], spec.temp_wien, spec.wl_flux_max)

    Tmap = np.log(T_cols / np.min(T_cols))
    color_map = color_map_temp()
    m = plt.cm.ScalarMappable(cmap=color_map)
    m.set_array(Tmap)
    facecolors = m.to_rgba(Tmap * 1.e-4)
    poly = PolyCollection(verts, facecolors=facecolors, linewidths=1.5)  # np.ones(len(t_data)))
    poly.set_alpha(0.5)
    ax.add_collection3d(poly, zs=t_data, zdir='y')

    # Create a color bar with 11 ticks
    cbar = plt.colorbar(m, shrink=0.85)
    ticks = np.linspace(min(Tmap), max(Tmap), 10)
    ticks_lbl = np.round(np.exp(ticks) * np.min(T_cols), -2)
    cbar.ax.set_yticklabels(ticks_lbl)
    cbar.ax.yaxis.set_label_text("Color temperature [K]")

    ax.set_xlim3d(x_lim)
    ax.set_ylim3d(min(t_data), max(t_data))
    ax.set_zlim3d(z_lim)

    # ax.set_xscale('log')
    ax.set_zscale('log')

    ax.xaxis.set_label_text('Wavelength [A]')
    ax.yaxis.set_label_text('Time [days]')
    ax.zaxis.set_label_text('Flux [a.u.]')

    return fig


def plot_fit_bands(model, series, set_bands, times):
    name = model.Name
    tt = model.get_tt().load()
    tt = tt[tt['time'] > min(times) - 1.]  # time cut  days
    Rph_spline = interpolate.splrep(tt['time'], tt['Rph'], s=0)
    distance = ps.phys.pc2cm(10.)  # pc for Absolute magnitude
    dic_results = {}  # dict((k, None) for k in names)
    # it = 0
    for it, time in enumerate(times):
        print("\nRun: %s t=%f [%d/%d]" % (name, time, it, len(times)))
        spec = series.get_spec_by_time(time)
        spec.cut_flux(max(spec.Flux) * .1e-6)  # cut flux

        star = ps.Star("%s: %f" % (name, time), spec=spec, is_flux_eq_luminosity=True)
        star.set_distance(distance)
        radius = interpolate.splev(time, Rph_spline)
        star.set_radius_ph(radius)
        star.set_redshift(0.)

        for bset in set_bands:
            Tcol, zeta = compute_tcolor(star, bset.split('-'))
            # save results
            star.set_Tcol(Tcol, bset)
            star.set_zeta(zeta, bset)
            print("\nRun: %s Tcol=%f zeta=%f " % (bset, Tcol, zeta))

        dic_results[it] = star
        # it += 1
    fig = plot_spec(dic_results, times, set_bands, is_planck=True, is_filter=True)
    return fig


def save_fit_wl(fsave, mname, time, Tdil, Wdil, wl_ab=None):
    # print
    # print("{:>10s}  {:>12s}  {:>12s}  {:>12s}  ".format("Time",  "Twien", "Tcol", "zeta", "Tdil", "Wdil"))
    if fsave.endswith('.hdf5'):
        import h5py
        with h5py.File(fsave, "w") as f:
            # data = np.array([time, Tdil, Wdil])
            data = np.zeros(len(time), dtype={'names': ['time', 'Tcolor', 'W'], 'formats': [np.float] * 3})
            data['time'] = time
            data['Tcolor'] = Tdil
            data['W'] = Wdil
            ds = f.create_dataset('bbfit', data=data)
            ds.attrs['model'] = mname
            if wl_ab is not None:
                ds.attrs['wl_ab'] = '-'.join(map(str, wl_ab))
    else:
        with open(fsave, "w+") as f:
            print("{:>8s}".format("Time") +
                  ' '.join("{:>12s}".format(s) for s in ("T_dil", "W_dil")), file=f)
            # ' '.join("{:>12s}".format(s) for s in ("T_wien", "Tcol", "W", "T_dil", "W_dil")), file=f)
            # print("{:>10s}  {:>12s}  {:>12s}  {:>12s}  ".format("Time",  "Twien", "Tcol","zeta"), file=f)
            for t, *p in zip(time, Tdil, Wdil):
                print("{:8.2f}".format(t) + ' '.join("{0:12.2f}".format(s) for s in p), file=f)
    print("The temperature has been saved to %s " % fsave)
    return None


def plot_fit_wl(model, series, wl_ab, times=None, fsave=None):
    tt = model.get_tt().load()
    series_cut = series.copy(wl_ab=wl_ab)

    time = series.Time
    radius = np.interp(time, tt['time'], tt['rbb'])

    if True:
        Tcol, Twien, zeta = [], [], []
        Tdil, Wdil = [], []
        for i, t in enumerate(time):
            sp = series_cut.get_spec(i)
            # sp = series_cut.get_spec_by_time(t)
            R = radius[i]
            star_cut = ps.Star("bb_cut", sp)
            star_cut.set_distance(R)
            sp_obs = ps.Spectrum('cut', star_cut.Freq, star_cut.FluxObs)
            Tc = sp_obs.T_color
            Tw = sp_obs.T_wien
            w = np.power(Tw / Tc, 4)
            Td, W_d = sp_obs.T_color_zeta()

            Tcol.append(Tc)
            Twien.append(Tw)
            zeta.append(w)
            Tdil.append(Td)
            Wdil.append(W_d)
    else:
        Tcol = series_cut.get_T_color()
        Twien = series_cut.get_T_wien()
        zeta = np.power(Twien / Tcol, 4)

        res = series_cut.get_T_color_zeta()
        Tdil, Wdil = res[:, 0], res[:, 1]

    print("{:>8}".format("Time") +
          ' '.join("{:>12s}".format(s) for s in ("T_dil", "W_dil")))
    for t, *p in zip(time, Tdil, Wdil):
        print("{:8.2f}".format(t) + ' '.join("{0:12.2f}".format(s) for s in p))
        # print("%10.2f  %12.6f  %12.6f  %12.6f  " % p)

    # save
    if fsave is not None:
        save_fit_wl(fsave, model.Name, time, Tdil, Wdil, wl_ab=wl_ab)
        return

    # plot
    fig, ax = plt.subplots()
    marker_style = dict(linestyle=':', color='cornflowerblue', markersize=10)

    ax.semilogy(time, Tcol, 'ks-', markerfacecolor='white', markersize=3, label="Tcolor")
    ax.semilogy(time, Twien, 'r:', label="Twien")
    ax.semilogy(time, Tdil, 'ys-', markersize=3, label="T_dil")
    ax.semilogy(tt['time'], tt['Tbb'], 'g:', label="tt Tbb")
    # ax.semilogy(tt['time'], tt['Teff'], 'y:', label="tt Teff")
    # ax.semilogy(results['time'], tt['T'], 'm:', label="tt Teff")
    ax.legend()
    # wl range

    if times is not None:
        for xc in times:
            ax.axvline(x=xc, color="grey", linestyle='--')

        fig = plot_spec_wl(times, series, tt, wl_ab)
    else:
        pass
    return fig


def plot_spec_wl(times, series, tt, wl_ab, **kwargs):
    font_size = kwargs.get('font_size', 12)
    nrow = np.math.ceil(len(times)/2.)
    ncol = 2
    fig = plt.figure(figsize=(12, nrow * 4))
    plt.matplotlib.rcParams.update({'font.size': font_size})

    series_cut = series.copy(wl_ab=wl_ab)
    # radiuses = np.interp(times, tt['time'], tt['Rph'], 0, 0)
    radiuses = np.interp(times, tt['time'], tt['rbb'], 0, 0)
    Tbbes = np.interp(times, tt['time'], tt['Tbb'], 0, 0)
    marker_style = dict(linestyle=':', markersize=5)

    for i, t in enumerate(times):
        ax = fig.add_subplot(nrow, ncol, i + 1)
        spec = series.get_spec_by_time(t)
        spec.cut_flux(max(spec.Flux) * 1e-6)  # cut flux
        R = radiuses[i]

        star_bb = ps.Star("bb", spec)
        star_bb.set_distance(R)

        spec_cut = series_cut.get_spec_by_time(t)
        star_cut = ps.Star("bb_cut", spec_cut)
        star_cut.set_distance(R)

        # spectrum
        ax.semilogy(star_bb.Wl * ps.phys.cm_to_angs, star_bb.FluxWlObs, label="Spec Ph")
        # Tcolor
        spec_obs = ps.Spectrum('wbb', star_cut.Freq, star_cut.FluxObs)
        Tcol = spec_obs.T_color
        T_wien = spec_obs.T_wien
        zeta = (T_wien / Tcol) ** 4
        wbb = ps.SpectrumDilutePlanck(spec.Freq, Tcol, zeta)
        ax.semilogy(wbb.Wl * ps.phys.cm_to_angs, wbb.FluxWl, label="Tcol={:.0f} W={:.2f}".format(Tcol, zeta),
                    marker='<', **marker_style)
        # diluted
        Tdil, W = spec_obs.T_color_zeta()
        dil = ps.SpectrumDilutePlanck(spec.Freq, Tdil, W)
        ax.semilogy(dil.Wl * ps.phys.cm_to_angs, dil.FluxWl, label="Tdil={:.0f} W={:.2f}".format(Tdil, W),
                    marker='>', **marker_style)
        # Tbb
        Tbb = Tbbes[i]
        bb = ps.SpectrumPlanck(spec.Freq, Tbb)
        ax.semilogy(bb.Wl * ps.phys.cm_to_angs, bb.FluxWl, label="Tbb={:.0f}".format(Tbb),
                    marker='d', **marker_style)

        # wl range
        if wl_ab is not None:
            for xc in wl_ab:
                plt.axvline(x=xc, color="grey", linestyle='--')
        ax.legend(loc="best", prop={'size': 9})
        ax.text(0.01, 0.05, "$t_d: {:.1f}$".format(t), horizontalalignment='left', transform=ax.transAxes,
                bbox=dict(facecolor='green', alpha=0.3))
        xlim = ax.get_xlim()
        ax.set_xlim(xlim[0], min(xlim[1], 2e4))

    # fig.subplots_adjust(wspace=0, hspace=0)
    # fig.subplots_adjust(wspace=0, hspace=0, left=0.07, right=0.96, top=0.97, bottom=0.06)
    plt.subplots_adjust(left=0.07, right=0.96, top=0.97, bottom=0.06)
    return fig


def plot_spec_t(series, wl_lim=None, moments=None):
    wl_lim = wl_lim or (1e1, 5e4)
    moments = moments or (1., 3., 5, 10., 15., 25., 35., 50., 70., 100., 120., 150., 200.)
    # moments = np.arange(0.5, 200., 3)
    pos = 0
    t_data = []
    spec_array = []
    for i, t in enumerate(series.Time):
        if t > moments[pos]:
            t_data.append(t)
            spec_array.append(series.get_spec(i).Flux)
            pos += 1

    y_data = series.Wl * ps.phys.cm_to_angs
    spec_array = np.array(spec_array)

    x_data, y_data = np.meshgrid(t_data, y_data)

    x = x_data.flatten()
    y = y_data.flatten()
    z = spec_array.flatten()

    # filters
    is_z = z > np.max(z) * 1.e-20
    is_y = (y > wl_lim[0]) & (y < wl_lim[1])

    is_good = is_y & is_z
    x = x[is_good]
    y = y[is_good]
    z = z[is_good]

    fig = plt.figure()

    # scatter3D
    if True:
        ax = Axes3D(fig)
        ax.scatter3D(x, y, z, c=z, cmap=plt.cm.viridis)  # jet)
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
        ax.scatter3D(x, y, z, c=z, cmap=plt.cm.viridis)
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
    sp = ps.SpectrumDilutePlanck(star.Freq, star.get_Tcol(bset), star.get_zeta(bset) ** 2)
    star_bb = ps.Star("bb", sp)
    star_bb.set_radius_ph(star.radius_ph)
    star_bb.set_distance(star.distance)
    return star_bb


def epsilon_mag(x, freq, mag, bands, radius, dist):
    temp_color, zeta = x
    sp = ps.SpectrumDilutePlanck(freq, temp_color, zeta ** 2)

    star = ps.Star("bb", sp)
    star.set_radius_ph(radius)
    star.set_distance(dist)
    mag_bb = {b: star.magAB(ps.band.band_by_name(b)) for b in bands}
    e = 0
    for b in bands:
        e += abs(mag[b] - mag_bb[b])
    return e


def compute_tcolor(star, bands):
    mags = {}
    for n in bands:
        b = ps.band.band_by_name(n)
        mags[n] = star.magAB(b)

    Tcol, zeta = fmin(epsilon_mag, x0=np.array([1.e4, 1]),
                      args=(star.Freq, mags, bands, star.radius_ph, star.distance), disp=False)
    return Tcol, zeta


def interval2float(v):
    a = str(v).split(':')
    if len(a) == 2 and len(a[0]) == 0:
        return 0., float(a[1])
    elif len(a) == 2 and len(a[1]) == 0:
        return float(a[0]), float("inf")
    elif len(a) == 2:
        return list(map(np.float, str(v).split(':')))
    elif len(a) == 1:
        return float(v[0])
    else:
        return None


def usage():
    bands = ps.band.band_get_names()
    print("Usage: show F(t,nu) from ph-file")
    print("  spec.py [params]")
    print("  -b <set_bands>: delimiter '_'. Default: B-V.\n"
          "     Available: " + '-'.join(sorted(bands)))
    print("  -f  force mode: rewrite tcolor-files even if it exists")
    print("  -i <model name>.  Example: cat_R450_M15_Ni007_E7")
    print("  -k  k-correction: z:Srest:Sobs. Example: 0.5:V:R")
    print("  -o  options: [fit, wl] - plot spectral fit")
    print("  -p <model path(directory)>, default: ./")
    print("  -s <file-name> without extension. Save plot to pdf-file. Default: spec_<file-name>.pdf")
    print("  -t  time interval [day]. Example: 5.:75.")
    print("  -x  wave length interval [A]. Example: 1.:25e3")
    print("  -w  write data to file. Example: flux to magAB, Tcolor and W [.hdf5]")
    print("  -h  print usage")


def write_magAB(series, d=10):
    """
    Write SED in AB mag
    :param series:
    :param d: distance [pc]
    :return:
    """
    print("times: {0} freqs: {1}  [{1}]".format(len(series.Time), len(series.get_spec(0).Freq), series.Name))
    print("{0}".format(' '.join(map(str, series.get_spec(0).Freq))))
    for i, t in enumerate(series.Time):
        spec = series.get_spec(i)
        # freq = s.Freq
        # flux = spec.Flux
        # mAB = rf.Flux2MagAB(flux)
        star = ps.Star("%s: %f" % (series.Name, t), spec=spec, is_flux_eq_luminosity=True)
        star.set_distance(d * ps.phys.pc)
        # print("{0}  {1} ".format(t, '  '.join(map(str, star.Flux))))
        print("{0}  {1} ".format(t, '  '.join(map(str, star.FluxAB))))
        # print("{0}  {1} ".format(t, '  '.join(map(str, flux))))
        # print("{0}  {1} ".format(t, '  '.join(map(str, mAB))))


def plot_kcorr(times, kcorr, figsize=(12, 12)):
    fig, ax = plt.subplots(figsize=figsize)

    ax.plot(times, kcorr, marker='o', ls='')
    ax.set_xlabel('Time')
    ax.set_ylabel('K-correction')

    return fig


def kcorr_save(fname, times, kcorr):
    with open(fname, "w") as f:
        print("{:>8s}  {:>8s}".format("Time", "kcorr"), file=f)
        for t, k in zip(times, kcorr):
            print("{:8.2f}  {:12.2f}".format(t, k), file=f)
    print("The k-corrections has been saved to %s " % fname)


def main():
    is_save_plot = False
    is_kcorr = False
    is_fit = False
    is_fit_wl = False
    is_write = False

    fsave = None
    fplot = None
    z_sn = 0.
    bn_rest = None
    bn_obs = None

    try:
        opts, args = getopt.getopt(sys.argv[1:], "b:fhsup:i:k:o:t:w:x:")
    except getopt.GetoptError as err:
        print(str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    name = ''
    path = os.getcwd()
    ps.Band.load_settings()
    # name = 'cat_R500_M15_Ni006_E12'

    if not name:
        if len(opts) == 0:
            usage()
            sys.exit(2)
        for opt, arg in opts:
            if opt == '-i':
                name = str(arg)
                break

    # set_bands = ['B-V']
    # set_bands = ['B-V', 'B-V-I']
    # set_bands = ['U-B', 'U-B-V', 'B-V']
    t_ab = None
    wl_ab = None
    set_bands = ['B-V', 'B-V-I', 'V-I']
    # set_bands = ['B-V', 'B-V-I', 'V-I', 'J-H-K']
    times = [5., 15., 30., 60., 90., 120.]

    for opt, arg in opts:
        if opt == '-b':
            set_bands = str(arg).split('_')
            for bset in set_bands:
                for b in bset.split('-'):
                    if not ps.band.band_is_exist(b):
                        print('No such band: ' + b)
                        sys.exit(2)
            continue
        if opt == '-w':
            is_write = True
            fsave = arg
            continue
        if opt == '-s':
            is_save_plot = True
            if len(arg) > 0:
                fplot = str(arg).strip()
            continue
        if opt == '-x':
            wl_ab = interval2float(arg)
            # wl_ab = [np.float(s) for s in (str(arg).split(':'))]
            continue
        if opt == '-t':
            t_ab = list(map(float, arg.split(':')))  # interval2float(arg)
            if len(t_ab) > 1:
                times = t_ab
            continue
        if opt == '-o':
            ops = str(arg).split(':')
            is_fit = "fit" in ops
            is_fit_wl = "wl" in ops
            continue
        if opt == '-k':
            ops = str(arg).split(':')
            if len(ops) == 3:
                z_sn = float(ops[0])
                bn_rest = ops[1].strip()
                bn_obs = ops[2].strip()
                is_kcorr = True
            else:
                raise ValueError('Args: {} should be string as "z:Srest:Sobs"'.format(arg))
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

    if not name:
        print("No model. Use key -i.")
        sys.exit(2)

    model = ps.Stella(name, path=path)
    series = model.get_ph(t_diff=1.05)

    if not model.is_ph:
        print("No ph-data for: " + str(model))
        return None

    if is_fit:
        if is_write:
            if fsave is None or len(fsave) == 0:
                fsave = "spec_%s" % name
            print("Save series to %s " % fsave)
            series_cut = series.copy(t_ab=t_ab, wl_ab=wl_ab)
            write_magAB(series_cut)
            sys.exit(2)

        if not model.is_tt:
            print("Error in fit-band: no tt-data for: " + str(model))
            sys.exit(2)
        series = model.get_ph(t_diff=1.05)
        series_cut = series.copy(t_ab=t_ab, wl_ab=wl_ab)
        fig = plot_fit_bands(model, series_cut, set_bands, times)
    elif is_kcorr:
        times, kcorr = [], []
        for t, k in ps.rf.rad_func.kcorrection(series, z_sn, bn_rest, bn_obs):
            times.append(t)
            kcorr.append(k)
        if is_write:
            if fsave is None or len(fsave) == 0 or fsave == '1':
                fsave = os.path.join(os.path.expanduser('~/'), "kcorr_%s" % name) + '.txt'
            kcorr_save(fsave, times, kcorr)
            sys.exit(3)
        else:
            fig = plot_kcorr(times, kcorr)
    elif is_fit_wl:
        if not model.is_tt:
            print("Error in fit-wave: no tt-data for: " + str(model))
            sys.exit(2)
        if is_write:
            if fsave is None or len(fsave) == 0 or fsave == '1':
                fsave = os.path.join(os.path.expanduser('~/'), "temp_%s" % name) + '.txt'
            series = series.copy(t_ab=t_ab)
            plot_fit_wl(model, series, wl_ab, times, fsave=fsave)  # just save data
            sys.exit(3)

        fig = plot_fit_wl(model, series, wl_ab, times)
    else:
        series = model.get_ph(t_diff=1.05)
        series_cut = series.copy(t_ab=t_ab, wl_ab=wl_ab)
        fig = plot_spec_poly(series_cut)
        print("Plot spectral F(t,nu): " + str(model))

    if fig is not None:
        if is_save_plot:
            if fplot is None or len(fplot) == 0:
                fplot = "spec_%s" % name
            d = os.path.expanduser('~/')
            fplot = os.path.join(d, os.path.splitext(fplot)[0]) + '.pdf'

            print("Save plot to %s " % fplot)
            fig.savefig(fplot, bbox_inches='tight')
        else:
            # plt.grid()
            plt.show()


if __name__ == '__main__':
    main()
