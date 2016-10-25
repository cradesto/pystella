import csv
import numpy as np
import os

import matplotlib.pyplot as plt

import pystella.rf.rad_func as rf
from pystella.model.stella import Stella
from pystella.rf import band
from pystella.rf import extinction
from pystella.rf.lc import LightCurve

__author__ = 'bakl'

lc_colors = band.bands_colors()

lc_lntypes = dict(U="-", B="-", V="-", R="-", I="-",
                  UVM2="-.", UVW1="-.", UVW2="-.",
                  F125W="--", F160W="-.", F140W="--", F105W="-.", F435W="--", F606W="-.", F814W="--",
                  u="--", g="--", r="--", i="--", z="--",
                  bol='-')

markers = {u'D': u'diamond', 6: u'caretup', u's': u'square', u'x': u'x',
           5: u'caretright', u'^': u'triangle_up', u'd': u'thin_diamond', u'h': u'hexagon1',
           u'+': u'plus', u'*': u'star', u'o': u'circle', u'p': u'pentagon', u'3': u'tri_left',
           u'H': u'hexagon2', u'v': u'triangle_down', u'8': u'octagon', u'<': u'triangle_left'}
markers = markers.keys()


def lbl(b, band_shift):
    shift = band_shift[b]
    l = b
    if shift == int(shift):
        shift = int(shift)
    if shift > 0:
        l += '+' + str(shift)
    elif shift < 0:
        l += '-' + str(abs(shift))
    return l


def plot_ubv_models(ax, models_dic, bands, band_shift, xlim=None, ylim=None,
                    colors=lc_colors, is_time_points=False):
    is_compute_x_lim = xlim is None
    is_compute_y_lim = ylim is None

    t_points = [0.2, 1, 2, 3, 4, 5, 10, 20, 40, 80, 150]

    lw = 2
    mi = 0
    ib = 0
    x_max = []
    y_mid = []
    lc_min = {}
    for mname, mdic in models_dic.iteritems():
        mi += 1
        for bname in bands:
            ib += 1
            x = mdic['time']
            y = mdic[bname] + band_shift[bname]
            bcolor = colors[bname]

            if len(models_dic) == 1:
                ax.plot(x, y, label='%s  %s' % (lbl(bname, band_shift), mname), color=bcolor, ls="-", linewidth=lw)
            else:
                ax.plot(x, y, marker=markers[mi % (len(markers) - 1)], label='%s  %s' % (lbl(bname, band_shift), mname),
                        markersize=3, color=bcolor, ls="--", linewidth=lw)

            if is_time_points:
                integers = [np.abs(x - t).argmin() for t in t_points]  # set time points
                for (X, Y) in zip(x[integers], y[integers]):
                    ax.annotate('{:.0f}'.format(X), xy=(X, Y), xytext=(-10, 20), ha='right',
                                textcoords='offset points', color=bcolor,
                                arrowprops=dict(arrowstyle='->', shrinkA=0))
            idx = np.argmin(y)
            lc_min[bname] = (x[idx], y[idx])
            if is_compute_x_lim:
                x_max.append(np.max(x))
            if is_compute_y_lim:
                y_mid.append(np.min(y))

    if is_compute_x_lim:
        xlim = [-10, np.max(x_max) + 10.]
    if is_compute_y_lim:
        ylim = [np.min(y_mid) + 7., np.min(y_mid) - 2.]

    ax.set_xlim(xlim)
    ax.invert_yaxis()
    ax.set_ylim(ylim)

    return lc_min


def plot_models_curves(ax, models_curves, band_shift=None, xlim=None, ylim=None, lc_types=None, colors=lc_colors,
                       lw=2.):
    is_compute_x_lim = xlim is None
    is_compute_y_lim = ylim is None

    if lc_types is None:
        lc_types = dict((name, '-') for name in models_curves.keys())  # solid line

    mi, ib = 0, 0
    x_max = []
    y_mid = []
    for mname, curves in models_curves.iteritems():
        mi += 1
        bands = curves.BandNames
        if band_shift is None:
            mshift = dict((bname, 0.) for bname in bands)  # no y-shift
        for bname in bands:
            ib += 1
            x = curves.TimeDef
            y = curves[bname].Mag + mshift[bname]
            ax.plot(x, y, label='%s  %s' % (lbl(bname, mshift), mname),
                    color=colors[bname], ls=lc_types[mname], linewidth=lw)
            if is_compute_x_lim:
                x_max.append(np.max(x))
            if is_compute_y_lim:
                y_mid.append(np.min(y))

    if is_compute_x_lim:
        xlim = [-10, np.max(x_max) + 10.]
    if is_compute_y_lim:
        ylim = [np.min(y_mid) + 7., np.min(y_mid) - 2.]

    ax.set_xlim(xlim)
    ax.invert_yaxis()
    ax.set_ylim(ylim)


def plot_models_curves_fixed_bands(ax, models_curves, bands, band_shift=None, xlim=None, ylim=None, lc_types=None,
                                   colors=lc_colors,
                                   lw=2.):
    is_compute_x_lim = xlim is None
    is_compute_y_lim = ylim is None

    if band_shift is None:
        band_shift = dict((bname, 0) for bname in bands)  # no y-shift

    if lc_types is None:
        lc_types = dict((name, '-') for name in models_curves.keys())  # solid line

    mi, ib = 0, 0
    x_max = []
    y_mid = []
    for mname, curves in models_curves.iteritems():
        mi += 1
        for bname in bands:
            ib += 1
            x = curves.TimeDef
            y = curves[bname].Mag
            ax.plot(x, y, label='%s  %s' % (lbl(bname, band_shift), mname),
                    color=colors[bname], ls=lc_types[mname], linewidth=lw)
            if is_compute_x_lim:
                x_max.append(np.max(x))
            if is_compute_y_lim:
                y_mid.append(np.min(y))

    if is_compute_x_lim:
        xlim = [-10, np.max(x_max) + 10.]
    if is_compute_y_lim:
        ylim = [np.min(y_mid) + 7., np.min(y_mid) - 2.]

    ax.set_xlim(xlim)
    ax.invert_yaxis()
    ax.set_ylim(ylim)


def plot_bands(dict_mags, bands, title='', fname='', distance=10., is_time_points=True):
    fig = plt.figure()
    ax = fig.add_axes((0.1, 0.3, 0.8, 0.65))
    # ax.set_title(''.join(bands) + ' filter response')

    # colors = band.bands_colors()
    # colors = dict(U="blue", B="cyan", V="black", R="red", I="magenta",
    #               J="blue", H="cyan", K="black",
    #               UVM2="green", UVW1="red", UVW2="blue",
    #               g="black", r="red", i="magenta", u="blue", z="magenta")
    # band_shift = dict(U=6.9, B=3.7, V=0, R=-2.4, I=-4.7,
    #                   UVM2=11.3, UVW1=10, UVW2=13.6,
    #                   u=3.5, g=2.5, r=-1.2, i=-3.7, z=-4.2)
    band_shift = dict((k, 0) for k, v in lc_colors.items())  # no y-shift

    t_points = [2, 5, 10, 20, 40, 80, 150]

    dm = 5 * np.log10(distance) - 5  # distance module
    # dm = 0
    xlims = [-10, 200]
    ylims = [-12, -19]
    ylims += dm
    is_auto_lim = True
    if is_auto_lim:
        ylims = [0, 0]

    x = dict_mags['time']
    is_first = True
    for n in bands:
        y = dict_mags[n]
        y += dm + band_shift[n]
        ax.plot(x, y, label=lbl(n, band_shift), color=lc_colors[n], ls=lc_lntypes[n], linewidth=2.0)
        # ax.plot(x, y, label=lbl(n, band_shift), color=colors[n], ls=lntypes[n], linewidth=2.0, marker='s')
        if is_time_points and is_first:
            is_first = False
            integers = [np.abs(x - t).argmin() for t in t_points]  # set time points
            for (X, Y) in zip(x[integers], y[integers]):
                plt.annotate('{:.0f}'.format(X), xy=(X, Y), xytext=(10, -30), ha='right',
                             textcoords='offset points',
                             arrowprops=dict(arrowstyle='->', shrinkA=0))

        if is_auto_lim:
            if ylims[0] < max(y[len(y) / 2:]) or ylims[0] == 0:
                ylims[0] = max(y[len(y) / 2:])
            if ylims[1] > min(y) or ylims[1] == 0:
                ylims[1] = min(y)

    ylims = np.add(ylims, [1, -1])
    ax.invert_yaxis()
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.legend()
    ax.set_ylabel('Magnitude')
    ax.set_xlabel('Time [days]')
    # ax.set_title(title)
    fig.text(0.17, 0.07, title, family='monospace')
    ax.grid()
    if fname != '':
        plt.savefig("ubv_%s.png" % fname, format='png')
    # ax.show()
    # plt.close()
    return ax


def curves_plot(curves, ax=None, xlim=None, ylim=None, title=None, fname='', **kwargs):
    is_line = 'ls' in kwargs or 'lt' not in kwargs
    ls = kwargs.pop('ls', {lc.Band.Name: '-' for lc in curves})
    lt = kwargs.pop('lt', {lc.Band.Name: 'o' for lc in curves})
    colors = kwargs.pop('colors', lc_colors)
    linewidth = kwargs.pop('linewidth', 2.0)
    markersize = kwargs.pop('markersize', 5)
    rect = kwargs.pop('rect', (0.1, 0.3, 0.8, 0.65))
    fontsize = kwargs.pop('fontsize', 18)

    is_new_fig = ax is None
    if is_new_fig:
        plt.matplotlib.rcParams.update({'font.size': 14})
        fig = plt.figure()
        # fig = plt.figure(num=None, figsize=(7, 11), dpi=100, facecolor='w', edgecolor='k')
        # ax = fig.add_axes()
        ax = fig.add_axes(rect)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                         ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fontsize)
            # ax = fig.add_axes((0.1, 0.3, 0.8, 0.65))

    # plt.title(''.join(bands) + ' filter response')
    is_xlim = False
    is_ylim = False
    if xlim is None:
        is_xlim = True
        xlim = [float('inf'), float('-inf')]
    if ylim is None:
        is_ylim = True
        ylim = [float('-inf'), float('inf')]

    for lc in curves:
        x = lc.Time
        y = lc.Mag
        bname = lc.Band.Name
        if is_line:
            ax.plot(x, y, label='%s' % bname, color=colors[bname], ls=ls[bname], linewidth=linewidth)
        else:
            ax.plot(x, y, label='%s' % bname, color=colors[bname], ls=None, marker=lt[bname], markersize=markersize)

        if is_xlim:
            xlim[0] = min(xlim[0], np.min(x))
            xlim[1] = max(xlim[1], np.max(x))
        if is_ylim:
            ylim[0] = max(ylim[0], np.max(y))
            ylim[1] = min(ylim[1], np.min(y))

    if is_ylim:
        ylim = [ylim[1] + 10, ylim[1] - 2]
    # ylim = np.add(ylim, [1, -1])
    ax.invert_yaxis()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.legend()
    ax.set_ylabel('Magnitude')
    ax.set_xlabel('Time [days]')
    ax.grid()
    if title is not None:
        plt.title(title)
        # ax.text(0.17, 0.07, title, family='monospace')
    if fname != '':
        # plt.savefig("ubv_%s.png" % fname, format='png')
        plt.savefig(fname)
    # if is_new_fig:
    #     plt.show()
    # plt.close()
    return ax


def compute_mag(name, path, bands, ext=None, z=0., distance=10., magnification=1., t_diff=1.05, is_show_info=True,
                is_save=False):
    """
        Compute magnitude in bands for the 'name' model.
    :param name: the name of a model and data files
    :param path: the directory with data-files
    :param bands: photometric bands
    :param ext: extinction
    :param z: redshift, default 0
    :param distance: distance to star in parsec, default 10 pc
    :param magnification: gravitational lensing magnification
    :param t_diff:  depression between time points
    :param is_show_info: flag to write some information, default True
    :param is_save: flag to save result in file, default False
    :return: dictionary with keys = bands, value = star's magnitudes
    """
    model = Stella(name, path=path)
    if is_show_info:
        print ''
        model.show_info()

    if not model.is_ph_data:
        model.show_info()
        print "Error: No data for: " + str(model)
        return None

    # serial_spec = model.read_serial_spectrum(t_diff=0.)
    serial_spec = model.read_series_spectrum(t_diff=t_diff)
    mags = serial_spec.mags_bands(bands, z=z, d=rf.pc_to_cm(distance), magnification=magnification)

    if mags is not None:
        fname = os.path.join(path, name + '.ubv')
        if is_save:
            mags_save(mags, bands, fname)
            print "Magnitudes have been saved to " + fname

    if is_show_info:
        # print the time of maximum LC
        tmin = 2.0
        t = mags['time']
        for n in bands:
            t_min = t[t > tmin][mags[n][t > tmin].argmin()]
            print "t_max(%s) = %f" % (n, t_min)

    if ext is not None:  # add extinction
        for n in bands:
            mags[n] = mags[n] + ext[n]

    return mags


def curves_save(curves, fname):
    """
       Save curves to CSV-format
    :param curves:
    :param fname:
    :return:
    """
    if curves.Length > 0:
        with open(fname, 'wb') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['{:^8s}'.format(x) for x in ['time'] + curves.BandNames])
            for i, (row) in enumerate(zip(curves.TimeDef, *[curves.get(b) for b in curves.BandNames])):
                # row = row[-1:] + row[:-1]  # make time first column
                writer.writerow(['{:8.3f}'.format(x) for x in row])
                # writer.writerow(['{:3.4e}'.format(x) for x in row])
    else:
        print "Nothing to save: curves.Length=%d" % curves.Length


def curves_compute(name, path, bands, z=0., distance=10., magnification=1.,
                   t_beg=0., t_end=None, is_show_info=False, is_save=False):
    """
        Compute magnitude in bands for the 'name' model.
    :param name: the name of a model and data files
    :param path: the directory with data-files
    :param bands: photometric bands
    :param z: redshift, default 0
    :param distance: distance to star in parsec, default 10 pc
    :param magnification: gravitational lensing magnification
    :param t_end:
    :param t_beg:
    :param is_show_info: flag to write some information, default True
    :param is_save: flag to save result in file, default False
    :return: dictionary with keys = bands, value = star's magnitudes
    """
    if len(bands) == 0:
        raise ValueError("You have not set any bands for model: " + str(name))

    model = Stella(name, path=path)
    if not model.is_ph_data:
        model.show_info()
        raise ValueError("Error: No spectral data for: " + str(model))

    if is_show_info:
        print ''
        model.show_info()

    # serial_spec = model.read_serial_spectrum(t_diff=0.)
    serial_spec = model.read_series_spectrum(t_diff=1.05, t_beg=t_beg, t_end=t_end)
    curves = serial_spec.flux_to_curves(bands, z=z, d=rf.pc_to_cm(distance), magnification=magnification)
    # curves = SetLightCurve(name)
    # for n in bands:
    #     b = band.band_by_name(n)
    #     lc = serial_spec.flux_to_curve(b, z=z, dl=rf.pc_to_cm(distance), magnification=magnification)
    #     # time = serial_spec.times * (1. + z)
    #     # lc = LightCurve(b, time, mags)
    #     curves.add(lc)

    if is_save:
        fname = os.path.join(path, name + '.ubv')
        curves_save(curves, fname)
        print "Magnitudes have been saved to " + fname

    if is_show_info:
        # print the time of maximum LC
        tmin = 2.0
        t = curves.TimeDef
        for n in bands:
            t_min = t[t > tmin][curves[n][t > tmin].argmin()]
            print "t_max(%s) = %f" % (n, t_min)

    return curves


def curves_reddening(curves, ebv, z=None, law=extinction.law_default):
    if ebv < 0.:
        raise ValueError("ebv should be > 0")

    if isinstance(curves, LightCurve):
        lc = curves
        bands = list(lc.Band.Name)
        if z is not None and z > 0.1:
            ext = extinction.reddening_law_z(ebv=ebv, bands=bands, z=z, law=law)
        else:
            ext = extinction.reddening_law(ebv=ebv, bands=bands, law=law)
        lc.mshift = lc.mshift + ext[bands[0]]
        return lc
    else:
        bands = tuple(curves.BandNames)
        if z is not None and z > 0.1:
            ext = extinction.reddening_law_z(ebv=ebv, bands=bands, z=z, law=law)
        else:
            ext = extinction.reddening_law(ebv=ebv, bands=bands, law=law)

        for n in bands:
            curves[n].mshift = curves[n].mshift + ext[n]
        return curves


def mags_save(dictionary, bands, fname):
    with open(fname, 'wb') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['{:^8s}'.format(x) for x in ['time'] + bands])
        for i, (row) in enumerate(zip(*[dictionary[k] for k in 'time'.split() + bands])):
            # row = row[-1:] + row[:-1]  # make time first column
            writer.writerow(['{:8.3f}'.format(x) for x in row])
            # writer.writerow(['{:3.4e}'.format(x) for x in row])
