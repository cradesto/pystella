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
from matplotlib import gridspec

import pystella.util.rf as rf
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
    print "Usage:"
    print "  plot_spec.py [params]"
    print "  -b <set_bands>: delimiter '_'. Default: B-V.\n" \
          "     Available: " + '-'.join(sorted(bands))
    print "  -i <model name>.  Example: cat_R450_M15_Ni007_E7"
    print "  -p <model path(directory)>, default: ./"
    print "  -s  silence mode: no info, no plot"
    print "  -f  force mode: rewrite tcolor-files even if it exists"
    print "  -t  plot time points"
    print "  -w  write magnitudes to file, default 'False'"
    print "  -h  print usage"


def main():
    is_silence = False
    is_fit = True

    try:
        opts, args = getopt.getopt(sys.argv[1:], "afhnsuwtp:i:b:")
    except getopt.GetoptError as err:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    name = ''
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
    times = [15., 30, 60, 120]

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
            continue
        if opt == '-n':
            is_plot_Tnu = True
            continue
        if opt == '-t':
            is_plot_time_points = True
            continue
        if opt == '-w':
            is_save = True
            continue
        if opt == '-f':
            is_force = True
            is_save = True
            continue
        if opt == '-a':
            is_fit = True
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

    if not model.is_spec_data:
        print "No ph-data for: " + str(model)
        return None

    if not model.is_tt_data:
        print "No tt-data for: " + str(model)
        return None

    serial_spec = model.read_serial_spectrum(t_diff=1.05)
    tt = model.read_tt_data()
    tt = tt[tt['time'] > min(times) - 1.]  # time cut  days

    Rph_spline = interpolate.splrep(tt['time'], tt['Rph'], s=0)

    distance = rf.pc_to_cm(10.)  # pc for Absolute magnitude
    dic_results = {}  # dict((k, None) for k in names)
    it = 0
    for time in times:
        print "\nRun: %s t=%f [%d/%d]" % (name, time, it, len(times))
        spec = serial_spec.get_spec_by_time(time)
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
