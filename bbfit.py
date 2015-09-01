#!/usr/bin/python
# -*- coding: utf-8 -*-
from matplotlib import gridspec

from pystella.rf.star import Star

__author__ = 'bakl'

import os
import sys
import getopt
from os.path import isfile, join, dirname
import csv

import numpy as np
import scipy.optimize as optim
from scipy import interpolate
import matplotlib.pyplot as plt

from pystella.rf import band, spectrum
from pystella.model.stella import Stella
import pystella.util.rf as rf

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))


def plot_zeta(models_dic, set_bands, title='', is_fit=False, is_time_points=False):
    colors = {"B-V": "blue", 'B-V-I': "cyan", 'V-I': "black"}
    lntypes = {"B-V": "-", 'B-V-I': "-.", 'V-I': "--"}
    markers = {u'D': u'diamond', 6: u'caretup', u's': u'square', u'x': u'x'
        , 5: u'caretright', u'^': u'triangle_up', u'd': u'thin_diamond', u'h': u'hexagon1'
        , u'+': u'plus', u'*': u'star', u'o': u'circle', u'p': u'pentagon', u'3': u'tri_left'
        , u'H': u'hexagon2', u'v': u'triangle_down', u'8': u'octagon', u'<': u'triangle_left'}
    markers = markers.keys()

    t_points = [1, 2, 3, 4,  5, 10, 20, 40, 80, 150]

    xlim = [0, 30000]
    ylim = [0, 3]
    is_auto_lim = False
    if is_auto_lim:
        xlim = [0, 0]
        ylim = [0, 0]

    # setup figure
    plt.matplotlib.rcParams.update({'font.size': 8})
    fig = plt.figure(num=None, figsize=(7, 7), dpi=100, facecolor='w', edgecolor='k')
    gs1 = gridspec.GridSpec(len(set_bands), 1)
    # gs1 = gridspec.GridSpec(2, 4, width_ratios=(8, 1, 8, 1))
    gs1.update(wspace=0.3, hspace=0.3, left=0.05, right=0.95)

    ax_cache = {}
    mi = 0
    for mname, mdic in models_dic.iteritems():
        mi += 1
        i = 0
        for bset in set_bands:
            x = mdic[bset]['Tcol']
            y = mdic[bset]['zeta']
            z = mdic[bset]['time']
            if bset in ax_cache:
                ax = ax_cache[bset]
            else:
                ax = fig.add_subplot(gs1[i, 0])
                ax_cache[bset] = ax
            ax.plot(x, y, marker=markers[mi % (len(markers) - 1)], label=mname, markersize=5, color=colors[bset],
                    ls="-.",
                    linewidth=1.5)

            if is_fit:  # dessart
                yd = zeta_fit_dessart(x, bset)
                ax.plot(x, yd, color='red', ls="-", linewidth=1.5)
            if is_auto_lim:
                xlim[0], xlim[1] = min(np.append(x, xlim[0])), max(np.append(x, xlim[1]))
                ylim[0], ylim[1] = min(np.append(x, ylim[0])), max(np.append(y, ylim[1]))
            if is_time_points:
                integers = [np.abs(z - t).argmin() for t in t_points]  # set time points
                for (X, Y, Z) in zip(x[integers], y[integers], z[integers]):
                    ax.annotate('{:.0f}'.format(Z), xy=(X, Y), xytext=(-10, 20), ha='right',
                                textcoords='offset points',
                                arrowprops=dict(arrowstyle='->', shrinkA=0))

            ax.legend(prop={'size': 5})
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_ylabel('Zeta')
            ax.set_xlabel('T_color')
            ax.set_title(bset)
            i += 1
            t_min = z[y.argmin()]
            print "t_min=%f in %s" % (t_min, bset)

    # plt.xlim(xlim)
    # ylim = np.add(ylim, [1, -1])
    # plt.ylim(ylim)
    #     plt.title('; '.join(set_bands) + ' filter response')
    plt.grid()
    plt.show()


def zeta_fit_dessart(Tcol, bset):
    """
    Zeta fit from Dessart, L., & Hillier, D. J. (2005). doi:10.1051/0004-6361:20053217
    :param Tcol:
    :param bset:
    :return:
    """
    a = {
        'B-V': [0.47188, -0.25399, 0.32630],
        'B-V-I': [0.63241, -0.38375, 0.28425],
        'V-I': [0.81662, -0.62896, 0.33852],
        'J-H-K': [0.10786, 1.12374, 0.]
    }
    zeta = 0
    i = 0
    for ai in a[bset]:
        zeta += ai * (1.e4 / Tcol) ** i
        i += 1
    return zeta


def epsilon(x, freq, mag, bands, radius, dist):
    temp_color, zeta = x
    sp = spectrum.SpectrumPlanck(freq, temp_color)
    sp.correct_zeta(zeta)

    star = Star("bb", sp)
    star.set_radius_ph(radius)
    star.set_distance(dist)
    mag_bb = {b: star.flux_to_mag(band.band_by_name(b)) for b in bands}
    e = 0
    for b in bands:
        e += abs(mag[b] - mag_bb[b])
    return e


def compute_Tcolor_zeta(mags, tt, bands, freq, dist):
    temp = list()
    zeta_radius = list()
    Rph_spline = interpolate.splrep(tt['time'], tt['Rph'], s=0)
    for nt in range(len(mags['time'])):
        t = mags['time'][nt]
        if t < min(tt['time']):
            continue
        if t > max(tt['time']):
            break
        mag = {b: mags[b][nt] for b in bands}
        radius = interpolate.splev(t, Rph_spline)
        tcolor, w = optim.fmin(epsilon, x0=np.array([1.e4, 1]), args=(freq, mag, bands, radius, dist), disp=0)
        temp.append(tcolor)
        zeta_radius.append(w)
    return temp, zeta_radius


def compute_tcolor(name, path, bands, t_cut=2.):
    model = Stella(name, path=path)

    if not model.is_spec_data:
        print "No ph-data for: " + str(model)
        return None

    if not model.is_tt_data:
        print "No tt-data for: " + str(model)
        return None

    # serial_spec = model.read_serial_spectrum(t_diff=0.)
    serial_spec = model.read_serial_spectrum(t_diff=1.05)

    distance = rf.pc_to_cm(10.)  # pc for Absolute magnitude
    l = list()
    for n in bands:
        b = band.band_by_name(n)
        l.append(serial_spec.flux_to_mags(b, dl=distance))

    mags = np.array(np.zeros(len(l[0])),
                    dtype=np.dtype({'names': ['time'] + bands, 'formats': [np.float64] * (len(bands) + 1)}))
    mags['time'] = serial_spec.times
    for n in bands:
        mags[n] = l.pop(0)
    # mags = np.array([serial_spec.times] + l,
    #                 dtype=np.dtype({'names': ['time']+bands, 'formats':  [np.float64] * (len(bands)+1)}))

    # read R_ph
    tt = model.read_tt_data()

    # time cut  days
    mags = mags[mags['time'] > t_cut]
    tt = tt[tt['time'] > t_cut]

    # fit mags by B(T_col) and get \zeta\theta & T_col
    Tcolors, zetaR = compute_Tcolor_zeta(mags, tt=tt, bands=bands, freq=serial_spec.freq, dist=distance)

    # show results
    res = np.array(np.zeros(len(Tcolors)),
                   dtype=np.dtype({'names': ['time', 'Tcol', 'zeta'], 'formats': [np.float64] * 3}))
    res['time'] = mags['time']
    res['Tcol'] = Tcolors
    res['zeta'] = zetaR

    return res


def cache_name(name, path, bands):
    fname = os.path.join(path, "%s.%s.zeta" % (name, bands))
    return fname


def cache_save(narr, fname):
    names = narr.dtype.names
    with open(fname, 'wb') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['{:^8s}'.format(x) for x in names])
        for i, (row) in enumerate(zip(*[narr[k] for k in names])):
            writer.writerow(['{:8.3f}'.format(x) for x in row])


def cache_load(fname):
    header = 'time Tcol zeta'
    names = map(str.strip, header.split())
    dtype = np.dtype({'names': names, 'formats': [np.float64] * len(names)})
    block = np.loadtxt(fname, skiprows=1, dtype=dtype)
    return block


def usage():
    bands = band.band_get_names().keys()
    print "Usage:"
    print "  bbfit.py [params]"
    print "  -b <set_bands>: delimiter '_'. Default: B-V-I_B-V_V-I.\n" \
          "     Available: " + '-'.join(sorted(bands))
    print "  -i <model name>.  Example: cat_R450_M15_Ni007_E7"
    print "  -d <model directory>, default: ./"
    print "  -e <model extension> is used to define model name, default: tt "
    print "  -s  silence mode: no info, no plot"
    print "  -f  force mode: rewrite tcolor-files even if it exists"
    print "  -a  plot Dessart&Hillier fit"
    print "  -w  write magnitudes to file, default 'False'"
    print "  -h  print usage"


def main(name='', path='./', is_force=False, is_save=False):
    is_silence = False
    is_fit = False
    model_ext = '.tt'

    try:
        opts, args = getopt.getopt(sys.argv[1:], "afswhd:e:i:b:")
    except getopt.GetoptError as err:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(opts) == 0:
        usage()
        sys.exit(2)

    if not name:
        for opt, arg in opts:
            if opt == '-i':
                path = ROOT_DIRECTORY
                name = str(arg)
                break

    set_bands = ['B-V-I', 'B-V', 'V-I']

    for opt, arg in opts:
        if opt == '-e':
            model_ext = '.' + arg
            continue
        if opt == '-b':
            set_bands = str(arg).split('_')
            for bset in set_bands:
                if not band.band_is_exist(bset):
                    print 'No such band: ' + bset
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
        if opt == '-a':
            is_fit = True
            continue
        if opt == '-d':
            path = str(arg)
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
                    dic[bset] = compute_tcolor(name, path, bset.split('-'), t_cut=0.9)
                    if is_save:
                        print "Save Tcolor & Zeta for %s in %s." % (bset, fname)
                        cache_save(dic[bset], fname=fname)

            dic_results[name] = dic
            print "Finish: %s" % name
        plot_zeta(dic_results, set_bands, is_fit=is_fit, is_time_points=True)
    else:
        print "There are no models in the directory: %s with extension: %s " % (path, model_ext)


if __name__ == '__main__':
    main()
    # main(name="cat_R1000_M15_Ni007_E15", path="~/Sn/Release/seb_git/res/tt",, is_force=True, is_save=True)
