#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import sys
import getopt
from os.path import dirname
import numpy as np

from matplotlib import gridspec
import matplotlib.pyplot as plt
from pystella.model.stella import Stella

from pystella.rf import band

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))

_colors = ["blue", "cyan", "brown", 'darkseagreen', 'tomato', 'olive', 'orange',
           'skyblue', 'darkviolet']
# colors = {"B-V": "blue", 'B-V-I': "cyan", 'V-I': "brown"}
lntypes = dict(U="-", B="-", V="-", R="-", I="-",
               UVM2="-.", UVW1="-.", UVW2="-.",
               u="--", g="--", r="--", i="--", z="--")
markers = {u'D': u'diamond', 6: u'caretup', u's': u'square', u'x': u'x',
           5: u'caretright', u'^': u'triangle_up', u'd': u'thin_diamond', u'h': u'hexagon1',
           u'+': u'plus', u'*': u'star', u'o': u'circle', u'p': u'pentagon', u'3': u'tri_left',
           u'H': u'hexagon2', u'v': u'triangle_down', u'8': u'octagon', u'<': u'triangle_left'}
markers = markers.keys()


def plot_all(models_dic, bands, title='', is_time_points=True):
    colors = dict(U="blue", B="cyan", V="black", R="red", I="magenta",
                  J="blue", H="cyan", K="black",
                  UVM2="green", UVW1="red", UVW2="blue",
                  g="black", r="red", i="magenta", u="blue", z="magenta")
    # band_shift = dict(U=6.9, B=3.7, V=0, R=-2.4, I=-4.7,
    #                   UVM2=11.3, UVW1=10, UVW2=13.6,
    #                   u=3.5, g=2.5, r=-1.2, i=-3.7, z=-4.2)
    band_shift = dict((k, 0) for k, v in colors.items())  # no y-shift

    t_points = [0.2, 1, 2, 3, 4, 5, 10, 20, 40, 80, 150]
    xlim = [-10, 210]
    ylim = [-8, -22]

    def lbl(b):
        shift = band_shift[b]
        l = b
        if shift > 0:
            l += '+' + str(shift)
        elif shift < 0:
            l += '-' + str(abs(shift))
        return l

    # setup figure
    plt.matplotlib.rcParams.update({'font.size': 14})
    fig = plt.figure(num=None, figsize=(7, 7), dpi=100, facecolor='w', edgecolor='k')
    gs1 = gridspec.GridSpec(1, 1)
    gs1.update(wspace=0.3, hspace=0.3, left=0.1, right=0.95)
    ax = fig.add_subplot(gs1[0, 0])
    lw = 1.

    mi = 0
    ib = 0
    for mname, mdic in models_dic.iteritems():
        mi += 1
        for bname in bands:
            ib += 1
            x = mdic['time']
            y = mdic[bname]
            bcolor = colors[bname]
            ax.plot(x, y, marker=markers[mi % (len(markers) - 1)], label='%s  %s' % (lbl(bname), mname),
                    markersize=4, color=bcolor, ls="-", linewidth=lw)
            if is_time_points:
                integers = [np.abs(x - t).argmin() for t in t_points]  # set time points
                for (X, Y) in zip(x[integers], y[integers]):
                    ax.annotate('{:.0f}'.format(X), xy=(X, Y), xytext=(-10, 20), ha='right',
                                textcoords='offset points', color=bcolor,
                                arrowprops=dict(arrowstyle='->', shrinkA=0))

    ax.legend(prop={'size': 8})
    ax.invert_yaxis()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.ylabel('Magnitude')
    plt.xlabel('Time [days]')
    # ax.set_title(bset)

    plt.title(title)
    plt.grid()
    plt.show()


def cache_load(fname, colnames, skiprows=1):
    dtype = np.dtype({'names': colnames, 'formats': [np.float64] * len(colnames)})
    block = np.loadtxt(fname, skiprows=skiprows, dtype=dtype)
    return block


def usage():
    bands = band.band_get_names().keys()
    print "Usage:"
    print "  compare_ubv.py [params]"
    print "  -b <set_bands>: delimiter '_'. Default: u-g-r-i.\n" \
          "     Available: " + '-'.join(sorted(bands))
    print "  -i <model name>.  Example: cat_R450_M15_Ni007_E7"
    print "  -p <model path(directory)>, default: ./"
    print "  -t  plot time points"
    print "  -w  write magnitudes to file, default 'False'"
    print "  -h  print usage"


def compare_ABzVSugri(mname, path, bands=['u', 'g', 'r', 'i'], is_plot_time_points=False):
    names = ['time'] + bands

    dic_results = {}
    fname = os.path.join(path, mname + '.ABz')
    dic_results['ABz'] = cache_load(fname, names, skiprows=4)

    fname = os.path.join(path, mname + '.ubv')
    dic_results['ubv'] = cache_load(fname, names)

    plot_all(dic_results, bands, title=mname, is_time_points=is_plot_time_points)


def compare_ttVSubv(mname, path, bands=['U', 'B', 'V', 'R', 'I'], t_cut=1., is_plot_time_points=False):
    dic_results = {}

    model = Stella(mname, path=path)
    tt = model.read_tt_data()

    # time cut  days
    mags = tt[tt['time'] > t_cut]
    n1 = mags.dtype.names

    def f(x):
        if len(x) == 2 and x.startswith('M'):
            return x.replace('M', '')
        else:
            return x

    n2 = map(lambda x: f(x), n1)
    mags.dtype.names = n2

    dic_results['tt'] = mags

    # ubv
    serial_spec = model.read_series_spectrum(t_diff=1.05)
    mags = serial_spec.old_compute_mags(bands)
    dic_results['ubv'] = mags

    plot_all(dic_results, bands, title=mname, is_time_points=is_plot_time_points)


def main():
    path = '/home/bakl/Sn/Release/seb_git/res/tt/tanaka'
    mname = 'cat_R500_M15_Ni008_E40'

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hp:i:b:")
    except getopt.GetoptError as err:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-i':
            mname = str(arg)
            break
        if opt == '-b':
            bands = str(arg).split('-')
            for b in bands:
                if not band.band_is_exist(b):
                    print 'No such band: ' + b
                    sys.exit(2)
            continue
        if opt == '-t':
            is_plot_time_points = True
            continue
        if opt == '-p':
            path = str(arg)
            if not (os.path.isdir(path) and os.path.exists(path)):
                print "No such directory: " + path
                sys.exit(2)
            continue
        elif opt == '-h':
            usage()
            sys.exit(2)

    compare_ABzVSugri(mname, path)

    path = '/home/bakl/Sn/Release/seb_git/res/tt'
    compare_ttVSubv(mname, path)


if __name__ == '__main__':
    main()
