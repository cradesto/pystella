#!/usr/bin/env python3
# #!/usr/bin/python3

import getopt
import numpy as np
import os
import sys
from os.path import isfile, join, dirname
from scipy import interpolate

import matplotlib.pyplot as plt
from matplotlib import gridspec

import pystella.rf.rad_func as rf
from pystella.model.stella import Stella
from pystella.rf import band

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))

_colors = ["blue", "cyan", "brown", 'darkseagreen', 'tomato', 'olive', 'orange',
           'skyblue', 'darkviolet']
markers = {u'D': u'diamond', 6: u'caretup', u's': u'square', u'x': u'x',
           5: u'caretright', u'^': u'triangle_up', u'd': u'thin_diamond', u'h': u'hexagon1',
           u'+': u'plus', u'*': u'star', u'o': u'circle', u'p': u'pentagon', u'3': u'tri_left',
           u'H': u'hexagon2', u'v': u'triangle_down', u'8': u'octagon', u'<': u'triangle_left'}
markers = markers.keys()


def plot_dmdt(models_dic, bands, is_time_points=True):
    t_points = [0.2, 1, 2, 3, 4, 5, 10, 20, 40, 80, 150]
    xlim = [0.1, 11]
    ylim = [-11, -21.]

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
            x = abs(1. / mdic['d'][bname])
            # x = np.append(x, x[len(x)-1])
            y = mdic['m'][bname]
            z = mdic['m']['time']
            bcolor = _colors[ib % (len(_colors) - 1)]
            ax.plot(x, y, marker=markers[mi % (len(markers) - 1)], label='%s dmdt %s' % (bname, mname),
                    markersize=4, color=bcolor, ls="", linewidth=lw)
            if is_time_points:
                integers = [np.abs(z - t).argmin() for t in t_points]  # set time points
                for (X, Y, Z) in zip(x[integers], y[integers], z[integers]):
                    ax.annotate('{:.0f}'.format(Z), xy=(X, Y), xytext=(-10, 20), ha='right',
                                textcoords='offset points', color=bcolor,
                                arrowprops=dict(arrowstyle='->', shrinkA=0))

            if min(y) < -18.:
                print(" Max: %s in %s band " % (mname, bname))

    ax.legend(prop={'size': 8})
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xscale('log')
    ax.set_xlabel(r'Rising timescale (day/mag)')
    ax.set_ylabel(r'Absolute magnitude')
    # ax.set_title(bset)

    #     plt.title('; '.join(set_bands) + ' filter response')
    plt.grid()
    plt.show()


def compute_mag(name, path, bands, z=0., distance=10., t_cut=0., t_up=400., tdiff=0.):
    """
        Compute magnitude in bands for the 'name' model.
    :param name: the name of a model and data files
    :param path: the directory with data-files
    :param bands: photometric bands
    :param z: redshift, default 0
    :param distance: distance to star in parsec, default 10 pc
    :param t_cut:
    :param t_up:
    :param tdiff:
    :return: dictionary with keys = bands, value = star's magnitudes
    """
    model = Stella(name, path=path)

    if not model.is_ph:
        print("No data for: " + str(model))
        return None

    # serial_spec = model.read_serial_spectrum(t_diff=0.)
    serial_spec = model.read_series_spectrum(t_diff=1.05)
    mags = serial_spec.mags_bands(bands, z=z, d=rf.pc_to_cm(distance))

    # t_cut
    time = mags['time']
    cut = (t_cut < time) & (time < t_up)
    mags['time'] = time[cut]
    for n in bands:
        mags[n] = mags[n][cut]

    time = mags['time']
    if tdiff > 0.:
        ts = np.arange(np.min(time), np.max(time), tdiff)
        for n in bands:
            tck = interpolate.splrep(time, mags[n])
            mags[n] = interpolate.splev(ts, tck)
        mags['time'] = ts

    return mags


def compute_dmdt(mags, bands, is_spline=True, s=0.):
    dic_dmdt = dict((k, None) for k in bands)
    t = mags['time']
    for n in bands:
        if is_spline:
            w = np.ones(len(t))
            tck = interpolate.splrep(t, mags[n], w=w, s=s)
            dmdt = interpolate.splev(t, tck, der=1)
        else:
            dt = np.diff(t)
            dmdt = np.diff(mags[n]) / dt
            dmdt = np.append(dmdt, dmdt[-1])
        dic_dmdt[n] = dmdt

    return dic_dmdt


def usage():
    print("Usage:")
    print("  dmdt.py [params]")
    print("  -b <set_bands>: delimiter '-'. Default: 'UVW1-U'.\n")
    print("  -i <model name>.  Example: cat_R450_M15_Ni007_E7")
    print("  -p <model path(directory)>, default: ./")
    print("  -t  plot time points")
    print("  -w  write magnitudes to file, default 'False'")
    print("  -h  print usage")

    band.print_bands()


def main(name='', path='./'):
    model_ext = '.ph'
    is_time_points = False
    z = 0
    distance = 10.  # pc

    band.Band.load_settings()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hwtp:i:b:")
    except getopt.GetoptError as err:
        print(str(err))  # will print something like "option -a not recognized"
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

    bands = ['UVW1', 'U']

    for opt, arg in opts:
        if opt == '-b':
            bands = str(arg).split('-')
            for b in bands:
                if not band.is_exist(b):
                    print('No such band: ' + b)
                    sys.exit(2)
        elif opt == '-p':
            path = os.path.expanduser(str(arg))
            if not (os.path.isdir(path) and os.path.exists(path)):
                print("No such directory: " + path)
                sys.exit(2)
        elif opt == '-t':
            is_time_points = True
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

    if len(names) > 0:
        dic_results = {}  # dict((k, None) for k in names)
        i = 0
        for name in names:
            i += 1
            mags = compute_mag(name, path, bands, z=z, distance=distance, t_cut=0.1, t_up=5.,
                               tdiff=0.5)
            dmdt = compute_dmdt(mags, bands, is_spline=True, s=0.)
            dic_results[name] = dict(m=mags, d=dmdt)
            print("Finish: %s [%d/%d]" % (name, i, len(names)))
        plot_dmdt(dic_results, bands, is_time_points=is_time_points)

    else:
        print("There are no models in the directory: %s with extension: %s " % (path, model_ext))


if __name__ == '__main__':
    main()
    # main(name="cat_R1000_M15_Ni007_E15", path="/home/bakl/Sn/Release/seb_git/res/tt",
    #      is_force=False, is_save=True, is_plot_time_points=True)
