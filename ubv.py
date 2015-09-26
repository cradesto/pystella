#!/usr/bin/python
# -*- coding: utf-8 -*-
from matplotlib import gridspec
from pystella.util import rf
import os
import sys
import getopt
from os.path import dirname
from os.path import isfile, join
import csv
import numpy as np

import matplotlib.pyplot as plt

from pystella.model.stella import Stella
from pystella.rf import band

__author__ = 'bakl'


ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))

markers = {u'D': u'diamond', 6: u'caretup', u's': u'square', u'x': u'x',
           5: u'caretright', u'^': u'triangle_up', u'd': u'thin_diamond', u'h': u'hexagon1',
           u'+': u'plus', u'*': u'star', u'o': u'circle', u'p': u'pentagon', u'3': u'tri_left',
           u'H': u'hexagon2', u'v': u'triangle_down', u'8': u'octagon', u'<': u'triangle_left'}
markers = markers.keys()


def plot_all(models_dic, bands, is_time_points=True):
    colors = dict(U="blue", B="cyan", V="black", R="red", I="magenta",
                  J="blue", H="cyan", K="black",
                  UVM2="green", UVW1="red", UVW2="blue",
                  g="black", r="red", i="magenta", u="blue", z="magenta")
    lntypes = dict(U="-", B="-", V="-", R="-", I="-",
                   UVM2="-.", UVW1="-.", UVW2="-.",
                   u="--", g="--", r="--", i="--", z="--")
    band_shift = dict(U=6.9, B=3.7, V=0, R=-2.4, I=-4.7,
                      UVM2=11.3, UVW1=10, UVW2=13.6,
                      u=3.5, g=2.5, r=-1.2, i=-3.7, z=-4.2)
    band_shift = dict((k, 0) for k, v in band_shift.items())  # no y-shift

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
            ax.plot(x, y,  marker=markers[mi % (len(markers) - 1)], label='%s  %s' % (bname, mname),
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

    #     plt.title('; '.join(set_bands) + ' filter response')
    plt.grid()
    plt.show()

def plot_bands(dict_mags, bands, title='', fname='', distance=10., is_time_points=True):
    plt.title(''.join(bands) + ' filter response')

    colors = dict(U="blue", B="cyan", V="black", R="red", I="magenta",
                  J="blue", H="cyan", K="black",
                  UVM2="green", UVW1="red", UVW2="blue",
                  g="black", r="red", i="magenta", u="blue", z="magenta")
    lntypes = dict(U="-", B="-", V="-", R="-", I="-",
                   UVM2="-.", UVW1="-.", UVW2="-.",
                   u="--", g="--", r="--", i="--", z="--")
    band_shift = dict(U=6.9, B=3.7, V=0, R=-2.4, I=-4.7,
                      UVM2=11.3, UVW1=10, UVW2=13.6,
                      u=3.5, g=2.5, r=-1.2, i=-3.7, z=-4.2)
    band_shift = dict((k, 0) for k, v in band_shift.items())  # no y-shift

    t_points = [2, 5, 10, 20, 40, 80, 150]

    dm = 5 * np.log10(distance) - 5  # distance module
    # dm = 0
    xlims = [-10, 200]
    ylims = [-12, -19]
    ylims += dm
    is_auto_lim = True
    if is_auto_lim:
        ylims = [0, 0]

    def lbl(b):
        shift = band_shift[b]
        l = b
        if shift > 0:
            l += '+' + str(shift)
        elif shift < 0:
            l += '-' + str(abs(shift))
        return l

    x = dict_mags['time']
    is_first = True
    for n in bands:
        y = dict_mags[n]
        y += dm + band_shift[n]
        plt.plot(x, y, label=lbl(n), color=colors[n], ls=lntypes[n], linewidth=2.0, marker='s')
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
    plt.gca().invert_yaxis()
    plt.xlim(xlims)
    plt.ylim(ylims)
    plt.legend()
    plt.ylabel('Magnitude')
    plt.xlabel('Time [days]')
    plt.title(title)
    plt.grid()
    if fname != '':
        plt.savefig("ubv_%s.png" % fname, format='png')
    plt.show()
    # plt.close()


def compute_mag(name, path, bands, z=0., distance=10., is_show_info=True, is_save=False):
    """
        Compute magnitude in bands for the 'name' model.
    :param name: the name of a model and data files
    :param path: the directory with data-files
    :param bands: photometric bands
    :param z: redshift, default 0
    :param distance: distance to star in parsec, default 10 pc
    :param is_show_info: flag to trun some information, default True
    :param is_save: flag to save result in file, default False
    :return: dictionary with keys = bands, value = star's magnitudes
    """
    model = Stella(name, path=path)
    if is_show_info:
        print ''
        model.show_info()

    if not model.is_spec_data:
        print "No data for: " + str(model)
        return None

    # serial_spec = model.read_serial_spectrum(t_diff=0.)
    serial_spec = model.read_serial_spectrum(t_diff=1.05)

    mags = dict((k, None) for k in bands)

    # z, distance = 0, 10.  # pc for Absolute magnitude
    for n in bands:
        b = band.band_by_name(n)
        mags[n] = serial_spec.flux_to_mags(b, z=z, dl=rf.pc_to_cm(distance))

    mags['time'] = serial_spec.times * (1. + z)

    if mags is not None:
        fname = os.path.join(path, name + '.ubv')
        if is_save:
            mags_save(mags, bands, fname)
            print "Magnitudes have been saved to " + fname

    if is_show_info:
        # print the time of maximum LC
        t = mags['time']
        idxes = [t > 2.]
        for n in bands:
            t_min = t[idxes][mags[n][idxes].argmin()]
            print "t_max(%s) = %f" % (n, t_min)

    return mags


def mags_save(dictionary, bands, fname):
    with open(fname, 'wb') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['{:^8s}'.format(x) for x in ['time'] + bands])
        for i, (row) in enumerate(zip(*[dictionary[k] for k in 'time'.split() + bands])):
            # row = row[-1:] + row[:-1]  # make time first column
            writer.writerow(['{:8.3f}'.format(x) for x in row])
            # writer.writerow(['{:3.4e}'.format(x) for x in row])


def usage():
    bands = band.band_get_names().keys()
    print "Usage:"
    print "  ubv.py [params]"
    print "  -b <bands>: string, default: U-B-V-R-I, for example U-B-V-R-I-u-g-i-r-z-UVW1-UVW2.\n" \
          "     Available: " + '-'.join(sorted(bands))
    print "  -i <model name>.  Example: cat_R450_M15_Ni007_E7"
    print "  -p <model directory>, default: ./"
    print "  -e <model extension> is used to define model name, default: tt "
    print "  -d <distance> [pc].  Default: 10 pc"
    print "  -z <redshift>.  Default: 0"
    print "  -s  silence mode: no info, no plot"
    print "  -t  plot time points"
    print "  -w  write magnitudes to file, default 'False'"
    print "  -h  print usage"



def main(name='', model_ext='.ph'):
    is_silence = False
    is_save_mags = False
    is_plot_time_points = False
    path = ''
    z = 0
    distance = 10.  # pc

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hswtd:p:e:i:b:z:")
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
                # if name == '':
                #     print 'Error: you should specify the name of model.'
                #     sys.exit(2)

    bands = ['U', 'B', 'V', 'R', "I"]
    # bands = ['U', 'B', 'V', 'R', "I", 'UVM2', "UVW1", "UVW2", 'g', "r", "i"]

    for opt, arg in opts:
        if opt == '-e':
            model_ext = '.' + arg
            continue
        if opt == '-b':
            bands = str(arg).split('-')
            for b in bands:
                if not band.band_is_exist(b):
                    print 'No such band: ' + b
                    sys.exit(2)
            continue
        if opt == '-s':
            is_silence = True
            continue
        if opt == '-w':
            is_save_mags = True
            continue
        if opt == '-t':
            is_plot_time_points = True
            continue
        if opt == '-z':
            z = float(arg)
            continue
        if opt == '-d':
            distance = float(arg)
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

    print "Plot magnitudes on z=%f at distance=%f" % (z, distance)

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
            mags = compute_mag(name, path, bands, z=z, distance=distance,
                               is_show_info=not is_silence, is_save=is_save_mags)
            dic_results[name] = mags
            print "Finish: %s [%d/%d]" % (name, i, len(names))

            if not is_silence:
                # z, distance = 0.145, 687.7e6  # pc for comparison with Maria
                plot_bands(mags, bands, title=name, fname='', is_time_points=is_plot_time_points)
        plot_all(dic_results, bands, is_time_points=is_plot_time_points)
    else:
        print "There are no models in the directory: %s with extension: %s " % (path, model_ext)


if __name__ == '__main__':
    main()
    # main(name="cat_R1000_M15_Ni007_E15")
