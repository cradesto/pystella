#!/usr/bin/python
# -*- coding: utf-8 -*-
__author__ = 'bakl'

import os
import sys, getopt
from os.path import dirname
import csv
import numpy as np
import matplotlib.pyplot as plt

from pystella.model.stella import Stella
from pystella.rf import band

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))


def plot_bands(dict_mags, bands):
    plt.title(''.join(bands) + ' filter response')

    colors = dict(U="blue", B="cyan", V="black", R="red", I="magenta", UVM2="green", UVW1="red", UVW2="blue",
                  g="black", r="red", i="magenta")
    lntypes = dict(U="-", B="-", V="-", R="-", I="-", UVM2="-.", UVW1="-.", UVW2="-.", g="--", r="--", i="--")
    band_shift = dict(U=6.9, B=3.7, V=0, R=-2.4, I=-4.7, UVM2=11.3, UVW1=10, UVW2=13.6, g=2.5, r=-1.2, i=-3.7)
    # band_shift = dict((k, 0) for k, v in band_shift.items())  # no y-shift

    dist = 24e6  # distance to SN 2013ab
    dm = 5 * np.log10(dist) - 5  # distance module
    # dm = 0
    lims = (-12, -23)
    lims = [39, 8]
    is_auto_lim = False

    def lbl(b):
        shift = band_shift[b]
        l = b
        if shift > 0:
            l += '+' + str(shift)
        elif shift < 0:
            l += '-' + str(abs(shift))
        return l

    x = dict_mags['time']
    for n in bands:
        y = dict_mags[n]
        y += dm + band_shift[n]
        plt.plot(x, y, label=lbl(n), color=colors[n], ls=lntypes[n], linewidth=2.0)
        if is_auto_lim:
            if lims[0] < max(y):
                lims[0] = max(y)
            if lims[1] > min(y):
                lims[1] = min(y)

    plt.gca().invert_yaxis()
    plt.xlim([-10, 200])
    plt.ylim(lims)
    plt.legend()
    plt.ylabel('Magnitude')
    plt.xlabel('Time [days]')
    plt.grid()
    plt.show()


def compute_mag(name, path, bands):
    model = Stella(name, path=path)
    model.show_info()

    if not model.is_any_data:
        print "No data for: " + str(model)
        return None

    serial_spec = model.read_serial_spectrum(t_diff=1.05)

    dict_mags = dict((k, None) for k in bands)
    dict_mags['time'] = serial_spec.times
    for n in bands:
        b = band.band_by_name(n)
        mags = serial_spec.flux_to_mags(b)
        if mags is not None:
            dict_mags[n] = mags

    return dict_mags


def mags_save(dictionary, bands, fname):
    with open(fname, 'wb') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['{:^8s}'.format(x) for x in ['time'] + bands])
        for i, (row) in enumerate(zip(*[dictionary[k] for k in 'time'.split() + bands])):
            # row = row[-1:] + row[:-1]  # make time first column
            writer.writerow(['{:8.3f}'.format(x) for x in row])
            # writer.writerow(['{:3.4e}'.format(x) for x in row])


def usage():
    print "Usage:"
    print "  sn.py [params]"
    print "  -b <bands>: string like U-B-V-R-I-g-r-i-UVM2-UVW1-UVW2, default: U-B-V-R-I"
    print "  -i <model name>.  Ex: cat_R1000_M15_Ni007_E15"
    print "  -d <model directory>, default: ./"
    print "  -h: print usage"


def main(name=''):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:i:b:")
    except getopt.GetoptError as err:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(opts) == 0:
        usage()
        sys.exit(2)

    # set prf from command line
    if not name and not ('-i' in opts):
        print 'You should specify model name'  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    bands = ['U', 'B', 'V', 'R', "I"]
    # bands = ['U', 'B', 'V', 'R', "I", 'UVM2', "UVW1", "UVW2", 'g', "r", "i"]

    path = ROOT_DIRECTORY
    for opt, arg in opts:
        if opt == '-b':
            bands = str(arg).split('-')
            for b in bands:
                if not band.band_is_exist(b):
                    print 'No such band: ' + b
                    sys.exit(2)
        if opt == '-d':
            path = str(arg)
            if not (os.path.isdir(path) and os.path.exists(path)):
                print "No such directory: " + path
                sys.exit(2)
        elif opt == '-h':
            usage()
            sys.exit(2)

    dict_mags = compute_mag(name, path, bands)

    if dict_mags is not None:
        mags_save(dict_mags, bands, name + '_' + ''.join(bands) + '.txt')
        plot_bands(dict_mags, bands)


if __name__ == '__main__':
    main(name="cat_R1000_M15_Ni007_E15")
