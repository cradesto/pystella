#!/usr/bin/python
# -*- coding: utf-8 -*-
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


def plot_zeta(res_dic, set_bands, title=''):
    colors = {"B-V": "blue", 'B-V-I': "cyan", 'V-I': "black"}
    lntypes = {"B-V": "-", 'B-V-I': "-.", 'V-I': "--"}

    xlim = [0, 25000]
    ylim = [0, 5]
    is_auto_lim = False
    if is_auto_lim:
        xlim = [0, 0]
        ylim = [0, 0]

    for n in set_bands:
        x = res_dic[n]['Tcol']
        y = res_dic[n]['zeta']
        plt.plot(x, y, label=n, color=colors[n], ls=lntypes[n], linewidth=2.0)
        if is_auto_lim:
            xlim[0], xlim[1] = min(np.append(x, xlim[0])), max(np.append(x, xlim[1]))
            ylim[0], ylim[1] = min(np.append(x, ylim[0])), max(np.append(y, ylim[1]))

    # plt.xlim(xlim)
    # ylim = np.add(ylim, [1, -1])
    # plt.ylim(ylim)
    plt.legend()
    plt.ylabel('Zeta')
    plt.xlabel('T_color')
    if title != '':
        plt.title(title)
    else:
        plt.title('; '.join(set_bands) + ' filter response')
    plt.grid()
    plt.show()


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
        tcolor, w = optim.fmin(epsilon, x0=np.array([1.e4, 1]), args=(freq, mag, bands, radius, dist))
        temp.append(tcolor)
        zeta_radius.append(w)
    return temp, zeta_radius


def compute_tcolor(name, path, bands, is_show_info=False, is_save=False):
    model = Stella(name, path=path)

    if not model.is_spec_data:
        print "No ph-data for: " + str(model)
        return None

    if not model.is_tt_data:
        print "No tt-data for: " + str(model)
        return None

    print "\nRun: %s for %s " % (name, join(bands))

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

    # time cut
    t_cut = 2.  # days
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

    print "\nFinish: %s for %d times" % (name, len(Tcolors))

    if is_save:
        fname = os.path.join(path, "%s.%s_%s.dat" % (name, "TcolorZeta", '-'.join(bands),))
        print "\nSave Tcolor & Zeta for %s in %s." % (bands, fname)
        res_save(res, fname=fname)

    return res


def res_save(narr, fname):
    names = narr.dtype.names
    with open(fname, 'wb') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['{:^8s}'.format(x) for x in names])
        for i, (row) in enumerate(zip(*[narr[k] for k in names])):
            writer.writerow(['{:8.3f}'.format(x) for x in row])


def usage():
    bands = band.band_get_names().keys()
    print "Usage:"
    print "  bbfit.py [params]"
    print "  -b <bands>: string, default: U-B-V-R-I, for example U-B-V-R-I-u-g-i-r-z-UVW1-UVW2.\n" \
          "     Available: " + '-'.join(sorted(bands))
    print "  -i <model name>.  Example: cat_R450_M15_Ni007_E7"
    print "  -d <model directory>, default: ./"
    print "  -e <model extension> is used to define model name, default: tt "
    print "  -s  silence mode: no info, no plot"
    print "  -w  write magnitudes to file, default 'False'"
    print "  -h  print usage"


def main(name='', model_ext='.tt'):
    is_silence = False
    is_save_mags = False
    path = ''

    try:
        opts, args = getopt.getopt(sys.argv[1:], "swhd:e:i:b:")
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

    set_bands = ['B', 'V']

    for opt, arg in opts:
        if opt == '-e':
            model_ext = '.' + arg
            continue
        if opt == '-b':
            set_bands = str(arg).split('-')
            for bset in set_bands:
                if not band.band_is_exist(bset):
                    print 'No such band: ' + bset
                    sys.exit(2)
            continue
        if opt == '-s':
            is_silence = True
            continue
        if opt == '-w':
            is_save_mags = True
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

    set_bands = ['B-V-I', 'B-V', 'V-I']

    if name != '':
        dic = dict((k, None) for k in set_bands)
        for bset in set_bands:
            dic[bset] = compute_tcolor(name, path, bset.split('-')
                                       , is_show_info=not is_silence, is_save=is_save_mags)
        plot_zeta(dic, set_bands)

    elif path != '':  # run for whole path
        names = []
        files = [f for f in os.listdir(path) if isfile(join(path, f)) and f.endswith(model_ext)]
        for f in files:
            names.append(os.path.splitext(f)[0])
        if len(names) > 0:
            for name in names:
                compute_tcolor(name, path, set_bands, is_show_info=not is_silence, is_save=is_save_mags)
        else:
            print "There are no models in the directory: %s with extension: %s " % (path, model_ext)


if __name__ == '__main__':
    main()
    # main(name="cat_R1000_M15_Ni007_E15")
