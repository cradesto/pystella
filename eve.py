#!/usr/bin/python
# -*- coding: utf-8 -*-


import getopt
import os
import sys
from os.path import dirname
from os.path import isfile, join

import matplotlib.pyplot as plt
from matplotlib import gridspec

import pystella.model.sn_eve  as sneve


__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))


def plot_all(models_vels, models_dic, bands, call=None, xlim=None, ylim=None,
             is_time_points=False, title='', fsave=None):
    colors = band.bands_colors()
    band_shift = dict((k, 0) for k, v in colors.items())  # no y-shift
    is_vel = models_vels is not None

    # setup figure
    plt.matplotlib.rcParams.update({'font.size': 14})
    fig = plt.figure(num=None, figsize=(7, 11), dpi=100, facecolor='w', edgecolor='k')

    if is_vel:
        gs1 = gridspec.GridSpec(4, 1)
        axUbv = fig.add_subplot(gs1[:-1, 0])
        axVel = fig.add_subplot(gs1[3, 0])
    else:
        gs1 = gridspec.GridSpec(1, 1)
        axUbv = fig.add_subplot(gs1[0, 0])
        axVel = None
    gs1.update(wspace=0.3, hspace=0.3, left=0.1, right=0.95)

    # plot the light curves
    lc.plot_ubv_models(axUbv, models_dic, bands, band_shift=band_shift, xlim=xlim, ylim=ylim,
                       is_time_points=is_time_points)

    # plot callback
    if call is not None:
        call.plot(axUbv, {'ax2': axVel})

    # finish plot
    axUbv.set_ylabel('Magnitude')
    axUbv.set_xlabel('Time [days]')

    axUbv.legend(prop={'size': 8}, loc=4)
    # ax.set_title(bset)
    if title:
        axUbv.set_title(title)

    # plot velocities
    if is_vel:
        vel.plot_vels_models(axVel, models_vels, xlim=axUbv.get_xlim())
        # vel.plot_vels_sn87a(axVel, z=1.49)
        axVel.legend(prop={'size': 8}, loc=4)

    plt.grid()

    plt.show()

    if fsave is not None:
        print "Save plot to %s " % fsave
        fig.savefig(fsave, bbox_inches='tight')
        # plt.savefig(fsave, format='pdf')


def usage():
    elements = sneve.eve_elements
    print "Usage:"
    print "  eve.py [params]"
    print "  -b <elements>: string, default: U-B-V-R-I, for example U-B-V-R-I-u-g-i-r-z-UVW1-UVW2.\n" \
          "     Available: " + '-'.join(elements)
    print "  -i <model name>.  Example: cat_R450_M15_Ni007"
    print "  -p <model directory>, default: ./"
    print "  -q  quiet mode: no info, no plot"
    print "  -s  save plot to pdf-file."
    print "  -w  write magnitudes to file, default 'False'"
    print "  -h  print usage"


def main(name='', model_ext='.ph'):
    is_quiet = False
    is_save_mags = False
    is_save_plot = False
    path = ''

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hqswp:i:b:")
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

    elements = sneve.eve_elements
    # bands = ['U', 'B', 'V', 'R', "I", 'UVM2', "UVW1", "UVW2", 'g', "r", "i"]

    for opt, arg in opts:
        if opt == '-e':
            e = float(arg)
            is_extinction = True
            continue
        if opt == '-b':
            elements = str(arg).split('-')
            for e in elements:
                if not sneve.eve_elements(e):
                    print 'No such element: ' + e
                    sys.exit(2)
            continue
        if opt == '-q':
            is_quiet = True
            continue
        if opt == '-s':
            is_save_plot = True
            continue
        if opt == '-w':
            is_save_mags = True
            continue
        if opt == '-d':
            distance = float(arg)
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

    print "Plot chemical : %s %s" % (path, name)

    eve = sneve.StellaEve(name, path=path)
    if eve.is_rho_data:
        eve.rho_load()
        eve.plot_chem(elements=elements)


if __name__ == '__main__':
    main()
