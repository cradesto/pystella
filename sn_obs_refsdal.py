#!/usr/bin/python
# -*- coding: utf-8 -*-


import getopt
import os
import sys
from os.path import dirname

import matplotlib.pyplot as plt
from matplotlib import gridspec

import plugin.plot_snrefsdal as sn_obs
import pystella.util.callback as cb
from pystella import velocity as vel
from pystella.rf import band
from pystella.rf import extinction
from pystella.rf import light_curve as lc
from pystella.util.phys_var import cosmology_D_by_z

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))


def plot_S4(models_dic, bands, call=None, xlim=None, ylim=None, title='', fsave=None):
    set_images = ['S1', 'S2', 'S3', 'S4']
    colors = band.bands_colors()
    band_shift = dict((k, 0) for k, v in colors.items())  # no y-shift

    # setup figure
    plt.matplotlib.rcParams.update({'font.size': 14})
    fig = plt.figure(num=len(set_images), figsize=(9, 9), dpi=100, facecolor='w', edgecolor='k')
    gs1 = gridspec.GridSpec(len(set_images) / 2 + len(set_images) % 2, 2)
    gs1.update(wspace=0., hspace=0., left=0.1, right=0.9)

    ax_cache = {}

    # create the grid of figures
    ib = 0
    for im in set_images:
        ib += 1
        icol = (ib - 1) % 2
        irow = (ib - 1) / 2
        ax = fig.add_subplot(gs1[irow, icol])
        ax_cache[ib] = ax

        # set axis
        if icol > 0:
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position("right")
        ax.set_ylabel('Obs. Magnitude (AB)')

        if irow == 1:
            ax.set_xlabel('Time [days]')

        if ib > 1:
            xlim = ax_cache[1].get_xlim()
            ylim = ax_cache[1].get_ylim()

        lc.plot_ubv_models(ax, {im: models_dic[im]}, bands, band_shift=band_shift, xlim=xlim, ylim=ylim)
        # plot callback
        if call is not None:
            call.plot(ax, {'image': im})

        ax.text(15, 23.5, im, bbox={'facecolor': 'blue', 'alpha': 0.2, 'pad': 10})

    # plt.legend(prop={'size': 8}, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    ax_cache[2].legend(prop={'size': 8}, loc='upper center', bbox_to_anchor=(0.02, 1.2),
                       ncol=4, fancybox=True, shadow=True)

    # plt.grid()
    plt.show()

    if fsave is not None:
        print "Save plot to %s " % fsave
        fig.savefig(fsave, bbox_inches='tight')


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
    bands = band.band_get_names().keys()
    print "Usage:"
    print "  ubv.py [params]"
    print "  -b <bands>: string, default: U-B-V-R-I, for example U-B-V-R-I-u-g-i-r-z-UVW1-UVW2.\n" \
          "     Available: " + '-'.join(sorted(bands))
    print "  -i <model name>.  Example: cat_R450_M15_Ni007_E7"
    print "  -p <model directory>, default: ./"
    print "  -e <extinction, E(B-V)> is used to define A_nu, default: 0 "
    print "  -c <callback> [plot_tolstov, plot_snrefsdal]. You can add parameters in format func:params"
    print "  -d <distance> [pc].  Default: 10 pc"
    print "  -m <magnification>.  Default: None, used for grav lens"
    print "  -z <redshift>.  Default: 0"
    print "  -q  quiet mode: no info, no plot"
    print "  -t  plot time points"
    print "  -s  save plot to pdf-file."
    print "  -v  plot model velocities."
    print "  -w  write magnitudes to file, default 'False'"
    print "  -h  print usage"


def lc_wrapper(param):
    a = param.split(':')
    func = a.pop(0)
    c = cb.CallBack(func, path=cb.plugin_path, args=a, load=1)
    return c


def main(name='', model_ext='.ph'):
    is_quiet = False
    is_save = False
    is_save_plot = False
    is_plot_time_points = False
    is_extinction = False
    is_vel = False
    path = ''
    z = 0
    e = 0.
    magnification = 1.
    distance = 10.  # pc
    callback = None

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hqswtc:d:p:e:i:b:m:vz:")
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
        if opt == '-b':
            bands = str(arg).split('-')
            for b in bands:
                if not band.band_is_exist(b):
                    print 'No such band: ' + b
                    sys.exit(2)
            continue
        if opt == '-c':
            callback = lc_wrapper(str(arg))
            continue
        if opt == '-d':
            distance = float(arg)
            continue
        if opt == '-e':
            e = float(arg)
            is_extinction = True
            continue
        if opt == '-h':
            usage()
            sys.exit(2)
        if opt == '-m':
            magnification = float(arg)
            continue
        if opt == '-p':
            path = os.path.expanduser(str(arg))
            if not (os.path.isdir(path) and os.path.exists(path)):
                print "No such directory: " + path
                sys.exit(2)
            continue
        if opt == '-q':
            is_quiet = True
            continue
        if opt == '-s':
            is_save_plot = True
            continue
        if opt == '-t':
            is_plot_time_points = True
            continue
        if opt == '-v':
            is_vel = True
            continue
        if opt == '-w':
            is_save = True
            continue
        if opt == '-z':
            z = float(arg)
            continue

    print "Plot magnitudes on z=%f at distance=%e [cosmology D(z)=%s Mpc]" % (z, distance, cosmology_D_by_z(z))

    if is_vel:
        run_ubv_vel(name, path, bands, e, z, distance, magnification, callback,
                    is_show_info=not is_quiet, is_save=is_save)
    else:
        run_S4(name, path, bands, e, z, distance, magnification, callback
               , is_save=is_save)


def run_S4(name, path, bands, e, z, distance, magnification, callback, is_save):
    if e > 0:
        if z > 1:
            ext = extinction.extinction_law_z(ebv=e, bands=bands, z=z)
        else:
            ext = extinction.extinction_law(ebv=e, bands=bands)
    else:
        ext = None

    sn_images = sn_obs.coef_magnification('oguri')
    if len(sn_images) > 0:
        models_mags = {}  # dict((k, None) for k in names)
        i = 0
        for im, mgf in sn_images.items():
            if im == 'SX':  # pass image
                continue
            i += 1
            mgf *= magnification
            mags = lc.compute_mag(name, path, bands, ext=ext, z=z, distance=distance, magnification=mgf,
                                  is_show_info=False, is_save=is_save)
            models_mags[im] = mags
            print "Finish image: %s [%d/%d]" % (im, i, len(sn_images))

        mgf = sn_images['S1'] * magnification
        if callback is not None:
            t = "ts=%s z=%4.2f D=%6.2e mu=%3.1f ebv=%4.2f" % (callback.arg_totext(0), z, distance, mgf, e)
        else:
            t = "z=%4.2f D=%6.2e mu=%3.1f ebv=%4.2f" % (z, distance, mgf, e)

        fsave = None
        if is_save:
            fsave = "ubv_%s" % name

            if ext is not None and ext > 0:
                fsave = "%s_e0%2d" % (fsave, int(ext * 100))  # bad formula for name

            d = os.path.expanduser('~/')
            # d = '/home/bakl/Sn/my/conf/2016/snrefsdal/img'
            fsave = os.path.join(d, fsave) + '.pdf'

        plot_S4(models_mags, bands, call=callback, title=t, fsave=fsave)
    else:
        print "There are no sn images"


def run_ubv_vel(name, path, bands, ext, z, distance, magnification, callback,
                is_vel, is_show_info, is_save):
    models_mags = {}
    models_vels = {}

    mags = lc.compute_mag(name, path, bands, ext=ext, z=z, distance=distance, magnification=magnification,
                          is_show_info=is_show_info, is_save=is_save)
    models_mags[name] = mags

    if is_show_info:
        lc.plot_bands(mags, bands, title=name)

    if is_vel:
        vels = vel.compute_vel(name, path, z=z)
        if vels is None:
            sys.exit("No data for: %s in %s" % (name, path))
        models_vels[name] = vels
        print "Finish velocity: %s " % name
    else:
        models_vels = None
        print "Finish mags: %s " % name

    if callback is not None:
        t = "ts=%s z=%4.2f D=%6.2e mu=%3.1f ebv=%4.2f" % (callback.arg_totext(0), z, distance, magnification, ext)
    else:
        t = "z=%4.2f D=%6.2e mu=%3.1f ebv=%4.2f" % (z, distance, magnification, ext)

    fsave = None
    if is_save:
        if is_vel:
            fsave = "ubv_vel_%s" % name
        else:
            fsave = "ubv_%s" % name

        if ext is not None and ext > 0:
            fsave = "%s_e0%2d" % (fsave, int(ext * 100))  # bad formula for name

        d = os.path.expanduser('~/')
        # d = '/home/bakl/Sn/my/conf/2016/snrefsdal/img'
        fsave = os.path.join(d, fsave) + '.pdf'

    plot_all(models_vels, models_mags, bands, call=callback, is_time_points=False, title=t, fsave=fsave)
    # else:
    #     plot_S4(models_mags, bands, call=callback, title=t, fsave=fsave)
    # plot_all(dic_results, bands,  xlim=(-10, 410), is_time_points=is_plot_time_points)
    # plot_all(dic_results, bands, xlim=(-10, 410), callback=callback, is_time_points=is_plot_time_points)
    # plot_all(dic_results, bands,  ylim=(40, 23),  is_time_points=is_plot_time_points)


if __name__ == '__main__':
    main()
    # main(name="cat_R1000_M15_Ni007_E15")
