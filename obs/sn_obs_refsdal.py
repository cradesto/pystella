#!/usr/bin/python
# -*- coding: utf-8 -*-

import getopt
import os
import sys
from os.path import dirname

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np

from pystella.rf import band
from pystella.rf import extinction
from pystella.rf import light_curve_func as lc
from pystella import velocity as vel
import pystella.util.callback as cb
from pystella.util.phys_var import cosmology_D_by_z
import plugin.plot_snrefsdal as sn_obs

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))


def plot_SX(models_dic, bands, call=None, xlim=None, ylim=None, title='', fsave=None):
    glens_models = sn_obs.coef_glens().keys()
    colors = band.bands_colors()
    band_shift = dict((k, 0) for k, v in colors.items())  # no y-shift

    # setup figure
    plt.matplotlib.rcParams.update({'font.size': 14})
    fig = plt.figure(num=len(glens_models), figsize=(12, 4), dpi=100, facecolor='w', edgecolor='k')
    gs1 = gridspec.GridSpec(1, len(glens_models))
    gs1.update(wspace=0., hspace=0., left=0.1, right=0.9)

    ax_cache = {}

    # create the grid of figures
    ib = 0
    im = 'SX'
    irow = 0
    for model in glens_models:
        ib += 1
        icol = ib - 1
        ax = fig.add_subplot(gs1[irow, icol])
        ax_cache[ib] = ax

        # set axis
        if icol == 0:
            # ax.yaxis.tick_right()
            # ax.yaxis.set_label_position("right")
            ax.set_ylabel('Obs. Magnitude (AB)')
        else:
            ax.set_yticklabels([])

        # if irow == 1:
        ax.set_xlabel('Time [days]')

        if ib > 1:
            xlim = ax_cache[1].get_xlim()
            ylim = ax_cache[1].get_ylim()

        lc.plot_ubv_models(ax, {im: models_dic[im]}, bands, band_shift=band_shift, xlim=xlim, ylim=ylim)
        # plot callback
        if call is not None:
            call.plot(ax, {'glens': model, 'image': im})

        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(start, end, 100))
        ax.text(15, 24.7, '%s: %s' % (im, model), bbox={'facecolor': 'blue', 'alpha': 0.2, 'pad': 10})
        ax.grid()

    # plt.legend(prop={'size': 8}, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    ax_cache[2].legend(prop={'size': 8}, loc='upper center', bbox_to_anchor=(0.02, 1.2), ncol=4)

    # plt.title(title)
    plt.show()

    if fsave is not None:
        print "Save plot to %s " % fsave
        fig.savefig(fsave, bbox_inches='tight')


def plot_BV(models_dic, bands, glens, call=None, xlim=None, title='', fsave=None):
    set_images = ['S1', 'S2', 'S3', 'S4']
    colors = band.bands_colors()
    band_shift = dict((k, 0) for k, v in colors.items())  # no y-shift
    ylim = [-5, 5]
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
        ax.set_ylabel('Color Index')

        if irow == 1:
            ax.set_xlabel('Time [days]')

        if ib > 1:
            xlim = ax_cache[1].get_xlim()
            ylim = ax_cache[1].get_ylim()

        # lc.plot_bv_models(ax, models_dic, bands, band_shift=band_shift, xlim=xlim, ylim=ylim)

        # plot callback
        if call is not None:
            call.plot(ax, {'bv': True, 'glens': glens, 'image': im})

        ax.text(400, 3.5, '%s' % im, bbox={'facecolor': 'blue', 'alpha': 0.2, 'pad': 10})
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    # plt.legend(prop={'size': 8}, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    ax_cache[2].legend(prop={'size': 8}, loc='upper center', bbox_to_anchor=(0.02, 1.1),
                       ncol=4, fancybox=True, shadow=True)

    # ax_cache[2].text(15, 23.5, title, bbox={'facecolor': 'blue', 'alpha': 0.2, 'pad': 10})

    # plt.grid()
    # plt.title(title)
    fig.suptitle(title, fontsize=11)
    plt.show()

    if fsave is not None:
        print "Save plot to %s " % fsave
        fig.savefig(fsave, bbox_inches='tight')


def old_plot_S4(models_dic, bands, glens, call=None, xlim=None, ylim=None, title='', fsave=None):
    # set_images = ['S1', 'S2', 'S3', 'SX']
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
            call.plot(ax, {'glens': glens, 'image': im})

        ax.text(15, 23.5, '%s: %s' % (im, glens), bbox={'facecolor': 'blue', 'alpha': 0.2, 'pad': 10})

    # plt.legend(prop={'size': 8}, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, , fancybox=True, shadow=True
    ax_cache[2].legend(prop={'size': 8}, loc='upper center', bbox_to_anchor=(0.02, 1.13), ncol=4)

    # ax_cache[2].text(15, 23.5, title, bbox={'facecolor': 'blue', 'alpha': 0.2, 'pad': 10})

    # plt.grid()
    # plt.title(title)
    fig.suptitle(title, fontsize=11)
    plt.show()

    if fsave is not None:
        print "Save plot to %s " % fsave
        fig.savefig(fsave, bbox_inches='tight')


def plot_S4_curves(models_curves, bands, glens, call=None, xlim=None, ylim=None, title='', fsave=None):
    # set_images = ['S1', 'S2', 'S3', 'SX']
    # set_images = models_curves.keys()
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

        lc.plot_models_curves(ax, {im: models_curves[im]}, bands, band_shift=band_shift, xlim=xlim, ylim=ylim)
        # plot callback
        if call is not None:
            call.plot(ax, {'glens': glens, 'image': im})

        ax.text(15, 23.5, '%s: %s' % (im, glens), bbox={'facecolor': 'blue', 'alpha': 0.2, 'pad': 10})

    # plt.legend(prop={'size': 8}, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, , fancybox=True, shadow=True
    ax_cache[2].legend(prop={'size': 8}, loc='upper center', bbox_to_anchor=(0.02, 1.13), ncol=4)

    # ax_cache[2].text(15, 23.5, title, bbox={'facecolor': 'blue', 'alpha': 0.2, 'pad': 10})

    # plt.grid()
    # plt.title(title)
    fig.suptitle(title, fontsize=11)
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
    lc_min = lc.plot_ubv_models(axUbv, models_dic, bands, band_shift=band_shift, xlim=xlim, ylim=ylim,
                                is_time_points=is_time_points)

    # show  times of spectral observations
    dt = 1.
    z = 1.49
    ts = np.array([-47, 16])  # see spectral dates in Kelly, 1512.09093
    ts *= 1. + z
    for bname, v in lc_min.items():
        if bname == 'F160W':
            for t in ts:
                axVel.axvspan(v[0] + t - dt, v[0] + t + dt, facecolor='g', alpha=0.5)
                print "Spectral obs: t=%8.1f, tmax=%8.1f" % (v[0] + t, v[0],)

    # plot callback
    if call is not None:
        call.plot(axUbv, {'ax2': axVel})

    # finish plot
    axUbv.set_ylabel('Magnitude')
    # axUbv.set_xlabel('Time [days]')

    axUbv.legend(prop={'size': 8}, loc=2, ncol=4)
    # ax.set_title(bset)
    if title:
        axUbv.set_title(title)

    # plot velocities
    if is_vel:
        vel.plot_vels_models(axVel, models_vels, xlim=axUbv.get_xlim())
        # vel.plot_vels_sn87a(axVel, z=1.49)
        axVel.legend(prop={'size': 8}, loc=1)
        # for xlim in zip(x-xerr, x+xerr):
        #     axVel.axvspan(xlim[0], xlim[1], facecolor='g', alpha=0.5)

    plt.grid()

    plt.show()

    if fsave is not None:
        print "Save plot to %s " % fsave
        fig.savefig(fsave, bbox_inches='tight')
        # plt.savefig(fsave, format='pdf')


def run_BV(name, path, bands, e, z, distance, magnification, callback, xlim, is_save):
    if e > 0:
        if z > 1:
            ext = extinction.extinction_law_z(ebv=e, bands=bands, z=z)
        else:
            ext = extinction.extinction_law(ebv=e, bands=bands)
    else:
        ext = None

    glens = callback.get_arg(2)
    if glens is None:
        glens = sn_obs.grav_lens_def

    mags = lc.compute_mag(name, path, bands, ext=ext, z=z, distance=distance, magnification=magnification,
                          is_show_info=False, is_save=is_save)
    models_mags = {name: mags}

    if callback is not None:
        t = "ts=%s z=%4.2f D=%6.2e mu=%3.1f ebv=%4.2f" % (callback.arg_totext(0), z, distance, magnification, e)
    else:
        t = "z=%4.2f D=%6.2e mu=%3.1f ebv=%4.2f" % (z, distance, magnification, e)

    fsave = None
    if is_save:
        fsave = "bv_%s" % name

        if ext is not None and ext > 0:
            fsave = "%s_e0%2d" % (fsave, int(e * 100))  # bad formula for name

        d = os.path.expanduser('~/')
        # d = '/home/bakl/Sn/my/conf/2016/snrefsdal/img'
        fsave = os.path.join(d, fsave) + '.pdf'

    plot_BV(models_mags, bands, glens, call=callback, xlim=xlim, title=t, fsave=fsave)


def old_run_S4(name, path, bands, e, z, distance, magnification, callback, xlim, is_save, is_glens):
    if e > 0:
        if z > 1:
            ext = extinction.extinction_law_z(ebv=e, bands=bands, z=z)
        else:
            ext = extinction.extinction_law(ebv=e, bands=bands)
    else:
        ext = None

    glens = callback.get_arg(2)
    if glens is None:
        glens = sn_obs.grav_lens_def

    sn_images = sn_obs.coef_magnification(glens)
    if len(sn_images) > 0:
        models_mags = {}  # dict((k, None) for k in names)
        i = 0
        for im, mgf in sn_images.items():
            # if im == 'SX':  # pass image
            #     continue
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
            fsave = "ubv_%s_%s" % (glens, name)

            if ext is not None and ext > 0:
                fsave = "%s_e0%2d" % (fsave, int(e * 100))  # bad formula for name

            d = os.path.expanduser('~/')
            # d = '/home/bakl/Sn/my/conf/2016/snrefsdal/img'
            fsave = os.path.join(d, fsave) + '.pdf'

        if is_glens:
            plot_SX(models_mags, bands, call=callback, xlim=xlim, title=t, fsave=fsave)
        else:
            if is_save:  # Don't save subtitle
                print t
                old_plot_S4(models_mags, bands, glens=glens, call=callback, xlim=xlim, fsave=fsave)
            else:
                old_plot_S4(models_mags, bands, glens=glens, call=callback, xlim=xlim, title=t, fsave=fsave)
    else:
        print "There are no sn images"


def run_S4_curves(name, path, bands, e, z, distance, magnification, callback, xlim, is_save, is_SX):
    if e > 0:
        if z > 1:
            ext = extinction.extinction_law_z(ebv=e, bands=bands, z=z)
        else:
            ext = extinction.extinction_law(ebv=e, bands=bands)
    else:
        ext = None

    glens = callback.get_arg(2)
    if glens is None:
        glens = sn_obs.grav_lens_def

    sn_images = sn_obs.coef_magnification(glens)
    if len(sn_images) > 0:
        models_curves = {}  # dict((k, None) for k in names)
        i = 0
        for im, mgf in sn_images.items():
            # if im == 'SX':  # pass image
            #     continue
            i += 1
            mgf *= magnification
            curves = lc.compute_curves(name, path, bands, ext=ext, z=z, distance=distance, magnification=mgf,
                                       is_show_info=False, is_save=is_save)
            models_curves[im] = curves
            print "Finish image: %s [%d/%d]" % (im, i, len(sn_images))

        mgf = sn_images['S1'] * magnification
        if callback is not None:
            t = "ts=%s z=%4.2f D=%6.2e mu=%3.1f ebv=%4.2f" % (callback.arg_totext(0), z, distance, mgf, e)
        else:
            t = "z=%4.2f D=%6.2e mu=%3.1f ebv=%4.2f" % (z, distance, mgf, e)

        fsave = None
        if is_save:
            fsave = "ubv_%s_%s" % (glens, name)

            if ext is not None and ext > 0:
                fsave = "%s_e0%2d" % (fsave, int(e * 100))  # bad formula for name

            d = os.path.expanduser('~/')
            # d = '/home/bakl/Sn/my/conf/2016/snrefsdal/img'
            fsave = os.path.join(d, fsave) + '.pdf'

        if is_SX:
            plot_SX(models_curves, bands, call=callback, xlim=xlim, title=t, fsave=fsave)
        else:
            if is_save:  # Don't save subtitle
                print t
                plot_S4_curves(models_curves, bands, glens=glens, call=callback, xlim=xlim, fsave=fsave)
            else:
                plot_S4_curves(models_curves, bands, glens=glens, call=callback, xlim=xlim, title=t, fsave=fsave)
    else:
        print "There are no sn images"


def run_ubv_vel(name, path, bands, e, z, distance, magnification, xlim, callback=None,
                is_vel=False, is_save=False):
    if e > 0:
        if z > 1:
            ext = extinction.extinction_law_z(ebv=e, bands=bands, z=z)
        else:
            ext = extinction.extinction_law(ebv=e, bands=bands)
    else:
        ext = None

    models_mags = {}
    models_vels = {}

    mags = lc.compute_mag(name, path, bands, ext=ext, z=z, distance=distance, magnification=magnification,
                          is_show_info=False, is_save=is_save)
    models_mags[name] = mags

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
        t = "ts=%s z=%4.2f D=%6.2e mu=%3.1f ebv=%4.2f" % (callback.arg_totext(0), z, distance, magnification, e)
    else:
        t = "z=%4.2f D=%6.2e mu=%3.1f ebv=%4.2f" % (z, distance, magnification, e)

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

    plot_all(models_vels, models_mags, bands, call=callback, xlim=xlim,
             is_time_points=False, title=t, fsave=fsave)


def usage():
    bands = band.band_get_names().keys()
    print "Usage: lens: "
    print "      grav. lens: ogu-a sha-a die-a "
    print "  ubv.py [params]"
    print "  -b <bands>: string, default: U-B-V-R-I, for example U-B-V-R-I-u-g-i-r-z-UVW1-UVW2.\n" \
          "     Available: " + '-'.join(sorted(bands))
    print "  -i <model name>.  Example: cat_R450_M15_Ni007_E7"
    print "  -p <model directory>, default: ./"
    print "  -e <extinction, E(B-V)> is used to define A_nu, default: 0 "
    print "  -c <callback> [plot_snrefsdal:-56950:1.49:ogu-a (%s)]." % ', '.join(sn_obs.coef_glens().keys())
    print "  -d <distance> [pc].  Default: 11e9 pc"
    print "  -m <magnification>.  Default: 15.4, used for grav lens"
    print "  -z <redshift>.  Default: 1.49"
    print "  -t  plot time points"
    print "  -s  save plot to pdf-file."
    print "  -o  options: [vel, gl, bv]  - plot model velocities,  plot SX with grav.lens, colors [B-V, ...]"
    print "  -w  write magnitudes to file, default 'False'"
    print "  -h  print usage"


#
# def lc_wrapper(param):
#     a = param.split(':')
#     func = a.pop(0)
#     c = cb.CallBack(func, path=cb.plugin_path, args=a, load=1)
#     return c


def main(name=''):
    is_save = False
    is_vel = False
    is_glens = False
    is_BV = False
    path = ''
    z = 1.49
    e = 0.
    gl = 'ogu-a'
    xlim = [-10., 450.]
    magnification = 15.4
    distance = 11e9  # pc
    jd_shift = -56950
    # callback = None
    callback = cb.lc_wrapper('plot_snrefsdal:%s:%s:%s' % (jd_shift, z, gl))

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hsc:d:o:p:e:i:b:m:z:")
    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit(2)

    if len(opts) == 0:
        usage()
        sys.exit(2)

    if not name:
        for opt, arg in opts:
            if opt == '-i':
                path = ROOT_DIRECTORY
                name = os.path.splitext(os.path.basename(str(arg)))[0]
                break
                # if name == '':
                #     print 'Error: you should specify the name of model.'
                #     sys.exit(2)

    # bands = 'F125W-F160W'.split('-')
    bands = 'F105W-F125W-F140W-F160W-F606W-F814W'.split('-')
    # bands = ['U', 'B', 'V', 'R', "I"]
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
            callback = cb.lc_wrapper(str(arg))
            continue
        if opt == '-d':
            distance = float(arg)
            continue
        if opt == '-e':
            e = float(arg)
            continue
        if opt == '-h':
            usage()
            sys.exit(2)
        if opt == '-m':
            magnification = float(arg)
            continue
        if opt == '-o':
            ops = str(arg).split(':')
            is_vel = "vel" in ops
            is_BV = "bv" in ops
            if is_BV:
                bands = str('F814W-F105W-F125W-F160W').split('-')
            is_glens = "gl" in ops
            if is_glens:
                bands = str('F125W-F160W').split('-')
            continue
        if opt == '-p':
            path = os.path.expanduser(str(arg))
            if not (os.path.isdir(path) and os.path.exists(path)):
                print "No such directory: " + path
                sys.exit(2)
            continue
        if opt == '-s':
            is_save = True
            continue
        if opt == '-z':
            z = float(arg)
            continue

    print "Plot magnitudes on z=%f at distance=%e [cosmology D(z)=%s Mpc]" % (z, distance, cosmology_D_by_z(z))

    if is_BV:
        run_BV(name, path, bands, e, z, distance, magnification, callback, xlim=xlim, is_save=is_save)
    elif is_vel:
        run_ubv_vel(name, path, bands, e, z, distance, magnification, xlim=xlim, callback=callback,
                    is_vel=is_vel, is_save=is_save)
    else:
        run_S4_curves(name, path, bands, e, z, distance, magnification, callback, xlim=xlim,
                      is_save=is_save, is_SX=is_glens)


if __name__ == '__main__':
    main()
    # main(name="cat_R1000_M15_Ni007_E15")
