#!/usr/bin/python3
# -*- coding: utf-8 -*-

import getopt
import os
import sys
from os.path import dirname

# matplotlib.use("Agg")
# matplotlib.rcParams['backend'] = "TkAgg"
# matplotlib.rcParams['backend'] = "Qt4Agg"
import math
import matplotlib.pyplot as plt
from matplotlib import gridspec

import pystella.util.callback as cb
from pystella import velocity as vel
from pystella.rf import band
from pystella.rf import light_curve_func as lcf
from pystella.rf import light_curve_plot as lcp
from pystella.util.path_misc import get_model_names
from pystella.util.phys_var import cosmology_D_by_z

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))


def plot_grid(models_dic, bnames, call=None, **kwargs):
    title = kwargs.get('title', '')
    # setup figure
    plt.matplotlib.rcParams.update({'font.size': kwargs.get('fontsize', 12)})
    fig, axs = plt.subplots(int(math.ceil(len(bnames)/2)), 2, sharex='col', sharey='row',
                            figsize=(8, 8))
    plt.subplots_adjust(wspace=0, hspace=0)

    for i, bname in enumerate(bnames):
        icol = i % 2
        irow = int(i / 2)
        ax = axs[irow, icol]
        lcp.plot_models_band(ax, models_dic, bname, **kwargs)

        # plot callback
        if call is not None:
            call.plot(ax, {'bnames': [bname], 'bcolors': {bname: 'black'}, 'markersize': 9})

        if icol == 0:
            ax.set_ylabel('Magnitude')
        if irow == int(len(bnames)/2)-1:
            ax.set_xlabel('Time [days]')

        ax.legend(prop={'size': 8}, loc=4)
        props = dict(facecolor='wheat')
        # props = dict(boxstyle='round', facecolor='white')
        # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(.95, .9, bname, horizontalalignment='right', transform=ax.transAxes, bbox=props)

        if kwargs.get('is_grid', False):
            ax.grid(linestyle=':')

    if title:
        plt.title(title)

    # fig.tight_layout()
    plt.show()
    return fig


def plot_all(models_vels, models_dic, bnames, d=10, call=None, **kwargs):
    # xlim = None, ylim = None,  is_time_points = False, title = '', bshift = None
    title = kwargs.get('title', '')
    # band_shift['UVW1'] = 3
    # band_shift['UVW2'] = 5
    # band_shift['i'] = -1
    is_vel = models_vels is not None

    # setup figure
    plt.matplotlib.rcParams.update({'font.size': 12})
    fig = plt.figure(figsize=(12, 12))
#    fig = plt.figure(num=None, figsize=(12, 12), dpi=100, facecolor='w', edgecolor='k')

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
    lcp.plot_ubv_models(axUbv, models_dic, bnames, **kwargs)
    # lcp.plot_ubv_models(axUbv, models_dic, bands, band_shift=band_shift, xlim=xlim, ylim=ylim,
    #                     is_time_points=is_time_points)

    # plot callback
    if call is not None:
        call.plot(axUbv, {'ax2': axVel})

    # finish plot
    axUbv.set_ylabel('Magnitude')
    axUbv.set_xlabel('Time [days]')
    axUbv.minorticks_on()

    axUbv.legend(prop={'size': 8}, loc=4)
    # ax.set_title(bset)
    if title:
        axUbv.set_title(title)

    # plot velocities
    if is_vel:
        vel.plot_vels_models(axVel, models_vels, xlim=axUbv.get_xlim())
        # vel.plot_vels_sn87a(axVel, z=1.49)
        axVel.legend(prop={'size': 8}, loc=4)

    # show right axes in absolute magnitudes
    if d > 10:
        dm = -5.*math.log10(d/10.)
        axUbvR = axUbv.twinx()
        axUbvR.set_ylim([x+dm for x in axUbv.get_ylim()])
        axUbvR.minorticks_on()
    axUbv.grid(linestyle=':')
    plt.show()
#    plt.show(block=True)
    return fig


def usage():
    print("Usage:")
    print("  ubv.py [params]")
    print("  -b <bands:shift>: string, default: U-B-V-R-I, for example U-B-V-R-I-u-g-i-r-z-UVW1-UVW2.\n"
          "     shift: to move lc along y-axe (minus is '_', for example -b R:2-V-I:_4-B:5 ")
    print("  -i <model name>.  Example: cat_R450_M15_Ni007_E7")
    print("  -p <model directory>, default: ./")
    print("  -e <extinction, E(B-V)> is used to define A_nu, default: 0 ")
    print("  -c <callback> [lcobs:fname:marker:dt:dm, popov[:R:M:E[FOE]:Mni]]. "
          "You can add parameters in format func:params")
    print("  -d <distance> [pc].  Default: 10 pc")
    print("  -g <single, grid, gridm, gridl> Select plot view.  single [default] = all models in one figure"
          ", grid = for each band separate figure.")
    print("  -m <magnification>.  Default: None, used for grav lens")
    print("  -q  turn off quiet mode: print info and additional plots")
    print("  -t  plot time points")
    print("  -s  <file-name> without extension. Save plot to pdf-file. Default: ubv_<file-name>.pdf")
    print("  -x  <xbeg:xend> - xlim, ex: 0:12. Default: None, used all days.")
    print("  -y  <ybeg:yend> - ylim, ex: 26:21. Default: None, used top-magnitude+-5.")
    print("  -v  plot model velocities.")
    print("  -w  write magnitudes to file, default 'False'")
    print("  -z <redshift>.  Default: 0")
    print("  --dt=<t_diff>  time difference between two spectra")
    print("  -l  write plot label")
    print("  -h  print usage")
    print("   --- ")
    band.print_bands()


def lc_wrapper(param, p=None):
    a = param.split(':')
    fname = a.pop(0)
    if p is None:
        if os.path.isfile(fname + '.py'):
            p, fname = os.path.split(fname)
        elif os.path.isfile(os.path.join(os.getcwd(), fname + '.py')):
            p = os.getcwd()
        else:
            p = cb.plugin_path
    print("Call: {} from {}".format(fname, p))
    c = cb.CallBack(fname, path=p, args=a, load=1)
    print("Call: %s from %s" % (c.Func, c.FuncFileFull))
    return c


def main(name='', model_ext='.ph'):
    is_quiet = True
    is_save_mags = False
    is_save_plot = False
    is_plot_time_points = False
    is_extinction = False
    is_vel = False
    view_opts = ('single', 'grid', 'gridl', 'gridm')
    opt_grid = view_opts[0]
    t_diff = 1.01

    label = None
    fsave = None
    # path = ''
    path = os.getcwd()
    z = 0.
    e = 0.
    magnification = 1.
    distance = None  # 10.  # pc
    callback = None
    xlim = None
    ylim = None
    # bshift = None
    bshift = None

    band.Band.load_settings()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hqwtc:d:p:e:g:i:b:l:m:vs:x:y:z:", longopts='dt')
    except getopt.GetoptError as err:
        print(str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(args) > 0:
        path, name = os.path.split(str(args[0]))
        path = os.path.expanduser(path)
        name = name.replace('.ph', '')
    elif len(opts) == 0:
        usage()
        sys.exit(2)

    bnames = ['U', 'B', 'V', 'R', "I"]
    # bands = ['U', 'B', 'V', 'R', "I", 'UVM2', "UVW1", "UVW2", 'g', "r", "i"]

    for opt, arg in opts:
        if opt == '-e':
            e = float(arg)
            is_extinction = True
            continue
        if opt == '-b':
            bnames = []
            for b in str(arg).split('-'):
                # extract band shift
                if ':' in b:
                    bname, shift = b.split(':')
                    bshift = {}
                    if '_' in shift:
                        bshift[bname] = -float(shift.replace('_', ''))
                    else:
                        bshift[bname] = float(shift)
                else:
                    bname = b
                if not band.band_is_exist(bname):
                    print('No such band: ' + bname)
                    sys.exit(2)
                bnames.append(bname)
            continue
        if opt == '-c':
            c = lc_wrapper(str(arg))
            if callback is not None:
                c = cb.CallBackArray((callback, c))
            callback = c
            continue
        if opt == '--dt':
            t_diff = float(arg)
            continue
        if opt == '-q':
            is_quiet = False
            continue
        if opt == '-g':
            opt_grid = str.strip(arg).lower()
            if opt_grid not in view_opts:
                print('No such view option: {0}. Can be '.format(opt_grid, '|'.join(view_opts)))
                sys.exit(2)
            continue
        if opt == '-s':
            is_save_plot = True
            fsave = str.strip(arg)
            continue
        if opt == '-w':
            is_save_mags = True
            continue
        if opt == '-t':
            is_plot_time_points = True
            continue
        if opt == '-v':
            is_vel = True
            continue
        if opt == '-m':
            magnification = float(arg)
            continue
        if opt == '-z':
            z = float(arg)
            continue
        if opt == '-d':
            distance = float(arg)
            continue
        if opt == '-l':
            label = str.strip(arg)
            continue
        if opt == '-x':
            xlim = map(float, str(arg).split(':'))
            continue
        if opt == '-y':
            ylim = map(float, str(arg).split(':'))
            continue
        if opt == '-p':
            path = os.path.expanduser(str(arg))
            if not (os.path.isdir(path) and os.path.exists(path)):
                print("No such directory: " + path)
                sys.exit(2)
            continue
        elif opt == '-h':
            usage()
            sys.exit(2)

    # Set model names
    names = []
    if not name:
        for opt, arg in opts:
            if opt == '-i':
                name = os.path.splitext(os.path.basename(str(arg)))[0]
                names.append(name)
    else:
        names.append(name)

    if len(names) == 0:  # run for all files in the path
        names = get_model_names(path, model_ext)

    # Set distance and redshift
    if distance is None:
        if z > 0:
            distance = cosmology_D_by_z(z)*1e6
            print("Plot magnitudes on z={0:F}. Use cosmology D(z)={1:E} pc".format(z, distance))
        else:
            distance = 10  # pc
    else:
        print("Plot magnitudes on z={0:F} at distance={1:E}".format(z, distance))
        if z > 0:
            print("  Cosmology D(z)={0:E} Mpc".format(cosmology_D_by_z(z)))

    # Run plotting
    if len(names) > 0:
        models_mags = {}  # dict((k, None) for k in names)
        models_vels = {}  # dict((k, None) for k in names)
        i = 0
        for name in names:
            i += 1
            # mags = lcf.compute_mag(name, path, bands, ext=ext, z=z, distance=distance, magnification=magnification,
            #                        t_diff=t_diff, is_show_info=not is_quiet, is_save=is_save_mags)
            curves = lcf.curves_compute(name, path, bnames, z=z, distance=distance,
                                        magnification=magnification, t_diff=t_diff)

            if is_extinction:
                lcf.curves_reddening(curves, ebv=e, z=z)
            # lcf.plot_curves(curves)
            # exit()
            # models_mags[name] = mags
            models_mags[name] = curves

            if is_save_mags:
                fname = os.path.join(path, name + '.ubv')
                lcf.curves_save(curves, fname)
                print("Magnitudes have been saved to " + fname)

            if not (is_save_mags or is_quiet):
                # z, distance = 0.145, 687.7e6  # pc for comparison with Maria
                # lcf.plot_bands(mags, bands, title=name, fname='', is_time_points=is_plot_time_points)
                lcp.plot_bands(curves, bnames, title=name, fname='', is_time_points=is_plot_time_points)

            if is_vel:
                vels = vel.compute_vel(name, path, z=z)
                if vels is None:
                    sys.exit("No data for: %s in %s" % (name, path))
                models_vels[name] = vels
                print("Finish velocity: %s [%d/%d]" % (name, i, len(names)))
            else:
                models_vels = None
                print("Finish mags: %s [%d/%d] in %s" % (name, i, len(names), path))

        if label is None:
            if callback is not None:
                label = "ts=%s z=%4.2f D=%6.2e mu=%3.1f ebv=%4.2f" % (
                    callback.arg_totext(0), z, distance, magnification, e)
            else:
                label = "z=%4.2f D=%6.2e mu=%3.1f ebv=%4.2f" % (z, distance, magnification, e)

        if not is_save_mags:
            if is_save_plot:
                if len(fsave) == 0:
                    if is_vel:
                        fsave = "ubv_vel_%s" % name
                    else:
                        fsave = "ubv_%s" % name

                if is_extinction and e > 0:
                    fsave = "%s_e0%2d" % (fsave, int(e * 100))  # bad formula for name

                d = os.path.expanduser('~/')
                # d = '/home/bakl/Sn/my/conf/2016/snrefsdal/img'

                fsave = os.path.join(d, os.path.splitext(fsave)[0]) + '.pdf'

            if opt_grid in view_opts[1:]:
                sep = opt_grid[:-1]
                if sep == 'd':
                    sep = 'l'  # line separator
                fig = plot_grid(models_mags, bnames, call=callback, xlim=xlim, ylim=ylim,
                                sep=sep, is_grid=False)
            else:
                fig = plot_all(models_vels, models_mags, bnames, d=distance, call=callback, xlim=xlim, ylim=ylim,
                               is_time_points=is_plot_time_points, title=label, bshift=bshift)
            if fsave is not None:
                print("Save plot to %s " % fsave)
                fig.savefig(fsave, bbox_inches='tight')
                # plt.savefig(fsave, format='pdf')

                # plot_all(dic_results, bands,  xlim=(-10, 410), is_time_points=is_plot_time_points)
            # plot_all(dic_results, bands, xlim=(-10, 410), callback=callback, is_time_points=is_plot_time_points)
            # plot_all(dic_results, bands,  ylim=(40, 23),  is_time_points=is_plot_time_points)
    else:
        print("There are no models in the directory: %s with extension: %s " % (path, model_ext))


if __name__ == '__main__':
    main()
    # main(name="cat_R1000_M15_Ni007_E15")
