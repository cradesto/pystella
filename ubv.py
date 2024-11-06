#!/usr/bin/env python3

from os.path import dirname, abspath

# matplotlib.use("Agg")
# matplotlib.rcParams['backend'] = "TkAgg"
# matplotlib.rcParams['backend'] = "Qt4Agg"
import math
import pystella as ps

import logging

mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.ERROR)

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(abspath(__file__)))


def plot_grid(models_dic, bnames, call=None, **kwargs):
    import matplotlib.pyplot as plt
    import numpy as np

    title = kwargs.get('title', False)
    fontsize = kwargs.get('fontsize', 12)
    figsize = kwargs.get('figsize', (8, 8))
    xtype = kwargs.get('xtype', 'lin')
    legloc = kwargs.get('legloc', 4)

    # setup figure
    # plt.matplotlib.rcParams.update({'font.size':})
    if len(bnames) == 1:
        fig, axs = plt.subplots(figsize=figsize)
        axs = np.array(axs)
    else:
        fig, axs = plt.subplots(math.ceil(len(bnames) / 2), 2, sharex='col', sharey='row', figsize=figsize)
    plt.subplots_adjust(wspace=0, hspace=0)

    for i, bname in enumerate(bnames):
        icol = i % 2
        irow = int(i / 2)
        # print(i, irow, icol, axs.shape[0])
        # ax = axs[irow, icol]
        ax = axs.ravel()[2 * irow + icol]
        ps.lcp.plot_models_band(ax, models_dic, bname, **kwargs)

        # plot callback
        if call is not None:
            call.plot(ax, {'bnames': [bname], 'bcolors': {bname: 'black'}, 'markersize': 9})

        if icol == 0:
            ax.set_ylabel('Magnitude')

        ax.legend(prop={'size': 8}, loc=legloc)
        props = dict(facecolor='wheat')
        # props = dict(boxstyle='round', facecolor='white')
        # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(.95, .9, bname, horizontalalignment='right', transform=ax.transAxes, bbox=props, fontsize=fontsize)

        if kwargs.get('is_grid', False):
            ax.grid(linestyle=':')

        # # Setup axes
        # for i, bname in enumerate(bnames):
        #     icol = i % 2
        #     irow = int(i / 2)
        #     ax = axs.ravel()[i * irow + icol]
        if irow == math.ceil(len(bnames) / 2) - 1:
            ax.set_xlabel('Time [days]')

        if xtype == 'log':
            ax.set_xscale('log')

    if title:
        axs.ravel()[0].set_title(title)

    # fig.tight_layout()
    # plt.show()
    return fig


def plot_all(models_vels, models_dic, bnames, d=10, call=None, **kwargs):
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    # xlim = None, ylim = None,  is_time_points = False, title = '', bshift = None
    title = kwargs.get('title', '')
    legend = kwargs.get('legend', 'box')
    is_axes_right = kwargs.get('is_axes_right', True)
    is_grid = kwargs.get('is_grid', True)
    legloc = kwargs.get('legloc', 4)
    fontsize = kwargs.get('fontsize', 12)
    bshift = kwargs.get('bshift', None)
    xtype = kwargs.get('xtype', 'lin')
    is_lum = kwargs.get('is_lum', False)
    # band_shift['UVW1'] = 3
    # band_shift['UVW2'] = 5
    # band_shift['i'] = -1
    is_vel = models_vels is not None

    # setup figure
    plt.matplotlib.rcParams.update({'font.size': fontsize})
    plt.matplotlib.rcParams['font.family'] = 'sans-serif'
    fig = plt.figure(figsize=(12, 12))
    #    fig = plt.figure(num=None, figsize=(12, 12), dpi=100, facecolor='w', edgecolor='k')

    if is_vel:
        gs1 = gridspec.GridSpec(4, 1)
        axUbv = fig.add_subplot(gs1[:-1, 0])
        # axVel = fig.add_subplot(gs1[3, 0])
        axVel = fig.add_subplot(gs1[3, 0], sharex=axUbv)
        # plt.setp(axUbv.get_xticklabels(), visible=False)
        axUbv.xaxis.tick_top()
        axUbv.tick_params(axis="x",direction="in", which='both')
        # axUbv.tick_params(direction="in")
    else:
        gs1 = gridspec.GridSpec(1, 1)
        axUbv = fig.add_subplot(gs1[0, 0])
        axVel = None
    gs1.update(wspace=0.3, hspace=0., left=0.1, right=0.95)

    # plot the light curves
    if is_lum:
        ps.lcp.plot_lum_models(axUbv, models_dic, bnames, **kwargs)
    else:
        ps.lcp.plot_ubv_models(axUbv, models_dic, bnames, **kwargs)

    # plot callback
    if call is not None:
        if is_vel:
            call.plot((axUbv, axVel), dic={'bnames': bnames, 'bshift': bshift})
        else:
            call.plot(axUbv, dic={'bnames': bnames, 'bshift': bshift})
    # finish plot
    if is_lum:
        axUbv.set_ylabel('Luminosity [erg/s]')  # Magnitude  Mag
        axUbv.set_yscale('log')
    else:
        axUbv.set_ylabel('Magnitude')  # Magnitude  Mag
    
    axUbv.set_xlabel('Time since explosion [days]')
    axUbv.minorticks_on()

    if legend == 'box':
        axUbv.legend(loc=legloc, frameon=False)

    if xtype == 'log':
        axUbv.set_xscale('log')

    # axUbv.legend(prop={'size': 8}, loc=legloc)
    if title:
        axUbv.set_title(title)

    # plot velocities
    if is_vel:
        ps.vel.plot_vels_models(axVel, models_vels, xlim=axUbv.get_xlim())
        # vel.plot_vels_sn87a(axVel, z=1.49)
        axVel.legend(loc=legloc)
        if xtype == 'log':
            axVel.set_xscale('log')
        # axVel.legend(prop={'size': 8}, loc=legloc)

    # show right axes in absolute magnitudes
    if is_axes_right and d > 10:
        dm = -5. * math.log10(d / 10.)
        axUbvR = axUbv.twinx()
        axUbvR.set_ylim([x + dm for x in axUbv.get_ylim()])
        axUbvR.minorticks_on()
    if is_grid:
        axUbv.grid(linestyle=':')
    #    plt.show(block=True)
    if is_vel:
        # gs1.tight_layout(fig, rect=[0, 0, 0, 0.])
        fig.tight_layout()
    return fig


def usage():
    print("Usage:")
    print("  ubv.py [params]")
    print("  -b <bands:shift>: string, default: U-B-V-R-I, for example U-B-V-R-I-u-g-i-r-z-UVW1-UVW2.\n"
          "     shift: to move lc along y-axe (minus is '_', for example -b R:2-V-I:_4-B:5 ")
    print("  -i <model-name OR pattern, like '*R450*'>.  Example: cat_R450_M15_Ni007_E7")
    print("  -p <model directory>, default: ./")
    print("  -e <extinction, E(B-V)> is used to define A_nu, default: 0 ")
    print("  -c <callback> [lcobs:fname:marker:dt:dm, velobs:fname:marker:dt:vnorm(to cm/s), "
          "popov[:R:M:E[FOE]:Mni], lcobssm as lcobs, but for sm-format data-files]. "
          "You can add parameters in format func:params")
    print("  -d <distance> [pc].  Default: 10 pc")
    print("  -g <single, grid, gridm, gridl> Select plot view.  single [default] = all models in one figure"
          ", grid = for each band separate figure.")
    print("  -o <is_axes_right>.  Default: empty string")
    print("  -m <magnification>.  Default: None, used for grav lens")
    print("  -q  turn off quiet mode: print info and additional plots")
    print("  -t  plot time points")
    print("  -s  <file-name> without extension. Save plot to pdf-file. Default: ubv_<file-name>.pdf")
    print("  -x  <xbeg:xend[:xtype]> - xlim, ex: -x 0:12. xtype is optional and can be 'lin'(default value) or 'log'."
          "if you use log scale for X you should set xbeg something like 0.1 or larger, e.g.  -x 0.1:200:log"
          " Default: None, used all days.")
    print("  -y  <ybeg:yend> - ylim, ex: 26:21. Default: None, used top-magnitude+-5.")
    print("  -v  <swd OR ttres[ttresold]> - plot model velocities computed from swd OR tt-res files[ttresold for old res "
          "format].")
    print("  -w  write magnitudes to out-file. Use '1' for the default name of out-file")
    print("  -z <redshift>.  Default: 0")
    print("  --dt=<t_diff>  time difference between two spectra")
    print("  --curve-old  - use old procedure")
    print("  --curve-tt  - take curves from tt-file: UBVRI+bol")
    print("  --legloc=<legloc>  legend position. It could be integer 0-10. Default: 0 "
            'Maybe: best 0, upper right 1, upper left 2, lower left 3, lower right 4, right 5, center left 6, center right 7, lower center 8, upper center 9, center 10')
    print("  --lum  plot luminosity at y-axe")
    print("  -l  write plot label")
    print("  -h  print usage")
    print("   --- ")
    ps.band.print_bands()




def main(name=None, model_ext='.ph'):
    import sys
    import os
    import getopt
    import fnmatch

    is_quiet = False
    is_save_mags = False
    is_save_plot = False
    is_plot_time_points = False
    is_extinction = False
    is_curve_old = False
    is_curve_tt = False
    is_axes_right = False
    is_grid = False
    is_lum = False

    vel_mode = None
    view_opts = ('single', 'grid', 'gridl', 'gridm')
    opt_grid = view_opts[0]
    t_diff = 1.01
    linestyles = ['-']

    label = None
    fsave = None
    fname = None
    # path = ''
    path = os.getcwd()
    z = 0.
    e = 0.
    magnification = 1.
    distance = None  # 10.  # pc
    callback = None
    xlim = None
    ylim = None
    xtype = 'lin'
    # bshift = None
    bshift = None
    legloc = 0
    level = logging.INFO
    
    ps.band.Band.load_settings()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hqtb:c:d:e:g:i:l:o:m:p:v:s:w:x:y:z:",
                                   ['dt=', 'curve-old', 'curve-tt', 'legloc=', 'lum'])
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

    bnames = ['U', 'B', 'V', 'R']
    # bands = ['U', 'B', 'V', 'R', "I", 'UVM2', "UVW1", "UVW2", 'g', "r", "i"]

    for opt, arg in opts:
        if opt == '-e':
            e = float(arg)
            is_extinction = True
            continue
        if opt == '-b':
            bnames = []
            bshift = {}
            for b in str(arg).split('-'):
                # extract band shift
                if ':' in b:
                    bname, shift = b.split(':')
                    if '_' in shift:
                        bshift[bname] = -float(shift.replace('_', ''))
                    else:
                        bshift[bname] = float(shift)
                else:
                    bname = b
                if not ps.band.is_exist(bname):
                    print('No such band: ' + bname)
                    sys.exit(2)
                bnames.append(bname)
            if len(bshift) == 0:
                bshift = None
            continue
        if opt == '-c':
            c = ps.cb.lc_wrapper(str(arg))
            if callback is not None:
                c = ps.cb.CallBackArray((callback, c))
            callback = c
            continue
        if opt == '--dt':
            t_diff = float(arg)
            continue
        if opt == '--legloc':
            legloc = int(arg)
            continue
        if opt == '-q':
            is_quiet = True
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
        if opt == '--curve-old':
            is_curve_old = True
            continue
        if opt == '--curve-tt':
            is_curve_tt = True
            continue
        if opt == '--lum':
            is_lum = True
            continue
        if opt == '-w':
            is_save_mags = True
            if arg != '1':
                fname = arg.strip()
            continue
        if opt == '-t':
            is_plot_time_points = True
            continue
        if opt == '-v':
            vel_mode = arg.strip()
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
            if 'log' in arg:
                xtype = 'log'
                s12 = arg.replace(xtype, '')
                s12 = s12.rstrip(':')
            else:
                s12 = arg
            # print(f's12= {s12}')
            xlim = ps.str2interval(s12, llim=0, rlim=float('inf'))
            continue
        if opt == '-y':
            ylim = ps.str2interval(arg, llim=-10, rlim=-22)
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
    is_set_model = False
    if name is None:
        for opt, arg in opts:
            if opt == '-i':
                nm = os.path.splitext(os.path.basename(str(arg)))[0]
                if os.path.exists(os.path.join(path, nm + model_ext)):
                    names.append(nm)
                elif is_curve_tt:
                    if os.path.exists(os.path.join(path, nm + '.tt')):
                        names.append(nm)
                else:
                    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))
                             and fnmatch.fnmatch(f, arg)]
                    for f in files:
                        nm = os.path.splitext(os.path.basename(f))[0]
                        names.append(nm)
                    names = list(set(names))
                    print('Input {} models: {}'.format(len(names), ' '.join(names)))
                is_set_model = True
    else:
        print(name)
        names.append(name)
        is_set_model = True

    if len(names) == 0 and not is_set_model:  # run for all files in the path
        names = ps.path.get_model_names(path, model_ext)

    # Set distance and redshift
    if distance is None:
        if z > 0:
            distance = ps.cosmology_D_by_z(z) * 1e6
            print("Plot magnitudes on z={0:F} with D(z)={1:E} pc (cosmological)".format(z, distance))
        else:
            distance = 10  # pc
    else:
        print("Plot magnitudes on z={0:F} at distance={1:E}".format(z, distance))
        if z > 0:
            print("  Cosmology D(z)={0:E} Mpc".format(ps.cosmology_D_by_z(z)))

    # Run models
    if len(names) > 0:
        models_mags = {}  # dict((k, None) for k in names)
        models_vels = {}  # dict((k, None) for k in names)
        for i, name in enumerate(names):
            mdl = ps.Stella(name, path=path)

            if is_curve_tt:  # tt
                print("The curves [UBVRI+bol] was taken from tt-file. IMPORTANT: distance: 10 pc, z=0, E(B-V) = 0")
                curves = mdl.get_tt().read_curves()
            elif is_curve_old:  # old
                print("Use old proc for Stella magnitudes")
                curves = ps.lcf.curves_compute(name, path, bnames, z=z, distance=distance,
                                               magnification=magnification, t_diff=t_diff)
                if is_extinction:
                    curves = ps.lcf.curves_reddening(curves, ebv=e, z=z)
            else:
                curves = mdl.curves(bnames, z=z, distance=distance, ebv=e, magnification=magnification,
                                    t_diff=t_diff)

            models_mags[name] = curves

            if vel_mode is not None:
                if vel_mode == 'swd':
                    vels = ps.vel.compute_vel_swd(name, path, z=z)
                elif vel_mode.startswith('ttres'):
                    vels = ps.vel.compute_vel_res_tt(name, path, z=z,
                                                     is_info=False, is_new_std='old' not in vel_mode.lower())
                else:
                    raise ValueError('This mode [{}] for velocity is not supported'.format(vel_mode))

                if vels is None:
                    sys.exit("Error: no data for: %s in %s" % (name, path))
                models_vels[name] = vels
                print("[%d/%d] Done mags & velocity for %s" % (i + 1, len(names), name))
            else:
                models_vels = None
                print("[%d/%d] Done mags for %s" % (i + 1, len(names), name))

        if label is None:
            if callback is not None:
                label = "ts=%s z=%4.2f D=%6.2e mu=%3.1f ebv=%4.2f" % (
                    callback.arg_totext(0), z, distance, magnification, e)
            else:
                label = "z=%4.2f D=%6.2e mu=%3.1f ebv=%4.2f" % (z, distance, magnification, e)

        # save curves to files
        if is_save_mags:
            for curves in models_mags.values():
                if fname is None:
                    fname = os.path.join(path, curves.Name)
                    if z > 0.:
                        fname = '{}_Z{:.2g}'.format(fname, z)
                    if distance > 10.:
                        fname = '{}_D{:.2e}'.format(fname, distance)
                    if e > 0:
                        fname = '{}_E{:0.2g}'.format(fname, e)
                    fname = '{}{}'.format(fname, '.ubv')
                if ps.lcf.curves_save(curves, fname):
                    print("Magnitudes of {} have been saved to {}".format(curves.Name, fname))
                else:
                    print("Error with Magnitudes saved to {}".format(curves.Name, fname))
        # plot
        elif not is_quiet:
            if opt_grid in view_opts[1:]:
                sep = opt_grid[:-1]
                if sep == 'd':
                    sep = 'l'  # line separator
                fig = plot_grid(models_mags, bnames, call=callback, xlim=xlim, xtype=xtype, legloc=legloc, ylim=ylim, title=label,
                                sep=sep, is_grid=False)
            else:
                # linestyles = ['--', '-.', '-', ':']
                fig = plot_all(models_vels, models_mags, bnames, d=distance, call=callback, xlim=xlim, xtype=xtype,
                               ylim=ylim, is_time_points=is_plot_time_points, title=label, bshift=bshift,
                               is_axes_right=is_axes_right, is_grid=is_grid, legloc=legloc, fontsize=18,
                               lines=linestyles, is_lum=is_lum)
                # lcp.setFigMarkersBW(fig)
                # lcp.setFigLinesBW(fig)

            if is_save_plot:
                if len(fsave) == 0:
                    if vel_mode is not None:
                        fsave = "ubv_vel_%s" % name
                    else:
                        fsave = "ubv_%s" % name

                if is_extinction and e > 0:
                    fsave = "%s_e0%2d" % (fsave, int(e * 100))  # bad formula for name

                fsave = os.path.expanduser(fsave)
                fsave = os.path.splitext(fsave)[0] + '.pdf'

                print("Save plot to %s " % os.path.abspath(fsave))
                fig.savefig(fsave, bbox_inches='tight', format='pdf')
            else:
                import matplotlib.pyplot as plt
                # plt.ion()
                plt.show()
                # plt.pause(0.0001)
                # print('')
                # input("===> Hit <return> to quit")

    else:
        print("There are no such models in the directory: %s with extension: %s " % (path, model_ext))


if __name__ == '__main__':
    main()
