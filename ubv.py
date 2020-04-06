#!/usr/bin/env python3

from os.path import dirname, abspath

# matplotlib.use("Agg")
# matplotlib.rcParams['backend'] = "TkAgg"
# matplotlib.rcParams['backend'] = "Qt4Agg"
import math
import pystella as ps

import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(abspath(__file__)))


def plot_grid(models_dic, bnames, call=None, **kwargs):
    import matplotlib.pyplot as plt

    title = kwargs.get('title', '')
    fontsize = kwargs.get('fontsize', 12)
    # setup figure
    # plt.matplotlib.rcParams.update({'font.size':})
    fig, axs = plt.subplots(math.ceil(len(bnames)/2), 2, sharex='col', sharey='row', figsize=(8, 8))
    plt.subplots_adjust(wspace=0, hspace=0)

    for i, bname in enumerate(bnames):
        icol = i % 2
        irow = int(i / 2)
        ax = axs[irow, icol]
        ps.lcp.plot_models_band(ax, models_dic, bname, **kwargs)

        # plot callback
        if call is not None:
            call.plot(ax, {'bnames': [bname], 'bcolors': {bname: 'black'}, 'markersize': 9})

        if icol == 0:
            ax.set_ylabel('Magnitude')

        ax.legend(prop={'size': 8}, loc=4)
        props = dict(facecolor='wheat')
        # props = dict(boxstyle='round', facecolor='white')
        # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(.95, .9, bname, horizontalalignment='right', transform=ax.transAxes, bbox=props, fontsize=fontsize)

        if kwargs.get('is_grid', False):
            ax.grid(linestyle=':')

    # Setup axes
    for ax in axs[-1, :]:
        ax.set_xlabel('Time [days]')

    if title:
        plt.title(title)

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
        axVel = fig.add_subplot(gs1[3, 0])
    else:
        gs1 = gridspec.GridSpec(1, 1)
        axUbv = fig.add_subplot(gs1[0, 0])
        axVel = None
    gs1.update(wspace=0.3, hspace=0.3, left=0.1, right=0.95)

    # plot the light curves
    ps.lcp.plot_ubv_models(axUbv, models_dic, bnames, **kwargs)
    # lcp.plot_ubv_models(axUbv, models_dic, bands, band_shift=band_shift, xlim=xlim, ylim=ylim,
    #                     is_time_points=is_time_points)

    # plot callback
    if call is not None:
        if is_vel:
            call.plot((axUbv, axVel), dic={'bnames': bnames, 'bshift': bshift})
        else:
            call.plot(axUbv, dic={'bnames': bnames, 'bshift': bshift})
    # finish plot
    axUbv.set_ylabel('Magnitude')  # Magnitude  Mag
    axUbv.set_xlabel('Time since explosion [days]')
    axUbv.minorticks_on()

    if legend == 'box':
        axUbv.legend(loc=legloc, frameon=False)

    # axUbv.legend(prop={'size': 8}, loc=legloc)
    if title:
        axUbv.set_title(title)

    # plot velocities
    if is_vel:
        ps.vel.plot_vels_models(axVel, models_vels, xlim=axUbv.get_xlim())
        # vel.plot_vels_sn87a(axVel, z=1.49)
        axVel.legend(loc=legloc)
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
    return fig


def usage():
    print("Usage:")
    print("  ubv.py [params]")
    print("  -b <bands:shift>: string, default: U-B-V-R-I, for example U-B-V-R-I-u-g-i-r-z-UVW1-UVW2.\n"
          "     shift: to move lc along y-axe (minus is '_', for example -b R:2-V-I:_4-B:5 ")
    print("  -i <model-name OR pattern, like '*R450*'>.  Example: cat_R450_M15_Ni007_E7")
    print("  -p <model directory>, default: ./")
    print("  -e <extinction, E(B-V)> is used to define A_nu, default: 0 ")
    print("  -c <callback> [lcobs:fname:marker:dt:dm, velobs:fname:marker:dt:vnorm(1e8), "
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
    print("  -x  <xbeg:xend> - xlim, ex: 0:12. Default: None, used all days.")
    print("  -y  <ybeg:yend> - ylim, ex: 26:21. Default: None, used top-magnitude+-5.")
    print("  -v  <swd OR ttres> - plot model velocities computed from swd OR tt-res files.")
    print("  -w  write magnitudes to out-file. Use '1' for the default name of out-file")
    print("  -z <redshift>.  Default: 0")
    print("  --dt=<t_diff>  time difference between two spectra")
    print("  --curve-old  - use old procedure")
    print("  --curve-tt  - take curves from tt-file: UBVRI+bol")
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
    # bshift = None
    bshift = None

    ps.band.Band.load_settings()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hqtb:c:d:e:g:i:l:o:m:p:v:s:w:x:y:z:",
                                   ['dt=', 'curve-old', 'curve-tt'])
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
                if not ps.band.band_is_exist(bname):
                    print('No such band: ' + bname)
                    sys.exit(2)
                bnames.append(bname)
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
            xlim = ps.str2interval(arg, llim=0, rlim=float('inf'))
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
                if os.path.exists(os.path.join(path, nm+model_ext)):
                    names.append(nm)
                else:
                    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))
                             and fnmatch.fnmatch(f, arg)]
                    for f in files:
                        nm = os.path.splitext(os.path.basename(f))[0]
                        names.append(nm)
                        print('input: {}'.format(nm))
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
                if vel_mode == 'ttres':
                    vels = ps.vel.compute_vel_res_tt(name, path, z=z)
                elif vel_mode == 'swd':
                    vels = ps.vel.compute_vel_swd(name, path, z=z)
                else:
                    raise ValueError('This mode [{}] for velocity is not supported'.format(vel_mode))

                if vels is None:
                    sys.exit("Error: no data for: %s in %s" % (name, path))
                models_vels[name] = vels
                print("[%d/%d] Done mags & velocity for %s" % (i+1, len(names), name))
            else:
                models_vels = None
                print("[%d/%d] Done mags for %s" % (i+1, len(names), name))

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
                fig = plot_grid(models_mags, bnames, call=callback, xlim=xlim, ylim=ylim,
                                sep=sep, is_grid=False)
            else:
                # linestyles = ['--', '-.', '-', ':']
                fig = plot_all(models_vels, models_mags, bnames, d=distance, call=callback, xlim=xlim, ylim=ylim,
                               is_time_points=is_plot_time_points, title=label, bshift=bshift,
                               is_axes_right=is_axes_right, is_grid=is_grid, legloc=1, fontsize=14,
                               lines=linestyles)
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

                d = os.path.expanduser('~/')
                fsave = os.path.join(d, os.path.splitext(fsave)[0]) + '.pdf'

                print("Save plot to %s " % fsave)
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
