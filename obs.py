#!/usr/bin/env python3
# #!/usr/bin/python3

import getopt
# import matplotlib
# matplotlib.use("Agg")
# matplotlib.rcParams['backend'] = "TkAgg"
# matplotlib.rcParams['backend'] = "Qt4Agg"
import math
import os
import sys
from os.path import dirname
import logging

import matplotlib.pyplot as plt
from matplotlib import gridspec

import pystella as ps
# import pystella.util.callback as cb
# from pystella.rf import band

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))

logger = logging.getLogger(__name__)
level = logging.INFO

def merge_obs(inp, type_el):
    result = None
    if isinstance(inp, (list, tuple)):
        for o in inp:
            if isinstance(o, type_el):
                result = type_el.Merge(result, o)
    else:
        if isinstance(inp, type_el):
            result = inp
    return result


def plot_grid(call, bnames, xlim=None, ylim=None, **kwargs):
    title = kwargs.get('title', '')
    markersize = kwargs.get('markersize', 9)
    # setup figure
    plt.matplotlib.rcParams.update({'font.size': kwargs.get('fontsize', 14)})

    nrows = math.ceil(len(bnames)/2)
    ncols = 1 if len(bnames) == 1 else 2
    # if len(bnames) > 1:
    fig, axs = plt.subplots(nrows, ncols, sharex='col', sharey='row', figsize=(8, 8))
    # else:
    #     fig, axs = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(8, 8))
    plt.subplots_adjust(wspace=0, hspace=0)

    res = []
    for i, bname in enumerate(bnames):
        icol = i % ncols
        irow = int(i / ncols)
        if nrows > 1:
            ax = axs[irow, icol]
        elif nrows == 1 and ncols > 1:
            ax = axs[icol]
        else:
            ax = axs

        if icol == 1:
            ax = ax.twinx()

        # plot callback
        if call is not None:
            r = call.plot(ax, {'bnames': [bname], 'bcolors': {bname: 'black'}, 'markersize': markersize})
            res.append(r)

        # if icol == 0:
        ax.set_ylabel('Magnitude')
        if irow == int(len(bnames) / 2) - 1:
            ax.set_xlabel('Time [days]')

        props = dict(facecolor='wheat')
        # props = dict(boxstyle='round', facecolor='white')
        # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(.95, .9, bname, horizontalalignment='right', transform=ax.transAxes, bbox=props)

        ax.legend(prop={'size': 8}, loc=4)
        # if icol == 0:
        ax.invert_yaxis()

        if xlim is not None:
            # xlim = ax.get_xlim()
            ax.set_xlim(xlim)
        if ylim is not None:
            # ylim = ax.get_ylim()
            ax.set_ylim(ylim)

        if kwargs.get('is_grid', False):
            ax.grid(linestyle=':')

    if title:
        plt.title(title)

    return fig, res


def plot_all(call, bnames, xlim=None, ylim=None, **kwargs):
    # xlim = None, ylim = None,  is_time_points = False, title = '', bshift = None
    title = kwargs.get('title', '')
    markersize = kwargs.get('markersize', 9)

    # setup figure
    plt.matplotlib.rcParams.update({'font.size': 12})
    fig = plt.figure(figsize=(12, 12))
    #    fig = plt.figure(num=None, figsize=(12, 12), dpi=100, facecolor='w', edgecolor='k')

    gs1 = gridspec.GridSpec(1, 1)
    axUbv = fig.add_subplot(gs1[0, 0])
    gs1.update(wspace=0.3, hspace=0.3, left=0.1, right=0.95)

    res = call.plot(axUbv, {'bnames': bnames, 'markersize': markersize})

    # finish plot
    axUbv.set_ylabel('Magnitude')
    axUbv.set_xlabel('Time [days]')
    axUbv.minorticks_on()
    if xlim is not None:
        # xlim = ax.get_xlim()
        axUbv.set_xlim(xlim)

    axUbv.invert_yaxis()
    if ylim is not None:
        # ylim = ax.get_ylim()
        axUbv.set_ylim(ylim)

    axUbv.legend(prop={'size': 8}, loc=4)
    # ax.set_title(bset)
    if title:
        axUbv.set_title(title)
    axUbv.grid(linestyle=':')
    return fig, res


def usage():
    print("Usage:")
    print("  obs.py [params]")
    print("  -b <bands:shift>: string, default: U-B-V-R-I.\n"
          "     shift: to move lc along y-axe (minus is '_', for example -b R:2-V-I:_4-B:5 ")
    print("  -c <callback> [lcobs:fname:marker:dt:dm, velobs:fname:marker:dt:vnorm(to cm/s), "
          "popov[:R:M:E[FOE]:Mni], lcobssm as lcobs, but for sm-format data-files]. "
          "You can add parameters in format func:params")
    print("  -g <single, grid, gridm, gridl> Select plot view.  single [default] = all models in one figure"
          ", grid = for each band separate figure.")
    print("  -s  without extension. Save plot to pdf-file. Default: if 1, fname =  'ubv_obs.pdf'")
    print("  -w  write magnitudes to out-file.")
    print("  -x  <xbeg:xend> - xlim, ex: 0:12. Default: None, used all days.")
    print("  -l  write plot label")
    print("  -h  print usage")
    print("   --- ")
    ps.band.print_bands()


def old_lc_wrapper(param, p=None):
    a = param.split(':')
    fname = a.pop(0)
    if p is None:
        if os.path.isfile(fname + '.py'):
            p, fname = os.path.split(fname)
        elif os.path.isfile(os.path.join(os.getcwd(), fname + '.py')):
            p = os.getcwd()
        else:
            p = ps.cb.plugin_path
    logger.info("Call: {} from {}".format(fname, p))
    c = ps.cb.CallBack(fname, path=p, args=a, load=1)
    logger.info("Call: %s from %s" % (c.Func, c.FuncFileFull))
    return c


def main():
    def is_level(lvl):
        levels = ['CRITICAL', 'FATAL','ERROR','WARN','WARNING','INFO','DEBUG','NOTSET']
        return lvl.upper() in levels
        
    view_opts = ('single', 'grid', 'gridl', 'gridm')
    opt_grid = view_opts[0]

    label = None
    fsave = None
    is_save_mags = False
    callback = None
    xlim = None
    ylim = None
    level = logging.INFO
    
    ps.band.Band.load_settings()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hqc:g:b:l:s:w:x:y:", '--log')
    except getopt.GetoptError as err:
        print(str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(opts) == 0:
        usage()
        sys.exit(2)

    bnames = None
    # bnames = ['U', 'B', 'V', 'R', "I"]
    # bands = ['U', 'B', 'V', 'R', "I", 'UVM2', "UVW1", "UVW2", 'g', "r", "i"]

    for opt, arg in opts:
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
                if not ps.band.is_exist(bname):
                    logger.error('No such band: ' + bname)
                    sys.exit(2)
                bnames.append(bname)
            continue
        if opt == '-c':
            c = ps.cb.lc_wrapper(str(arg))
            if callback is not None:
                c = ps.cb.CallBackArray((callback, c))
            callback = c
            continue
        if opt == '-g':
            opt_grid = str.strip(arg).lower()
            if opt_grid not in view_opts:
                logger.error('No such view option: {0}. Can be '.format(opt_grid, '|'.join(view_opts)))
                sys.exit(2)
            continue
        if opt == '-l':
            label = str.strip(arg)
            continue
        if opt == '--log':
            if is_level(args.log):
                level = logging.getLevelName(args.log.upper())
            else:
                level = logging.INFO
                print(f"ERROR: Bad value for log: {level}. See help.")
            continue
        if opt == '-s':
            fsave = str.strip(arg)
            continue
        if opt == '-w':
            fsave = str.strip(arg)
            is_save_mags = True
            continue
        if opt == '-x':
            xlim = list(float(x) for x in str(arg).split(':'))
            continue
        if opt == '-y':
            ylim = list(float(x) for x in str(arg).split(':'))
            continue
        elif opt == '-h':
            usage()
            sys.exit(2)
    
    logger.setLevel(level)
    logging.basicConfig(level=level)
    for nm in ['velobs', 'lcobs']:
        log = logging.getLogger(nm).setLevel(level)

    if callback is None:
        logger.error('No  obs data. You my use lcobs or other callbacks.')
        usage()
        sys.exit(2)

    if opt_grid in view_opts[1:]:
        sep = opt_grid[:-1]
        if sep == 'd':
            sep = 'l'  # line separator
        fig, res = plot_grid(callback, bnames, xlim=xlim, ylim=ylim, sep=sep, is_grid=False)
    else:
        fig, res = plot_all(callback, bnames, xlim=xlim, ylim=ylim, title=label)

    if is_save_mags:
        if len(res) == 0:
            logger.error('No data to write.')
        else:
            curves = merge_obs(res, ps.SetLightCurve)
            fsave = os.path.expanduser(fsave)
            if ps.lcf.curves_save(curves, fsave):
                print("Magnitudes of {} have been saved to {}".format(curves.Name, fsave))
            else:
                logger.error("Error with Magnitudes saved to {}".format(curves.Name, fsave))

    else:
        plt.show()
    # plt.show(block=False)

    if not is_save_mags and fsave is not None:
        if fsave == 1:
            fsave = "ubv_obs"
        fsave = os.path.expanduser(os.path.splitext(fsave)[0] + '.pdf')
        print("Save plot to %s " % fsave)
        fig.savefig(fsave, bbox_inches='tight')


if __name__ == '__main__':
    main()
