#!/usr/bin/env python3

import logging
import os
import sys
from itertools import cycle

import numpy as np
import pystella as ps




if True:
    try:
        sys.path.append(os.path.expanduser('~/Sn/Release/python/smplotlib/src'))
        import smplotlib
    except ImportError as ex:
        print('  Probably, you should install module: {}'.format('smplotlib'))
        
try:
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
except ImportError as ex:
    # import traceback
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fn = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fn, exc_tb.tb_lineno, ex)
    print('  Probably, you should install module: {}'.format('matplotlib'))
    #    print(ex)
    plt = None
    gridspec = None
    pass



__author__ = 'bakl'

# ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))
logging.basicConfig()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.ERROR)


def uph_write(dictionary, fname, sep='\t'):
    """
    Save dict to file. Keys are column's names, values are column data
    :param dictionary:
    :type fname: the name of the file
    :param sep: string separator
    """
    import csv

    col = ('time', 'zone', 'R', 'V')
    with open(fname, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=sep, quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['# {:^8s}'.format(x) for x in col])
        for i, (row) in enumerate(zip(*[dictionary[c] for c in col])):
            writer.writerow(['{:8.3e}'.format(x) for x in row])


def plot_uph(uph, vnorm=1.e8, label='', lw=2, fontsize=18, ls="-", color='blue'):
    import matplotlib.pyplot as plt

    # setup plot
    plt.matplotlib.rcParams.update({'font.size': fontsize})
    fig = plt.figure(figsize=(8, 8))

    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('Time [day]')
    ax.set_ylabel(r'Photospheric Velocity [$\times 10^{0}$ km/s]'.format(int(np.log10(vnorm / 1e5))))
    #    ax.set_title(fname_uph)

    x = uph['time']
    y = uph['V'] / vnorm
    ax.plot(x, y, label=label, color=color, ls=ls, linewidth=lw)
    ax.legend()
    return fig


def make_cartoon(swd, times, vnorm, rnorm, lumnorm, is_legend, fout=None):
    import subprocess
    import matplotlib.pyplot as plt

    time = np.ma.masked_outside(swd.Times, times[0], times[-1])
    # time = np.exp(np.linspace(np.log(times[0]), np.log(times[-1]), 50))
    for i, t in enumerate(time.compressed()):
        fig = plot_shock_details(swd, times=[t], vnorm=vnorm, rnorm=rnorm,
                                        lumnorm=lumnorm, is_legend=is_legend)
        fsave = os.path.expanduser("img{0}{1:04d}.png".format(swd.Name, i))
        print("Save plot to {0} at t={1}".format(fsave, t))
        fig.savefig(fsave, bbox_inches='tight')
        plt.close(fig)
        plt.clf()
    if fout is None:
        fout = 'video_{0}.mp4'.format(swd.Name)
    print("Convert images to movie: {0}".format(fout))
    subprocess.call("convert -delay 1x2 -quality 100 img*.png {0}".format(fout), shell=True)
    print("Done")
    # os.subprocess.call("ffmpeg -f image2 -r 1/5 -i img%04d.png -vcodec mpeg4 -y {}".format(fout), shell=False)


def get_parser(times='41:75:150:300'):
    import argparse
    from argparse import RawTextHelpFormatter

    parser = argparse.ArgumentParser(description='Process Supremna Shock Wave Details.',
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--input', action='append', nargs=1,
                        metavar='Model name', help='Key -i can be used multiple times.')
    # parser.add_argument('-i', '--input',
    #                     required=False,
    #                     dest="input",
    #                     help="Model name, example: cat_R450_M15_Ni007")

    parser.add_argument('-p', '--path',
                        required=False,
                        type=str,
                        default='./',
                        dest="path",
                        help="Model directory")
    parser.add_argument('--lumnorm',
                        required=False,
                        type=float,
                        default=1e40,
                        dest="lumnorm",
                        help="Luminously normalization, example: 1e40")
    parser.add_argument('--rnorm',
                        required=False,
                        default='lgr',  # 1e14,
                        dest="rnorm",
                        help="Radius normalization, example: 'm' or 'sun' or 1e13")
    parser.add_argument('--tnorm',
                        nargs='?',
                        required=False,
                        const=None,
                        type=float,
                        # default=1e3,
                        dest="tnorm",
                        help="Temperature normalization, example: 1e3")
    parser.add_argument('--vnorm',
                        required=False,
                        type=float,
                        default=1e8,
                        dest="vnorm",
                        help="Velocity normalization, example: 1e8")
    parser.add_argument('-t', '--time',
                        required=False,
                        type=str,
                        default=times,
                        dest="times",
                        help="Plot shock wave snap for selected time moments. Default: {}".format(times))
    parser.add_argument('--ylim_par',
                        required=False,
                        type=str,
                        default='0.1:14.5',
                        dest="ylim_par",
                        help="Ylim for the parameter axes. Default: 0.1:14.5")
    parser.add_argument('-c', action='store_const', dest='constant_value',
                        const='value-to-store',
                        help='Store a constant value')
    parser.add_argument('--uph', action='store_const', dest='is_uph',
                        const=True,
                        help='To compute the photospheric velocity')
    parser.add_argument('--mult', action='store_const', dest='is_mult',
                        const=True,
                        help='To make cartoon for evolution. Use: '
                             + "ffmpeg -framerate 4 -i img%%04d.png -c:v libx264 -r 30 out.mp4 "
                             + r'convert -quality 100 img*.png out.mp4 ')
    # .format("ffmpeg -framerate 4 -i img%%04d.png -c:v libx264 -r 30 out.mp4 "))
    parser.add_argument('--chem',
                        required=False,
                        type=str,
                        default=None,
                        dest="frho",
                        help='Set RhoFile to plot the chemical composition via +RhoFile+ListElements. '
                             'Example: --chem RhoFile.rho+H:He:O:Fe')
    parser.add_argument('-s', '--save',
                        action='store_const',
                        const=True,
                        dest="is_save",
                        help="To save the result plot to pdf-file. Format: swd_[name]_t[times].pdf.")
    parser.add_argument('-w', '--write',
                        action='store_const',
                        const=True,
                        dest="is_write",
                        help="To write the data to txt-file.")
    return parser


def ticks_on(ax, minor=3, major=6):
    ax.minorticks_on()
    ax.tick_params(direction='in', which='minor', length=minor)
    ax.tick_params(direction='in', which='major', length=major)
    return ax


def plot_shock_details(swd, times, **kwargs):
    from pystella.model import supr_swd

    is_legend = kwargs.get('is_legend', False)
    # ylim_par = kwargs.get('ylim_par', (0.001, 11))
    font_size = kwargs.get('font_size', 12)
    # is_grid = kwargs.get('is_grid', False)
    is_adjust = kwargs.get('is_adjust', True)
    is_axes = kwargs.get('is_axes', False)
    dic_axes = kwargs.get('dic_axes', None)
    is_ax_old = False
    xlim = None
    ylim_rho = None
    nrow = len(times)
    ncol = 2

    axes1 = []
    if dic_axes is None:
        dic_axes = {'r': [], 'm': []}
        fig = plt.figure(figsize=(12, nrow * 4))
        # plt.minorticks_on()
        # fig = plt.figure(num=None, figsize=(12, len(times) * 4), dpi=100, facecolor='w', edgecolor='k')
        # gs1 = gridspec.GridSpec(len(times), 2)
    else:
        fig = (dic_axes['r'][0]['rho']).get_figure()
        is_ax_old = True
        is_adjust = False
        # is_legend = False
        kwargs['is_day'] = False

    plt.matplotlib.rcParams.update({'font.size': font_size})
    # plot radius column
    for i, t in enumerate(times):
        if is_ax_old:
            axrho, axpar = dic_axes['r'][i]['rho'], dic_axes['r'][i]['par']
        else:
            axrho = fig.add_subplot(nrow, ncol, ncol * i + 1, label='radius {}'.format(i))
            axpar = None
        axes1.append(axrho)
        legmask = supr_swd.LEGEND_MASK_None
        if is_legend and i == 0:
            legmask = supr_swd.LEGEND_MASK_Rho
        # plot swd(radius)
        b = swd.block_nearest(t)
        axrho, axpar = supr_swd.plot_swd((axrho, axpar), b, name=swd.Name, is_xlabel=(i == len(times) - 1),
                                       legmask=legmask, is_yrlabel=False, text_posy=0.88,
                                       **kwargs)
        if not is_ax_old:
            x = axrho.get_xlim()
            if xlim is None:
                xlim = x
            else:
                xlim = (min(x[0], xlim[0]), max(x[1], xlim[1]))
            y = axrho.get_ylim()
            if ylim_rho is None:
                ylim_rho = y
            else:
                ylim_rho = (min(y[0], ylim_rho[0]), max(y[1], ylim_rho[1]))
            # axpar.tick_params(direction='in', which='both', length=4)
            ticks_on(axrho)
            ticks_on(axpar)
            dic_axes['r'].append({'itime': i, 't': t, 'rho': axrho, 'par': axpar})

    if 'rnorm' in kwargs:
        kwargs.pop('rnorm')
    axes2 = []
    # Plot mass column
    for i, t in enumerate(times):
        if is_ax_old:
            # ax2 = dic_axes['m'][i]['rho']
            axrho, axpar = dic_axes['m'][i]['rho'], dic_axes['m'][i]['par']
        else:
            axrho = fig.add_subplot(nrow, ncol, ncol * i + 2, label='mass {}'.format(i))
            axrho.tick_params(direction='in', which='minor', length=3)
            axpar = None
        axes2.append(axrho)
        legmask = supr_swd.LEGEND_MASK_None

        if is_legend and i == 0:
            legmask = supr_swd.LEGEND_MASK_Vars
        b = swd.block_nearest(t)
        axrho, axpar = supr_swd.plot_swd((axrho, axpar), b, name=swd.Name, is_xlabel=(i == len(times) - 1),
                                       rnorm='m', legmask=legmask, is_yllabel=False, text_posy=0.88,
                                       **kwargs)
        if not is_ax_old:
            dic_axes['m'].append({'itime': i, 't': t, 'rho': axrho, 'par': axpar})
            ticks_on(axrho)
            axpar.tick_params(direction='in', which='major', length=5)
            ticks_on(axpar)

    # Set limits
    for i, ax in enumerate(axes1):
        ax.set_xlim(xlim)
        ax.set_ylim(ylim_rho)
        # remove labels between subplots
        if not i == len(times) - 1:
            plt.setp(ax.get_xticklabels(), visible=False)

    for i, ax2 in enumerate(axes2):
        ax2.set_ylim(ylim_rho)
        # remove labels between subplots
        if not i == len(times) - 1:
            plt.setp(ax2.get_xticklabels(), visible=False)

    if is_adjust:
        fig.subplots_adjust(wspace=0., hspace=0.)

    # print(len(axes1), len(axes2))
    if is_axes:
        return fig, dic_axes
    return fig


def main():
    import sys
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        plt = None
    # try:
    #     import seaborn as sns
    #     # sns.set()
    #     sns.set_style("ticks")
    #
    # except ImportError:
    #     sns = None
    fsave = None
    dic_axes = None
    ylim_par = None
    is_legend = True
    ls_cycle = cycle(ps.linestyles_extend)
    marker_cycle = cycle(ps.lcp.markers)

    parser = get_parser()
    args, unknownargs = parser.parse_known_args()

    if args.ylim_par:
        ylim_par = ps.str2interval(args.ylim_par, llim=0, rlim=9.9, sep=':')

    if args.path:
        pathDef = os.path.expanduser(args.path)
    else:
        pathDef = os.getcwd()
    # Set model names
    names = []
    if args.input:
        for nm in args.input:
            names.append(nm[0])  # remove extension
    else:
        if len(unknownargs) > 0:
            names.append(unknownargs[0])

    if len(names) == 0:
        # logger.error(" No data. Use key '-i' ")
        parser.print_help()
        sys.exit(2)

    times = list(map(float, args.times.split(':')))

    for i, nm in enumerate(names):
        path, name = os.path.split(nm)
        if len(path) == 0:
            path = pathDef
        name = name.replace('.swd', '')  # remove extension

        print("Run swd-model %s %s for %s moments" % (path, name, args.times))
        supr = ps.Supremna(name, path=path)
        swd = supr.get_swd().load()

        if args.is_uph:
            logger.info(' Compute and print uph')
            duph = swd.params_ph()
            # print uph
            print(duph.keys())
            for row in zip(*duph.values()):
                print(['{:12f}'.format(x) for x in row])
            # save uph
            if args.is_write:
                fsave = os.path.join(path, "{0}.uph".format(name))
                print("Save uph to {0}".format(fsave))
                uph_write(duph, fsave)
            else:
                fig = plot_uph(duph, vnorm=args.vnorm, label=name)
                if args.is_save:
                    fsave = os.path.expanduser("~/uph_{0}.pdf".format(name))
        elif args.is_mult:
            make_cartoon(swd, times, vnorm=args.vnorm, rnorm=args.rnorm,
                         lumnorm=args.lumnorm, is_legend=is_legend)
        else:
            # ls = next(ls_cycle) # skip solid
            fig, dic_axes = plot_shock_details(swd, times=times,
                                                      vnorm=args.vnorm, rnorm=args.rnorm, tnorm=args.tnorm,
                                                      lumnorm=args.lumnorm, is_legend=is_legend, is_axes=True,
                                                      ylim_par=ylim_par,
                                                      dic_axes=dic_axes, ls=next(ls_cycle))
            if args.frho is not None:
                ps.lcp.plot_swd_chem(dic_axes, args.frho, Supremna.Path)

            if args.is_save:
                fsave = os.path.expanduser("~/rswd_{0}_t{1}.pdf".format(name, str.replace(args.times, ':', '-')))
        #  Save the figure or show
    if fsave is not None:
        print("Save plot to {0}".format(fsave))
        fig.savefig(fsave, bbox_inches='tight')
    else:
        # plt.ion()
        plt.show()
        # plt.pause(0.0001)
        # print('')
        # input("===> Hit <return> to quit")


if __name__ == '__main__':
    main()
