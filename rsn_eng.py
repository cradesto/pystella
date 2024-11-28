#!/usr/bin/env python3

import logging
import os
import sys
from itertools import cycle

import numpy as np
import pystella as ps
        
try:
    import matplotlib.pyplot as plt
except ImportError as ex:
    # import traceback
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fn = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fn, exc_tb.tb_lineno, ex)
    print('  Probably, you should install module: {}'.format('matplotlib'))
    #    print(ex)
    plt = None
    pass   

__author__ = 'bakl'

# ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))
logging.basicConfig()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.ERROR)

def make_cartoon(eng, times, engnorm, is_legend, fout=None, **kwargs):
    import subprocess
    import matplotlib.pyplot as plt

    time = np.ma.masked_outside(eng.Times, times[0], times[-1])
    # time = np.exp(np.linspace(np.log(times[0]), np.log(times[-1]), 50))
    for i, t in enumerate(time.compressed()):
        fig = plot_eng_details(eng, times=[t], engnorm=engnorm, is_legend=is_legend, **kwargs)
        fsave = os.path.expanduser("img{0}{1:04d}.png".format(eng.Name, i))
        print("Save plot to {0} at t={1}".format(fsave, t))
        fig.savefig(fsave, bbox_inches='tight')
        plt.close(fig)
        plt.clf()
    if fout is None:
        fout = 'video_{0}.mp4'.format(eng.Name)
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
    parser.add_argument('--engnorm',
                        required=False,
                        type=float,
                        default=1.,
                        dest="engnorm",
                        help="Energy normalization, example: 1.")
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
    parser.add_argument('--mult', action='store_const', dest='is_mult',
                        const=True,
                        help='To make cartoon for evolution. Use: '
                             + "ffmpeg -framerate 4 -i img%%04d.png -c:v libx264 -r 30 out.mp4 "
                             + r'convert -quality 100 img*.png out.mp4 ')
    # .format("ffmpeg -framerate 4 -i img%%04d.png -c:v libx264 -r 30 out.mp4 "))
    parser.add_argument('-s', '--save',
                        action='store_const',
                        const=True,
                        dest="is_save",
                        help="To save the result plot to pdf-file. Format: eng_[name]_t[times].pdf.")
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


def plot_eng_details(eng, times, **kwargs):
    from pystella.model import supr_eng

    is_legend = kwargs.get('is_legend', False)
    # ylim_par = kwargs.get('ylim_par', (0.001, 11))
    font_size = kwargs.get('font_size', 12)
    # is_grid = kwargs.get('is_grid', False)
    is_adjust = kwargs.get('is_adjust', True)
    is_axes = kwargs.get('is_axes', False)

    is_ax_old = False
    xlim = None
    ylim_eng = None
    nrow = len(times)
    ncol = 1

    axes1 = []
    fig = plt.figure(figsize=(12, nrow * 4))

    plt.matplotlib.rcParams.update({'font.size': font_size})
    # plot radius column
    for i, t in enumerate(times):
        axeng = fig.add_subplot(nrow, ncol, ncol * i + 1, label='radius {}'.format(i))
        axes1.append(axeng)
        legmask = supr_eng.LEGEND_MASK_None
        if is_legend and i == 0:
            legmask = supr_eng.LEGEND_MASK_Eng
        # plot eng(radius)
        b = eng.block_nearest(t)
        axeng = supr_eng.plot_eng(axeng, b, name=eng.Name, is_xlabel=(i == len(times) - 1),
                                       legmask=legmask, is_yrlabel=False, text_posy=0.88, **kwargs)
        if not is_ax_old:
            x = axeng.get_xlim()
            if xlim is None:
                xlim = x
            else:
                xlim = (min(x[0], xlim[0]), max(x[1], xlim[1]))
            y = axeng.get_ylim()
            if ylim_eng is None:
                ylim_eng = y
            else:
                ylim_eng = (min(y[0], ylim_eng[0]), max(y[1], ylim_eng[1]))
            ticks_on(axeng)   

    # Set limits
    for i, ax in enumerate(axes1):
        ax.set_xlim(xlim)
        ax.set_ylim(ylim_eng)
        ax.set_yscale('log')
        # remove labels between subplots
        if not i == len(times) - 1:
            plt.setp(ax.get_xticklabels(), visible=False)

    if is_adjust:
        fig.subplots_adjust(wspace=0., hspace=0.)

    if is_axes:
        return fig
    return fig


def main():
    # import sys
    # try:
    #     import matplotlib.pyplot as plt
    # except ImportError:
    #     plt = None

    # if True:
    #     try:
    #         sys.path.append(os.path.expanduser('~/Sn/Release/python/smplotlib/src'))
    #         import smplotlib
    #     except ImportError as ex:
    #         print('  Probably, you should install module: {}'.format('smplotlib'))
     
    # try:
    #     import seaborn as sns
    #     # sns.set()
    #     sns.set_style("ticks")
    #
    # except ImportError:
    #     sns = None
    fsave = None
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
        name = name.replace('.eng', '')  # remove extension

        print("Run reng-model %s %s for %s moments" % (path, name, args.times))
        supr = ps.Supremna(name, path=path)
        eng = supr.get_eng().load()

        if args.is_mult:
            make_cartoon(eng, times, engnorm=args.engnorm, is_legend=is_legend)
        else:
            # ls = next(ls_cycle) # skip solid
            fig = plot_eng_details(eng, times=times, engnorm=args.engnorm, is_legend=is_legend, is_axes=True,
                                                    ylim_par=ylim_par, ls=next(ls_cycle))
            if args.is_save:
                fsave = os.path.expanduser("~/reng_{0}_t{1}.pdf".format(name, str.replace(args.times, ':', '-')))
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
