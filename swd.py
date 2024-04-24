#!/usr/bin/env python3

import logging
import os
from itertools import cycle

import numpy as np
import pystella as ps

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


def make_cartoon(swd, times, vnorm, axeX, lumnorm, is_legend, fout=None):
    import subprocess
    import matplotlib.pyplot as plt

    time = np.ma.masked_outside(swd.Times, times[0], times[-1])
    print("Create images from {}  to {} days. Use -t to set other time ranges ".format(times[0], times[-1]))

    # time = np.exp(np.linspace(np.log(times[0]), np.log(times[-1]), 50))
    for i, t in enumerate(time.compressed()):
        fig = ps.lcp.plot_shock_details(swd, times=[t], vnorm=vnorm, axeX=axeX,
                                        lumnorm=lumnorm, is_legend=is_legend)
        fsave = os.path.expanduser("img{0}.{1:04d}.png".format(swd.Name, i))
        print("Save plot to {0} at t={1}".format(fsave, t))
        fig.savefig(fsave, bbox_inches='tight')
        plt.close(fig)
        plt.clf()
    if fout is None:
        fout = 'video_{0}.mp4'.format(swd.Name)
    fout = os.path.abspath(os.path.expanduser(fout))
    print("Convert images to movie: {0}".format(fout))
    # subprocess.call("convert -delay 1x2 -quality 100 img*.png {0}".format(fout), shell=True)
    cmd = " ".join( ("ffmpeg -y -framerate 1  -i img{0}.%04d.png ".format(swd.Name),
           '-r 30 ',
        #    '-vf -c:v libx264 -r 30 -pix_fmt yuv420p ',
           " {}".format(fout))
           )
    print("cmd:  {}".format(cmd))
    subprocess.call(cmd, shell=True)
    print("Done")
    # os.subprocess.call("ffmpeg -f image2 -r 1/5 -i img%04d.png -vcodec mpeg4 -y {}".format(fout), shell=False)


def get_parser(times='1:4:15:65', bnames='U:B:V:R', tau_ph=2. / 3):
    import argparse
    from argparse import RawTextHelpFormatter

    parser = argparse.ArgumentParser(description='Process Stella Shock Wave Details.',
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-b', '--band',
                        nargs='?',
                        required=False,
                        # default=bnames,
                        const=bnames,
                        type=str,
                        dest="bnames",
                        help="-b <bands>: string. If set only -b BNAMES is {}".format(bnames))
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
    parser.add_argument('--x',
                        required=False,
                        default='lgr',  # 1e14,
                        dest="axeX",
                        help="Radius normalization, example: 'lgr, 'm', 'z' or 'sun' or 1e13")
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
                        default='0.001:11',
                        dest="ylim_par",
                        help="Ylim for the parameter axes. Default: 0:9.9")
    parser.add_argument('-c', action='store_const', dest='constant_value',
                        const='value-to-store',
                        help='Store a constant value')
    parser.add_argument('--uph', action='store_const', dest='is_uph',
                        const=True,
                        help='To compute the photospheric velocity. If one uses -w, then the velocities is saved in a file.')
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
    parser.add_argument('--tau',
                        nargs='?',
                        required=False,
                        const=tau_ph,
                        type=float,
                        dest="tau",
                        help="To show the optical depth at tau. It is needed the tau-file at the same directory as swd."
                             "  Default: tau={:.3f}".
                        format(tau_ph))
    # parser.add_argument('--tau',
    #                     required=False,
    #                     type=str,
    #                     default=None,
    #                     dest="ftau",
    #                     help='Set TauFile to plot the photospheric data via +RhoFile+Bands. Default band: B '
    #                          'Example: -i TauFile+U:B:V')
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
        stella = ps.Stella(name, path=path)
        swd = stella.get_swd().load()

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
            make_cartoon(swd, times, vnorm=args.vnorm, axeX=args.axeX,
                         lumnorm=args.lumnorm, is_legend=is_legend)
            return
        else:
            # ls = next(ls_cycle) # skip solid
            fig, dic_axes = ps.lcp.plot_shock_details(swd, times=times,
                                                      vnorm=args.vnorm, axeX=args.axeX, tnorm=args.tnorm,
                                                      lumnorm=args.lumnorm, is_legend=is_legend, is_axes=True,
                                                      ylim_par=ylim_par,
                                                      dic_axes=dic_axes, ls=next(ls_cycle))
            if args.frho is not None:
                ps.lcp.plot_swd_chem(dic_axes, args.frho, stella.Path)

            if args.tau:
                ps.Band.load_settings()
                # Set band names
                bnames = ('B',)
                if args.bnames:
                    ps.Band.load_settings()
                    bnames = []
                    for bname in args.bnames.split(':'):
                        if not ps.band.is_exist(bname):
                            print('No such band: ' + bname)
                            parser.print_help()
                            sys.exit(2)
                        bnames.append(bname)

                ps.lcp.plot_swd_tau(dic_axes, stella, times=times, bnames=bnames, tau_ph=args.tau,
                                    is_obs_time=False, marker=next(marker_cycle),
                                    vnorm=args.vnorm, tnorm=args.tnorm)

            if args.is_save:
                fsave = os.path.expanduser("~/swd_{0}_t{1}.pdf".format(name, str.replace(args.times, ':', '-')))
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
