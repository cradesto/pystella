#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
from os.path import dirname, join, abspath
import logging

import matplotlib

# matplotlib.use("Agg")
from matplotlib import gridspec

from pystella.model import sn_swd
from pystella.model.stella import Stella
from pystella.rf import light_curve_plot as lcp

matplotlib.rcParams['backend'] = "Qt4Agg"
import matplotlib.pyplot as plt
import pystella.model.sn_eve as sneve

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))
logging.basicConfig()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def usage():
    print("Usage:")
    print("  swd.py [params]")
    print("  -t <bands:shift>: string, default: 1:10:50.")
    print("  -i <model name>.  Example: cat_R450_M15_Ni007_E7")
    print("  -p <model directory>, default: ./")
    print("  -h  print usage")
    print("   --- ")


def main():
    parser = argparse.ArgumentParser(description='Process Stella Shock Wave Details.')

    parser.add_argument('-i', '--input',
                        required=False,
                        dest="name",
                        help="Model name, example: cat_R450_M15_Ni007")

    parser.add_argument('-p', '--path',
                        required=False,
                        type=str,
                        default='./',
                        dest="path",
                        help="Model directory")

    parser.add_argument('--rnorm',
                        required=False,
                        default='m',
                        dest="rnorm",
                        help="Radius normalization, example: 'm' or 'sun' or 1e13")

    parser.add_argument('--vnorm',
                        required=False,
                        type=float,
                        default=1e8,
                        dest="vnorm",
                        help="Velocity normalization, example: 1e8")

    parser.add_argument('--lumnorm',
                        required=False,
                        type=float,
                        default=1e40,
                        dest="lumnorm",
                        help="Luminously normalization, example: 1e40")

    d = '1:4:15:65'
    parser.add_argument('-t', '--time',
                        required=False,
                        type=str,
                        default=d,
                        dest="times",
                        help="Plot shock wave snap for selected time moments. Default: {0}".format('2:10:50'))

    args, unknownargs = parser.parse_known_args()

    name = None
    if len(unknownargs) > 0:
        path, name = os.path.split(unknownargs[0])
        path = os.path.expanduser(path)
        name = name.replace('.swd', '')
    else:
        if args.path:
            path = args.path
        else:
            path = ROOT_DIRECTORY
        if args.name is not None:
            name = os.path.splitext(os.path.basename(args.name))[0]

    if name is None:
        # logger.error(" No data. Use key '-i' ")
        parser.print_help()
        sys.exit(2)

    times = list(map(float, args.times.split(':')))
    print("Run swd-model %s in %s for %s moments" % (name, path, args.times))
    stella = Stella(name, path=path)
    swd = stella.get_swd().load()

    lcp.plot_shock_details(swd, times=times, vnorm=args.vnorm, rnorm=args.rnorm,
                           lumnorm=args.lumnorm, is_legend=False)

    # times = [1., 2., 3.]
    #
    # fig = plt.figure(num=None, figsize=(8, 3*len(times)), dpi=100, facecolor='w', edgecolor='k')
    # gs1 = gridspec.GridSpec(len(times), 1)
    # plt.matplotlib.rcParams.update({'font.size': 14})
    #
    # i = 0
    # for t in times:
    #     ax = fig.add_subplot(gs1[i, 0])
    #     b = swd.block_nearest(t)
    #     sn_swd.plot_swd(ax, b, is_xlabel=i == len(times) - 1, rnorm=args.rnorm, vnorm=args.rnorm, lumnorm=args.lumnorm, is_legend=False)
    #     # ax.grid()
    #     i += 1
    #
    # plt.show()


if __name__ == '__main__':
    main()
