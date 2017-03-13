#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import logging
import os
import sys
from os.path import dirname

import matplotlib

from pystella.model.stella import Stella
from pystella.rf import light_curve_plot as lcp

# matplotlib.rcParams['backend'] = "TkAgg"
# matplotlib.rcParams['backend'] = "Qt4Agg"
import matplotlib.pyplot as plt

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))
logging.basicConfig()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


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
                        default='lgr',  # 1e14,
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

    parser.add_argument('-c', action='store_const', dest='constant_value',
                        const='value-to-store',
                        help='Store a constant value')

    parser.add_argument('-s', '--save',
                        action='store_const',
                        const=True,
                        dest="is_save",
                        help="To save plot to pdf-file. Default: False".format('2:10:50'))

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
    print("Run swd-model %s %s for %s moments" % (path, name, args.times))
    stella = Stella(name, path=path)
    swd = stella.get_swd().load()

    fig = lcp.plot_shock_details(swd, times=times, vnorm=args.vnorm, rnorm=args.rnorm,
                                 lumnorm=args.lumnorm, is_legend=True)

    plt.show()
    if args.is_save:
        fsave = os.path.expanduser("~/swd_{0}_t{1}.pdf".format(name, str.replace(args.times, ':', '-')))
        print("Save plot to {0}".format(fsave))
        fig.savefig(fsave, bbox_inches='tight')


if __name__ == '__main__':
    main()
