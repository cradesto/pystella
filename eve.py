#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
from os.path import dirname
import logging

import matplotlib

# matplotlib.use("Agg")
matplotlib.rcParams['backend'] = "Qt4Agg"
import matplotlib.pyplot as plt
import pystella.model.sn_eve as sneve

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))
logging.basicConfig()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def main():
    parser = argparse.ArgumentParser(description='Process PreSN configuration.')

    parser.add_argument('-r', '--rho',  nargs="?",
                        required=False,
                        const=True,
                        dest="rho",
                        metavar="<r OR m>",
                        help="Plot Rho-figure instead Composition-figure  ")

    parser.add_argument('-x',
                        required=False,
                        dest="x",
                        metavar="<r OR m OR lgR>",
                        help="setup abscissa: Rho(r) or Rho(m)")

    parser.add_argument('-s', '--save',
                        required=False,
                        type=bool,
                        default=False,
                        dest="is_save_plot",
                        help="save plot to pdf-file")

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

    parser.add_argument('-e', '--elements',
                        required=False,
                        type=str,
                        default='H-He-C-O-Si-Fe-Ni-Ni56',
                        dest="elements",
                        help="Elements directory. \n   Available: {0}".format('-'.join(sneve.eve_elements)))

    args, unknownargs = parser.parse_known_args()

    name = None
    if len(unknownargs) > 0:
        path, name = os.path.split(unknownargs[0])
        path = os.path.expanduser(path)
        name = name.replace('.rho', '')
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

    print("Run eve-model %s in %s" % (name, path))
    rho_file = os.path.join(path, name + '.rho')
    eve = sneve.load_rho(rho_file)
    if args.rho:
        ax = eve.plot_rho(x=args.x, is_save=args.is_save_plot)
    else:
        # print "Plot eve-model %s" % name
        elements = args.elements.split('-')
        for e in elements:
            if e not in sneve.eve_elements:
                logger.error('No such element: ' + e)
                sys.exit(2)
        ax = eve.plot_chem(elements=elements, x=args.x, ylim=(1e-8, 1.), is_save=args.is_save_plot)
    ax.grid()
    plt.show()


if __name__ == '__main__':
    main()
