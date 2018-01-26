#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import logging
import os
import sys
from itertools import cycle
from os.path import dirname

import matplotlib.lines as mlines
import matplotlib.pyplot as plt

import pystella.model.sn_eve as sneve
import pystella.rf.light_curve_plot as lcp

# matplotlib.use("Agg")
# import matplotlib
# matplotlib.rcParams['backend'] = "TkAgg"
# matplotlib.rcParams['backend'] = "Qt5Agg"
# matplotlib.rcParams['backend'] = "Qt4Agg"

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))
logging.basicConfig()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

markers = {u'x': u'x', u'd': u'thin_diamond',
           u'+': u'plus', u'*': u'star', u'o': u'circle', u'v': u'triangle_down', u'<': u'triangle_left'}
markers_style = list(markers.keys())
lines_style = lcp.linestyles


def get_parser():
    parser = argparse.ArgumentParser(description='Process PreSN configuration.')

    parser.add_argument('-r', '--rho', nargs="?",
                        required=False,
                        const=True,
                        dest="rho",
                        metavar="<r OR m>",
                        help="Plot Rho-figure")

    parser.add_argument('-x',
                        required=False,
                        dest="x",
                        default='m',
                        metavar="<m OR r OR lgR>",
                        help="Setup abscissa: Rho(r) or Rho(m)")

    parser.add_argument('-s', '--save',
                        required=False,
                        type=bool,
                        default=False,
                        dest="is_save_plot",
                        help="save plot to pdf-file, default: False")

    # parser.add_argument('-c', '--chem',
    #                     required=False,
    #                     type=bool,
    #                     default=True,
    #                     dest="is_chem",
    #                     help="Show chemical composition, default: True")
    parser.add_argument('--chem', dest='is_chem', action='store_true', help="Show chemical composition [default].")
    parser.add_argument('--no-chem', dest='is_chem', action='store_false', help="Not show chemical composition")
    parser.set_defaults(is_chem=True)

    parser.add_argument('-i', '--input', action='append', nargs=1,
                        metavar='model name', help='Key -i can be used multiple times')

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
    parser.add_argument('-w', '--write',
                        action='store_const',
                        const=True,
                        dest="is_write",
                        help="To write the data to txt-file.")
    return parser


def main():
    parser = get_parser()
    args, unknownargs = parser.parse_known_args()
    markersize = 4

    if args.path:
        pathDef = os.path.expanduser(args.path)
    else:
        pathDef = os.getcwd()

    # if args.elements:
    elements = args.elements.split('-')
    for e in elements:
        if e not in sneve.eve_elements:
            logger.error('No such element: ' + e)
            sys.exit(2)

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

    if len(names) > 1:  # special case
        markers_cycler = cycle(markers_style)
        lines_cycler = cycle(lines_style)
    else:
        markers_cycler = cycle([None])
        lines_cycler = cycle(['-'])

    ax = None
    ax2 = None
    handles_nm = []
    for nm in names:
        print("Run eve-model %s" % nm)
        path, name = os.path.split(nm)
        if len(path) == 0:
            path = pathDef
        name = name.replace('.rho', '')  # remove extension
        # print("Run eve-model %s in %s" % (name, path))

        rho_file = os.path.join(path, name + '.rho')
        eve = sneve.load_rho(rho_file)

        marker = next(markers_cycler)
        ls = next(lines_cycler)

        if args.is_write:
            fname = os.path.join(path, name)
            f = fname+'.eve.abn'
            if eve.write_abn(f):
                print(" abn have been saved to {}".format(f))
            else:
                print("Error with abn saving to {}".format(f))
            if eve.write_hyd(f):
                print(" hyd have been saved to {}".format(f))
            else:
                print("Error with hyd saving to {}".format(f))

        if args.is_chem:
            # print "Plot eve-model %s" % name
            ax = eve.plot_chem(elements=elements, ax=ax, x=args.x, ylim=(1e-8, 1.), marker=marker,
                               markersize=markersize)

        if args.rho:
            if args.is_chem:
                if ax2 is None:
                    ax2 = ax.twinx()
                    ax2.set_ylabel(r'$\rho, [g/cm^3]$ ')
            else:
                ax2 = ax
            ax = eve.plot_rho(x=args.x, ax=ax2, ls=ls)
        else:
            ls = 'None'

        handle = mlines.Line2D([], [], color='black', marker=marker,
                               markersize=markersize, label=name, linestyle=ls)
        handles_nm.append(handle)

    # plt.legend()
    # ax.legend(prop={'size': 9}, loc=3, ncol=4, fancybox=False, frameon=True)

    if len(names) > 1:
        if ax2 is None:
            ax2 = ax.twinx()
        ax2.legend(handles=handles_nm, loc=4, fancybox=False, frameon=False)
    # plt.gca().add_artist(lgd)

    # plt.grid()
    plt.show()

    if args.is_save_plot:
        if args.rho:
            fsave = os.path.join(os.path.expanduser('~/'), 'rho_%s.pdf' % names[0])
        else:
            fsave = os.path.join(os.path.expanduser('~/'), 'chem_%s.pdf' % names[0])
        logger.info(" Save plot to %s " % fsave)
        ax.get_figure().savefig(fsave, bbox_inches='tight', format='pdf')


if __name__ == '__main__':
    main()
