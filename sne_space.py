#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import logging
import os
import sys

import matplotlib.pyplot as plt

import pystella as ps
# from pystella.model.SneSpace import SneSpace
# from pystella.rf import light_curve_func as lcf
# from pystella.rf import light_curve_plot as lcp

__author__ = 'bakl'

logging.basicConfig()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_parser():
    parser = argparse.ArgumentParser(description='Plot SN using json-file from https://sne.')
    parser.add_argument('-i', '--input',
                        required=False,
                        dest="fname",
                        help="File name with json-data")
    parser.add_argument('--lumnorm',
                        required=False,
                        type=float,
                        default=1e40,
                        dest="lumnorm",
                        help="Luminously normalization, example: 1e40")
    parser.add_argument('-s', '--save',
                        action='store_const',
                        const=True,
                        dest="is_save",
                        help="To save plot to pdf-file. Default: False".format('2:10:50'))
    parser.add_argument('-w', '--write',
                        action='store_const',
                        const=True,
                        dest="is_write",
                        help="To write the data to txt-file.")
    return parser


def main():
    # fname = "~/Sn/Release/svn_kepler/stella/branches/lucy/run/res/sncurve/sn1999em/SN1999em.json"
    fname_saved = None
    parser = get_parser()
    args, unknownargs = parser.parse_known_args()

    if len(unknownargs) > 0:
        fname = unknownargs[0]
    else:
        if args.fname:
            fname = args.fname
        else:
            # logger.error(" No data. Use key '-i' ")
            parser.print_help()
            sys.exit(2)
    fname = os.path.expanduser(fname)
    if not os.path.isfile(fname):
        logger.error(" No such file: {}".format(fname))
        parser.print_help()
        sys.exit(2)

    sn = ps.SneSpace()
    print("Load {0}".format(fname))
    sn.load(fname)
    print("{} is loaded.".format(sn.Name))

    curves = sn.to_curves()
    print("Curves {} is computed. There are {} bands.".format(sn.Name, '-'.join(curves.BandNames)))
    if args.is_write:
        if fname_saved is None:
            path = os.getcwd()
            fname_saved = os.path.join(path, curves.Name)
            fname_saved = '{}_{}.ubv'.format(fname_saved, '-'.join(sorted(curves.BandNames)))
        if ps.light_curve_func.curves_save(curves, fname_saved):
            print("Magnitudes of {} have been saved to {}".format(curves.Name, fname_saved))
        else:
            print("Error: magnitudes have not been saved [{}]".format(fname_saved))
        sys.exit(3)

    # plotting
    ax = ps.light_curve_plot.curves_plot(curves, title='Photometry', is_line=False)

    # Spectra
    if False:
        serial_spec = sn.to_spectra()
        print("Spectra {} is computed: {} times points.".format(sn.Name, serial_spec.Length))
        # curves_spectra = serial_spec.flux_to_curves(['V'], d=sn.lumdist*phys.pc)
        mu = 1.
        # mu = (10./sn.lumdist)**2
        curves_spectra = serial_spec.flux_to_curves(curves.BandNames, d=None, magnification=mu)
        ps.light_curve_plot.curves_plot(curves_spectra, title='From spectra', is_line=False)
        # lcp.curves_plot(curves_spectra, ax=ax, is_line=True)

    plt.show()


if __name__ == '__main__':
    main()
