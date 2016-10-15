#!/usr/bin/python
# -*- coding: utf-8 -*-


import getopt
import os
import sys
from os.path import dirname
import logging

import pystella.model.sn_eve as sneve

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))
logging.basicConfig()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def usage():
    elements = sneve.eve_elements
    print "Usage:"
    print "  eve.py [params]"
    print "  -b <elements>: string, default: 'H-He-C-O-Si-Fe-Ni-Ni56'.\n" \
          "     Available: " + '-'.join(elements)
    print "  -i <model name>.  Example: cat_R450_M15_Ni007"
    print "  -p <model directory>, default: ./"
    print "  -s  save plot to pdf-file."
    print "  -h  print usage"
    print "  "
    print "  ipython: sneve.StellaEve(name=name, path=path).rho_load().plot_chem(ylim=[-6, 0.]) \n"


def main(name=''):
    is_save_plot = False
    path = ''

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hsp:i:b:")
    except getopt.GetoptError as err:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(opts) == 0:
        usage()
        sys.exit(2)

    if not name:
        for opt, arg in opts:
            if opt == '-i':
                path = ROOT_DIRECTORY
                name = os.path.splitext(os.path.basename(str(arg)))[0]
                break

    # elements = sneve.eve_elements
    elements = ['H', 'He', 'C', 'O', 'Si', 'Fe', 'Ni', 'Ni56']

    for opt, arg in opts:
        if opt == '-b':
            elements = str(arg).split('-')
            for e in elements:
                if e not in sneve.eve_elements:
                    logger.error('No such element: ' + e)
                    sys.exit(2)
            continue
        if opt == '-s':
            is_save_plot = True
            continue
        if opt == '-p':
            path = os.path.expanduser(str(arg))
            if not (os.path.isdir(path) and os.path.exists(path)):
                logger.error("No such directory: " + path)
                sys.exit(2)
            continue
        elif opt == '-h':
            usage()
            sys.exit(2)

    print "Run eve-model %s in %s" % (name, path)

    eve = sneve.StellaEve(name, path=path)
    if eve.is_rho_data:
        eve.load()
        eve.plot_chem(elements=elements, ylim=[-6, 0.], is_save=is_save_plot)
    else:
        eve.show_info()
        logger.error("No chemical data for " + eve.rho_file)

    # sneve.StellaEve(name=name, path=path).rho_load().plot_chem(ylim=[-6, 0.])

if __name__ == '__main__':
    main()
