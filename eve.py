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


def usage():
    elements = sneve.eve_elements
    print "Usage:"
    print "  eve.py [params]"
    print "  -b <elements>: string, default: U-B-V-R-I, for example U-B-V-R-I-u-g-i-r-z-UVW1-UVW2.\n" \
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
                name = str(arg)
                break

    elements = sneve.eve_elements
    for opt, arg in opts:
        if opt == '-b':
            elements = str(arg).split('-')
            for e in elements:
                if not sneve.eve_elements(e):
                    print 'No such element: ' + e
                    sys.exit(2)
            continue
        if opt == '-s':
            is_save_plot = True
            continue
        if opt == '-p':
            path = os.path.expanduser(str(arg))
            if not (os.path.isdir(path) and os.path.exists(path)):
                print "No such directory: " + path
                sys.exit(2)
            continue
        elif opt == '-h':
            usage()
            sys.exit(2)

    print "Plot chemical : %s %s" % (path, name)

    eve = sneve.StellaEve(name, path=path)
    if eve.is_rho_data:
        eve.rho_load()
        eve.plot_chem(elements=elements, ylim=[-6, 0.], is_save=is_save_plot)

    # sneve.StellaEve(name=name, path=path).rho_load().plot_chem(ylim=[-6, 0.])

if __name__ == '__main__':
    main()
