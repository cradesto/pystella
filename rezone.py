#!/usr/bin/python
# -*- coding: utf-8 -*-
from scipy.optimize import fmin
from scipy import interpolate
import os
import sys
import getopt
from os.path import isfile, join, dirname
import csv
import numpy as np

from matplotlib import gridspec
import matplotlib.pyplot as plt

from pystella.rf import band, spectrum
from pystella.rf.star import Star
from pystella.model.stella import Stella
import pystella.util.rf as rf
from pystella.util.phys_var import phys

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))


def usage():
    bands = band.band_get_names().keys()
    print "\n Create hyd- abn-files from res-file."
    print "Usage:"
    print "  rezone.py [params]"
    print "  -i <model name>.  Example: cat_R450_M15_Ni007_E7"
    print "  -n  new zon number, default: 100"
    print "  -p <model path(directory)>, default: ./"
    print "  -s  silence mode: no info, no plot"
    print "  -t  time moment, default: 10"
    print "  -h  print usage"


def main(name=False):
    is_silence = False
    nzon = 100
    t = 1  # days
    ext_file = '.res'

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hsp:n:i:t:")
    except getopt.GetoptError as err:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if not name:
        if len(opts) == 0:
            usage()
            sys.exit(2)
        for opt, arg in opts:
            if opt == '-i':
                path = ROOT_DIRECTORY
                name = str(arg)
                break

    for opt, arg in opts:
        if opt == '-n':
            nzon = int(arg)
            continue
        if opt == '-s':
            is_silence = True
            continue
        if opt == '-t':
            t = float(arg)
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

    model = Stella(name, path=path)

    if not model.is_res_data:
        print "There are no  %s in the directory: %s " % (name, path)

    res = model.get_res()
    block = res.read_at_time(time=t)

    print("Len(block) = %i " % len(block))
    if not is_silence:
        plt.plot(block['M'], block['V8'])
        plt.show()
    # write_data(res, path, fname=name)



if __name__ == '__main__':
    main()
