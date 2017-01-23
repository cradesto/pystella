#!/usr/bin/python
# -*- coding: utf-8 -*-
import fnmatch
import getopt
import os
import sys
from os import listdir
from os.path import isfile, join, dirname
from types import FunctionType

import pandas as pd
import numpy as np

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))
sys.path.append(ROOT_DIRECTORY)

from pystella.model.stella import Stella
from pystella.model.sn_res import StellaResInfo

__author__ = 'bakl'


class PrintDF:
    @staticmethod
    def name(df, columns=None):
        for n in df['name']:
            print n

    @staticmethod
    def tex(df, columns, lend=''):
        frm_header = " \\mbox{{Name}} " + " & \\mbox{{{:s}}} " * (len(columns)-1)
        print frm_header.format(*columns[1:])
        frm = " \\mbox{{{:s}}} " + "&  {:6.2f} " * (len(columns) - 1) + " \\\ {:s} "
        for i, r in df.iterrows():
            # print(r['name'])
            s = frm.format(*([r[n] for n in columns] + [lend]))
            print s

    @staticmethod
    def display(df, columns):
        # print "| %40s |  %7.2f |  %6.2f | %6.2f | %6.2f |" % (self.name, self.R, self.M, self.E, self.Mni)
        frm_header = "{:>40s}" + " | {:>14s}" * (len(columns)-1)
        print frm_header.format(*(['Name']+columns[1:]))
        frm = "{:>40s}" + " | {:>14.5f}" * (len(columns) - 1)
        for i, r in df.iterrows():
            # print(r['name'])
            s = frm.format(*([r[n] for n in columns]))
            print s

    @staticmethod
    def csv(df, columns=None):
        df.to_csv(sys.stdout)

    @staticmethod
    def methods():
        # return dir(PrintDF)
        return [func for func in dir(PrintDF) if callable(getattr(PrintDF, func)) and func != 'methods']
        # return [x for x, y in PrintDF.__dict__.items() if type(y) == FunctionType]


def usage():
    print "Print info from stella res-files"
    print "Usage:"
    print "  %s [params]" % __file__
    print "  -p <directory with models>"
    print "  -f <format: %s >" % PrintDF.methods()
    print "  -m <mask: reg expression>"


def main():
    dir_target = './'
    ext_def = '.res'
    mask = None
    fformat = 'display'

    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:f:m:p:")
    except getopt.GetoptError as err:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(opts) == 0:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-p':
            dir_target = os.path.expanduser(str(arg))
            if not (os.path.isdir(dir_target) and os.path.exists(dir_target)):
                print "No such directory: " + dir_target
                sys.exit(2)
            continue
        if opt == '-f':
            frm = str(arg).strip().lower()
            if frm in PrintDF.methods():
                fformat = frm
            else:
                print "No such format: " + frm
                print "Use: %s " % PrintDF.methods()
                sys.exit(2)
        if opt == '-m':
            mask = str(arg).strip()
        elif opt == '-h':
            usage()
            sys.exit(2)

    files = [f for f in listdir(dir_target) if isfile(join(dir_target, f)) and f.endswith(ext_def)]
    if len(files) == 0:
        print "No res-files in  directory: " + dir_target
        sys.exit(2)

    columns = ['name'] + list(StellaResInfo.Params)
    index = range(len(files))
    # df = pd.DataFrame(index=index, columns=columns)
    df = pd.DataFrame(columns=(columns))

    # get info
    i = 1
    for f in files:
        # print 'Read: %s' % f
        fname, ext = os.path.splitext(f)
        if mask is not None:
            if not fnmatch.fnmatch(fname, mask):
                continue
        stella = Stella(fname, path=dir_target)
        info = stella.get_res().Info

        df.loc[i] = [info.name] + [info.get(k) for k in StellaResInfo.Params]
        i += 1

    # sort info
    df.sort(list(StellaResInfo.Params), ascending=False)
    # df.sort([StellaResInfo.sRinit, StellaResInfo.sMtot, StellaResInfo.sEburst, StellaResInfo.sMni], ascending=False)

    # print info
    invert_op = getattr(PrintDF, fformat, None)
    if callable(invert_op):
        invert_op(df, columns)


if __name__ == '__main__':
    main()
