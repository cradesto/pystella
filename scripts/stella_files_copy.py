#!/usr/bin/python
# -*- coding: utf-8 -*-

import getopt
import os
import sys
from os import listdir
from os.path import isfile, join
from shutil import copyfile

__author__ = 'bakl'


def usage():
    print "Copy files with the same names as *.tt"
    print "Usage:"
    print "  %s [params]" % __file__
    print "  -p <directory with models>"
    print "  -t <target directory>, default ./"
    print "  -e <extension to copy>, default: swd"
    print "  -o <pattern extension>, default: tt"


def main():
    dir_target = './'
    ext_def = '.tt'
    ext = '.swd'
    d = '/home/bakl/Sn/Release/seb_git/run/strad/'

    try:
        opts, args = getopt.getopt(sys.argv[1:], "e:h:p:o:t:")
    except getopt.GetoptError as err:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(opts) == 0:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-e':
            ext = str(arg)
            continue
        if opt == '-o':
            ext_def = str(arg)
            continue
        if opt == '-p':
            d = os.path.expanduser(str(arg))
            if not (os.path.isdir(d) and os.path.exists(d)):
                print "No such directory: " + d
                sys.exit(2)
            continue
        if opt == '-t':
            dir_target = os.path.expanduser(str(arg))
            if not (os.path.isdir(dir_target) and os.path.exists(dir_target)):
                print "No target directory: " + dir_target
                sys.exit(2)
            continue
        elif opt == '-h':
            usage()
            sys.exit(2)

    files = [f for f in listdir(dir_target) if isfile(join(dir_target, f)) and f.endswith(ext_def)]
    for f in files:
        s = f.replace(ext_def, ext)
        print 'Copy: %s to %s' % (join(d, s), join(dir_target, s))
        copyfile(join(d, s), os.path.abspath(join(dir_target, s)))


if __name__ == '__main__':
    main()
