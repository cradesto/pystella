#!/usr/bin/python
# -*- coding: utf-8 -*-

import getopt
import os
import sys
from os import listdir
from os.path import isfile
from shutil import copyfile

__author__ = 'bakl'

stl_exts = ('tt', 'res', 'ph', 'swd')


def simplecopy(dir_src, dir_target, ext_ptn, ext):
    enew = '.' + ext.replace('.', '')
    eold = '.' + ext_ptn.replace('.', '')
    files = [f for f in listdir(dir_src) if isfile(os.path.join(dir_src, f)) and f.endswith(eold)]
    for f in files:
        filename, file_extension = os.path.splitext(f)
        s = filename+enew
        print 'Copy: %s to %s' % (os.path.join(dir_src, s), os.path.join(dir_target, s))
        copyfile(os.path.join(dir_src, s), os.path.abspath(os.path.join(dir_target, s)))


def usage():
    print "Copy files with the same names as *.tt [see key: -o]"
    print "Usage:"
    print "  %s [params]" % __file__
    print "  -p <source directory with models>"
    print "  -t <target directory>, default ./"
    print "  -e <extension to copy>, default: swd or %s" % '-'.join(stl_exts)
    print "  -m <cat>. if mode is 'cat', then copy %s-files to target.  Default: simple copy" % '-'.join(stl_exts)
    print "  -o <pattern extension>. Default: tt"


def catcopy(dir_src, dir_target, ext='tt', exts=stl_exts):
    if not (os.path.isdir(dir_src) and os.path.exists(dir_src)):
        print "No source directory: " + dir_src
        sys.exit(2)
    if not (os.path.isdir(dir_target) and os.path.exists(dir_target)):
        print "No target directory: " + dir_target
        sys.exit(2)

    for e in exts:
        simplecopy(dir_src, dir_target, ext, e)


def main():
    dir_target = './'
    ext_ptn = 'tt'
    exts = stl_exts
    ext = 'swd'
    dir_src = '/home/bakl/Sn/Release/seb_git/run/strad/'
    mode = None

    try:
        opts, args = getopt.getopt(sys.argv[1:], "e:h:m:o:p:t:")
    except getopt.GetoptError as err:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(opts) == 0:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-e':
            ext = str.strip(arg)
            exts = map(str.strip, str(arg).split('-'))
            continue
        if opt == '-o':
            ext_ptn = str.strip(arg)
            continue
        if opt == '-m':
            mode = str.strip(arg)
            continue
        if opt == '-p':
            dir_src = os.path.expanduser(str.strip(arg))
            if not (os.path.isdir(dir_src) and os.path.exists(dir_src)):
                print "No such directory: " + dir_src
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

    if mode == 'cat':
        catcopy(dir_src, dir_target, ext_ptn, exts)
    else:
        simplecopy(dir_src, dir_target, ext_ptn, ext)


if __name__ == '__main__':
    main()
