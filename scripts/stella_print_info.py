#!/usr/bin/python
# -*- coding: utf-8 -*-
import fnmatch
import getopt
import os
import sys
from os import listdir
from os.path import isfile, join, dirname

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))
sys.path.append(ROOT_DIRECTORY)

from pystella.model.stella import Stella

__author__ = 'bakl'


def usage():
    print "Print info from stella res-files"
    print "Usage:"
    print "  %s [params]" % __file__
    print "  -p <directory with models>"
    print "  -f <format: tex, simple>"
    print "  -m <mask: reg expression>"


def main():
    dir_target = './'
    ext_def = '.res'
    mask = None
    is_tex = False
    is_simple = False

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
            if frm == 'tex':
                is_tex = True
            if frm == 'simple':
                is_simple = True
        if opt == '-m':
                mask = str(arg).strip()
        elif opt == '-h':
            usage()
            sys.exit(2)

    files = [f for f in listdir(dir_target) if isfile(join(dir_target, f)) and f.endswith(ext_def)]
    if len(files) == 0:
        print "No res-files in  directory: " + dir_target
        sys.exit(2)

    for f in files:
        # print 'Read: %s' % f
        fname, ext = os.path.splitext(f)
        if mask is not None:
            if not fnmatch.fnmatch(fname, mask):
                continue
        stella = Stella(fname, path=dir_target)
        info = stella.get_res().Info
        if is_tex:
            info.print_tex(lend=' \hline')
        elif is_simple:
            print "%s" % fname
        else:
            info.show()


if __name__ == '__main__':
    main()
