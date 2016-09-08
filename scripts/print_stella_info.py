#!/usr/bin/python
# -*- coding: utf-8 -*-

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
    print "Print info from stella files"
    print "Usage:"
    print "  %s [params]" % __file__
    print "  -p <directory with models>"
    print "  -f <format: tex>"


def main():
    dir_target = './'
    ext_def = '.res'
    is_tex = False

    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:f:p:")
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
            if str(arg) == 'tex':
                is_tex = True
        elif opt == '-h':
            usage()
            sys.exit(2)

    files = [f for f in listdir(dir_target) if isfile(join(dir_target, f)) and f.endswith(ext_def)]
    if len(files) == 0:
        print "No res-files in  directory: " + dir_target
        sys.exit(2)

    for f in files:
        # print 'Read: %s' % f
        name, ext = os.path.splitext(f)
        stella = Stella(name, path=dir_target)
        info = stella.get_res().Info
        if is_tex:
            info.print_tex(lend=' \hline')
        else:
            info.show()


if __name__ == '__main__':
    main()
