#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import logging
import os
import sys
from os.path import dirname

from pystella.model.stella import Stella, stella_extensions, stella_ls
from pystella.util.string_misc import print_table, list_to_table

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))
logging.basicConfig()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_parser():
    parser = argparse.ArgumentParser(description='Show stella models info.')
    parser.add_argument('-p', '--path',
                        required=False,
                        type=str,
                        default='./',
                        dest="path",
                        help="The directory of stella models")
    parser.add_argument('-m', '--match',
                        required=False,
                        default='*',
                        dest="pattern",
                        help="File pattern")
    parser.add_argument('-l', '--long',
                        action='store_const',
                        const=True,
                        dest="is_long",
                        help="use a long info format")
    return parser


def main():
    parser = get_parser()
    args, unknownargs = parser.parse_known_args()
    if len(unknownargs) > 0:
        path = os.path.expanduser(unknownargs[0])
    else:
        path = args.path

    if path is None:
        parser.print_help()
        sys.exit(2)

    models = stella_ls(path, pattern=args.pattern)

    if len(models) == 0:
        print("No models in {0}".format(path))
        return None

    if args.is_long:
        print("| %40s |  %7s |  %6s | %6s |  %s" % ('Name', 'R', 'M', 'E', 'comment'))
        print("| %40s |  %7s |  %6s | %6s |  %s" % ('-'*40, '-'*7, '-'*6, '-'*6, '-'*7))
        for mdl, exts in models.items():
            stella = Stella(mdl, path=path)
            if stella.is_tt_data:
                info = stella.get_tt().Info
                info.show(comment=' {0}'.format(exts))
            else:
                print("{0} No tt.".format(stella))
    else:
        print('Models: {0} in {1}'.format(len(models), path))
        tbl = list_to_table(models.keys())
        print_table(tbl)
        # print("{0}".format(' '.join(models.keys())))


if __name__ == '__main__':
    main()