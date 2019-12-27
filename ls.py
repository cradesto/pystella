#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from os.path import dirname

import pystella as ps

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
    parser.add_argument('--pars',
                        required=False,
                        default=':'.join(('R', 'M', 'Mni', 'E')),
                        dest="pars",
                        help="List of parameters. Default: " + ':'.join(('R', 'M', 'Mni', 'E')))
    parser.add_argument('--sep',
                        required=False,
                        default=' |',
                        dest="sep",
                        help="String separator Default: |")
    parser.add_argument('--summary',
                        action='store_const',
                        const=True,
                        dest="is_summary",
                        help="show summary for a directory")
    parser.add_argument('-v', '--verbose',
                        action='store_const',
                        const=True,
                        dest="is_verbose",
                        help="show verbose information")
    return parser


def main():
    # pars = ('R', 'M', 'Mni', 'E')
    parser = get_parser()
    args, unknownargs = parser.parse_known_args()

    if len(unknownargs) > 0 and len(unknownargs[0]) > 0:
        path = os.path.expanduser(unknownargs[0])
        # print("1 path: {} | {} | {}".format(path, len(unknownargs), '-'.join(unknownargs)))
    else:
        path = args.path

    if path is None:
        parser.print_help()
        sys.exit(2)

    models = ps.ls.short(path, pattern=args.pattern)

    if len(models) == 0:
        print("No tt-files in {0}".format(path))
        return None

    pars = args.pars.split(':')
    sep = args.sep

    if args.is_summary:
        ps.ls.summary(models.keys(), path, pars, sep=sep, is_verbose=args.is_verbose)
    elif args.is_long:
        ps.ls.long(models, path, pars, sep=sep, is_verbose=args.is_verbose)
    else:
        print('Models: {0} in {1}'.format(len(models), path))
        tbl = list_to_table(models.keys())
        print_table(tbl)
        # print("{0}".format(' '.join(models.keys())))


if __name__ == '__main__':
    main()
