#!/usr/bin/env python3
# #!/usr/bin/python3

import argparse
import logging
import os
import sys
from os.path import dirname

import pystella as ps

# from pystella.model.stella import Stella, stella_ls
from pystella.util.string_misc import print_table, list_to_table

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))
logging.basicConfig()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def ls_long(models, path='.', pars=('R', 'M', 'Mni', 'E'), sep=' & ', is_verbose=False):
    """
    Print model list with parameters
        TODO  filter by parameters
    :param is_verbose:
    :param sep:
    :param models:  list of model names
    :param path: working directory. Default: ./
    :param pars: the list of printed parameters
    :return:
    """

    is_dict = type(models) is dict
    if is_dict:
        mnanes = list(models.keys())
    else:
        mnanes = models

    # print header
    s1 = '{:>40s} '.format('Name')
    s2 = '{:40s} '.format('-' * 40)
    for k in pars:
        s1 += '{}  {:7s}'.format(sep, k)
        s2 += '{}  {:7s}'.format(sep, '-' * 7)
    s1 += '{}  {}'.format(sep, 'comment')
    s2 += '{}  {}'.format(sep, '-' * 7)
    print(s1)
    print(s2)
    # print("| %40s |  %7s |  %6s | %6s |  %s" % ('Name', 'R', 'M', 'E', 'comment'))
    # print("| %40s |  %7s |  %6s | %6s |  %s" % ('-' * 40, '-' * 7, '-' * 6, '-' * 6, '-' * 7))
    for mdl in mnanes:
        stella = ps.Stella(mdl, path=path)
        exts = models[mdl] if is_dict else ''
        if stella.is_tt:
            info = stella.get_tt().Info
            try:
                s = '{:>40s} '.format(info.Name)
                for k in pars:
                    v = getattr(info, k)
                    s += '{}  {:7.3f}'.format(sep, v)
                s += '{}  {}'.format(sep, exts)
                print(s)
            except KeyError as ex:
                print("| %40s |  %7s |  %6s | %6s | %s  | KeyError: %s" % (info.Name, '', '', '', exts, ex))
            except:
                print("| %40s |  %7s |  %6s | %6s | %s  | Unexpected error: %s"
                      % (info.Name, '', '', '', exts, sys.exc_info()[0]))
        else:
            if is_verbose:
                print("| %40s |  %26s |  %s" % (stella.name, ' ', exts))


def summary(models, path='.', pars=('R', 'M', 'Mni', 'E'), sep='  ', is_verbose=False):
    """
    Print list main param
    :param models: list of model names
    :param path: working directory. Default: ./
    :param pars: the parameters for summary. Default: pars=('R', 'M', 'Mni', 'E')
    :param sep: string separation. Default: space
    :param is_verbose:
    :return:
    """
    data = {k: [] for k in pars}
    for mdl in models:
        stella = ps.Stella(mdl, path=path)
        if stella.is_tt:
            info = stella.get_tt().Info
            try:
                for k in data.keys():
                    v = getattr(info, k)
                    if v not in data[k]:
                        data[k].append(v)
            except KeyError as ex:
                print(" KeyError: %s. %s" % (info.Name, ex))
        else:
            if is_verbose:
                print("{0} No tt.".format(stella))
    if len(list(data.values())[0]) > 0:
        print('Summary in {}'.format(path))
        for k in data.keys():
            print("{:5s} {}".format(k, sep.join(map(str, sorted(data[k])))))
    else:
        if is_verbose:
            print('No tt-files in {}'.format(path))


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
                        default=' | ',
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

    models = ps.ls(path, pattern=args.pattern)

    if len(models) == 0:
        print("No tt-files in {0}".format(path))
        return None

    pars = args.pars.split(':')
    sep = args.sep

    if args.is_summary:
        summary(models.keys(), path, pars, sep=sep, is_verbose=args.is_verbose)
    elif args.is_long:
        ls_long(models, path, pars, sep=sep, is_verbose=args.is_verbose)
    else:
        print('Models: {0} in {1}'.format(len(models), path))
        tbl = list_to_table(models.keys())
        print_table(tbl)
        # print("{0}".format(' '.join(models.keys())))


if __name__ == '__main__':
    main()
