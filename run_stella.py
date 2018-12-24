#!/usr/bin/env python3

import io
import os
# from os.path import dirname
import argparse
import logging
import ConfigParser

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# import numpy as np

# ROOT_DIRECTORY = dirname(os.path.abspath(__file__))


class Config(object):
    def __init__(self, fname, isload=True):
        self._fconfig = fname
        self._config = None
        if isload:
            self.load()

    def load(self):
            with open(self.ConfigFile) as f:
                sample_config = f.read()
            self._config = ConfigParser.RawConfigParser(allow_no_value=True)
            self._config.readfp(io.BytesIO(sample_config))


    @property
    def ConfigFile(self):
        return self._fconfig

    @property
    def Root(self):
        return self._path_root

    @property
    def Eve(self):
        return os.path.join(self.Root, 'eve')

    @property
    def Eve1(self):
        return os.path.join(self.Eve, 'eve.1')

    @property
    def Strad(self):
        return os.path.join(self.Root, 'strad')

    @property
    def Strad1(self):
        return os.path.join(self.Strad, 'strad.1')

    @property
    def Vladsf(self):
        return os.path.join(self.Root, 'vladsf')

    @property
    def Ronfict1(self):
        return os.path.join(self.Vladsf, 'ronfict.1')


def get_parser():
    parser = argparse.ArgumentParser(description='Run STELLA modelling.')
    # print(" Observational data could be loaded with plugin, ex: -c lcobs:filedata:tshift:mshift")

    parser.add_argument('-i', '--input',
                        # nargs='+',
                        required=False,
                        dest="input",
                        help="config file, example: cat_R450_M15_Ni007.config")
    parser.add_argument('-so', '--show-only',
                        action='store_const',
                        default=False,
                        const=True,
                        dest="is_show_only",
                        help="Just show config, nothing will be done.")
    return parser


def main():
    parser = get_parser()
    args, unknownargs = parser.parse_known_args()

    if args.is_show_only:
        print(">>>>>   JUST SHOW!   <<<<<<")

    fconfig = None
    if len(unknownargs) > 0:
        fconfig = unknownargs[0]
    else:
        if args.input:
            fconfig = os.path.expanduser(args.input)

    if fconfig is None:
        logger.error('  No any config data.')
        parser.print_help()
    else:
        logger.info(' Read config: {}'.format(fconfig))


if __name__ == '__main__':
    main()
