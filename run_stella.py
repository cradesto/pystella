#!/usr/bin/env python3


import os
import argparse
import logging
import numpy as np

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class Runner(object):
    def __init__(self, config):
        self._config = config

    @property
    def Model(self):
        return self._config.Model

    @property
    def Cfg(self):
        return self._config

    def add_line(self, fname, nline=2, pattern=None):
        """
        Add a control line to the file named fname.

        :param fname:  file with type *.1
        :param nline:  position new line. Default: 2
        :param pattern: string to set the place to rename, as "111001", if it's None [default] rename all substrings.
        :return:
        """
        if nline < 1:
            raise ValueError('nline should be more 0: {}'.format(nline))

        nline -= 1  # shift to zero as array indexes

        # read the file
        with open(fname, "r") as f:
            lines = f.readlines()

        # take pattern as nline
        subs = lines[nline].split()
        if pattern is None:
            pattern = np.ones(len(subs))
        else:
            pattern = pattern.replace(' ', '')
            pattern = [int(p) for p in pattern]

        for i, p in enumerate(pattern):
            if bool(p):
                filename = subs[i]
                if '.' in filename:
                    (prefix, sep1, suffix) = filename.rpartition('.')
                else:
                    (prefix, sep1, suffix) = filename, '', ''
                (predir, sep2, nm) = prefix.rpartition('/')
                subs[i] = '{}{}{}{}{}'.format(predir, sep2, self.Model, sep1, suffix)
        newline = '  '.join(subs)

        # add new line at nline position
        lines.insert(nline, newline + "\n")
        # write to the file
        with open(fname, "w") as fout:
            for l in lines:
                fout.write(l)
            # fout.writelines(lines)

    def run(self):

        # create first line in stella files
        for sec in self.Cfg.Sections:
            opts = self.Cfg.Options(sec)
            pattern = None
            if 'pattern' in opts:
                pattern = self.Cfg.get(sec, 'pattern')

            if 'file_add_line' in opts:
                fname = os.path.join(self.Cfg.Root, self.Cfg.get(sec, 'file_add_line'))

                self.add_line(fname, pattern=pattern)
                logger.info(' {}: Added new line to {} with pattern= {}'.format(sec, fname, pattern))

        # # runner.add_line(cfg.Eve, pattern=None)
        # #  EVE Section
        # # pattern = '10101001'
        # logger.info(' Added new line to {} with pattern= {}'.format(self.Cfg.Eve1, self.Cfg.EvePattern))
        #
        # # VLADSF Section
        # # pattern = '1101'
        # self.add_line(self.Cfg.Vladsf1, pattern=self.Cfg.VladsfPattern)
        # logger.info(' Added new line to {} with pattern= {}'.format(self.Cfg.Vladsf1, self.Cfg.VladsfPattern))
        #
        # # STRAD Section
        # # pattern = '11111'
        # self.add_line(self.Cfg.Strad1, pattern=self.Cfg.StradPattern)
        # logger.info(' Added new line to {} with pattern= {}'.format(self.Cfg.Strad1, self.Cfg.StradPattern))


class Config(object):
    def __init__(self, fname, isload=True):
        self._fconfig = fname
        self._config = None
        if isload:
            self.load()

    def load(self):
        from configparser import ConfigParser
        # parser = ConfigParser()
        # parser.optionxform = str
        # parser.read(self.ConfigFile)
        # with open(self.ConfigFile) as f:
        #     sample_config = f.read()
        # self._config = ConfigParser(allow_no_value=True)
        # print('Load {}'.format(self.ConfigFile))
        # self._config.read_file(self.ConfigFile)
        parser = ConfigParser()
        parser.optionxform = str
        parser.read(self.ConfigFile)
        self._config = parser

    @property
    def ConfigFile(self):
        return self._fconfig

    @property
    def Config(self):
        return self._config

    @property
    def Sections(self):
        return self._config.sections()

    def Options(self, sec):
        # a = self.Config.items(sec)
        # if len(a)
        return self.Config.options(sec)

    def get(self, sec, name):
        return self.Config.get(sec, name)

    @property
    def Root(self):
        return self.Config.get('DEFAULT', 'root')

    @property
    def Model(self):
        return self.Config.get('DEFAULT', 'mname')

    @property
    def Eve(self):
        return os.path.join(self.Root, 'eve')

    @property
    def Eve1(self):
        fname = self.Config.get('EVE', 'file')
        return os.path.join(self.Eve, fname)

    @property
    def EvePattern(self):
        return self.Config.get('EVE', 'pattern')

    @property
    def Strad(self):
        return os.path.join(self.Root, 'strad')

    @property
    def Strad1(self):
        fname = self.Config.get('STRAD', 'file')
        return os.path.join(self.Strad, fname)

    @property
    def StradPattern(self):
        return self.Config.get('STRAD', 'pattern')

    @property
    def Vladsf(self):
        return os.path.join(self.Root, 'vladsf')

    @property
    def Vladsf1(self):
        fname = self.Config.get('VLADSF', 'file')
        return os.path.join(self.Vladsf, fname)

    @property
    def VladsfPattern(self):
        return self.Config.get('VLADSF', 'pattern')

    @property
    def as_dict(self):
        return {section: dict(self._config[section]) for section in self._config.sections()}

    def print(self):
        print(self.as_dict)


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
    import sys

    parser = get_parser()
    args, unknownargs = parser.parse_known_args()

    if args.is_show_only:
        print(">>>>>   JUST SHOW!   <<<<<<")

    fconfig = None
    if len(unknownargs) > 0:
        fconfig = unknownargs[0]
    else:
        if args.input:
            fconfig = args.input

    if fconfig is None:
        logger.error('  No any config data.')
        parser.print_help()
        sys.exit(2)

    fconfig = os.path.expanduser(fconfig)
    logger.info(' Read config: {}'.format(fconfig))
    cfg = Config(fconfig)
    # cfg.print()

    runner = Runner(cfg)
    logger.info(' Start runner for the model: {} in {} '.format(runner.Model, cfg.Root))

    runner.run()


if __name__ == '__main__':
    main()
