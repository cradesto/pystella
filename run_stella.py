#!/usr/bin/env python3


import os
import argparse
import logging
import numpy as np
import shutil

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class Runner(object):
    def __init__(self, config):
        self._config = config

    @property
    def Cfg(self):
        return self._config

    @staticmethod
    def add_line(fname, mname, nline=2, pattern=None):
        """
        Add a control line to the file named fname.

        :param fname:  file with type *.1
        :param mname:  model name
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
                subs[i] = '{}{}{}{}{}'.format(predir, sep2, mname, sep1, suffix)
        newline = '  '.join(subs)

        # add new line at nline position
        lines.insert(nline, newline + "\n")
        # write to the file
        with open(fname, "w") as fout:
            for l in lines:
                fout.write(l)
            # fout.writelines(lines)

    @staticmethod
    def cp_sample(fin, fout, mode=1):
        if os.path.exists(fout):
            if mode == 2:
                raise ValueError('The target file: {}  is exist.'.format(fout))
            if mode == 1:
                logger.warning('The target file: {} was exist and left unchanged.'.format(fout))
                return

            logger.warning('The target file: {}  was rewritten.'.format(fout))
        logger.info(' Copy  {} to {}'.format(fin, fout))

        try:
            shutil.copy2(fin, fout)
        except shutil.SameFileError:
            logger.info('    {} and {} are the same file'.format(fin, fout))
            pass

    def run(self, mname):
        mode_sample = 1
        # create first line in stella files
        for sec in self.Cfg.Sections:
            opts = self.Cfg.Options(sec)
            pattern = None
            if 'pattern' in opts:
                pattern = self.Cfg.get(sec, 'pattern')

            if 'file_add_line' in opts:
                fname = os.path.join(self.Cfg.Root, self.Cfg.get(sec, 'file_add_line'))

                self.add_line(fname, mname, pattern=pattern)
                logger.info(' {}: Added new line to {} with pattern= {}'.format(sec, fname, pattern))

            if 'sample' in opts:
                fname = os.path.join(self.Cfg.Root, self.Cfg.get(sec, 'sample'))
                extension = os.path.splitext(fname)[1]
                fout = os.path.join(os.path.dirname(fname), '{}{}'.format(mname, extension))
                if 'mode_sample' in opts:
                    mode_sample = int(self.Cfg.get(sec, 'mode_sample'))
                self.cp_sample(fname, fout, mode=mode_sample)

        # todo Make run command


class Config(object):
    def __init__(self, fname, isload=True):
        self._fconfig = fname
        self._parser = None
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
        l = parser.read(self.ConfigFile)
        if len(l) > 0:
            self._parser = parser
            return True
        else:
            raise ValueError('Problem with reading the config file {}'.format(self.ConfigFile))

    @property
    def ConfigFile(self):
        return self._fconfig

    @property
    def Parser(self):
        return self._parser

    @property
    def Sections(self):
        return self._parser.sections()

    def Options(self, sec):
        # a = self.Config.items(sec)
        # if len(a)
        return self.Parser.options(sec)

    def get(self, sec, name):
        return self.Parser.get(sec, name)

    @property
    def Root(self):
        d = os.path.dirname(self._fconfig)
        if len(d) == 0:
            d = './'
        return d
        # return self.Config.get('DEFAULT', 'root')
        # return self.Config.get('DEFAULT', 'root')

    # @property
    # def Model(self):
    #     return self.Parser.get('DEFAULT', 'mname')

    @property
    def Eve(self):
        return os.path.join(self.Root, 'eve')

    @property
    def Eve1(self):
        fname = self.Parser.get('EVE', 'file')
        return os.path.join(self.Eve, fname)

    @property
    def EvePattern(self):
        return self.Parser.get('EVE', 'pattern')

    @property
    def Strad(self):
        return os.path.join(self.Root, 'strad')

    @property
    def Strad1(self):
        fname = self.Parser.get('STRAD', 'file')
        return os.path.join(self.Strad, fname)

    @property
    def StradPattern(self):
        return self.Parser.get('STRAD', 'pattern')

    @property
    def Vladsf(self):
        return os.path.join(self.Root, 'vladsf')

    @property
    def Vladsf1(self):
        fname = self.Parser.get('VLADSF', 'file')
        return os.path.join(self.Vladsf, fname)

    @property
    def VladsfPattern(self):
        return self.Parser.get('VLADSF', 'pattern')

    @property
    def as_dict(self):
        return {section: dict(self._parser[section]) for section in self._parser.sections()}

    def print(self):
        print(self.as_dict)


def get_parser():
    parser = argparse.ArgumentParser(description='Run STELLA modelling.')
    # print(" Observational data could be loaded with plugin, ex: -c lcobs:filedata:tshift:mshift")

    parser.add_argument('-i', '--input',
                        nargs='+',
                        required=True,
                        dest="input",
                        help="Model name, example: cat_R450_M15_Ni007 OR set of names: model1  model2")
    parser.add_argument('-r', '--run_config',
                        required=True,
                        dest="run_config",
                        help="config file, example: run.config")
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

    if len(unknownargs) > 0:
        logger.error('  Undefined strings in the command line')
        parser.print_help()
        sys.exit(2)
    else:
        if args.run_config:
            fconfig = os.path.expanduser(args.run_config)
        else:
            logger.error('  No any config data.')
            parser.print_help()
            sys.exit(2)

    models = args.input
    logger.info(' Run {} for  config: {}'.format(models, fconfig))
    cfg = Config(fconfig)
    # cfg.print()

    runner = Runner(cfg)

    for mname in models:
        logger.info(' Start runner for the model: {} in {} '.format(mname, cfg.Root))
        runner.run(mname)


if __name__ == '__main__':
    main()
