#!/usr/bin/env python3
import os
import sys
import argparse
import logging

from io import IOBase
from sys import stdout
from select import select
from threading import Thread
from time import sleep

from io import StringIO

import shutil
from datetime import datetime

import numpy as np

logging.basicConfig(filename=datetime.now().strftime('run_%Y%m%d_%H_%M.log'), level=logging.DEBUG)

# handler = logging.root.handlers.pop()
# assert logging.root.handlers == [], "root logging handlers aren't empty"
# handler.stream.close()
# handler.stream = stdout
# logging.root.addHandler(handler)
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# Log to screen
console_logger = logging.StreamHandler(sys.stdout)
logger.addHandler(console_logger)


# -------------------------------
class StreamLogger(IOBase, logging.Handler):
    _run = None

    def __init__(self, logger_obj, level):
        super(StreamLogger, self).__init__()
        self.logger_obj = logger_obj
        self.level = level
        self.pipe = os.pipe()
        self.thread = Thread(target=self._flusher)
        self.thread.start()

    def __call__(self):
        return self

    def _flusher(self):
        self._run = True
        buf = b''
        while self._run:
            for fh in select([self.pipe[0]], [], [], 1)[0]:
                buf += os.read(fh, 1024)
                while b'\n' in buf:
                    data, buf = buf.split(b'\n', 1)
                    self.write(data.decode())
        self._run = None

    def write(self, data):
        return self.logger_obj.log(self.level, data)

    emit = write

    def fileno(self):
        return self.pipe[1]

    def close(self):
        if self._run:
            self._run = False
            while self._run is not None:
                sleep(1)
            for pipe in self.pipe:
                os.close(pipe)
            self.thread.join(1)


class LevelRangeFilter(logging.Filter):
    def __init__(self, min_level, max_level, name=''):
        super(LevelRangeFilter, self).__init__(name)
        self._min_level = min_level
        self._max_level = max_level

    def filter(self, record):
        return super(LevelRangeFilter, self).filter(record) and (
            self._min_level is None or self._min_level <= record.levelno) and (
                   self._max_level is None or record.levelno < self._max_level)


class Runner(object):
    def __init__(self, config):
        self._config = config

    @property
    def Cfg(self):
        return self._config

    @staticmethod
    def add_line(fname, mname, nline=2, pattern=None, is_demo=False):
        """
        Add a control line to the file named fname.

        :param is_demo:
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
        if not is_demo:
            with open(fname, "w") as fout:
                for l in lines:
                    fout.write(l)
        else:
            logger.info('Demo: to {} add line \n {} '.format(fname, newline))

    @staticmethod
    def cp_sample(fin, fout, mode=1, is_demo=False):
        if os.path.exists(fout):
            if mode == 2:
                raise ValueError('The target file: {}  is exist.'.format(fout))
            if mode == 1:
                logger.warning('The target file: {} was exist and left unchanged.'.format(fout))
                return

            logger.warning('The target file: {}  was rewritten.'.format(fout))
        logger.info(' Copy  {} to {}'.format(fin, fout))

        try:
            if not is_demo:
                shutil.copy2(fin, fout)
            else:
                logger.info('Demo:  copy {}  to {}'.format(fin, fout))
        except shutil.SameFileError:
            logger.info('    {} and {} are the same file'.format(fin, fout))
            pass

    @staticmethod
    def eval_cmd_log2(cmd, path, is_demo=False, **kwargs):
        import subprocess
        # # subprocess.check_output(['ls','-l']) #all that is technically needed...
        # print        subprocess.check_output(['ls', '-l'])
        logger.info(' Run cmd: {}   in {}'.format(cmd, path))
        stdin = kwargs.get('stdin', None)
        stdin_flag = None
        if not is_demo:
            if not stdin is None:
                stdin_flag = subprocess.PIPE

            try:
                stderr_stream = logging.StreamHandler(StringIO())
                stderr_stream.addFilter(LevelRangeFilter(logging.ERROR, None))
                logger.addHandler(stderr_stream)
                logger.setLevel(logging.ERROR)

                stdout_stream = logging.StreamHandler(StringIO())
                logger.addFilter(LevelRangeFilter(logging.INFO, logging.ERROR))
                logger.addHandler(stdout_stream)
                logger.setLevel(logging.INFO)

                with StreamLogger(logger, logging.INFO) as out, StreamLogger(logger, logging.ERROR) as err:
                    proc = subprocess.Popen(cmd, cwd=path, shell=True, stdout=out, stderr=err)

                # print
                # 'stderr_tee =', stderr_tee.stream.getvalue()
                # print
                # 'stdout_tee =', stdout_tee.stream.getvalue()
            finally:
                for handler in logger.handlers:
                    logger.removeHandler(handler)
                    handler.stream.close()
                    handler.close()
                # stderr_tee.stream.close()
                # stdout_tee.stream.close()

            return proc.returncode, False, False  #, stdout, stderr
        else:
            returncode, stdout, stderr = 0, False, False
            logger.info('Demo: system cmd: {} : '.format(cmd))
            return returncode, stdout, stderr

    @staticmethod
    def eval_cmd_log(cmd, path, is_demo=False, **kwargs):
        import subprocess
        # # subprocess.check_output(['ls','-l']) #all that is technically needed...
        # print        subprocess.check_output(['ls', '-l'])
        # logger.info(' Run cmd: {}   in {}'.format(cmd, path))
        stdin = kwargs.get('stdin', None)
        stdin_flag = None
        if not is_demo:
            logger.info(' Run cmd: {}   in {}'.format(cmd, path))

            # External script output
            logger.info(
                subprocess.check_output(cmd, cwd=path, shell=True)
            )
            return 0, False, False
        else:
            returncode, stdout, stderr = 0, '', ''
            logger.info('Demo: system cmd: {} : '.format(cmd))
            return returncode, stdout, stderr

    @staticmethod
    def eval_cmd(cmd, path, is_demo=False, **kwargs):
        import subprocess
        # # subprocess.check_output(['ls','-l']) #all that is technically needed...
        # print        subprocess.check_output(['ls', '-l'])
        logger.info(' Run cmd: {}   in {}'.format(cmd, path))
        stdin = kwargs.get('stdin', None)
        stdin_flag = None
        if not is_demo:
            if not stdin is None:
                stdin_flag = subprocess.PIPE

            proc = subprocess.Popen(
                cmd, cwd=path, shell=True,
                stdin=stdin_flag,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate(stdin)
            return proc.returncode, stdout, stderr
        else:
            returncode, stdout, stderr = 0, '', ''
            logger.info('Demo: system cmd: {} : '.format(cmd))
            return returncode, stdout, stderr

    def run(self, mname, is_sys=False, is_demo=False):
        mode_sample = 1
        # create first line in stella files
        for sec in self.Cfg.Sections:
            opts = self.Cfg.Options(sec)
            pattern = None
            
            if 'mode_sample' in opts:
                mode_sample = int(self.Cfg.get(sec, 'mode_sample'))

            if 'dir' in opts:
                path = os.path.join(self.Cfg.Root, self.Cfg.get(sec, 'dir'))
            else:
                path = self.Cfg.Root
            path = os.path.realpath(path)

            # print("{}: mode_sample= {}".format(sec, mode_sample));
                
            if 'pattern' in opts:
                pattern = self.Cfg.get(sec, 'pattern')

            if 'file_add_line' in opts:
                fname = os.path.join(self.Cfg.Root, self.Cfg.get(sec, 'file_add_line'))

                self.add_line(fname, mname, pattern=pattern, is_demo=is_demo)
                logger.info(' {}: Added new line to {} with pattern= {}'.format(sec, fname, pattern))

            if 'sample' in opts:
                fname = os.path.join(self.Cfg.Root, self.Cfg.get(sec, 'sample'))
                extension = os.path.splitext(fname)[1]
                fout = os.path.join(os.path.dirname(fname), '{}{}'.format(mname, extension))
                self.cp_sample(fname, fout, mode=mode_sample, is_demo=is_demo)

            if mode_sample == -1:
                cmd = './delstel.pl  {};'.format(mname)
                return_code, stdout, stderr = self.eval_cmd(cmd, path, is_demo=is_demo)

                
            if is_sys and 'cmd' in opts:
                cmd = self.Cfg.get(sec, 'cmd')
                return_code, stdout, stderr = self.eval_cmd(cmd, path, is_demo=is_demo)
                # return_code, stdout, stderr = self.eval_cmd_log(cmd, path, is_demo=is_demo)

                if stdout:
                    for line in stdout.decode('utf8').strip().split("\n"):
                        logger.info('stdout: {}'.format(line))

                if stderr:
                    for line in stderr.decode('utf8').strip().split("\n"):
                        logger.error('stderr: {}'.format(line))
                        # logger.error(line)


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
    parser.add_argument('--run',
                        action='store_const',
                        default=False,
                        const=True,
                        dest="is_sys",
                        help="Run with system commands.")

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
        runner.run(mname, is_sys=args.is_sys, is_demo=args.is_show_only)


if __name__ == '__main__':
    main()
