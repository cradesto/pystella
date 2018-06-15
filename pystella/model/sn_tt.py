import logging
import numpy as np
import os
import re

from pystella.rf import band
from pystella.rf.lc import SetLightCurve, LightCurve

__author__ = 'bakl'


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class StellaTt:
    def __init__(self, name, path='./'):
        """Creates a StellaTt model instance.  Required parameters:  name."""
        self.name = name
        self.path = path  # path to model files

    def __str__(self):
        return "%s, path: %s" % (self.name, self.path)

    def __repr__(self):
        return "%s, path: %s" % (self.name, self.path)
        # return "%s" % self.name

    @property
    def Info(self):
        return StellaTtInfo(self.name, self.path)
        # return "%s" % self.name

    def read(self, ext='tt', line_header=80):
        """
        Read tt-data
        Columns: time Tbb rbb Teff Rlast_sc Rph Mbol MU MB MV MI MR   Mbolavg  gdepos
        :type ext: extension of tt-file, default: 'tt'
        :param line_header: skip_rows - 1
        :return: data in np.array
        """

        fname = os.path.join(self.path, self.name + '.' + ext)
        header = ''
        i = 0
        with open(fname, "r") as f:
            for line in f:
                i += 1
                if i < line_header:
                    continue
                if line.lstrip().startswith("time"):
                    header = line
                    line_header = i
                    break
        # time Tbb rbb Teff Rlast_sc R(tau2/3) Mbol MU MB MV MI MR Mbolavg  gdepos
        # time Tbb rbb Teff Rlast_sc R(tau2/3) Mbol MU MB MV MI MR   Mbolavg  gdepos
        if header != '':
            names = map(str.strip, header.split())
            names = [w.replace('R(tau2/3)', 'Rph') for w in names]
            dtype = np.dtype({'names': names, 'formats': [np.float64] * len(names)})
            block = np.loadtxt(fname, skiprows=line_header + 1, dtype=dtype)
            return block
        else:
            return None

    def read_curves(self):
        block = self.read()
        header = 'Mbol MU MB MV MI MR'.split()
        curves = SetLightCurve(self.name)
        time = block['time']
        for col in header:
            b = band.band_by_name(col.replace('M', ''))
            mag = block[col]
            lc = LightCurve(b, time, mag)
            curves.add(lc)
        return curves

    def read_curves_gri(self):
        block = self.read(ext='gri', line_header=1)
        # header = 'L_bol    Mu        MB    Mg         Mr         Mi'.split()
        # header = 'MB    MV'.split()
        header = 'L_bol    Mu        MB   MV    Mg    Mr Mi  J  H  K '.split()
        # header = 'MB        MV '.split()
        # header = 'L_bol   L_ubvgri      Mu        MB        MV       Mg         Mr         Mi'.split()
        curves = SetLightCurve(self.name)
        time = block['time']
        for col in header:
            b = band.band_by_name(col.replace('M', '').replace('L_', ''))
            mag = block[col]
            lc = LightCurve(b, time, mag)
            curves.add(lc)
        return curves


class StellaTtInfo:
    def __init__(self, name, path='./'):
        """Creates a StellaTtInfo model instance.  Required parameters:  name."""
        self._dict = None
        self._name = name
        self._path = path  # path to model files
        self.parse()

    def parse(self, header_end=29):
        fname = os.path.join(self._path, self._name + ".tt")
        self._dict = {'fname': fname}

        # prepare pattern
        pattern = re.compile(r"(.*?)\s*=\s*([+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)")
        # run pattern
        i = 0
        with open(fname) as f:
            for line in f:
                i += 1
                if i == header_end:
                    break
                if '= ' not in line:
                    continue
                res = pattern.findall(line)
                if len(res) > 0:
                    for k, v in res:
                        logger.debug("key: %s  v: %f " % (k, float(v)))
                        self._dict[str.strip(k)] = float(v)

        return self

    @property
    def Name(self):
        return self._name

    @property
    def Path(self):
        return self._path

    @property
    def Data(self):
        return self._dict

    @property
    def R(self):
        return self._dict['RADIUS(SOLAR)']

    @property
    def M(self):
        return self._dict['MASS(SOLAR)']

    @property
    def E(self):
        return self._dict['Ebstht'] / 10.  # to FOE

    @property
    def Mni(self):
        return self._dict['XMNI']

    @property
    def TcurA(self):
        return self._dict['TcurA']

    @property
    def TcurB(self):
        return self._dict['TcurB']

    @property
    def Rce(self):
        return self._dict['Rce']

    def show(self, comment=''):
        # print "INFO %s" % self.name
        # print " %40s: R = %7.2f M = %6.2f E = %6.2f " % (self.name, self.R, self.M, self.E)
        print("| %40s |  %7.2f |  %6.2f | %6.2f |  %s" % (self._name, self.R, self.M, self.E, comment
                                                          ))

    def print_tex(self, o=None, lend=''):
        # print "INFO %s" % self.name
        # print " %40s: R = %7.2f M = %6.2f E = %6.2f " % (self.name, self.R, self.M, self.E)
        s = " \\mbox{%s} &  %7.2f &  %6.2f & %6.2f \\\ %s " % (self._name, self.R, self.M, self.E, lend)
        if o is not None and o:
            return s
        else:
            print(o)
