import os
import numpy as np
import itertools
import re

__author__ = 'bakl'


class StellaRes:
    def __init__(self, name, path='./'):
        """Creates a StellaRes model instance.  Required parameters:  name."""
        self.name = name
        self.path = path  # path to model files
        self.nzon = None
        self.block = None

    def __str__(self):
        return "%s, path: %s" % (self.name, self.path)

    def __repr__(self):
        return "%s, path: %s" % (self.name, self.path)
        # return "%s" % self.name

    @property
    def Info(self):
        return StellaResInfo(self.name, self.path)
        # return "%s" % self.name

    def read_at_time(self, time):
        t, s, e = self.find_block(time)
        return self.read_res_block(s, e)

    def read_res_block(self, start, end):
        """
        Read one part of res-data
        :param start: line number to start
        :param end: line number to finish
        :return:
        """
        fname = os.path.join(self.path, self.name + ".res")
        nrows = end - start + 1
        colstr = "ZON M R14 V8 T5 Trad5 lgDm6 lgP7  lgQv lgQRT XHI ENG LUM CAPPA ZON1 n_bar n_e Fe II III"
        cols = map(str.strip, colstr.split())
        # FORMAT(1X,I3,1P,E9.2,0P,F12.7,F8.4,F10.3,F8.3,4F7.2,1P,4E10.2,I5, 11e10.2);

        # delimiter = (4, 9, 12, 8, 10, 8, 7, 7, 7, 7, 10, 10, 10, 10, 5, 10, 10, 10, 10, 10)
        col_w = "i4 f9 f12 f8 f10 f8 f7 f7 f7 f7 e10 e10 e10 e10 i5 e10 e10 e10 e10 e10"
        delimiter = (int(re.sub(r'\D', '', w)) for w in map(str.strip, col_w.split()))
        # dt = np.dtype({'names': cols, 'formats': map(str.strip, col_w.split())})
        dt = np.dtype({'names': cols, 'formats': [np.float] * len(cols)})
        # tmp = [x for x in map(str.strip, col_w.split())]
        # dt = np.dtype(','.join(tmp))
        # dt = np.dtype({'names': cols, 'formats': [np.float64] * len(cols)})
        with open(fname) as f_in:
            b = np.genfromtxt(itertools.islice(f_in, start - 1, end),
                              names=cols, dtype=dt, delimiter=delimiter)

        # cols = ["z", "m", "r", "rho", "v", "T", "Trad"]
        # dt = np.dtype({'names': cols, 'formats': [np.float64] * len(cols)})
        # c = np.vstack((block['ZON'],  # z
        #               block['M'],  # m
        #               block['R14'] * 1e14,  # r
        #               10 ** (block['lgDm6'] - 6.),  # rho
        #               block['V8'] * 1e8,
        #               block['T5'] * 1e5,
        #               block['Trad5'] * 1e5))
        # b = np.array(c, dtype=dt)
        m = b['M']
        dm = -1. * m[m < 0.]
        if len(dm) > 0:
            mtot = m
            mtot[m < 0.] = 0.
            dm[:len(dm) - 2] = dm[:len(dm) - 2] - dm[1:len(dm) - 1]
            idx = len(m) - len(dm) - 1
            dm = np.cumsum(dm)
            mtot[idx + 1:] = m[idx] + dm[:]
            b['M'] = mtot

        # make mass
        return b

    def find_block(self, time):
        t = 0.
        start = end = 0
        is_block = False
        fname = os.path.join(self.path, self.name + '.res')
        i = 0
        with open(fname, "r") as f:
            for line in f:
                i += 1
                if is_block and line.lstrip().startswith("%B"):
                    end = i - 1  # %B - next line
                    break
                if not is_block and line.lstrip().startswith("OBS.TIME"):
                    header = line
                    names = map(str.strip, header.split())
                    t = float(names[1])
                    if t > time:
                        is_block = True
                        start = i + 2  # data after OBS.TIME line

        if start < end:
            return t, start, end
        else:
            return None


class StellaResInfo:
    def __init__(self, name, path='./'):
        """Creates a StellaRes model instance.  Required parameters:  name."""
        self._dict = None
        self.name = name
        self.path = path  # path to model files
        self.nzon = None
        self.parse()

    def parse(self, header_end=28):
        fname = os.path.join(self.path, self.name + ".res")
        self._dict = {'fname': fname}

        # prepare pattern
        res_header = """
             MASS(SOLAR)= 6.000       RADIUS(SOLAR)=    3.200
             EXPLOSION ENERGY(10**50 ERG)= 2.50000E-01

                  <=====                                 =====>


  INPUT PARAMETERS
 EPS   =         0.00300          Rce   =     1.00000E-02 SOL.Rad.
 HMIN  =     1.00000E-25 S        AMht  =     2.60000E+00 SOL.MASS
 HMAX  =     5.00000E+04 S        Tburst=     1.00000E+04 S
 THRMAT=     1.00000E-30          Ebstht=     2.50000E-01 1e50 ergs
 METH  =               3          CONV  =               F
 JSTART=               0          EDTM  =               T
 MAXORD=               4          CHNCND=               T
 NSTA  =           -4500          FitTau=     3.00000E+02
 NSTB  =              -1          TauTol=     1.30000E+00
 NOUT  =             100          IOUT  =              -1
 TcurA =         0.00000          Rvis   =        1.00000
 TcurB =        60.00000          BQ    =     1.00000E+00
 NSTMAX=        14400000          DRT   =     1.00000E+00
 XMNI  =     0.00000E+00 SOL.MASS NRT   =               1
 XNifor=     0.00000E+00
 MNicor=     1.16292E-01 SOL.MASS SCAT  =               T
        """

        # pattern = filter(lambda x: len(x) > 0, res_header.splitlines())
        # patternDigit = map(lambda s:
        #                    # re.findall(r"[-+]?\d*\.\d+|\d+", s)
        #                    re.sub(r"[-+]?\d*\.\d+|\d+", '()' s)
        #                    , pattern)
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
                        # print "key: %s  v: %f " % (k, float(v))
                        self._dict[str.strip(k)] = float(v)
                    #
                    # p = r"(.*?)\s*=\s+([-+]?\d*\.\d+|\d+)"
                    # for line in pattern:
                    #     res = re.findall(p, line)
                    #     if len(res) > 0:
                    #         for k, v in res:
                    #             print "key: %s  v: %f " % (k, float(v))


                                # exp = re.compile(".*RADIAT.*=(.*)TOTAL.*ENERGY.*=(.*)")
        # tmp = re.match(exp, buffer)
        # self.rad_e = float(tmp.groups()[0])
        # self.tot_e = float(tmp.groups()[1])

        return self

    @property
    def R(self):
        return self._dict['RADIUS(SOLAR)']

    @property
    def M(self):
        return self._dict['MASS(SOLAR)']

    @property
    def E(self):
        return self._dict['Ebstht']

    def show(self):
        # print "INFO %s" % self.name
        print "%40s: R = %7.2f M = %6.2f E = %6.2f" % (self.name, self.R, self.M, self.E)
