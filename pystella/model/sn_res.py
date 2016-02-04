import os
import numpy as np
import itertools
import re
from reportlab.graphics.charts.axes import _AxisG

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
        dt = np.dtype({'names': cols, 'formats': map(str.strip, col_w.split())})
        # dt = np.dtype({'names': cols, 'formats': [np.float64] * len(cols)})
        # tmp = [x for x in map(str.strip, col_w.split())]
        # dt = np.dtype(','.join(tmp))
        # dt = np.dtype({'names': cols, 'formats': [np.float64] * len(cols)})
        with open(fname) as f_in:
            b = np.genfromtxt(itertools.islice(f_in, start-1, end),
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
            dm[:len(dm)-2] = dm[:len(dm)-2]-dm[1:len(dm)-1]
            idx = len(m)-len(dm) - 1
            dm = np.cumsum(dm)
            mtot[idx+1:] = m[idx] + dm[:]
            b['M'] = mtot

        #  make mass
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
