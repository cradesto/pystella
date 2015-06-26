import os
from os.path import dirname

__author__ = 'bakl'

import sys
import numpy as np
import string


class Band:
    def __init__(self, name=None, fname=None, load=1):
        """Creates a band instance.  Required parameters:  name and file.)."""
        self.name = name
        self.file = fname  # location of the filter response
        self.wl = None  # wavelength of response
        self.resp = None  # response
        if fname is not None and load == 1:
            self.read()

    def __str__(self):
        return "%s" % self.name

    def __repr__(self):
        return "%s" % self.name

    def read(self):
        """Reads in the response for file and updates several member functions."""
        if self.file is not None:
            f = open(self.file)
            try:
                lines = f.readlines()
                self.wl = np.array([float(string.split(line)[0]) for line in lines if line[0] != "#"])
                self.resp = np.array([float(string.split(line)[1]) for line in lines if line[0] != "#"])
                f.close()
            except Exception:
                print"Error in band file: %s.  Exception:  %s" % (self.file, sys.exc_info()[0])
            finally:
                f.close()

    def wave_range(self):
        if self.wl is not None:
            return np.min(self.wl), np.max(self.wl)
        else:
            return None

ROOT_DIRECTORY = dirname(dirname(dirname(os.path.abspath(__file__))))
  #  d = dirname(dirname(os.path.abspath(__file__)))


def get_full_path(fname):
    return os.path.join(ROOT_DIRECTORY, fname)


def band_get_names():
    bands = dict(U="kait_U.dat", B="kait_B.dat", V="kait_V.dat", R="kait_R.dat", I="kait_I.dat")
    d = os.path.join(ROOT_DIRECTORY, "data/bands/KAIT")
    for k, v in bands.items():
        bands[k] = os.path.join(d, v)
    return bands


def band_is_exist(name):
    b = band_get_names()
    return name in b


def band_by_name(name):
    """
        Get band by name, for example "U"
    :param name: for example "U" or "B"
    :return: class Band with waves and respons
    """
    if band_is_exist(name):
        bands = band_get_names()
        b = Band(name=name, fname=bands[name], load=1)
        return b
    else:
        None