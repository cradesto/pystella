import os
import sys
import numpy as np
import string
from os.path import dirname
from pystella.util.phys_var import phys
from pystella.util.arr_dict import merge_dicts

__author__ = 'bakl'


class Band(object):
    Cache = dict()

    def __init__(self, name=None, fname=None, load=1):
        """Creates a band instance.  Required parameters:  name and file."""
        self.name = name
        self.file = fname  # location of the filter response
        self.__freq = None  # frequencies of response [Hz]
        self.__wl = None  # wavelength of response [cm]
        self.resp = None  # response
        self.zp = None  # zero points
        self.file_zp = 'filters.dat'  # file with zero points data
        if fname is not None and load == 1:
            self.read()

    @property
    def freq(self):
        return self.__freq

    @property
    def wl(self):
        return self.__wl

    @wl.setter
    def wl(self, wl):
        self.__wl = wl
        self.__freq = phys.c / self.__wl

    def __str__(self):
        return "%s" % self.name

    def __repr__(self):
        return "%s" % self.name

    def read(self):
        """Reads the waves and response from file. Read zero point"""
        if self.file is not None:
            f = open(self.file)
            try:
                lines = f.readlines()
                wl = np.array([float(string.split(line)[0]) for line in lines if line[0] != "#"])
                self.resp = np.array([float(string.split(line)[1]) for line in lines if line[0] != "#"])
                self.wl = wl * phys.angs_to_cm
                f.close()
            # except Exception:
            #     print"Error in rf file: %s.  Exception:  %s" % (self.file, sys.exc_info()[0])
            finally:
                f.close()

            # read zero point
            ptn_file = os.path.basename(self.file)
            fname = os.path.join(dirname(os.path.abspath(self.file)), self.file_zp)

            self.zp = read_zero_point(fname, ptn_file)

    @property
    def wave_range(self):
        if self.wl is not None:
            return np.min(self.wl), np.max(self.wl)
        else:
            return None

    def summary(self, out=sys.stdout):
        """Get a quick summary of the band

        Args:
         out (str or open file): where to write the summary

        Returns:
             None
        """
        print >> out, '-' * 80
        print >> out, "Band: ", self.name
        print >> out, "wave length = %.3f, %3.f " % (min(self.wl)*1e8, max(self.wl)*1e8)
        print >> out, ""


def read_zero_point(fname, ptn_file):
    f = open(fname)
    try:
        lines = f.readlines()
        for line in lines:
            if line[0] == "#":
                continue
            if string.split(line)[1] == ptn_file:
                zp = float(string.split(line)[2])
                return zp
    # except Exception:
    #     print"Error in zero points data-file: %s.  Exception:  %s" % (fname, sys.exc_info()[0])
    #     return np.nan
    finally:
        f.close()


ROOT_DIRECTORY = dirname(dirname(dirname(os.path.abspath(__file__))))


def get_full_path(fname):
    return os.path.join(ROOT_DIRECTORY, fname)


def bands_dict_Bessell():
    bands = dict(U="U-bessell.dat", B="B-bessell.dat", V="V-bessell.dat", R="R-bessell.dat", I="I-bessell.dat")
    d = os.path.join(ROOT_DIRECTORY, "data/bands/Bessell")
    for k, v in bands.items():
        bands[k] = os.path.join(d, v)

    return bands

def bands_dict_KAIT():
    bands = dict(U="kait_U.dat", B="kait_B.dat", V="kait_V.dat", R="kait_R.dat", I="kait_I.dat")
    d = os.path.join(ROOT_DIRECTORY, "data/bands/KAIT")
    for k, v in bands.items():
        bands[k] = os.path.join(d, v)

    return bands

def bands_dict_Persson():
    bands = dict(J="jfilter", H="hfilter", K="kfilter")
    d = os.path.join(ROOT_DIRECTORY, "data/bands/Persson")
    for k, v in bands.items():
        bands[k] = os.path.join(d, v)

    return bands

def bands_dict_USNO():
    # bands = dict(V="usno_g.res", I="usno_i.res", R="usno_r.res", U="usno_u.res", B="usno_z.res") # for comparison
    bands = dict(g="usno_g.res", i="usno_i.res", r="usno_r.res", u="usno_u.res", z="usno_z.res")
    d = os.path.join(ROOT_DIRECTORY, "data/bands/USNO40")
    for k, v in bands.items():
        bands[k] = os.path.join(d, v)
    return bands


def bands_dict_SDSS():
    bands = dict(g="sdss_g.dat", i="sdss_i.dat", r="sdss_r.dat", u="sdss_u.dat", z="sdss_z.dat")
    d = os.path.join(ROOT_DIRECTORY, "data/bands/SDSS")
    for k, v in bands.items():
        bands[k] = os.path.join(d, v)
    return bands


def bands_dict_PS1():
    """
        The Pan-STARRS1 Photometric System
        see http://ipp.ifa.hawaii.edu/ps1.filters/

    :return: band-pass filters
    """
    bands = dict(gps1="ps1_g.dat", ips1="ps1_i.dat", rps1="ps1_r.dat", zps1="ps1_z.dat",
                 y="ps1_y.dat", w="ps1_w.dat")
    d = os.path.join(ROOT_DIRECTORY, "data/bands/PS1")
    for k, v in bands.items():
        bands[k] = os.path.join(d, v)
    return bands


def bands_dict_SWIFT():
    bands4 = dict(UVM2="photonUVM2.dat", UVW1="photonUVW1.dat", UVW2="photonUVW2.dat", U_UVOT="photonU_UVOT.dat"
                  , B_UVOT="photonB_UVOT.dat", V_UVOT="photonV_UVOT.dat")
    d = os.path.join(ROOT_DIRECTORY, "data/bands/SWIFTUVOT")
    for k, v in bands4.items():
        bands4[k] = os.path.join(d, v)
    return bands4


def band_get_names():
    # STANDARD
    # bands0 = bands_dict_STANDARD()

    # KAIT
    # bands1 = bands_dict_KAIT()
    # Bessell
    bands1 = bands_dict_Bessell()

    # HJK
    bandsJHK = bands_dict_Persson()

    # USNO40
    bands2 = bands_dict_USNO()

    # SDSS
    bands3 = bands_dict_SDSS()

    # SWIFT
    bands4 = bands_dict_SWIFT()

    # The Pan-STARRS1 Photometric System
    bands5 = bands_dict_PS1()

    return merge_dicts(bands1, bandsJHK, bands3, bands4, bands5)


def band_is_exist(name):
    b = band_get_names()
    return name in b


def band_by_name(name):
    """
        Get rf by name, for example "U"
    :param name: for example "U" or "B"
    :return: class Band with waves and respons
    """
    if name in Band.Cache:
        return Band.Cache[name]

    if band_is_exist(name):
        bands = band_get_names()
        b = Band(name=name, fname=bands[name], load=1)
        Band.Cache[name] = b
        return b
    else:
        return None
