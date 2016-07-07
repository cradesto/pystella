import os
import sys
import numpy as np
import string
from os.path import dirname
from pystella.util.phys_var import phys
from pystella.util.arr_dict import merge_dicts

__author__ = 'bakl'

# see bands: http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse&gname=GALEX


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

    @property
    def Name(self):
        return self.name

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
            # except Exception:
            #     print"Error in band file: %s.  Exception:  %s" % (self.file, sys.exc_info()[0])
            #     sys.exit("Error while parse band-file: %s in %s" % self.file)
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
        print >> out, "wave length = %.3f, %3.f " % (min(self.wl) * 1e8, max(self.wl) * 1e8)
        print >> out, ""


class BandUni(Band):
    def __init__(self, name='Uniform', wlrange=(1e1, 5e4), length=100):
        """Creates a band with uniform responce.
        :param name:  default 'Uniform'.
        :param wlrange:  the wavelength range, default (1e1, 5e4) [A]
        :param length: numbers of bins default 100
        """
        super(BandUni, self).__init__(name)
        # self.name = name
        wl = np.exp(np.linspace(np.log(wlrange[0]) * phys.angs_to_cm
                                , np.log(wlrange[1]) * phys.angs_to_cm
                                , length))
        self.wl = wl  # wavelength of response [cm]
        self.resp = 1.  # response


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


def bands_colors():
    colors = dict(U="blue", B="cyan", V="black", R="red", I="magenta",
                  J="green", H="cyan", K="black",
                  UVM2="skyblue", UVW1="orange", UVW2="blue",
                  F105W="magenta", F435W="skyblue",  F606W="cyan", F125W="g", F140W="orange", F160W="r", F814W="blue",
                  g="green", r="red", i="magenta", u="blue", z="chocolate",
                  y='olive', w='tomato')
    # for Subaru HCS: colors
    for b in list('grizy'):
        colors['HSC' + b] = colors[b]

    return colors


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


def bands_dict_HST():
    d = os.path.join(ROOT_DIRECTORY, "data/bands/HST")
    # bands = dict(F125W="hst_wfc3_ir_f125w.dat", F160W="hst_wfc3_ir_f160w.dat")
    fname = os.path.join(d, 'filters.dat')
    bands = {}
    with open(fname, "r") as f:
        for line in f:
            names = map(str.strip, line.split())
            bands[names[0]] = os.path.join(d, names[1])
    return bands


def bands_dict_PS1():
    """
        The Pan-STARRS1 Photometric System
        see http://ipp.ifa.hawaii.edu/ps1.filters/

    :return: band-pass filters
    """
    bands = dict(PS1g="ps1_g.dat", PS1i="ps1_i.dat", PS1r="ps1_r.dat", PS1z="ps1_z.dat",
                 y="ps1_y.dat", w="ps1_w.dat")
    d = os.path.join(ROOT_DIRECTORY, "data/bands/PS1")
    for k, v in bands.items():
        bands[k] = os.path.join(d, v)
    return bands


def bands_dict_SWIFT():
    bands4 = dict(UVM2="photonUVM2.dat", UVW1="photonUVW1.dat", UVW2="photonUVW2.dat",
                  UVOTU="photonU_UVOT.dat", UVOTB="photonB_UVOT.dat", UVOTV="photonV_UVOT.dat")
    d = os.path.join(ROOT_DIRECTORY, "data/bands/SWIFTUVOT")
    for k, v in bands4.items():
        bands4[k] = os.path.join(d, v)
    return bands4


def bands_dict_SubaruHSC():
    bands4 = dict(HSCg="HSC-g.txt", HSCr="HSC-r.txt", HSCi="HSC-i.txt",
                  HSCz="HSC-z.txt", HSCy="HSC-Y.txt")
    d = os.path.join(ROOT_DIRECTORY, "data/bands/SubaruHSC")
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
    # bands2 = bands_dict_USNO()

    # SDSS
    bands3 = bands_dict_SDSS()

    # SWIFT
    bands4 = bands_dict_SWIFT()

    # The Pan-STARRS1 Photometric System
    bands5 = bands_dict_PS1()

    # The HSC filters
    bands6 = bands_dict_SubaruHSC()

    # The HSC filters
    bands7 = bands_dict_HST()

    return merge_dicts(bands1, bandsJHK, bands3, bands4, bands5, bands6, bands7)


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
