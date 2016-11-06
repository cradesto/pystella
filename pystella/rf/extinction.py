import numpy as np
import os

__author__ = 'bakl'

ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
bands_alias = dict(U="Landolt_U", B="Landolt_B", V="Landolt_V", R="Landolt_R", I="Landolt_I",
                   J="UKIRT_J", H="UKIRT_H", K="UKIRT_K",  # UVM2="skyblue", UVW1="orange", UVW2="blue",
                   F105W="WFC3_F105W", F435W="ACS_F435W", F606W="ACS_F606W", F125W="WFC3_F125W", F140W="WFC3_F140W",
                   F160W="WFC3_F160W", F814W="ACS_F814W",
                   u="SDSS_u", g="SDSS_g", r="SDSS_r", i="SDSS_i", z="SDSS_z",
                   y='PS1_y', w='PS1_w')

law_default = 'Rv3.1'

EXT_DATA = None
EXT_LAWS = None


def read_data():
    laws = ('Rv2.1', 'Rv3.1', 'Rv4.1', 'Rv5.1')
    d = os.path.join(ROOT_DIRECTORY, "data/extinction")
    fname = os.path.join(d, 'extinction_schlafly.csv')
    data = np.loadtxt(fname, comments='#', skiprows=1,
                          dtype={'names': ('band', 'lambdaeff', 'Rv2.1', 'Rv3.1', 'Rv4.1', 'Rv5.1'),
                                 'formats': ('|S12', np.float, np.float, np.float, np.float, np.float)}, )
    return data, laws


def read_cache_data():
    global EXT_DATA, EXT_LAWS
    if EXT_DATA is None or EXT_LAWS is None:
        EXT_DATA, EXT_LAWS = read_data()
    return EXT_DATA, EXT_LAWS


def reddening_law(ebv, bands, law=law_default, is_info=True):
    """Compute  total extinction A(band)
    :param ebv: color excess E_{B-V}
    :param bands: list of bands
    :param law: reddening law 'Rv2.1', 'Rv3.1' (default), 'Rv4.1', 'Rv5.1'
    :param is_info:
    :return: dictionary of extinction for the bands
    """
    # see http://iopscience.iop.org/article/10.1088/0004-637X/737/2/103/meta
    if is_info:
        print "Reddening law [%s] Ebv=%6.3f " % (law, ebv)

    ext_all, laws = read_cache_data()
    if law not in laws:
        raise ValueError("The law should be in %s" % ' '.join(laws))

    ext_dic = dict(zip(ext_all['band'], ext_all[law]))

    if isinstance(bands, (str, unicode)):  # 1 bname
        return ext_dic[bands_alias[bands]] * ebv

    ext = {}
    for b in bands:
        ext[b] = ext_dic[bands_alias[b]] * ebv

    return ext


def reddening_law_z(ebv, bands, z, law=law_default, is_info=True):
    # see http://iopscience.iop.org/article/10.1088/0004-637X/737/2/103/meta

    # ext_all = dict((k, 0) for k, v in colors.items())   # it's just mock
    ext_all, laws = read_cache_data()
    lmb_dic = dict(zip(ext_all['band'], ext_all['lambdaeff']))
    # ext_dic = dict(zip(ext_all['band'], ext_all['Rv3.1']))

    def calc_Rv(bname):
        alias = bands_alias[bname]
        lmb_rest = lmb_dic[alias] / (1. + z)
        # find nearest band
        idx = np.abs(ext_all['lambdaeff'] - lmb_rest).argmin()
        rest_band = ext_all['band'][idx]
        rest_lmbeff = ext_all['lambdaeff'][idx]
        Rv = ext_all[law][idx]
        if is_info:
            print "Obs: %12s rest: %12s | obs.lambda=%8.1f rest.lambda=%8.1f nearest.rest.band.lambda=%8.1f Rv=%6.1f" % \
                  (alias, rest_band, lmb_dic[alias], lmb_rest, rest_lmbeff, Rv)
        return Rv

    ext = {}
    if is_info:
        print "Reddening law [%s] Ebv=%6.3f for z=%6.3f" % (law, ebv, z)

    if isinstance(bands, (str, unicode)):  # 1 bname
        return calc_Rv(bands) * ebv

    for b in bands:
        Rv = calc_Rv(b)
        ext[b] = Rv * ebv
    return ext


    return Rv
