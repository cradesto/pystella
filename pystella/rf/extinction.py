import numpy as np
import os

__author__ = 'bakl'

ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def ext_read_data():
    d = os.path.join(ROOT_DIRECTORY, "data/extinction")
    fname = os.path.join(d, 'extinction_schlafly.csv')
    ext_data = np.loadtxt(fname, comments='#', skiprows=1,
                          dtype={'names': ('band', 'lambdaeff', 'Rv2.1', 'Rv3.1', 'Rv4.1', 'Rv5.1'),
          'formats': ('|S12', np.float, np.float, np.float, np.float, np.float)},)
    res = dict(zip(ext_data['band'], ext_data['Rv3.1']))
    return res  # ext_data[:, [0, col31]]


def extinction_law(ebv, bands):
    # see http://iopscience.iop.org/article/10.1088/0004-637X/737/2/103/meta
    alias = dict(U="Landolt_U", B="Landolt_B", V="Landolt_V", R="Landolt_R", I="Landolt_I",
                 J="UKIRT_J", H="UKIRT_H", K="UKIRT_K",  # UVM2="skyblue", UVW1="orange", UVW2="blue",
                 F105W="WFC3_F105W", F435W="ACS_F435W",  F606W="ACS_F606W", F125W="WFC3_F125W", F140W="WFC3_F140W",
                 F160W="WFC3_F160W", F814W="ACS_F814W",
                 u="SDSS_u", g="SDSS_g", r="SDSS_r", i="SDSS_i", z="SDSS_z",
                 y='PS1_y', w='PS1_w')
    
    # ext_all = dict((k, 0) for k, v in colors.items())   # todo it's just mock
    ext_all = ext_read_data()
    ext = {}
    for b in bands:
        ext[b] = ext_all[alias[b]]*ebv
    return ext
