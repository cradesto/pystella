import numpy as np
import logging

from pystella.rf import band
from pystella.rf.lc import SetLightCurve, LightCurve

__author__ = 'bakl'

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def read_table_header_float(fname, header=None, skip=0):
    if header is None:
        i = 0
        with open(fname, "r") as f:
            for line in f:
                i += 1
                if i <= skip:
                    continue
                header = line
                break
    names = map(str.strip, header.split())
    dt = np.dtype({'names': names, 'formats': [np.float64] * len(names)})
    block = np.loadtxt(fname, skiprows=skip+1, dtype=dt)
    return block


def table2curves(tbl, bands, name):
    time = tbl['time']
    curves = SetLightCurve(name)
    for bname in bands:
        b = band.band_by_name(bname)
        mag = tbl[bname]
        lc = LightCurve(b, time, mag)
        curves.add(lc)
        return curves
