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
    names = [s for s in header.split()]
    dt = np.dtype({'names': names, 'formats': [np.float64] * len(names)})
    block = np.loadtxt(fname, skiprows=skip+1, dtype=dt, comments='#')
    return block


def table2curves(name, tbl, bands=None, colt=('time', 'JD', 'MJD')):
    # time = None
    for nm in colt:
        if nm in tbl.dtype.names:
            time = tbl[nm]
            break
    else:
        raise ValueError("THe table should contain a column with name in {0}", ', '.join(colt))

    curves = SetLightCurve(name)

    if bands is None:
        bands = [n for n in tbl.dtype.names if band.band_is_exist(n)]

    for bname in bands:
        b = band.band_by_name(bname)
        mag = tbl[bname]
        # filter
        mask = np.where(mag != 0)  # filter out values not equal 0
        t = time[mask]
        m = mag[mask]
        for err_name in ('err'+bname, 'err_'+bname):
            if err_name in tbl.dtype.names:
                err = tbl[err_name]
                e = err[mask]
                lc = LightCurve(b, t, m, e)
                break
        else:
            lc = LightCurve(b, t, m)
        curves.add(lc)
    return curves
