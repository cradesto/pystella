import numpy as np
import logging
import collections

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


def read_obs_table_header(fname, header=None, skip=0, colt=('time', 'JD', 'MJD'), include_names=None,
                          is_out=False, comments='#'):
    """
    Load tabular data from file.
    :param fname: The name of the file with data
    :param header: str, optional
        The string is used to build data-type of the resulting array.
        Default: None.
    :param skip:  int, optional
        Skip the first rows before header.
        Default: 0.
    :param colt: list, optional.
        Possible names for time column.
        Default:  ('time', 'JD', 'MJD')
    :param include_names:  list or None
        Which columns to read.
        Default None, results in all columns being read.
        Example: ['B','V','R']
        The columns with errors, like 'err'+use_names, also will be read.
    :param is_out:  bool, optional
        If True the skipped, header and first file-rows are printed.
        Default: False
    :param comments:  str or sequence, optional
        Default: '#'.
    :return:  ndarray - data is read from the file
    """
    lskip = 0
    if header is None:
        with open(fname, "r") as f:
            i = 0
            for line in f:
                i += 1
                if is_out:
                    print(line.strip())
                if i <= skip:
                    continue
                if line.strip().startswith(comments):
                    lskip += 1
                    continue
                header = line
                break
            else:
                raise ValueError('Could not get header. Check the file: {}. Probably skip [{}] is too large.'
                                 .format(fname, skip))
            # print first lines
            if is_out:
                line = f.readline().strip()
                print(line)
    cols_names = header.split()
    cols = {i: nm for i, nm in enumerate(cols_names)}

    is_time = False
    cols_data = {}
    for k, v in cols.items():
        # time
        if not is_time and v in colt:
            cols_data[k] = v
            is_time = True
        # data
        elif include_names is None or v in include_names:
            cols_data[k] = v
            # band error
            for err_name in ('err' + v, 'err_' + v):
                for i, bn in cols.items():
                    if err_name == bn:
                        cols_data[i] = err_name

    od = collections.OrderedDict(sorted(cols_data.items()))
    usecols = list(od.keys())
    names = list(od.values())

    dt = np.dtype({'names': names, 'formats': [np.float64] * len(names)})
    block = np.loadtxt(fname, skiprows=max(lskip, skip)+1, dtype=dt, comments=comments, usecols=usecols)
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
