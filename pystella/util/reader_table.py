import numpy as np
import logging
import collections

from pystella.rf import band
from pystella.rf.lc import SetLightCurve, LightCurve

__author__ = 'bakl'

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

err_prefix = ('er', 'er_', 'err', 'err_')


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
    dt = np.dtype({'names': names, 'formats': [float] * len(names)})
    block = np.loadtxt(fname, skiprows=skip+1, dtype=dt, comments='#')
    return block


def read_obs_table_header(fname, header=None, skip=0, colt=('time', 'JD', 'MJD'),
                          include_names=None, include_patterns=None,
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
    :param include_patterns:  list or None
        Which columns to read as the pattern of regular expression.
        Default None, results in all columns being read.
        Example: ['Vel\d+','Vel.*']
        The columns with errors, like 'err'+use_names, also will be read.
    :param is_out:  bool, optional
        If True the skipped, header and first file-rows are printed.
        Default: False
    :param comments:  str or sequence, optional
        Default: '#'.
    :return:  ndarray - data is read from the file
    """
    lskip = 0
    if isinstance(colt, str):
        colt = [colt]

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
                    if not line.strip().startswith('###'):
                        lskip += 1
                        continue
                    else:
                        line = line.replace('###', '')
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
    cols_used = {}
    cols_data = {}

    def check_col_nm(nm, names, patterns):
        import re
        if names is None and patterns is None:
            return True

        if names is not None and v in names:
            return True

        if patterns is not None:
            for n in patterns:
                if re.match(n, nm):
                    return True
        return False

    for k, v in list(cols.items()):
        # time
        if not is_time and v.lower() in map(str.lower, colt):
            cols_used[k] = v
            is_time = True
        # data
        elif check_col_nm(v, include_names, include_patterns):
            cols_used[k] = v
            cols_data[k] = v
            #  error
            for err_name in (es+v for es in err_prefix):
                for i, bn in list(cols.items()):
                    if err_name.upper() == bn.upper():
                        cols_used[i] = bn
                        cols.pop(i)
                        break

    od = collections.OrderedDict(sorted(cols_used.items()))
    usecols = list(od.keys())
    names = list(od.values())

    dt = np.dtype({'names': names, 'formats': [float] * len(names)})
    block = np.loadtxt(fname, skiprows=max(lskip, skip)+1, dtype=dt, comments=comments, usecols=usecols)
    return block, cols_data


def table2curves(name, tbl, bands=None, colt=('time', 'JD', 'MJD'), is_filter_zero=True):
    # time = None
    for nm in colt:
        if nm in tbl.dtype.names:
            time = tbl[nm]
            break
    else:
        raise ValueError("THe table should contain a column with name in [{0}]".format(', '.join(colt)))

    curves = SetLightCurve(name)

    if bands is None:
        bands = [n for n in tbl.dtype.names if band.is_exist(n)]

    for bname in bands:
        b = band.band_by_name(bname)
        mag = tbl[bname]
        mask = ~np.isnan(mag)
        # filter
        if is_filter_zero:
            mask = np.logical_and(mask, mag != 0)  # filter out values not equal 0
        mask = np.logical_and(mask, mag < 99)

        t = time[mask]
        m = mag[mask]
        for err_name in (prefix+bname for prefix in err_prefix):
            if err_name in tbl.dtype.names:
                err = tbl[err_name]
                e = err[mask]
                lc = LightCurve(b, t, m, e)
                break
        else:
            lc = LightCurve(b, t, m)
        curves.add(lc)
    return curves


def curves2table(curves):
    def add(a, vals):
        # a = np.append(a, [lc.Mag], axis=0)
        # # a[:, :-1] = lc.Mag
        combined = np.vstack((a, vals))
        return combined
    # a = np.array(curves.TimeCommon, dtype=[('time', float, (len(curves.TimeCommon)))])
    # a = np.empty([0, len(curves.TimeCommon)])
    # a = np.array([0,100])
    # a = np.append(a, [curves.TimeCommon], axis=0)
    a = np.array(curves.TimeCommon)
    names = ['time']
    for lc in curves:
        a = add(a, lc.Mag)
        names.append(lc.Band.Name)
        if lc.IsErr:
            a = add(a, lc.MagErr)
            names.append('err'+lc.Band.Name)
    # dt = {'names': names, 'formats': [float] * len(names)}
    dt = list(zip(names, [float] * len(names)))
    a = a.T
    a.dtype = np.dtype(dt)
    # a.dtype = np.dtype(dt)
    # a.dtype.names = names

    return a
