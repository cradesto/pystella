import csv

import numpy as np

__author__ = 'bakl'


def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")


def cache_save(tbl, fname):
    names = tbl.dtype.names
    with open(fname, 'wb') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['{:^8s}'.format(x) for x in names])
        for i, (row) in enumerate(zip(*[tbl[k] for k in names])):
            writer.writerow(['{:8.3f}'.format(x) for x in row])


def cache_load(fname):
    with open(fname, 'r') as f:
        header = f.readline()
    # header = 'time Tcol zeta Tnu Teff W'
    names = map(str.strip, header.split())
    dtype = np.dtype({'names': names, 'formats': [np.float64] * len(names)})
    tbl = np.loadtxt(fname, skiprows=1, dtype=dtype)
    return tbl
