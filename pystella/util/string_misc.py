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


def list_to_table(l, col=4):
    table = []
    row = []
    for i, x in enumerate(l):
        row.append(x)
        if len(row) % col == 0:
            table.append(row)
            row = []
    if len(row) > 0:
        table.append(row)
    return table


def print_table(tbl):
    # col_width = [max(len(x) for x in col) for col in zip(*tbl)]
    col_width = []
    for line in tbl:
        for i, x in enumerate(line):
            if i >= len(col_width):
                col_width.append(len(x))
            else:
                col_width[i] = max(col_width[i], len(x))
    for line in tbl:
        print("  ".join("{:{}}".format(x, col_width[i])
                        for i, x in enumerate(line)))
        # print("| " + " | ".join("{:{}}".format(x, col_width[i])
        #                         for i, x in enumerate(line)) + " |")
