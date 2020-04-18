import csv

import numpy as np

__author__ = 'bakl'


def str2bool(v):
    return v.lower() in ("yes", "true", "t", 'y', "1")


def str2interval(arg, llim=0., rlim=float('inf'), sep=':'):
    if sep in arg:
        if str.endswith(arg, sep):
            xlim = [float(arg[:-1]), rlim]
        elif str.startswith(arg, sep):
            xlim = [llim, float(arg[1:])]
        else:
            xlim = list(map(float, arg.split(sep)))
    else:
        xlim = [llim, float(arg)]
    return xlim


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


def list_to_table(lst, col=4):
    table = []
    row = []
    for i, x in enumerate(lst):
        row.append(x)
        if len(row) % col == 0:
            table.append(row)
            row = []
    if len(row) > 0:
        table.append(row)
    return table


def parse_float_table(fname, start, nzon):
    from itertools import islice
    b = []
    with open(fname, "rb") as f:
        for i, line in enumerate(islice(f, start - 1, start - 1 + nzon)):
            items = [float(v) for v in line.split()]
            b.append(items)
    # print(start, start-1+nzon, nzon)
    return np.array(b)


def print_table(tbl):
    # col_width = [max(len(x) for x in col) for col in zip(*tbl)]
    col_width = []
    # find the size of longest value
    for line in tbl:
        for i, x in enumerate(line):
            if i >= len(col_width):
                col_width.append(len(x))
            else:
                col_width[i] = max(col_width[i], len(x))
    # print
    for line in tbl:
        print("  ".join("{:{}}".format(x, col_width[i])
                        for i, x in enumerate(line)))
        # print("| " + " | ".join("{:{}}".format(x, col_width[i])
        #                         for i, x in enumerate(line)) + " |")
