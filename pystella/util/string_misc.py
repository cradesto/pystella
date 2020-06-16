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


def df2tex(a, pfx_err='err', pfx_errpm=('errp', 'errm'), fmt='.3f', fmt_err='.4f', math=' $ ', sep=' & ',
           str_end=' \\\ \n'):
    """
    Convert pandas dataframe to latex table.
    Usage:
    arr = pd2np(df_fit_images)
    s = df2tex(arr, pfx_err=None, fmt_err='.2f', sep='  &  ', strend='   \\\ \n')
    print(s)
    """
    names = a.dtype.names
    if not isinstance(fmt, dict):
        fmt = {nm: fmt for nm in names}
    if not isinstance(fmt_err, dict):
        fmt_err = {nm: fmt_err for nm in names}

    is_pm = pfx_err is not None
    is_up_down = not is_pm and pfx_errpm is not None

    var_names_errp, var_names_errm, var_names_errpm = {}, {}, {}
    if is_up_down:
        var_names = [x for x in names if not pfx_errpm[0] in x and not pfx_errpm[1] in x]
        var_names_errp = {xx: x for xx in var_names for x in names if xx in x and pfx_errpm[0] in x}
        var_names_errm = {xx: x for xx in var_names for x in names if xx in x and pfx_errpm[1] in x}
    elif is_pm:
        var_names = [x for x in names if not pfx_err in x]
        var_names_errpm = {xx: x for xx in var_names for x in names if xx in x and pfx_err in x}
    else:
        var_names = names
    #     print(var_names_errm)
    s_res = sep.join([nm.replace('_', '\_') for nm in var_names]) + str_end

    for i, row in enumerate(a):
        srow = ''
        for j, nm in enumerate(var_names):
            scol = ''
            x = row[nm]
            try:
                scol += '{:{fmt}}'.format(x, fmt=fmt[nm])
            except ValueError:
                scol += '{}'.format(x.replace('_', '\_'))
            # print uncertanties
            if is_up_down and nm in var_names_errp:
                #                 print(nm, var_names_errp[nm], var_names_errm[nm])
                ep, em = row[var_names_errp[nm]], row[var_names_errm[nm]]
                #                 try:
                scol += '^{{{ep:{fe}}}}_{{{em:{fe}}}}'.format(ep=ep, em=em, fe=fmt_err[nm])
            #                 except ValueError:
            # #                     print(f'  err: {nm+pfx_errpm[0]}   {nm+pfx_errpm[1]}  {row[nm+pfx_errpm[0]]}')
            #                     pass
            elif is_pm:
                e = row[var_names_errpm[nm]]
                scol += '\pm {e:{fe}}'.format(e=e, fe=fmt_err[nm])
            else:
                pass

            if math is not None:
                scol = '{math}{}{math}'.format(scol, math=math)
            #                 print(scol)

            if j < len(var_names) - 1:
                scol += '{sep}'.format(sep=sep)

            srow += scol

        s_res += '{}{}'.format(srow, str_end)
    return s_res


def pd2np(df):
    """
    Convert pandas dataframe to numpy array with column names
    """
    arr_ip = [tuple(i) for i in df.to_numpy()]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr
