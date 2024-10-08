import csv

import numpy as np

__author__ = 'bakl'


def str2bool(v):
    return v.lower() in ("yes", "true", "t", 'y', "1")


def str2interval(arg, llim=0., rlim=float('inf'), sep=':'):
    if arg is None or len(arg) == 0:
        # xlim = (llim, rlim)
        return None

    if sep in arg:
        if str.endswith(arg, sep):
            xlim = (float(arg[:-1]), rlim)
        elif str.startswith(arg, sep):
            xlim = (llim, float(arg[1:]))
        else:
            xlim = tuple(map(float, arg.split(sep)))
    else:
        xlim = (llim, float(arg))
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
    dtype = np.dtype({'names': names, 'formats': [float] * len(names)})
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


#
#
# def df2tex(a, pfx_err='err', pfx_errpm=('errp', 'errm'), fmt='.3f', fmt_err='.4f', math=' $ ', sep=' & ',
#            str_end=' \\\ \n'):
#     """
#     Convert pandas dataframe to latex table.
#     Usage:
#     arr = pd2np(df_fit_images)
#     s = df2tex(arr, pfx_err=None, fmt_err='.2f', sep='  &  ', strend='   \\\ \n')
#     print(s)
#     """
#     names = a.dtype.names
#     if not isinstance(fmt, dict):
#         fmt = {nm: fmt for nm in names}
#     if not isinstance(fmt_err, dict):
#         fmt_err = {nm: fmt_err for nm in names}
#
#     is_pm = pfx_err is not None
#     is_up_down = not is_pm and pfx_errpm is not None
#
#     var_names_errp, var_names_errm, var_names_errpm = {}, {}, {}
#     if is_up_down:
#         var_names = [x for x in names if not pfx_errpm[0] in x and not pfx_errpm[1] in x]
#         var_names_errp = {xx: x for xx in var_names for x in names if xx in x and pfx_errpm[0] in x}
#         var_names_errm = {xx: x for xx in var_names for x in names if xx in x and pfx_errpm[1] in x}
#     elif is_pm:
#         var_names = [x for x in names if not pfx_err in x]
#         var_names_errpm = {xx: x for xx in var_names for x in names if xx in x and pfx_err in x}
#     else:
#         var_names = names
#     #     print(var_names_errm)
#     s_res = sep.join([nm.replace('_', '\_') for nm in var_names]) + str_end
#
#     for i, row in enumerate(a):
#         srow = ''
#         for j, nm in enumerate(var_names):
#             scol = ''
#             x = row[nm]
#             try:
#                 scol += '{:{fmt}}'.format(x, fmt=fmt[nm])
#             except ValueError:
#                 scol += '{}'.format(x.replace('_', '\_'))
#             # print uncertanties
#             if is_up_down and nm in var_names_errp:
#                 #                 print(nm, var_names_errp[nm], var_names_errm[nm])
#                 ep, em = row[var_names_errp[nm]], row[var_names_errm[nm]]
#                 #                 try:
#                 scol += '^{{{ep:{fe}}}}_{{{em:{fe}}}}'.format(ep=ep, em=em, fe=fmt_err[nm])
#             #                 except ValueError:
#             # #                     print(f'  err: {nm+pfx_errpm[0]}   {nm+pfx_errpm[1]}  {row[nm+pfx_errpm[0]]}')
#             #                     pass
#             elif is_pm:
#                 e = row[var_names_errpm[nm]]
#                 scol += '\pm {e:{fe}}'.format(e=e, fe=fmt_err[nm])
#             else:
#                 pass
#
#             if math is not None:
#                 scol = '{math}{}{math}'.format(scol, math=math)
#             #                 print(scol)
#
#             if j < len(var_names) - 1:
#                 scol += '{sep}'.format(sep=sep)
#
#             srow += scol
#
#         s_res += '{}{}'.format(srow, str_end)
#     return s_res


def pd2np(df):
    """
    Convert pandas dataframe to numpy array with column names
    :param df: dataframe
    :return:  numpy array
    """
    arr_ip = [tuple(i) for i in df.to_numpy()]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr


def np2tex(a, sfx_err='err', sfx_errpm=('errp', 'errm'), fmt='.3f', fmt_err='.4f', math=' $ ', sep=' & ',
           str_end=r' \\', is_info=True):
    """
        Convert numpy array  to latex table.
        Add uncertainties if  exist the column with uncertainties: col+sfx_err or col+sfx_errpm
        Usage:
        arr = pd2np(df) # if df is pandas dataframe
        s = np2tex(arr, sfx_err=None, fmt_err='.2f', sep='  &  ', str_end='   \\')
        print(s)

        :type a: numpy array with dtype.names
        :param sfx_err: str, suffix for uncertainties: COL /pm COL+sfx_err. Default: sfx_err='err'
        :param sfx_errpm: tuple of 2 strings, suffix for uncertainties: COL^{COL+sfx_errpm[0]}_{COL+sfx_errpm[1]}.
                          Default: sfx_errpm=('errp', 'errm')
        :param fmt: dict or str. If fmt is type = dict, it should contains formats for all columns: COLs and COL_ERRs.
                    If fmt is string, when fmt is used for COLs and fmt_err is used for COL_ERRs
        :type fmt_err: str. Default:  '.4f'
        :type math: str.  Default:  ' $ '
        :type sep: str.  Default:  ' & '
        :type str_end: str.  Default:  r' \\'
        :type is_info: bool.  Print names and fmt. Default: True
        :return: string array with text of latex table
        """
    """
    Convert pandas dataframe to latex table.
    
    """
    names = a.dtype.names
    if not isinstance(fmt, dict):
        fmt = {nm: fmt for nm in names}
        if not isinstance(fmt_err, dict):
            fmt_err = {nm: fmt_err for nm in names}
    else:
        fmt_err = fmt.copy()

    if is_info:
        print(names)
        print('fmt_err: ', fmt_err)

    var_names_errp, var_names_errm, var_names_err = {}, {}, {}

    if sfx_err is not None:
        var_names_err = {xx: x for xx in names for x in names if
                         not xx.endswith(sfx_err) and xx in x and x.endswith(sfx_err)}
    if sfx_errpm is not None:
        var_names_errp = {xx: x for xx in names for x in names if
                          not xx.endswith(sfx_errpm[0]) and xx in x and x.endswith(sfx_errpm[0])}
        var_names_errm = {xx: x for xx in names for x in names if
                          not xx.endswith(sfx_errpm[1]) and xx in x and x.endswith(sfx_errpm[1])}

    # Columns without uncertainties
    z = list(var_names_err.values()) + list(var_names_errp.values()) + list(var_names_errm.values())
    #     print('z', z)
    var_names = []
    for nm in names:
        if nm in z:
            continue
        var_names.append(nm)
    #     print('var_names', var_names)

    s_res = [sep.join([nm.replace('_', r'\_') for nm in var_names]) + str_end]

    for i, row in enumerate(a):
        srow = ''
        for j, nm in enumerate(var_names):
            scol = ''
            x = row[nm]
            # if nm == 'name':
            #     x = f'M{i + 1}'
            try:
                scol += '{:{fmt}}'.format(x, fmt=fmt[nm])
            except ValueError:
                scol += '{}'.format(x.replace('_', r'\_'))
            # print uncertainties
            if nm in var_names_errp and nm in var_names_errm:
                #                 print(nm, var_names_errp[nm], var_names_errm[nm])
                ep, em = row[var_names_errp[nm]], row[var_names_errm[nm]]
                scol += '^{{+{ep:{fe}}}}_{{-{em:{fe}}}}'.format(ep=ep, em=em, fe=fmt_err[var_names_errp[nm]])
            #                 except ValueError:
            # #                     print(f'  err: {nm+sfx_errpm[0]}   {nm+sfx_errpm[1]}  {row[nm+sfx_errpm[0]]}')
            #                     pass
            elif nm in var_names_err:
                e = row[var_names_err[nm]]
                scol += r'\pm {e:{fe}}'.format(e=e, fe=fmt_err[var_names_err[nm]])
            else:
                pass

            if math is not None:
                scol = '{math}{}{math}'.format(scol, math=math)

            if j < len(var_names) - 1:
                scol += '{sep}'.format(sep=sep)

            srow += scol

        s_res.append('{}{}'.format(srow, str_end))
    return s_res
