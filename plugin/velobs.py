##
#  Plot from file
##

import os

import pystella.util.reader_table  as rtbl
from pystella.velocity import VelocityCurve, SetVelocityCurve


def plot(ax, dic=None, colt=('time', 'JD', 'MJD')):
    """
    Plot points from dat-files. Format fname:marker:jd_shift:mshift
    Header should be:  time  U [errU]  B [errB]  ...
    :param ax:
    :param dic:
    :param colt:
    :return:
    """
    if dic is None:
        dic = {}
    cname = dic.get('cname', 'Vel')
    fname = dic.get('fname', None)
    tshift = dic.get('tshift', 0.)
    mshift = dic.get('mshift', 0.)
    marker = dic.get('marker', 'o')
    markersize = dic.get('markersize', 9)
    color = dic.get('color', 'blue')

    arg = dic.get('args', [])
    if len(arg) > 0:
        fname = arg.pop(0)
        fname = os.path.expanduser(fname)
    if len(arg) > 0:
        marker = arg.pop(0)
    if len(arg) > 0:
        tshift = float(arg.pop(0))
    if len(arg) > 0:
        mshift = float(arg.pop(0))

    print("Plot {0} [{1}]  jd_shift={2}  mshift={3}".format(fname, marker, tshift, mshift))

    # read data
    tbl, cols_data = rtbl.read_obs_table_header(fname, is_out=True)
    vel_o = tbl2vel(tbl, cname, colt, mshift)
    vel_o.tshift = tshift

    x = vel_o.Time
    y = vel_o.Vel

    if vel_o.IsErr:
        yyerr = abs(vel_o.Err)
        ax.errorbar(x, y, label='{0} {1}'.format(vel_o.Name, fname), yerr=yyerr, fmt=marker,
                    color=color, ls='')
    else:
        ax.plot(x, y, label='{0} {1}'.format(vel_o.Name, fname), color=color, ls='',
                marker=marker, markersize=markersize)


def load(dic=None, colt=('time', 'JD', 'MJD')):
    """
    Load points from dat-files.
    :return SetLightCurves:
    """
    if dic is None:
        dic = {}
    is_debug = dic.get('is_debug', False)
    cnames = ('Vel', 'vHeI', 'vFeII', 'vHa', 'vHb', 'vHg')
    # cnames = ('Vel', 'vHeI', 'vFeII', 'vHa', 'vHb', 'vHg', 'Vel*', 'vel*', )
    cpatterns = ('Vel*', 'vel*')

    fname = None
    tshift = 0.
    mshift = 1.
    arg = []

    if 'args' in dic:
        arg = dic['args']

    if len(arg) > 0:
        fname = arg.pop(0)
        fname = os.path.expanduser(fname)
    if len(arg) > 0:
        s = arg.pop(0)
        if s.isnumeric():
            tshift = float(s)
        elif len(arg) > 0:
            tshift = float(arg.pop(0))
    if len(arg) > 0:
        mshift = float(arg.pop(0))

    print("Load {0} tshift={1}  mshift={2}".format(fname, tshift, mshift))

    # read data
    # tbl = read_table_header_float(fname)
    tbl, cols_data = rtbl.read_obs_table_header(fname, colt=colt, include_names=cnames,
                                     include_patterns=cpatterns, is_out=is_debug)

    vels_o = SetVelocityCurve("Vel-{}".format(fname))
    for i, cname in cols_data.items():
        vel_o = tbl2vel(tbl, cname, colt, mshift)
        vels_o.add(vel_o)

    vels_o.set_tshift(tshift)
    return vels_o


def tbl2vel(tbl, cname, colt, mshift):
    for nm in colt:
        if nm in tbl.dtype.names:
            time = tbl[nm]
            break
    else:
        raise ValueError("The table should contain a column with name in [{0}]".format(', '.join(colt)))
    values = tbl[cname] * mshift
    # filter
    mask = values > 0
    t = time[mask]
    v = values[mask]

    for err_name in (prefix + cname for prefix in rtbl.err_prefix):
        if err_name in tbl.dtype.names:
            err = tbl[err_name]
            e = err[mask]
            vel_o = VelocityCurve(cname, t, v, e)
            break
    else:
        vel_o = VelocityCurve(cname, t, v)

    return vel_o
