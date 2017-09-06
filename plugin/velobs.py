##
#  Plot from file
##

import os

from pystella.util.reader_table import read_obs_table_header
from pystella.velocity import VelocityCurve


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
    tbl = read_obs_table_header(fname, is_out=True)
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
    cname = dic.get('cname', 'Vel')

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
    tbl = read_obs_table_header(fname, colt=colt, include_names=[cname], is_out=is_debug)

    vel_o = tbl2vel(tbl, cname, colt, mshift)
    vel_o.tshift = tshift
    return vel_o


def tbl2vel(tbl, cname, colt, mshift):
    for nm in colt:
        if nm in tbl.dtype.names:
            time = tbl[nm]
            break
    else:
        raise ValueError("THe table should contain a column with name in {0}", ', '.join(colt))
    values = tbl[cname] * mshift
    vel_o = VelocityCurve('Vel', time, values)
    return vel_o
