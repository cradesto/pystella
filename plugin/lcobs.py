##
#  Plot from file
##

import os

from pystella.rf import band
from pystella.util.reader_table import read_table_header_float, table2curves


def plot(ax, dic=None):
    """
    Plot points from dat-files. Format fname:marker:jd_shift:mshift
    Header should be:  time  U [errU]  B [errB]  ...
    :param ax:
    :param dic:
    :return:
    """
    marker_size = 8
    colors = band.bands_colors()
    fname = None
    jd_shift = 0.
    mshift = 0.
    marker = 'o'
    arg = []
    if dic is not None and 'args' in dic:
        arg = dic['args']

    if len(arg) > 0:
        fname = arg.pop(0)
        fname = os.path.expanduser(fname)
    if len(arg) > 0:
        marker = arg.pop(0)
    if len(arg) > 0:
        jd_shift = float(arg.pop(0))
    if len(arg) > 0:
        mshift = float(arg.pop(0))

    print("Plot {0} [{1}]  jd_shift={2}  mshift={3}".format(fname, marker, jd_shift, mshift))

    # read data
    tbl = read_table_header_float(fname)
    curves = table2curves(os.path.basename(fname), tbl)

    # plot data
    for lc in curves:
        bname = lc.Band.Name
        x = lc.Time + jd_shift
        y = lc.Mag + mshift
        if lc.IsErr:
            yyerr = abs(lc.MagErr)
            ax.errorbar(x, y, label='{0} {1}'.format(bname, fname), yerr=yyerr, fmt=marker,
                        color=colors[bname], ls='')
        else:
            ax.plot(x, y, label='{0} {1}'.format(bname, fname), color=colors[bname], ls='',
                    marker=marker, markersize=marker_size)


def load(dic=None):
    """
    Load points from dat-files.
    :return SetLightCurves:
    """
    fname = None
    jd_shift = 0.
    mshift = 0.
    arg = []
    if dic is not None and 'args' in dic:
        arg = dic['args']

    if len(arg) > 0:
        fname = arg.pop(0)
        fname = os.path.expanduser(fname)
    if len(arg) > 0:
        jd_shift = float(arg.pop(0))
    if len(arg) > 0:
        mshift = float(arg.pop(0))

    print("Load {0} jd_shift={1}  mshift={2}".format(fname, jd_shift, mshift))

    # read data
    tbl = read_table_header_float(fname)
    curves = table2curves(os.path.basename(fname), tbl)
    curves.tshift = jd_shift
    curves.mshift = mshift
    return curves

