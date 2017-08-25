##
#  Plot from file
##

import os

from pystella.rf import band
from pystella.rf.lc import SetLightCurve, LightCurve
from pystella.util.reader_table import table2curves, read_obs_table_header


def plot(ax, dic=None, mag_lim=30.):
    """
    Plot points from dat-files. Format fname:marker:jd_shift:mshift
    Header should be:  time  U [errU]  B [errB]  ...
    :param ax:
    :param dic:
    :param mag_lim:
    :return:
    """
    # colors = band.bands_colors()
    if dic is None:
        dic = {}
    fname = dic.get('fname', None)
    jd_shift = dic.get('jd_shift', 0.)
    mshift = dic.get('mshift', 0.)
    marker = dic.get('marker', 'o')
    markersize = dic.get('markersize', 9)
    bnames = dic.get('bnames', None)
    bcolors = dic.get('bcolors', band.bands_colors())

    arg = dic.get('args', [])
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
    tbl = read_obs_table_header(fname)
    # tbl = read_table_header_float(fname)
    curves = table2curves(os.path.basename(fname), tbl)

    # filter bands
    if bnames is not None:
        bands = curves.BandNames
        for b in bands:
            if b not in bnames:
                curves.rm(b)
    # plot data
    for lc in curves:
        bname = lc.Band.Name
        is_good = lc.Mag < (mag_lim-mshift)
        x = lc.Time[is_good] + jd_shift
        y = lc.Mag[is_good] + mshift
        if lc.IsErr:
            yyerr = abs(lc.MagErr[is_good])
            ax.errorbar(x, y, label='{0} {1}'.format(bname, fname), yerr=yyerr, fmt=marker,
                        color=bcolors[bname], ls='')
        else:
            ax.plot(x, y, label='{0} {1}'.format(bname, fname), color=bcolors[bname], ls='',
                    marker=marker, markersize=markersize)


def load(dic=None):
    """
    Load points from dat-files.
    :return SetLightCurves:
    """
    fname = None
    jd_shift = 0.
    mshift = 0.
    mag_lim = 30.
    arg = []
    if dic is not None:
        is_debug = dic.get('is_debug', False)
        mag_lim = dic.get('mag_lim', 30.)
        if 'args' in dic:
            arg = dic['args']

    if len(arg) > 0:
        fname = arg.pop(0)
        fname = os.path.expanduser(fname)
    if len(arg) > 0:
        s = arg.pop(0)
        if s.isnumeric():
            jd_shift = float(s)
        elif len(arg) > 0:
            jd_shift = float(arg.pop(0))
    if len(arg) > 0:
        mshift = float(arg.pop(0))

    print("Load {0} jd_shift={1}  mshift={2}".format(fname, jd_shift, mshift))

    # read data
    # tbl = read_table_header_float(fname)
    tbl = read_obs_table_header(fname,is_out=is_debug)
    curves = table2curves(os.path.basename(fname), tbl)

    # remove bad data
    res_curves = SetLightCurve(curves.Name)
    for lc_orig in curves:
        is_good = lc_orig.Mag < mag_lim
        t = lc_orig.Time[is_good]
        m = lc_orig.Mag[is_good]
        e = None
        if lc_orig.IsErr:
            e = lc_orig.MagErr[is_good]
        lc = LightCurve(lc_orig.Band, t, m, e)
        res_curves.add(lc)

    res_curves.set_tshift(jd_shift)
    res_curves.set_mshift(mshift)
    return res_curves
