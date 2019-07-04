##
#  Same as lcobs except load data files.
# Comments = !
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

    if isinstance(ax, (list, tuple)):
        ax = ax[0]

    if dic is None:
        dic = {}
    fname = dic.get('fname', None)
    jd_shift = dic.get('jd_shift', 0.)
    mshift = dic.get('mshift', 0.)
    marker = dic.get('marker', 'o')
    markersize = dic.get('markersize', 5)
    bnames = dic.get('bnames', None)
    bcolors = dic.get('bcolors', band.bands_colors())
    comments = dic.get('comments', '#')

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

    # read header
    with open(fname) as f:
        a = f.readline().split()
        header = ' '.join(a[1:])  # to remove #
    # read data
    tbl, cols_data = read_obs_table_header(fname, header=header, include_names=band.band_get_names_alias(),
                                           is_out=True, comments=comments)
    # tbl = read_table_header_float(fname)
    curves = table2curves(os.path.basename(fname), tbl)

    # filter bands
    if bnames is not None:
        bands = curves.BandNames
        for b in bands:
            if b not in bnames:
                curves.pop(b)
    # plot data
    for lc in curves:
        bname = lc.Band.Name
        is_good = lc.Mag < (mag_lim - mshift)
        x = lc.Time[is_good] + jd_shift
        y = lc.Mag[is_good] + mshift
        if lc.IsErr:
            yyerr = abs(lc.Err[is_good])
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
    from pystella.util.math import is_number

    fname = None
    tshift = 0.
    mshift = 0.
    mag_lim = 99.
    arg = []
    is_debug = False
    comments = '#'

    if dic is not None:
        is_debug = dic.get('is_debug', False)
        mag_lim = dic.get('mag_lim', mag_lim)
        if 'args' in dic:
            arg = dic['args']

    if len(arg) > 0:
        fname = arg.pop(0)
        fname = os.path.expanduser(fname)
    if len(arg) > 0:
        s = arg.pop(0)
        if is_number(s):
            tshift = float(s)
        elif len(arg) > 0:
            tshift = float(arg.pop(0))
    if len(arg) > 0:
        mshift = float(arg.pop(0))

    print("Load {0} tshift={1}  mshift={2}".format(fname, tshift, mshift))

    # read data
    # tbl = read_table_header_float(fname)
    tbl, cols_data = read_obs_table_header(fname, include_names=band.band_get_names_alias(),
                                           is_out=is_debug, comments=comments)
    curves = table2curves(os.path.basename(fname), tbl)

    # remove bad data
    res_curves = SetLightCurve(curves.Name)
    for lc_orig in curves:
        is_good = lc_orig.Mag < mag_lim
        t = lc_orig.Time[is_good]
        m = lc_orig.Mag[is_good]
        e = None
        if lc_orig.IsErr:
            e = lc_orig.Err[is_good]
        lc = LightCurve(lc_orig.Band, t, m, e)
        res_curves.add(lc)

    res_curves.set_tshift(tshift)
    res_curves.set_mshift(mshift)
    return res_curves


def load_curves(fname, skiprows=1):
    from pystella import curves_read_mix
    return curves_read_mix(fname, skiprows)
