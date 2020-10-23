##
#  Plot from file
##

import os

from pystella.rf import band
from pystella.rf.lc import SetLightCurve, LightCurve
from pystella.util.reader_table import table2curves, read_obs_table_header


def lbl(b, band_shift, length=0):
    shift = band_shift[b]
    if shift == int(shift):
        shift = int(shift)
    # if shift == 0:
    #     return b

    s = b
    if shift > 0:
        s += '+'
        s += str(abs(shift))
    elif shift < 0:
        s += '-'
        s += str(abs(shift))

    if length > 0:
        s = ("{0:<" + str(length) + "s}").format(s)
    return s


def lbl_length(bshifts):
    return max((len(lbl(b, bshifts)) for b in bshifts.keys()))


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
    bcolors = dic.get('bcolors', band.colors())
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

    # read data
    tbl, cols_data = read_obs_table_header(fname, include_names=band.band_get_names_alias(),
                                           is_out=True, comments=comments)
    # tbl = read_table_header_float(fname)
    curves = table2curves(os.path.basename(fname), tbl)

    # filter bands
    if bnames is None:
        bnames = curves.BandNames
    else:
        bands = curves.BandNames
        for b in bands:
            if b not in bnames:
                curves.pop(b)
    # band shifts
    bshift = dic.get('bshift', None)
    if bshift is None:
        bshift = {b: 0. for b in bnames}
    else:
        bshiftzero = {b: 0. for b in bnames if b not in bshift}
        bshift.update(bshiftzero)

    lbl_len = lbl_length(bshift)
    # print(bshift, lbl_len)

    # plot data
    for lc in curves:
        bname = lc.Band.Name
        is_good = lc.Mag < (mag_lim - mshift)
        x = lc.Time[is_good] + jd_shift
        y = lc.Mag[is_good] + mshift + bshift[bname]
        if lc.IsErr:
            yyerr = abs(lc.Err[is_good])
            ax.errorbar(x, y, label='{0} {1}'.format(lbl(bname, bshift, lbl_len), fname), yerr=yyerr, fmt=marker,
                        color=bcolors[bname], ls='')
        else:
            ax.plot(x, y, label='{0} {1}'.format(lbl(bname, bshift, lbl_len), fname), color=bcolors[bname], ls='',
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
    tbl, cols_data = read_obs_table_header(fname, include_names=band.band_get_names_alias(), is_out=is_debug)
    curves = table2curves(os.path.basename(fname), tbl)

    # remove bad data
    res_curves = SetLightCurve(curves.Name)
    for lc_orig in curves:
        is_good = lc_orig.Mag < mag_lim
        t = lc_orig.Time[is_good] + tshift
        m = lc_orig.Mag[is_good] + mshift
        e = None
        if lc_orig.IsErr:
            e = lc_orig.Err[is_good]
        lc = LightCurve(lc_orig.Band, t, m, e)
        res_curves.add(lc)

    # res_curves.set_tshift(tshift)
    # res_curves.set_mshift(mshift)
    return res_curves
