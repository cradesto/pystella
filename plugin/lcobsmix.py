##
#  Plot from file
##

import os

import numpy as np

from pystella.rf import band
from pystella.rf.lc import SetLightCurve, LightCurve


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
    # tbl, cols_data = read_obs_table_header(fname, include_names=band.band_get_names_alias(), is_out=True)
    # # tbl = read_table_header_float(fname)
    # curves = table2curves(os.path.basename(fname), tbl)
    curves = load(dic={'arg': [fname]})

    # filter bands
    if bnames is not None:
        bands = curves.BandNames
        for b in bands:
            if b not in bnames:
                curves.rm(b)
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
    Reader data-file with mix bname data, like:
    >>% head photometry.txt
    jd filter mag mage
    2457059.6228778586 V 17.493766309999998 0.0592200135089
    2457059.6244578934 V 17.539956019999998
    0.0542402986717 2457059.6261980557 g 17.782871193345898
    0.0454000142503 2457059.6287036575 g 17.7782469177482 0.0395424488201

    :return SetLightCurves:
    """
    fname = None
    tshift = 0.
    mshift = 0.
    mag_lim = 30.
    skiprows = 1.
    arg = []

    if dic is not None:
        mag_lim = dic.get('mag_lim', 30.)
        skiprows = dic.get('skiprows', 1)
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
    dtype = [('time', '<f4'), ('b', 'S1'), ('mag', '<f4'), ('err', '<f4')]
    lc_data = np.loadtxt(fname, skiprows=skiprows, dtype=dtype, comments='#')  # jd filter mag mage
    b_tot = lc_data['b']
    bnames = np.unique(b_tot)

    curves = SetLightCurve()
    for bname in bnames:
        if band.band_is_exist(bname):
            # filter of the current band
            is_good = np.array(map(lambda x: x == bname, b_tot), dtype=bool)
            t = lc_data['time'][is_good]
            m = lc_data['mag'][is_good]
            e = lc_data['err'][is_good]
            # add light curve
            b = band.band_by_name(bname)
            lc = LightCurve(b, t, m, e)
            curves.add(lc)
        else:
            print('Could read the light curve. There is no band: {}. '
                  'You may try to add it to dir data/bands'.format(bname))

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
