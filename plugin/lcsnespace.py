##
#  Plot from file
##

import os
import numpy as np
import pystella as ps


def parse_arg(arg):
    fname = None
    tshift = 0.
    mshift = 0.
    marker = 'o'

    if len(arg) > 0:
        fname = os.path.expanduser(arg.pop(0))
    if len(arg) > 0:
        marker = arg.pop(0)
    if len(arg) > 0:
        tshift = float(arg.pop(0))
    if len(arg) > 0:
        mshift = float(arg.pop(0))

    return fname, marker, tshift, mshift


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
    tshift = dic.get('tshift', 0.)
    mshift = dic.get('mshift', 0.)
    marker = dic.get('marker', 'o')
    markersize = dic.get('markersize', 9)
    bnames = dic.get('bnames', None)
    bcolors = dic.get('bcolors', ps.band.colors())

    arg = dic.get('args', [])
    fname, marker, tshift, mshift = parse_arg(arg)

    print("Plot {0} [{1}]  jd_shift={2}  mshift={3}".format(fname, marker, tshift, mshift))

    curves = load(dic={'args': [fname, marker, str(tshift), str(mshift)]})

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
        x = lc.Time[is_good] + tshift
        y = lc.Mag[is_good] + mshift
        if lc.IsErr: # todo Make plot with errors
            yyerr = np.abs(lc.Err[is_good])
            # is_e = np.isnan(yyerr)
            ax.errorbar(x, y, label='{0} {1}'.format(bname, fname), yerr=yyerr, fmt=marker,
                        color=bcolors[bname], ls='')
            # ax.plot(x, y, label='{0} {1}'.format(bname, fname), color=bcolors[bname], ls='',
            #         marker=marker, markersize=markersize)
        else:
            ax.plot(x, y, label='{0} {1}'.format(bname, fname), color=bcolors[bname], ls='',
                    marker=marker, markersize=markersize)


def load(dic=None):
    """
    Load points from Sne Space file.json

    :return SetLightCurves:
    """
    fname, marker, tshift, mshift = parse_arg(dic['args'])
    # print("Loading SneSpace: {0} tshift={1}  mshift={2}".format(fname, tshift, mshift))

    obs = ps.SneSpace()
    if not obs.load(fname):
        print("No obs data from {}".format(fname))
        return None

    curves = obs.to_curves()
    timeMin = curves.TimeMin
    d = obs.comovingdist
    md = ps.rf.distance_modulus(d)
    print("Obs data loaded from {}. BandTimeMin= {} comovingdist= {} [MD={:.2f}]".format(fname, timeMin, d, md))

    if 'mag_lim' in dic:
        mag_lim = dic.get('mag_lim', 30.)
        # remove bad data
        res_curves = ps.SetLightCurve(curves.Name)
        for lc_orig in curves:
            is_good = lc_orig.Mag < mag_lim
            t = lc_orig.Time[is_good]
            m = lc_orig.Mag[is_good]
            e = None
            if lc_orig.IsErr:
                e = lc_orig.Err[is_good]
            lc = ps.LightCurve(lc_orig.Band, t, m, e)
            res_curves.add(lc)

        res_curves.set_tshift(tshift)
        res_curves.set_mshift(mshift)
    else:
        res_curves = curves
    return res_curves
