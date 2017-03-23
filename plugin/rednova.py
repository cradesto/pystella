##
#  Data of Red Nova
##

import numpy as np
import os

from pystella.rf import band
from pystella.rf.lc import SetLightCurve, LightCurve

sn_path = os.path.expanduser('~/Sn/Release/svn_kepler/stella/branches/lucy/run/res/sncurve/rednovaM31')


def plot(ax, dic=None):
    arg = []
    if dic is not None and 'args' in dic:
        arg = dic['args']

    jd_shift = 0.
    if len(arg) > 0:
        jd_shift = float(arg.pop(0))

    print("Plot Red Nova  jd_shift=%f  path: %s" % (jd_shift, sn_path))
    plot_ubv(ax=ax, path=sn_path, jd_shift=jd_shift)


def plot_ubv(ax, path, jd_shift=0., mshift=0.):
    path = os.path.expanduser(path)

    colors = band.bands_colors()
    curves = read_curves_master(path)
    for lc in curves:
        x = lc.Time + jd_shift
        y = lc.Mag  # todo + mshift
        bcolor = colors[lc.Band.Name]
        ax.plot(x, y, label='%s Master' % lc.Band.Name,
                ls=":", color=bcolor, markersize=7, marker="o")
        ax.errorbar(x, y, yerr=lc.MagErr, color='gray', fmt='.', zorder=1)

    print("jd_shift=%f mshift=%f " % (jd_shift, mshift))

    curves = read_curves_kurtenkov(path)
    for lc in curves:
        x = lc.Time + jd_shift
        y = lc.Mag + mshift
        bcolor = colors[lc.Band.Name]
        ax.plot(x, y, label='%s Kurtenkov' % lc.Band.Name,
                ls=":", color=bcolor, markersize=7, marker="*")
        ax.errorbar(x, y, yerr=lc.MagErr, color='gray', fmt='.', zorder=1)


def read_curves_master_abs_mag(path=sn_path):
    jd = 2457036
    header = 'V  I  R'
    bnames = map(str.strip, header.split())
    curves = SetLightCurve('Red Nova')
    for i, n in enumerate(bnames):
        b = band.band_by_name(n)
        time = np.loadtxt(os.path.join(path, n + '_jd.txt')) + jd
        mags = np.loadtxt(os.path.join(path, n + '_mag.txt'))
        errs = np.loadtxt(os.path.join(path, n + '_err.txt'))

        # filter bad values
        # is_good = mags < 30.
        # t = time[is_good]
        # mags = mags[is_good]
        # errs = errs[is_good]
        # add
        lc = LightCurve(b, time, mags, errs=errs)
        curves.add(lc)

    return curves


def read_curves_master(path=sn_path):
    header = 'V  I  R'
    bnames = map(str.strip, header.split())
    curves = SetLightCurve('Red Nova')
    for i, n in enumerate(bnames):
        b = band.band_by_name(n)
        d = np.loadtxt(os.path.join(path, n + '.txt'),
                       dtype=[('JD', '<f4'), ('mag', '<f4'), ('err', '<f4')])
        time = d['JD']
        mags = d['mag']
        errs = d['err']
        lc = LightCurve(b, time, mags, errs=errs)
        curves.add(lc)

    return curves


def read_curves_kurtenkov(path=sn_path):
    jd = 2457000
    lc_data = np.loadtxt(os.path.join(path, 'lrn_aa26564-15_p5.csv'), skiprows=2, usecols=(0, 1, 2, 3),
                         dtype=[('JD', '<f4'), ('b', '|S1'), ('mag', '<f4'), ('err', '<f4')])

    # mshift = 24.43  # Distance Module to M31
    curves = SetLightCurve('Red Nova')

    bnames = np.unique(lc_data['b'])
    for i, n in enumerate(bnames):
        b = band.band_by_name(n.decode("utf-8"))
        d = lc_data[lc_data['b'] == n, ]
        time = d['JD'] + jd
        mags = d['mag']
        errs = d['err']

        # filter bad values
        # is_good = mags < 30.
        # t = time[is_good]
        # mags = mags[is_good]
        # errs = errs[is_good]
        # add
        lc = LightCurve(b, time, mags, errs=errs)
        # lc.tshift = -tshift
        # lc.mshift = -mshift
        curves.add(lc)

    return curves
