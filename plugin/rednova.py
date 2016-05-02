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

    if len(arg) > 0:
        jd_shift = float(arg.pop(0))
    else:
        jd_shift = 0.

    print "Plot Red Nova  jd_shift=%f  path: %s" % (jd_shift, sn_path)
    plot_ubv(ax=ax, path=sn_path, jd_shift=jd_shift)


def plot_ubv(ax, path, jd_shift=0., mshift=0.):
    colors = band.bands_colors()
    curves = read_curves(path)
    for lc in curves:
        x = lc.Time + jd_shift
        y = lc.Mag + mshift
        bcolor = colors[lc.Band.Name]
        ax.plot(x, y, label='%s SN Red Nova' % lc.Band.Name,
                ls=".", color=bcolor, markersize=8, marker="o")


def read_curves(path=sn_path):
    header = 'V  I  R'
    bnames = map(str.strip, header.split())
    curves = SetLightCurve('Red Nova')
    for i, n in enumerate(bnames):
        b = band.band_by_name(n)
        time = np.loadtxt(os.path.join(path, n+'_jd.txt'))
        mags = np.loadtxt(os.path.join(path, n+'_mag.txt'))
        errs = np.loadtxt(os.path.join(path, n+'_err.txt'))

        # filter bad values
        # is_good = mags < 30.
        # t = time[is_good]
        # mags = mags[is_good]
        # errs = errs[is_good]
        # add
        lc = LightCurve(b, time, mags, errs=errs)
        curves.add(lc)

    return curves
