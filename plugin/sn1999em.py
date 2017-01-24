##
#  Data of SN 1999em
##

import numpy as np
import os

from pystella.rf import band
from pystella.rf.lc import SetLightCurve, LightCurve

sn_path = os.path.expanduser('~/Sn/Release/svn_kepler/stella/branches/lucy/run/res/sncurve/sn1999em')


def plot(ax, dic=None):
    jd_shift = 0.  # -2451482.  # for date 30/10/1999
    m_shift = 0.  # DM = 29.4  d=7.8 Mpc
    arg = []
    if dic is not None and 'args' in dic:
        arg = dic['args']

    if len(arg) > 0:
        jd_shift = float(arg.pop(0))
    if len(arg) > 0:
        m_shift = float(arg.pop(0))

    print("Plot Sn 1987A  jd_shift=%f  mshift=%f path: %s" % (jd_shift, m_shift, sn_path))
    plot_ubv(ax=ax, path=sn_path, jd_shift=jd_shift, mshift=m_shift)


def plot_ubv(ax, path, jd_shift=0., mshift=0.):
    colors = band.bands_colors()
    curves = read_curves(path)
    for lc in curves:
        x = lc.Time + jd_shift
        y = lc.Mag + mshift
        bcolor = colors[lc.Band.Name]
        ax.plot(x, y, label='%s SN 1999em' % lc.Band.Name,
                ls="", color=bcolor, markersize=8, marker="o")


def read_curves(path=sn_path):
    header = 'U B  V  I  R'
    cols = map(str.strip, header.split())
    lc_data = np.loadtxt(os.path.join(path, '1999emubvir.dat'), skiprows=1, usecols=range(1, 12))
    time = lc_data[:, 0]
    curves = SetLightCurve('Sn99em')
    for i, n in enumerate(cols):
        b = band.band_by_name(n)
        mags = lc_data[:, 2*i+1]
        errs = lc_data[:, 2*i+2]
        # filter bad values
        is_good = mags < 30.
        t = time[is_good]
        mags = mags[is_good]
        errs = errs[is_good]
        # add
        lc = LightCurve(b, t, mags, errs=errs)
        curves.add(lc)

    return curves
