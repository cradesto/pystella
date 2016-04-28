##
#  Callbacks
##

import numpy as np
import os

from pystella.rf import band
from pystella.rf.lc import SetLightCurve, LightCurve

sn_path = os.path.expanduser('~/Sn/Release/svn_kepler/stella/branches/lucy/run/res/sncurve/sn1987a')


def plot(ax, dic=None):
    arg = []
    if dic is not None and 'args' in dic:
        arg = dic['args']

    if len(arg) > 0:
        jd_shift = float(arg.pop(0))
    else:
        jd_shift = -2446850.

    print "Plot Sn 1987A  jd_shift=%f  path: %s" % (jd_shift, sn_path)
    plot_ubv(ax=ax, path=sn_path, jd_shift=jd_shift)


def plot_vels_sn87a(ax, path, z=0):
    print "Plot the velocities of Sn 87A "

    jd_shift = 2446850  # moment of explosion SN 1987A, Hamuy 1988, doi:10.1086/114613

    # Blanco's data from plot
    fs = {'Halpha':   'Halpha_blanco.csv',   'Hbeta': 'Hbeta_blanco.csv',
          'Hgamma':   'Hgamma_blanco.csv',   'NaID': 'NaID_blanco.csv',
          'FeII5018': 'FeII5018_blanco.csv', 'FeII5169': 'FeII5169_blanco.csv'
          }
    fs = dict((k, os.path.join(path, v)) for k, v in fs.items())
    elcolors = {'Halpha': "black", 'Hbeta': "cyan", 'Hgamma': "orange", 'NaID': "red", 'FeII5018': "orange",
                'FeII5169': "magenta"}
    elmarkers = {'Halpha': u's', 'Hbeta': u'x', 'Hgamma': u'd', 'NaID': u'+', 'FeII5018': u'D', 'FeII5169': u'o'}

    for el, fname in fs.items():
        data = np.loadtxt(fname, comments='#')
        x = data[:, 0] - jd_shift
        x *= 1. + z  # redshift
        y = data[:, 1]
        ax.plot(x, y, label='%s, SN 87A' % el, ls=".", color=elcolors[el], markersize=6, marker=elmarkers[el])


def plot_ubv(ax, path, jd_shift=0., mshift=0.):
    colors = band.bands_colors()
    curves = read(path)
    for lc in curves:
        x = lc.Time + jd_shift
        y = lc.Mag + mshift
        bcolor = colors[lc.Band.Name]
        ax.plot(x, y, label='%s SN 1987A' % lc.Band.Name,
                ls=".", color=bcolor, markersize=6, marker=".")


def read(path=sn_path):
    fs = {'U': 'pav_U.csv', 'B': 'pav_B.csv', 'V': 'pav_V.csv', 'R': 'pav_R.csv', 'I': 'pav_I.csv'}
    fs = dict((k, os.path.join(path, v)) for k, v in fs.items())

    curves = SetLightCurve('Sn87A')
    for n, fname in fs.items():
        data = np.loadtxt(fname, comments='#', skiprows=1, delimiter="\t", usecols=(0, 1))
        b = band.band_by_name(n)
        t = data[:, 0]
        mags = data[:, 1]
        lc = LightCurve(b, t, mags)
        curves.add(lc)

    return curves
