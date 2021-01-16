##
#  Data of Red Nova
##

import numpy as np
import os

from pystella.rf import band
from pystella.rf.lc import SetLightCurve, LightCurve
from pystella.util.reader_table import read_table_header_float, table2curves

sn_path = os.path.expanduser('~/Sn/Release/svn_kepler/stella/branches/lucy/run/res/sncurve/rednovaM31')


def plot(ax, dic=None):
    arg = []
    bnames = None
    if dic is not None and 'args' in dic:
        arg = dic['args']
        # bnames = dic.get('bnames', None)
        # bcolors = dic.get('bcolors', None)

    jd_shift = 0.
    if len(arg) > 0:
        jd_shift = float(arg.pop(0))

    print("Plot Red Nova  jd_shift=%f  path: %s" % (jd_shift, sn_path))
    plot_ubv(ax=ax, path=sn_path, jd_shift=jd_shift, **dic)


def plot_ubv(ax, path, jd_shift=0., mshift=0., **kwargs):
    markersize = kwargs.pop('markersize', 7)
    is_mas = kwargs.pop('is_mas', True)
    is_kur = kwargs.pop('is_kur', True)
    is_wil = kwargs.pop('is_wil', True)
    bnames_fix = kwargs.pop('bnames', None)
    colors = kwargs.pop('bcolors', band.colors())
    path = os.path.expanduser(path)

    print("Plot RedNova: jd_shift=%f mshift=%f " % (jd_shift, mshift))

    def plot_lc(cs, band_names, lbl, fmt):
        if band_names is None:
            band_names = cs.BandNames
        for lc in cs:
            if lc.Band.Name in band_names:
                x = lc.Time + jd_shift
                y = lc.Mag + mshift
                bcolor = colors[lc.Band.Name]
                # ax.plot(x, y, label='%s Master' % lc.Band.Name,
                #         ls="", color=bcolor, markersize=markersize, marker="o")
                ax.errorbar(x, y, yerr=lc.Err, label=lbl % lc.Band.Name,
                            color=bcolor, fmt=fmt, fillstyle='none', ecolor=bcolor,
                            mec=bcolor, markersize=markersize, zorder=1)

    if is_mas:
        curves_mst = read_curves_master(path)
        plot_lc(curves_mst, bnames_fix, lbl='%s Master', fmt='o')

    if is_kur:
        curves_kur = read_curves_kurtenkov(path)
        plot_lc(curves_kur, bnames_fix, lbl='%s Kurtenkov', fmt='.')
        # for lc in curves:
        #     x = lc.Time + jd_shift
        #     y = lc.Mag + mshift
        #     bcolor = colors[lc.Band.Name]
        #     # ax.plot(x, y, label='%s Kurtenkov' % lc.Band.Name,
        #     #         ls="", color=bcolor, markersize=markersize, marker="x")
        #     ax.errorbar(x, y, yerr=lc.MagErr,  label='%s Kurtenkov' % lc.Band.Name,
        #                 color=bcolor, fmt='.')

    if is_wil:            
        curves_wil = read_curves_williams(path)
        plot_lc(curves_wil, bnames_fix, lbl='%s Williams', fmt='d')
        # for lc in curves_wil:
        #     x = lc.Time + jd_shift
        #     y = lc.Mag + mshift
        #     bcolor = colors[lc.Band.Name]
        #     # ax.plot(x, y, label='%s Williams' % lc.Band.Name,
        #     #         ls="", color=bcolor, markersize=markersize, marker="+")
        #     ax.errorbar(x, y, yerr=lc.MagErr, label='%s Williams' % lc.Band.Name,
        #                 color=bcolor,  fmt='d', fillstyle='none')


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


def read_curves_williams(path=sn_path):
    # jd0 = 2457000.
    tbl = read_table_header_float(os.path.join(path, 'williams_lnrm31.txt'))
    # tbl['time'] -= jd0
    curves = table2curves('LRN Williams', tbl)
    return curves
