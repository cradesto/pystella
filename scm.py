#!/usr/bin/python3
# -*- coding: utf-8 -*-
import getopt
import numpy as np
import os
import sys
from os.path import isfile, join, dirname
from scipy import interpolate
from scipy.optimize import fmin

import matplotlib.pyplot as plt
from matplotlib import gridspec

from pystella import velocity
from pystella.rf import band
from pystella.rf import light_curve_func as lcf
from pystella.util.phys_var import phys

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))

_colors = ["blue", "cyan", "brown", 'darkseagreen', 'tomato', 'olive', 'orange',
           'skyblue', 'darkviolet']
# colors = {"B-V": "blue", 'B-V-I': "cyan", 'V-I': "brown"}
lntypes = {"B-V": "-", 'B-V-I': "-.", 'V-I': "--",
           "u-g": "-", 'g-r-i': "-.", 'r-i': "--"}
markers = {u'D': u'diamond', 6: u'caretup', u's': u'square', u'x': u'x',
           5: u'caretright', u'^': u'triangle_up', u'd': u'thin_diamond', u'h': u'hexagon1',
           u'+': u'plus', u'*': u'star', u'o': u'circle', u'p': u'pentagon', u'3': u'tri_left',
           u'H': u'hexagon2', u'v': u'triangle_down', u'8': u'octagon', u'<': u'triangle_left'}
markers = markers.keys()


def plot_scm(models_data, bands, z, xlim=None, is_fit=False):
    ylim = (1.5, 4.)

    # setup figure
    plt.matplotlib.rcParams.update({'font.size': 14})
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    fig = plt.figure(num=len(bands), figsize=(9, 9), dpi=100, facecolor='w', edgecolor='k')
    gs1 = gridspec.GridSpec(len(bands) // 2 + len(bands) % 2, 2)
    # gs1 = gridspec.GridSpec(2, 4, width_ratios=(8, 1, 8, 1))
    gs1.update(wspace=0., hspace=0., left=0.1, right=0.9)

    ax_cache = {}

    # create the grid of figures
    ib = 0
    for bname in bands:
        ib += 1
        icol = (ib - 1) % 2
        irow = (ib - 1) / 2
        ax = fig.add_subplot(gs1[irow, icol])
        ax_cache[bname] = ax
        # ax.legend(prop={'size': 6})
        # x
        if icol > 0:
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position("right")
        ax.set_ylabel(r'$V_{ph}$ [1000 km/s]')

        if irow == 0:
            ax.set_xlabel(r'$M_{%s}$' % bname)

        if xlim is None:
            ax.set_xlim(max(models_data[bname]) + 2., min(models_data[bname]) - 2.)
        else:
            ax.set_xlim(xlim)
        # xstart, xend = 0, 20000.
        # ax.xaxis.set_ticks(np.arange(5000, xend, (xend - xstart) / 4.))
        # ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
        # y
        ax.set_ylim(ylim)
        # ax.set_yticks(np.arange(0.5, ylim[1], 0.5))
        # ax.text(.5, .9, bset, horizontalalignment='center', transform=ax.transAxes)
        # ax.set_title(bset)
        ax.set_title(bname, x=0.5, y=0.9)

    # plot data
    mi = 0
    vel = models_data['v'] / 1e8  # convert to 1000 km/c
    for bname in bands:
        ax = ax_cache[bname]
        x = models_data[bname]
        bcolor = _colors[ib % (len(_colors) - 1)]
        ax.plot(x, vel, marker=markers[mi % (len(markers) - 1)], label='Models',
                markersize=5, color=bcolor, ls="", linewidth=1.5)
        # выбросы todo сделать относительно фита

    if is_fit:
        yh = np.linspace(ylim[0], ylim[1], num=50)
        # hamuy
        for bname in bands:
            ax = ax_cache[bname]
            xx = scm_fit(yh, Av=0, bname=bname, z=z, src='hamuy')
            if yh is not None:
                ax.plot(xx, yh, color="darkviolet", ls="--", linewidth=2.5, label='Hamuy')
        # kasen
        for bname in bands:
            ax = ax_cache[bname]
            xx = scm_fit(yh, Av=0, bname=bname, z=z, src='hamuy')
            if yh is not None:
                ax.plot(xx, yh, color="red", ls="--", linewidth=2.5, label='Kasen')
        # nugent
        yn = np.linspace(ylim[0], ylim[1], num=len(models_data['V']))
        for bname in bands:
            ax = ax_cache[bname]
            V_I = models_data['V'] - models_data['I']  # todo make for any bands
            xn = scm_fit(yn, Av=V_I, src='nugent')
            if yn is not None:
                ax.plot(xn, yn, label='Nugent', color="blue", ls="--")
                # ax.plot(xn, yn, marker='x', label='Nugent', markersize=5, color="blue", ls="")
                # ax.plot(xx, yy, color="orange", ls="--", linewidth=2.5, label='Nugent')
        # bakl
        for bname in bands:
            ax = ax_cache[bname]
            mag = models_data[bname]
            Av = 0.
            a = scm_fit_bakl(mag, vel, z=z, Av=Av)
            mag_bakl = hamuy_fit(a, vel, z, Av)
            if yn is not None:
                ax.plot(mag_bakl, vel, label='Baklanov', color="orange", ls="--", lw=2)
                # ax.plot(mag_bakl, vel, marker='o', label='Baklanov', markersize=5, color="blue", ls="")
            print("Bakl fit for %s: %4.2f, %4.2f" % (bname, a[0], a[1]))
            # ax.plot(xx, yy, color="orange", ls="--", linewidth=2.5, label='Nugent')

    # legend
    for bname in bands:
        ax_cache[bname].legend(prop={'size': 6})

    # plt.title('; '.join(set_bands) + ' filter response')
    # plt.grid()
    plt.show()
    return fig


def scm_fit(v, Av=0, bname=None, z=0.003, src='hamuy'):
    """
    Magnitude V from Hamuy & Pinto 2002, Nugent 2006
    :param bname:  band name
    :param v:  expansion velocity [1000 km/c]
    :param Av:  host and Galaxy extinction; Av  => V-I for Nugent method
    :param z:  redshift in the cosmic microwave background frame
    :param src: fit source
    :return: magnitude
    """
    coef = {'hamuy': {'V': [6.504, 1.294], 'I': [5.820, 1.797]},
            'nugent': {'alf': 6.69, 'V_I_0': 0.53, 'RI': 1.36, 'MI0': -17.49}
            }
    if src == 'hamuy':
        a = coef[src][bname]
        mag = hamuy_fit(a, v, z, Av)
        return mag
    if src == 'nugent':
        a = coef[src]
        mag = -1. * a['alf'] * np.log10(v / 5.) - a['RI'] * (Av - a['V_I_0']) + a['MI0']
        return mag
    if src == 'kasen':
        t = 50.  # day
        mag = -17.4 - 6.9 * np.log10(v / 5.) + 3.1 * np.log10(t / 50)
        return mag
    return None


def hamuy_fit(a, v, z, Av):
    # v /= 5.  # convert to units of 5000 km/c as Hamuy
    res = Av - a[0] * np.log10(v / 5.) + 5 * np.log10(phys.c / 1e5 * z) - a[1]
    return res


def scm_fit_bakl(mag, vel, z, Av=0.):
    a_init = [5, 0.]

    def cost(a, m, v):
        m_fit = hamuy_fit(a, v, z, Av)
        e = np.sum((m - m_fit) ** 2) / len(m)
        return e

    res = fmin(cost, x0=a_init, args=(mag, vel), disp=0)
    return res


def usage():
    print("Usage:")
    print("  %s [params]" % __file__)
    print("  -b <set_bands>: delimiter '_'. Default: B-V-I_B-V_V-I")
    print("  -i <model name>.  Example: cat_R450_M15_Ni007_E7")
    print("  -p <model path(directory)>, default: ./")
    print("  -f  force mode: rewrite zeta-files even if it exists")
    print("  -o  options: <fit:fitb:time:ubv:Tnu> - fit E&D: fit bakl:  show time points: plot UBV")
    print("  -s  save plot to file, default 'False'")
    print("  -h  print usage")
    print("   --- ")
    band.print_bands()


def extract_time(t, times, mags):
    if t < times[0] or t > times[-1]:
        raise ValueError("Time (%f) should be in range [%f, %f]" % (t, times[0], times[-1]))

    if len(times) > 4:
        y_spline = interpolate.splrep(times, mags, s=0)
        res = interpolate.splev(t, y_spline)
    else:
        res = np.interp(t, times, mags, 0, 0)
    return res


def run_scm(bands, distance, names, path, t50, t_beg, t_end, z):
    res = np.array(np.zeros(len(names)),
                   dtype=np.dtype({'names': ['v'] + bands,
                                   'formats': [np.float] * (1 + len(bands))}))
    im = 0
    for name in names:
        vels = velocity.compute_vel(name, path, z=z, t_beg=t_beg, t_end=t_end)
        if vels is None:
            print("No enough data for %s " % name)
            continue

        curves = lcf.curves_compute(name, path, bands, z=z, distance=distance,
                                                 t_beg=t_beg, t_end=t_end)
        v = extract_time(t50, vels['time'], vels['vel'])
        res[im]['v'] = v
        # dic_results[name]['v'] = v
        for bname in bands:
            m = extract_time(t50, curves.get(bname).Time, curves.get(bname).Mag)
            res[im][bname] = m
        print("Run: %s [%d/%d]" % (name, im + 1, len(names)))
        im += 1
    res = res[res[:]['v'] > 0]
    return res


def main(name='', path='./'):
    model_ext = '.tt'
    ubv_args = ''
    is_save = False

    band.Band.load_settings()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "fhstp:i:b:o:")
    except getopt.GetoptError as err:
        print(str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if not name:
        if len(opts) == 0:
            usage()
            sys.exit(2)
        for opt, arg in opts:
            if opt == '-i':
                path = ROOT_DIRECTORY
                name = str(arg)
                break

    bands = ['V', 'I']
    # set_bands = ['B-V', 'B-V-I', 'V-I', 'J-H-K']
    # set_bands = ['U-B-V-I', 'U-B-V-R-I', 'U-B', 'V-R', 'B-V-I', 'B-V', 'V-I']

    for opt, arg in opts:
        if opt == '-e':
            model_ext = '.' + arg
            continue
        if opt == '-b':
            bands = str(arg).split('-')
            for b in bands:
                if not band.band_is_exist(b):
                    print('No such band: ' + b)
                    sys.exit(2)
            continue
        if opt == '-o':
            ubv_args += " %s %s ".format(opt, arg)
            continue
        if opt == '-s':
            is_save = True
            ubv_args += opt + ' '
            continue
        if opt == '-f':
            is_save = True
            continue
        if opt == '-p':
            path = os.path.expanduser(str(arg))
            if not (os.path.isdir(path) and os.path.exists(path)):
                print("No such directory: " + path)
                sys.exit(2)
            continue
        elif opt == '-h':
            usage()
            sys.exit(2)

    names = []
    if name != '':
        names.append(name)
    else:  # run for all files in the path
        files = [f for f in os.listdir(path) if isfile(join(path, f)) and f.endswith(model_ext)]
        for f in files:
            names.append(os.path.splitext(f)[0])

    distance = 10.  # pc for Absolute magnitude
    # distance = 10e6  # pc for Absolute magnitude
    z = phys.H0 * (distance / 1e6) / (phys.c / 1e5)  # convert D to Mpc, c to km/c
    t50 = 50.
    t_beg = max(0., t50 - 10.)
    t_end = t50 + 10.

    if len(names) > 0:
        res = run_scm(bands, distance, names, path, t50, t_beg, t_end, z)
        if len(res) > 0:
            fig = plot_scm(res, bands, z, is_fit=True)
            if is_save:
                fsave = os.path.join(os.path.expanduser('~/'), 'scm_' + '_'.join(bands) + '.pdf')
                print("Save plot in %s" % fsave)
                fig.savefig(fsave, bbox_inches='tight', format='pdf')
    else:
        print("There are no models in the directory: %s with extension: %s " % (path, model_ext))


if __name__ == '__main__':
    main()
