#!/usr/bin/env python3
# #!/usr/bin/python3

import argparse
import os
import sys
from os.path import dirname

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from scipy import interpolate
from scipy.optimize import fmin

from pystella import velocity
from pystella.rf import band
from pystella.rf import light_curve_func as lcf
from pystella.util.path_misc import get_model_names
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
markers = list(markers.keys())


def plot_scm(models_data, bands, z, is_fit=True, xlim=None, ylim=None):
    # setup figure
    plt.matplotlib.rcParams.update({'font.size': 14})
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    fig = plt.figure(num=len(bands), figsize=(9, 9), dpi=100, facecolor='w', edgecolor='k')
    gs1 = gridspec.GridSpec(len(bands) // 2 + len(bands) % 2, 2)
    gs1.update(wspace=0., hspace=0., left=0.1, right=0.9)

    ax_cache = {}

    # create the grid of figures
    for ib, bname in enumerate(bands):
        icol = ib % 2
        irow = ib // 2
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
        if ylim is not None:
            ax.set_ylim(ylim)
        # ax.set_yticks(np.arange(0.5, ylim[1], 0.5))
        # ax.text(.5, .9, bset, horizontalalignment='center', transform=ax.transAxes)
        # ax.set_title(bset)
        ax.set_title(bname, x=0.5, y=0.9)

    # plot data
    mi = 0
    vel = models_data['v'] / 1e8  # convert to 1000 km/c
    for ib, bname in enumerate(bands):
        ax = ax_cache[bname]
        x = models_data[bname]
        bcolor = _colors[ib % (len(_colors) - 1)]
        ax.plot(x, vel, marker=markers[mi % (len(markers) - 1)], label='Models',
                markersize=5, color=bcolor, ls="", linewidth=1.5)
        # выбросы todo сделать относительно фита
        # bakl
        if is_fit:
            ax = ax_cache[bname]
            mag = models_data[bname]
            Av = 0.
            a = scm_fit_bakl(mag, vel, z=z, Av=Av)
            mag_bakl = hamuy_fit(a, vel, z, Av)
            # sort
            sorti = np.argsort(mag_bakl)
            xb = mag_bakl[sorti]
            yb = vel[sorti]
            # if mag_bakl is not None:
            ax.plot(xb, yb, label='Baklanov', color="orange", ls="-", lw=2)
            # ax.plot(mag_bakl, vel, marker='o', label='Baklanov', markersize=5, color="blue", ls="")
            print("Bakl fit for %s: %4.2f, %4.2f" % (bname, a[0], a[1]))

    # plt.title('; '.join(set_bands) + ' filter response')
    # plt.grid()
    return fig, ax_cache


def plot_scm_fit(ax_cache, bands, models_data, z):
    # hamuy
    for bname in bands:
        ax = ax_cache[bname]
        ylim = ax.get_ylim()
        yh = np.linspace(ylim[0], ylim[1], num=50)
        xx = scm_fit(yh, Av=0, bname=bname, z=z, src='hamuy')
        if yh is not None:
            ax.plot(xx, yh, color="darkviolet", ls="--", linewidth=2.5, label='Hamuy')
    # kasen
    for bname in bands:
        ax = ax_cache[bname]
        yh = np.linspace(ylim[0], ylim[1], num=50)
        xx = scm_fit(yh, Av=0, bname=bname, z=z, src='kasen')
        if yh is not None:
            ax.plot(xx, yh, color="red", ls="-.", linewidth=2.5, label='Kasen')
    # nugent
    bname = 'I'
    if bname in bands:
        ax = ax_cache[bname]
        ylim = ax.get_ylim()
        yn = np.linspace(ylim[0], ylim[1], num=len(models_data['V']))
        V_I = models_data['V'] - models_data['I']  # todo make for any bands
        xn = scm_fit(yn, Av=V_I, src='nugent')
        xnZhl = scm_fit(yn, Av=V_I, src='nugentZhl')
        if yn is not None:
            ax.plot(xn, yn, label='Nugent', color="blue", ls="--")
            ax.plot(xnZhl, yn, label='Nugent, whole sample', color="blue", ls="-.")
            # ax.plot(xn, yn, marker='x', label='Nugent', markersize=5, color="blue", ls="")
            # ax.plot(xx, yy, color="orange", ls="--", linewidth=2.5, label='Nugent')


def scm_fit(v, Av=0, bname=None, z=None, src='hamuy'):
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
            'nugent': {'alf': 6.69, 'V_I_0': 0.53, 'RI': 1.36, 'MI0': -17.49},
            'nugentZhl': {'alf': 5.81, 'V_I_0': 0.53, 'RI': 1.36, 'MI0': -17.52},
            # 'nugent': {'alf': 6.69, 'V_I_0': 0.53, 'RI': 0., 'MI0': -17.49}
            }
    if src == 'hamuy':
        a = coef[src][bname]
        mag = hamuy_fit(a, v, z, Av)
        return mag
    if src in ['nugent', 'nugentZhl']:
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


def extract_time(t, times, mags):
    if t < times[0] or t > times[-1]:
        raise ValueError("Time (%f) should be in range [%f, %f]" % (t, times[0], times[-1]))

    if len(times) > 4:
        y_spline = interpolate.splrep(times, mags, s=0)
        res = interpolate.splev(t, y_spline)
    else:
        res = np.interp(t, times, mags, 0, 0)
    return res


def run_scm(bands, distance, names, path, t50, t_beg, t_end, z, method='tt'):
    res = np.array(np.zeros(len(names)),
                   dtype=np.dtype({'names': ['v'] + bands,
                                   'formats': [np.float] * (1 + len(bands))}))
    for im, name in enumerate(names):
        vels = None
        if method=='tt':
            try:
                vels = velocity.compute_vel_res_tt(name, path, z=z, t_beg=t_beg, t_end=t_end)
            except ValueError as ex:
                print(ex)
        elif method == 'swd':
            try:
                vels = velocity.compute_vel_swd(name, path, z=z)
            except ValueError as ex:
                print(ex)
        else:
            raise ValueError('The method [{}] for velocity computation is not supported.'.format(method))

        if vels is None:
            print("No enough data for %s " % name)
            continue

        v = extract_time(t50, vels['time'], vels['vel'])
        res[im]['v'] = v
        # dic_results[name]['v'] = v
        curves = lcf.curves_compute(name, path, bands, z=z, distance=distance,
                                    t_beg=t_beg, t_end=t_end)
        for bname in bands:
            m = extract_time(t50, curves.get(bname).Time, curves.get(bname).Mag)
            res[im][bname] = m
        print("Run: %s [%d/%d]" % (name, im + 1, len(names)))
    res = res[res[:]['v'] > 0]
    return res


def get_parser():
    parser = argparse.ArgumentParser(description='Standard Candle Method.')
    print(" Compute and plot the velocity-magnitude diagram for STELLA models")
    parser.add_argument('-b', '--band',
                        required=False,
                        default='V-I',
                        dest="bnames",
                        help="-b <bands>: string, default: V-I, for example B-R-V-I")
    parser.add_argument('-i', '--input',
                        required=False,
                        dest="input",
                        help="Model name, example: cat_R450_M15_Ni007")
    parser.add_argument('-p', '--path',
                        required=False,
                        type=str,
                        default='./',
                        dest="path",
                        help="Model directory")
    parser.add_argument('-v', '--vel',
                        required=False,
                        type=str,
                        default='tt',
                        dest="method_vel",
                        help="tt OR swd - the way to get photospheric velocity, default: tt")
    parser.add_argument('-s', '--save',
                        required=False,
                        type=str,
                        default=None,
                        dest="save_file",
                        help="To save the result plot to pdf-file.")
    parser.add_argument('-t',
                        required=False,
                        type=float,
                        default=50,
                        dest="tplateau",
                        help="Time moment at plateau stage [day].  Default: 50 ")

    return parser


def main(name=None):
    model_ext = '.ph'

    band.Band.load_settings()

    parser = get_parser()
    args, unknownargs = parser.parse_known_args()

    # Set model names
    names = []

    if len(unknownargs) > 0:
        path, name = os.path.split(unknownargs[0])
        path = os.path.expanduser(path)
        name = name.replace(model_ext, '')
    else:
        if args.path:
            path = os.path.expanduser(args.path)
        else:
            path = os.getcwd()
        if args.input is not None:
            name = os.path.splitext(os.path.basename(args.input))[0]  # remove extension

    # if len(unknownargs) == 0:
    #     parser.print_help()
    #     sys.exit(2)

    if name is None:
        names = get_model_names(path, model_ext)  # run for all files in the path
    else:
        names.append(name)

    # Set band names
    bnames = args.bnames.split('-')
    # check bnames
    for bname in bnames:
        if not band.is_exist(bname):
            print('No such band: ' + bname)
            parser.print_help()
            sys.exit(2)

    distance = 10.  # pc for Absolute magnitude
    # distance = 10e6  # pc for Absolute magnitude
    z = phys.H0 * (distance / 1e6) / (phys.c / 1e5)  # convert D to Mpc, c to km/c
    # z = 0.
    t50 = args.tplateau
    t_beg = max(0., t50 - 10.)
    t_end = t50 + 10.

    if len(names) > 0:
        print("Parameters: vel-method {}, time of plateau = {:f}".format(args.method_vel, t50))
        res = run_scm(bnames, distance, names, path, t50, t_beg, t_end, z, method=args.method_vel)
        if len(res) > 0:
            fig, ax_cache = plot_scm(res, bnames, z)
            # fig, ax_cache = plot_scm(res, bands, z, ylim=(1.5, 4.))

            if all(x in bnames for x in ('V', 'I')):  # fit works only for 'V', 'I'
                # if len(bands) == 2 and all(x in bands for x in ('V', 'I')):  # fit works only for 'V', 'I'
                plot_scm_fit(ax_cache, ('V', 'I'), res, z)

            # legend
            for bname in bnames:
                ax_cache[bname].legend(prop={'size': 6})

            plt.legend(prop={'size': 8})
            plt.show()

            if args.save_file:
                if args.save_file == '1':
                    fsave = os.path.join(os.path.expanduser('~/'), 'scm_' + '_'.join(bnames) + '.pdf')
                else:
                    fsave = args.save_file
                print("Save plot in %s" % fsave)
                fig.savefig(fsave, bbox_inches='tight', format='pdf')

    else:
        print("There are no models in the directory: %s with extension: %s " % (path, model_ext))


if __name__ == '__main__':
    main()
