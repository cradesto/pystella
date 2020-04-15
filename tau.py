#!/usr/bin/env python3

import argparse
import logging

import matplotlib.pyplot as plt
import numpy as np
# from matplotlib import cm
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap

import pystella as ps

mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

__author__ = 'bakl'


def color_map_temp():
    cdict = {'blue': ((0.0, 0.0, 0.0),
                      (0.25, 0.0, 0.0),
                      (0.25, 0.0, 0.0),
                      (0.5, 0.3, 0.3),
                      (0.75, 0.75, 1.0),
                      (1.0, 0.9, 1.0)),

             'green': ((0.0, 0.0, 0.0),
                       (0.25, 0., 0.0),
                       (0.5, 0.9, 0.9),
                       (0.75, 0., 0.0),
                       (1.0, 0.0, 0.0)),

             'red': ((0.0, 0.0, 0.4),
                     (0.15, 0.9, 0.9),
                     (0.25, 0.8, 0.8),
                     (0.5, 0.5, 0.5),
                     (0.75, 0.0, 0.0),
                     (1.0, 0.0, 0.0))
             }
    return LinearSegmentedColormap('UWR', cdict, 256)


def plot_tau(tau, moments=None, fcut=1.e-20, is_info=False):
    # moments = moments or (1., 3., 5, 10., 15., 25., 35., 50., 70., 100., 120., 150., 200.)
    # moments = moments or np.arange(0., 200., 3)
    moments = moments or np.exp(np.linspace(np.log(0.5), np.log(400.), 40))

    # init graph
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    pos = 0
    t_data = []
    for i, t in enumerate(tau.Time):
        if t > moments[pos]:
            t_data.append(t)
            pos += 1

    verts = []
    T_cols = []
    x_lim = [float("inf"), float("-inf")]
    z_lim = [float("inf"), float("-inf")]
    for i, t in enumerate(t_data):
        spec = tau.get_block_by_time(t_data[i])
        spec.cut_flux(fcut * max(spec.Flux))  # cut flux
        ys = spec.Flux
        # ys = np.log10(ys)
        wl = spec.Wl * ps.phys.cm_to_angs
        verts.append(list(zip(wl, ys)))

        x_lim[0] = min(x_lim[0], np.min(wl))
        x_lim[1] = max(x_lim[1], np.max(wl))

        z_lim[0] = min(z_lim[0], np.min(ys))
        z_lim[1] = max(z_lim[1], np.max(ys))

        T_cols.append(spec.T_color)
        if is_info:
            print("time: %f  T_color=%f" % (t, spec.T_color))
            # print "time: %f  T_wien=%f wl_max=%f" % (t_data[i], spec.temp_wien, spec.wl_flux_max)

    Tmap = np.log(T_cols / np.min(T_cols))
    color_map = color_map_temp()
    m = plt.cm.ScalarMappable(cmap=color_map)
    m.set_array(Tmap)
    facecolors = m.to_rgba(Tmap * 1.e-4)
    poly = PolyCollection(verts, facecolors=facecolors, linewidths=1.5)  # np.ones(len(t_data)))
    poly.set_alpha(0.5)
    ax.add_collection3d(poly, zs=t_data, zdir='y')

    # Create a color bar with 11 ticks
    cbar = plt.colorbar(m, shrink=0.85)
    ticks = np.linspace(min(Tmap), max(Tmap), 10)
    ticks_lbl = np.round(np.exp(ticks) * np.min(T_cols), -2)
    cbar.ax.set_yticklabels(ticks_lbl)
    cbar.ax.yaxis.set_label_text("Color temperature [K]")

    ax.set_xlim3d(x_lim)
    ax.set_ylim3d(min(t_data), max(t_data))
    ax.set_zlim3d(z_lim)

    # ax.set_xscale('log')
    ax.set_zscale('log')

    ax.xaxis.set_label_text('Wavelength [A]')
    ax.yaxis.set_label_text('Time [days]')
    ax.zaxis.set_label_text('Tau')

    return fig


def plot_spec_wl(times, series, tt, wl_ab, **kwargs):
    font_size = kwargs.get('font_size', 12)
    nrow = np.math.ceil(len(times)/2.)
    ncol = 2
    fig = plt.figure(figsize=(12, nrow * 4))
    plt.matplotlib.rcParams.update({'font.size': font_size})

    series_cut = series.copy(wl_ab=wl_ab)
    # radiuses = np.interp(times, tt['time'], tt['Rph'], 0, 0)
    radiuses = np.interp(times, tt['time'], tt['rbb'], 0, 0)
    Tbbes = np.interp(times, tt['time'], tt['Tbb'], 0, 0)
    marker_style = dict(linestyle=':', markersize=5)

    for i, t in enumerate(times):
        ax = fig.add_subplot(nrow, ncol, i + 1)
        spec = series.get_spec_by_time(t)
        spec.cut_flux(max(spec.Flux) * 1e-6)  # cut flux
        R = radiuses[i]

        star_bb = ps.Star("bb", spec)
        star_bb.set_distance(R)

        spec_cut = series_cut.get_spec_by_time(t)
        star_cut = ps.Star("bb_cut", spec_cut)
        star_cut.set_distance(R)

        # spectrum
        ax.semilogy(star_bb.Wl * ps.phys.cm_to_angs, star_bb.FluxWlObs, label="Spec Ph")
        # Tcolor
        spec_obs = ps.Spectrum('wbb', star_cut.Freq, star_cut.FluxObs)
        Tcol = spec_obs.T_color
        T_wien = spec_obs.T_wien
        zeta = (T_wien / Tcol) ** 4
        wbb = ps.SpectrumDilutePlanck(spec.Freq, Tcol, zeta)
        ax.semilogy(wbb.Wl * ps.phys.cm_to_angs, wbb.FluxWl, label="Tcol={:.0f} W={:.2f}".format(Tcol, zeta),
                    marker='<', **marker_style)
        # diluted
        Tdil, W = spec_obs.T_color_zeta()
        dil = ps.SpectrumDilutePlanck(spec.Freq, Tdil, W)
        ax.semilogy(dil.Wl * ps.phys.cm_to_angs, dil.FluxWl, label="Tdil={:.0f} W={:.2f}".format(Tdil, W),
                    marker='>', **marker_style)
        # Tbb
        Tbb = Tbbes[i]
        bb = ps.SpectrumPlanck(spec.Freq, Tbb)
        ax.semilogy(bb.Wl * ps.phys.cm_to_angs, bb.FluxWl, label="Tbb={:.0f}".format(Tbb),
                    marker='d', **marker_style)

        # wl range
        if wl_ab is not None:
            for xc in wl_ab:
                plt.axvline(x=xc, color="grey", linestyle='--')
        ax.legend(loc="best", prop={'size': 9})
        ax.text(0.01, 0.05, "$t_d: {:.1f}$".format(t), horizontalalignment='left', transform=ax.transAxes,
                bbox=dict(facecolor='green', alpha=0.3))
        xlim = ax.get_xlim()
        ax.set_xlim(xlim[0], min(xlim[1], 2e4))

    # fig.subplots_adjust(wspace=0, hspace=0)
    # fig.subplots_adjust(wspace=0, hspace=0, left=0.07, right=0.96, top=0.97, bottom=0.06)
    plt.subplots_adjust(left=0.07, right=0.96, top=0.97, bottom=0.06)
    return fig


def get_parser():
    parser = argparse.ArgumentParser(description='Standard Candle Method.')
    print(" Plot the tau-wave diagram for STELLA models")
    parser.add_argument('-b', '--band',
                        required=False,
                        default='V-I',
                        dest="bnames",
                        help="-b <bands>: string, default: V-I, for example B-R-V-I")
    parser.add_argument('-i', '--input',
                        required=True,
                        dest="input",
                        help="Model name, example: cat_R450_M15_Ni007")
    parser.add_argument('-p', '--path',
                        required=False,
                        type=str,
                        default=False,
                        dest="path",
                        help="Model directory")
    parser.add_argument('-t',
                        required=False,
                        type=str,
                        default=None,
                        dest="dt",
                        help="Time interval [day]. Example: 5.:75. Default: all times in the tau-file")
    parser.add_argument('-x',
                        required=False,
                        type=str,
                        default=None,
                        dest="dw",
                        help="wave length interval [A]. Example: 1.:25e3. Default: all waves in the tau-file")

    return parser

#
# def usage():
#     bands = ps.band.band_get_names()
#     print("Usage: show tau(t,nu) from tau-file")
#     print("  tau.py [params]")
#     print("  -b <set_bands>: delimiter '_'. Default: B-V.\n"
#           "     Available: " + '-'.join(sorted(bands)))
#     print("  -i <model name>.  Example: cat_R450_M15_Ni007_E7")
#     # print("  -o  options: [fit, wl] - plot spectral fit")
#     print("  -p <model path(directory)>, default: ./")
#     # print("  -s <file-name> without extension. Save plot to pdf-file. Default: spec_<file-name>.pdf")
#     print("  -t  time interval [day]. Example: 5.:75.")
#     print("  -x  wave length interval [A]. Example: 1.:25e3")
#     # print("  -w  write data to file. Example: flux to magAB, Tcolor and W [.hdf5]")
#     print("  -h  print usage")


def main():
    import os
    import sys

    ps.Band.load_settings()

    model_ext = '.tau'

    parser = get_parser()
    args, unknownargs = parser.parse_known_args()

    path = os.getcwd()
    if args.path:
        path = os.path.expanduser(path)

    # Set model names
    fname = None
    if args.input:
        for arg in args.input:
            fname = args.input.strip()
            fname = fname.replace(model_ext, '')

    if fname is None:
        parser.print_help()
        sys.exit(2)

    model = ps.Stella(fname, path=path)

    if not model.is_tau:
        print("No tau-data for: " + str(model))
        return None

    tau = model.get_tau().load()
    print('Loaded Tau from {} \n    nzone= {} ntimes= {} nfreq= {}'.format(tau.FName, tau.Nzon, tau.Ntimes, tau.NFreq))
    print(tau.Wl2angs)
    # fig = plot_tau(tau)
    # print("Plot  Tau(t,nu): " + str(model))
    # plt.show()


if __name__ == '__main__':
    main()
