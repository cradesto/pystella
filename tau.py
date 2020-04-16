#!/usr/bin/env python3

import argparse
import logging

import numpy as np
import pystella as ps

mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

__author__ = 'bakl'


def plot_tau_phot_moments(tau, par, moments, tau_ph=2./3., xlim=None):
    """
    Plot photosphere as  Func(nu). Maybe: R, V, V8, T
    :param tau:
    :param par: photosphere parameter. Maybe: R, V, V8, T
    :param moments: time moments
    :param tau_ph:  the photosphere location
    :param xlim:   wave length interval [A]
    :param is_info:
    :return: figure
    """
    import matplotlib.pyplot as plt

    # moments = moments or np.exp(np.linspace(np.log(tlim[0]), np.log(tlim[1]), 40))

    if isinstance(par, str):
        par = [par]

    fig, axs = plt.subplots(len(par)+1, figsize=(12, 12), sharex=True, gridspec_kw={'hspace': 0})

    # Setup
    ax = axs[0]
    ax.set_ylabel('Tau_ph')
    ax.set_title(tau.Name)
    # ax.xaxis.set_ticks_position('top')
    ax.xaxis.tick_top()

    for i, p in enumerate(par, 1):
        ax = axs[i]
        ax.set_ylabel(p+'_ph')
        if i < len(axs)-1:
            ax.set_xlabel('')
        else:
            ax.set_xlabel('Wavelength [A]')

    # Plot data
    for j, time in enumerate(moments):
        b = tau.block_nearest(time)
        wl = b.Wl2angs
        n = 2 if b.Time >= 10. else 4  # label format
        idxs = b.ph_indexes(tau_ph=tau_ph)
        if len(wl) != len(idxs):
            raise ValueError("Error in photo. indexes: len(wl)= {}  len(idxs)= {}".format(len(wl), len(idxs)))

        # Plot Tau
        ax = axs[0]
        lbl = "t= {:.{}f}".format(b.Time, n)
        ll = ax.semilogx(wl, idxs, label=lbl)
        color = ll[0].get_color()

        for i, p in enumerate(par, 1):
            ax = axs[i]
            is_log = p.startswith('log')
            p = p.replace('log', '')
            arr = getattr(b, p)
            x = wl
            y = [arr[idx] for idx in idxs]

            if is_log:
                ll = ax.loglog(x, y, color=color)
            else:
                ll = ax.semilogx(x, y, color=color)
            # color = ll[0].get_color()
            # axT.loglog(b.R, b.T, label="t={:.2f}".format(time), color=color)

    axs[0].legend(frameon=False)
    # Post
    for i, ax in enumerate(axs):
        if xlim is not None:
            ax.set_xlim(xlim)

    fig.tight_layout()
    return fig


def plot_tau_moments(tau, moments=None, xlim=None):
    import matplotlib.pyplot as plt

    moments = moments or np.exp(np.linspace(np.log(0.5), np.log(400.), 40))

    fig, (axV, axT) = plt.subplots(2, figsize=(12, 12), sharex=True, gridspec_kw={'hspace': 0})
    axV.set_title(tau.Name)
    axV.set_xlabel('')
    axV.set_ylabel('Velocity [1000 km/s]')

    axT.set_xlabel('Radius [cm]')
    axT.set_ylabel('Temperature [K]')

    for i, time in enumerate(moments):
        b = tau.block_nearest(time)
        n = 2 if b.Time >= 10. else 4  # label format
        p = axV.semilogx(b.R, b.V8, label="t= {:.{}f}".format(b.Time, n))
        color = p[0].get_color()
        axT.loglog(b.R, b.T, label="t={:.2f}".format(time), color=color)

    axV.legend(frameon=False)
    # axT.legend(frameon=False)

    fig.tight_layout()
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
    parser.add_argument('-ph', '--phot',
                        required=False,
                        type=str,
                        default=False,
                        dest="phot",
                        help='Plot photosphere parameter. Maybe: R, V, V8, T. '
                             'Yoy may use prefix log, e.g. logT or logV8')
    parser.add_argument('-s', '--save',
                        action='store_const',
                        const=True,
                        dest="is_save",
                        help="To save the result plot to pdf-file. Format: swd_[name]_t[times].pdf.")
    d = '0.1:0.5:1:4:15:65:99'
    parser.add_argument('-t', '--time',
                        required=False,
                        type=str,
                        default=d,
                        dest="times",
                        help="Plot tau snap for selected time moments. Default: {0}".format(d))
    parser.add_argument('-x', '--xlim',
                        required=False,
                        type=str,
                        default=None,
                        dest="xlim",
                        help="wave length interval [A]. Example: 1.:25e3. Default: all waves in the tau-file")

    return parser


def str2float(s):
    return list(map(float, s.split(':')))


def main():
    import os
    import sys
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        plt = None

    ps.Band.load_settings()

    model_ext = '.tau'

    parser = get_parser()
    args, unknownargs = parser.parse_known_args()

    path = os.getcwd()
    if args.path:
        path = os.path.expanduser(args.path)

    # Set model names
    fname = None
    if args.input:
        fname = args.input.strip()
        fname = fname.replace(model_ext, '')

    if fname is None:
        parser.print_help()
        sys.exit(2)

    model = ps.Stella(fname, path=path)

    if not model.is_tau:
        print("No tau-data for: " + str(model))
        return None

    xlim = None
    fsave = None
    times = str2float(args.times)
    if args.xlim is not None:
        xlim = str2float(args.xlim)
        print("   xlim: ", xlim)

    tau = model.get_tau().load(is_info=True)
    print("Times: {:.3e} - {:3e} days".format(min(tau.Times), max(tau.Times)))
    # print(tau.Wl2angs)
    # tau = b.Tau
    # print(tau.shape)

    ###
    # Plot
    if args.phot:
        fig = plot_tau_phot_moments(tau, par=args.phot.split(':'), moments=times, xlim=xlim)
        fsave = os.path.expanduser("~/tau_{}_{}.pdf".format(fname, args.phot))
    else:
        fig = plot_tau_moments(tau, moments=times, xlim=xlim)

    if args.is_save:
        if fsave is None:
            fsave = os.path.expanduser("~/tau_{0}.pdf".format(fname))
        print("Save plot to {0}".format(fsave))
        fig.savefig(fsave, bbox_inches='tight')
    else:
        plt.show()


if __name__ == '__main__':
    main()
