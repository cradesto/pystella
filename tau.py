#!/usr/bin/env python3

import argparse
import logging

import numpy as np
import pystella as ps
from pystella.model.sn_tau import StellaTauDetail

mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

__author__ = 'bakl'

# todo Show filters
# todo show valuse for filters
# todo compute SED = 4 pi R^2 sig T^4


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
        n = int(2 - np.log10(max(1e-03, abs(b.Time))))  # if b.Time >= 10. else 4  # label format
        p = axV.semilogx(b.R, b.V8, label="t= {:.{}f}".format(b.Time, n))
        color = p[0].get_color()
        axT.loglog(b.R, b.T, label="t={:.2f}".format(time), color=color)

    axV.legend(frameon=False)

    if xlim is not None:
        axT.set_xlim(xlim)
        axV.set_xlim(xlim)

    fig.tight_layout()
    return fig


def plot_bands(ax, bnames, amp=30, alpha=0.5):
    """Plot the filter responses"""
    color_dic = ps.band.bands_colors()
    res = {}
    for bname in bnames:
        b = ps.band.band_by_name(bname)
        wl = b.wl * ps.phys.cm_to_angs
        ax.plot(wl, b.resp_wl*amp, color_dic[bname], alpha=alpha)

        wl_eff = b.wl_eff_angs
        ax.axvline(x=wl_eff, ymin=0., ymax=0.99, linestyle='--', color=color_dic[bname], alpha=alpha)
        ax.text(wl_eff, 10, bname, fontsize=12)
        ax.text(wl_eff*.95, 3, "{:.0f}".format(wl_eff), fontsize=6)
        res[bname] = (wl_eff, color_dic[bname])
    return res


def plot_tau_phot(tau_data, pars, tau_ph, xlim=None, title='', bnames=None):
    """
    Plot photosphere as  Func(nu). Maybe: R, V, V8, T
    :param pars: the parameters of photosphere
    :param tau_data: the data at the  optical depth tau_ph
    :param tau_ph:  the photosphere location
    :param xlim:   wave length interval [A]
    :param title: the plot title
    :param bnames: array of filter names to show the filter responses
    :return: figure
    """
    import matplotlib.pyplot as plt

    def fr2wv(nu):
        return ps.phys.c / nu * ps.phys.cm_to_angs

    fig, axs = plt.subplots(len(pars)+1, figsize=(12, 12), sharex=True, gridspec_kw={'hspace': 0})

    # Setup
    ax = axs[0]
    ax.set_ylabel(r'Zone ($\tau_{{ph}}= {:.2f}$)'.format(tau_ph))
    ax.set_title(title)
    ax.xaxis.set_ticks_position('top')
    # ax.xaxis.tick_top()
    # ax.tick_params(axis="x", direction="in", pad=-22)
    # ax.tick_params(direction='in')

    for i, p in enumerate(pars, 1):
        ax = axs[i]
        ax.set_ylabel(r'{}$_{{ph}}$'.format(p))
        if i < len(axs)-1:
            ax.set_xlabel('')
            ax.tick_params(which='both', top=False, bottom=False)
        else:
            ax.set_xlabel('Wavelength [A]')

    # Plot Zone_ph
    colors = []
    for j, (t, freq, y) in enumerate(tau_data[StellaTauDetail.col_zon]):
        axzon = axs[0]
        n = int(3 - np.log10(max(1e-03, abs(t))))  # label format
        lbl = "t= {:.{}f} d".format(t, n)

        ll = axzon.semilogx(fr2wv(freq), y, label=lbl)
        color = ll[0].get_color()
        colors.append(color)

    bnames_waves = None
    if bnames is not None:
        ylim = axzon.get_ylim()
        bnames_waves = plot_bands(axzon, bnames, amp=ylim[1]*0.25, alpha=0.5)

    # Plot other params
    for i, p in enumerate(pars, 1):
        is_log = p.startswith('log')
        p_data = p.replace('log', '') if is_log else p
        ax = axs[i]
        for j, (t, freq, y) in enumerate(tau_data[p_data]):
            x = fr2wv(freq)
            if is_log:
                ax.loglog(x, y, color=colors[j])
            else:
                ax.semilogx(x, y, color=colors[j])

        if bnames_waves is not None:
            for bn,  (wl, col) in bnames_waves.items():
                ax.axvline(x=wl, ymin=0., ymax=0.99, linestyle='--', color=col, alpha=0.5)

    # Post-plotting
    for i, ax in enumerate(axs):
        ax.tick_params(which='both', left=True, right=True, direction="in")
        # ax.grid(axis="x", color="grey", alpha=.5, linewidth=1, linestyle=":")

        if xlim is not None:
            ax.set_xlim(xlim)

    axs[0].legend(frameon=False)

    fig.tight_layout()
    return fig


def get_parser(times='0.1:1:10:25:65', bnames='U:B:V:R'):
    parser = argparse.ArgumentParser(description='Standard Candle Method.')
    print(" Plot the tau-wave diagram for STELLA models")
    parser.add_argument('-b', '--band',
                        nargs='?',
                        required=False,
                        # default=bnames,
                        const=bnames,
                        type=str,
                        dest="bnames",
                        help="-b <bands>: string. If set only -b BNAMES is {}".format(bnames))
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
                        help='Plot photosphere parameter. Maybe: R, V, V8, T. Example: -ph R:V8:T ' 
                             'You may use prefix log, e.g. logT or logV8')
    parser.add_argument('-s', '--save',
                        action='store_const',
                        const=True,
                        dest="is_save",
                        help="To save the result plot to pdf-file. Format: tau_[name]_t[times].pdf.")
    parser.add_argument('-t', '--time',
                        required=False,
                        type=str,
                        default=times,
                        dest="times",
                        help="Plot tau snap for selected time moments. Default: {0}".format(times))
    parser.add_argument('--tau_ph',
                        required=False,
                        type=float,
                        default=2./3.,
                        dest="tau_ph",
                        help="The optical depth at the photosphere. Default: 2/3")
    parser.add_argument('-x', '--xlim',
                        required=False,
                        type=str,
                        default=None,
                        dest="xlim",
                        help="wave length interval [A]. Example: 1.:25e3. Default: all waves in the tau-file")
    parser.add_argument('-w', '--write',
                        required=False,
                        type=str,
                        default=None,
                        dest="write_prefix",
                        help="The prefix of file + -ParamName.dat")

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

    fig = None
    xlim = None
    fplot = None
    print('\n Arguments')
    times = str2float(args.times)
    print(' The time moments: ', args.times)
    print(' The optical depth ', args.tau_ph)
    if args.phot:
        print(' The photospheric parameters ', args.phot)
    if args.xlim is not None:
        xlim = str2float(args.xlim)
        print(" xlim: ", xlim)
    # Set band names
    bnames = ('B',)
    ps.Band.load_settings()
    if args.bnames:
        bnames = []
        for bname in args.bnames.split('-'):
            if not ps.band.band_is_exist(bname):
                print('No such band: ' + bname)
                parser.print_help()
                sys.exit(2)
            bnames.append(bname)

    tau = model.get_tau().load(is_info=False)
    print('\n Loaded data from {}'.format(tau.FName))
    print('Model has Nzone= {} Ntimes= {}'.format(tau.Nzon, tau.Ntimes))
    print("The model time interval: {:.3e} - {:3e} days".format(min(tau.Times), max(tau.Times)))
    print("The bnames are  {}".format(', '.join(bnames)))
    # print(tau.Wl2angs)
    # tau = b.Tau
    # print(tau.shape)

    ###
    # Plot
    if args.phot:
        pars = args.phot.split(':')
        if isinstance(pars, str):
            pars = [pars]
        pars_data = [p.replace('log', '') for p in pars]
        tau_data = tau.params_ph(pars=pars_data, moments=times, tau_ph=args.tau_ph)

        if args.write_prefix:
            fwrite = os.path.expanduser(args.write_prefix)
            tau.data_save(fwrite, tau_data, pars_data)
        else:
            # Print parameters
            print('\nPhotospheric parameters:')
            for ii, p in enumerate(pars_data):
                print('{:9s} {}'.format('t_real', ' '.join([f'{p}_{b:10s}' for b in bnames])))
                for i, (t, freq, y) in enumerate(tau_data[p]):
                    s = '{:9.4f} '.format(t)
                    for bname in bnames:
                        b = ps.band.band_by_name(bname)
                        fr_eff = b.freq_eff
                        idx = (np.abs(freq - fr_eff)).argmin()
                        s += ' {:10e}'.format( y[idx])
                    print(s)
            # Plot
            fig = plot_tau_phot(tau_data, pars,  tau_ph=args.tau_ph, xlim=xlim, title=tau.Name, bnames=bnames)
            fplot = os.path.expanduser("~/tau_{}_{}.pdf".format(fname, str.replace(args.phot, ':', '-')))
    else:
        fig = plot_tau_moments(tau, moments=times, xlim=xlim)

    if args.is_save:
        if fplot is None:
            fplot = os.path.expanduser("~/tau_{0}_t{1}.pdf".format(fname, str.replace(args.times, ':', '-')))
        print("Save plot to {0}".format(fplot))
        fig.savefig(fplot, bbox_inches='tight')
    else:
        plt.show()


if __name__ == '__main__':
    main()
