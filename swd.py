#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import csv
import logging
import os
import sys
from os.path import dirname

import numpy as np

from pystella.model.stella import Stella
from pystella.rf import light_curve_plot as lcp

# matplotlib.rcParams['backend'] = "TkAgg"
# matplotlib.rcParams['backend'] = "Qt4Agg"
import matplotlib.pyplot as plt

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))
logging.basicConfig()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def uph_compute_tau(swd):
    """Compute tau for each moment of time
    Return: 2d-array[i,k], where i - time index, k - zone index.
    """
    taus = np.zeros((swd.Ntimes, swd.Nzon))
    for i, time in enumerate(swd.Times):
        s = swd[i]
        tau = 0.
        for k in range(swd.Nzon-1, 1, -1):
            tau += s.Cappa[k] * s.Rho[k] * (s.R[k] - s.R[k-1])
            taus[i, k] = tau
    return taus


def uph_find_uph(swd, taus, tau_ph=2. / 3.):
    col = ('time', 'R', 'V', 'zone')
    res = {k: np.zeros(swd.Ntimes) for k in col}
    # res = {'time': [], 'R': [], 'V': []}
    for i, time in enumerate(swd.Times):
        s = swd[i]
        kph = 1  # todo check the default value
        for k in range(swd.Nzon-1, 1, -1):
            if taus[i][k] >= tau_ph:
                kph = k
                break
        res['time'][i] = time
        res['zone'][i] = kph
        res['R'][i] = s.R[kph]
        res['V'][i] = s.Vel[kph] / 1e8  # to 1000 km/s

    return res


def uph_save(dictionary, fname, sep='\t'):
    """
    Save dict to file. Keys are column's names, values are column data
    :param dictionary:
    :type fname: the name of the file
    :param sep: string separator
    """
    col = ('time', 'zone', 'R', 'V')
    with open(fname, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=sep, quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['# {:^8s}'.format(x) for x in col])
        for i, (row) in enumerate(zip(*[dictionary[c] for c in col])):
            writer.writerow(['{:8.3e}'.format(x) for x in row])


def plot_uph(uph, label='', lw=2):
    # setup plot
    plt.matplotlib.rcParams.update({'font.size': 12})
    fig = plt.figure(figsize=(8, 8))

    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('Time [day]')
    ax.set_ylabel(r'Photospheric Velocity [$\times 10^3$ km/s]')
    #    ax.set_title(fname_uph)

    x = uph['time']
    y = uph['V']
    ax.plot(x, y, label=label, color='blue', ls="-", linewidth=lw)
    ax.legend()
    plt.show()
    return fig


def main():
    parser = argparse.ArgumentParser(description='Process Stella Shock Wave Details.')

    parser.add_argument('-i', '--input',
                        required=False,
                        dest="name",
                        help="Model name, example: cat_R450_M15_Ni007")

    parser.add_argument('-p', '--path',
                        required=False,
                        type=str,
                        default='./',
                        dest="path",
                        help="Model directory")

    parser.add_argument('--rnorm',
                        required=False,
                        default='lgr',  # 1e14,
                        dest="rnorm",
                        help="Radius normalization, example: 'm' or 'sun' or 1e13")

    parser.add_argument('--vnorm',
                        required=False,
                        type=float,
                        default=1e8,
                        dest="vnorm",
                        help="Velocity normalization, example: 1e8")

    parser.add_argument('--lumnorm',
                        required=False,
                        type=float,
                        default=1e40,
                        dest="lumnorm",
                        help="Luminously normalization, example: 1e40")

    d = '1:4:15:65'
    parser.add_argument('-t', '--time',
                        required=False,
                        type=str,
                        default=d,
                        dest="times",
                        help="Plot shock wave snap for selected time moments. Default: {0}".format('2:10:50'))

    parser.add_argument('-c', action='store_const', dest='constant_value',
                        const='value-to-store',
                        help='Store a constant value')

    parser.add_argument('--uph', action='store_const', dest='is_uph',
                        const=True,
                        help='To compute the photospheric velocity')

    parser.add_argument('-s', '--save',
                        action='store_const',
                        const=True,
                        dest="is_save",
                        help="To save plot to pdf-file. Default: False".format('2:10:50'))

    parser.add_argument('-w', '--write',
                        action='store_const',
                        const=True,
                        dest="is_write",
                        help="To write the data to txt-file.")

    args, unknownargs = parser.parse_known_args()

    name = None
    if len(unknownargs) > 0:
        path, name = os.path.split(unknownargs[0])
        path = os.path.expanduser(path)
        name = name.replace('.swd', '')
    else:
        if args.path:
            path = args.path
        else:
            path = ROOT_DIRECTORY
        if args.name is not None:
            name = os.path.splitext(os.path.basename(args.name))[0]

    if name is None:
        # logger.error(" No data. Use key '-i' ")
        parser.print_help()
        sys.exit(2)

    times = list(map(float, args.times.split(':')))
    print("Run swd-model %s %s for %s moments" % (path, name, args.times))
    stella = Stella(name, path=path)
    swd = stella.get_swd().load()

    if args.is_uph:
        logger.info(' Compute and print uph')
        taus = uph_compute_tau(swd)
        duph = uph_find_uph(swd, taus)
        # print uph
        print(duph.keys())
        for row in zip(*duph.values()):
            print(['{:12f}'.format(x) for x in row])
        # save uph
        if args.is_write:
            fsave = os.path.join(path, "{0}.uph".format(name))
            print("Save uph to {0}".format(fsave))
            uph_save(duph, fsave)
        else:
            fig = plot_uph(duph, label=name)
            if args.is_save:
                fsave = os.path.expanduser("~/uph_{0}.pdf".format(name))
                print("Save plot to {0}".format(fsave))
                fig.savefig(fsave, bbox_inches='tight')

    else:
        fig = lcp.plot_shock_details(swd, times=times, vnorm=args.vnorm, rnorm=args.rnorm,
                                     lumnorm=args.lumnorm, is_legend=True)

        plt.show()
        if args.is_save:
            fsave = os.path.expanduser("~/swd_{0}_t{1}.pdf".format(name, str.replace(args.times, ':', '-')))
            print("Save plot to {0}".format(fsave))
            fig.savefig(fsave, bbox_inches='tight')


if __name__ == '__main__':
    main()
