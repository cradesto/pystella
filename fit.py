#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import logging
import os
import sys

import numpy as np

import pystella.util.callback as cb
from pystella.fit.fit_mcmc import FitLcMcmc
from pystella.fit.fit_mpfit import FitMPFit
from pystella.rf import band
from pystella.rf import light_curve_func as lcf

# from pystella.rf import light_curve_plot as lcp

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def old_lc_wrapper(param, p=None):
    a = param.split(':')
    fname = a.pop(0)
    if p is None:
        if os.path.isfile(fname + '.py'):
            p, fname = os.path.split(fname)
        elif os.path.isfile(os.path.join(os.getcwd(), fname + '.py')):
            p = os.getcwd()
        else:
            p = cb.plugin_path
    # print("Call: {} from {}".format(fname, p))
    c = cb.CallBack(fname, path=p, args=a, method='load', load=1)
    print("Call: %s from %s" % (c.Func, c.FuncFileFull))
    return c


def rel_errors(mu, sig, func, num=10000):
    x_norm = []
    for x, s in zip(mu, sig):
        x_norm.append(np.random.normal(x, s, num))
    # x_norm = np.random.normal(mu,sig, num)
    f_dist = func(x_norm)
    return np.mean(f_dist), np.std(f_dist)


def m_mu(x):
    return 10 ** (-x * 0.4)


def get_parser():
    parser = argparse.ArgumentParser(description='Process light curve fitting.')
    print("  -b <bands:shift>: string, default: U-B-V-R-I, for example U-B-V-R-I-u-g-i-r-z-UVW1-UVW2.\n"
          "     shift: to move lc along y-axe (minus is '_', for example -b R:2-V-I:_4-B:5 ")

    parser.add_argument('-b', '--band',
                        required=False,
                        dest="bnames",
                        help="-b <bands>: string, default: U-B-V-R-I, for example g-i-r-UVW2")
    parser.add_argument('-c', '--call',
                        required=True,
                        type=str,
                        dest='call',
                        help='Call observational data')
    parser.add_argument('-d',
                        required=False,
                        type=float,
                        default=10,
                        dest="distance",
                        help="Distance to the model [pc].  Default: 10 pc")
    parser.add_argument('-z',
                        required=False,
                        type=float,
                        default=0,
                        dest="redshift",
                        help="Redshift for the model .  Default: 0")
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
    parser.add_argument('-t', '--time',
                        required=False,
                        type=str,
                        default=None,
                        dest="times",
                        help="The range of fitting in model LC. Default: None (all points). Format: {0}".format('2:50'))
    parser.add_argument('-w', '--write',
                        action='store_const',
                        const=True,
                        dest="is_write",
                        help="To write the data to txt-file.")
    return parser


def plot_curves(curves, curves_o):
    from pystella.rf import light_curve_plot as lcp
    # matplotlib.rcParams['backend'] = "TkAgg"
    # matplotlib.rcParams['backend'] = "Qt4Agg"
    from matplotlib import pyplot as plt

    ax = lcp.curves_plot(curves)

    lt = {lc.Band.Name: 'o' for lc in curves_o}
    lcp.curves_plot(curves_o, ax, lt=lt, xlim=(-10, 300))
    plt.show()


def main():
    is_legend = True
    is_debug = True
    t_diff = 1.01

    callback = None
    times = None
    z = 0.
    distance = 10  # pc
    bnames = ['U', 'B', 'V', 'R', "I"]

    band.Band.load_settings()

    parser = get_parser()
    args, unknownargs = parser.parse_known_args()

    name = None
    if len(unknownargs) > 0:
        path, name = os.path.split(unknownargs[0])
        path = os.path.expanduser(path)
        name = name.replace('.ph', '')
    else:
        if args.path:
            path = os.path.expanduser(args.path)
        else:
            path = os.getcwd()
        if args.name is not None:
            name = os.path.splitext(os.path.basename(args.name))[0]

    if name is None:
        parser.print_help()
        sys.exit(2)

    if args.bnames:
        bnames = []
        for bname in str(args.bnames).split('-'):
            if not band.band_is_exist(bname):
                print('No such band: ' + bname)
                parser.print_help()
                sys.exit(2)
            bnames.append(bname)

    if args.distance:
        distance = float(args.distance)

    if args.call:
        callback = cb.lc_wrapper(str(args.call), method='load')
    else:
        print('No obs data. Use key: -c: ')
        parser.print_help()
        sys.exit(2)

    # Get observations
    curves_o = callback.load()

    if len(bnames) == 0:
        bnames = curves_o.BandNames

    # Get models
    if args.times:
        times = list(map(float, args.times.split(':')))

    if times is not None:
        print("Run fitting for model %s %s for %s moments" % (path, name, args.times))
        curves_mdl = lcf.curves_compute(name, path, bnames, z=args.redshift, distance=distance,
                                        t_beg=times[0], t_end=times[1], t_diff=t_diff)
    else:
        print("Run fitting for model %s %s " % (path, name))
        curves_mdl = lcf.curves_compute(name, path, bnames, z=z, distance=distance, t_diff=t_diff)

    # curves_o = curves_o.tmin
    # res = fit_curves_bayesian_2d(curves_o, curves, is_debug=is_debug, is_info=True)

    # fitter = FitLcMcmc()
    # tshift, tsigma = fitter.fit(curves_o.get('V'), curves.get('V'))
    # fitter = FitLcMcmc()
    fitter = FitMPFit(is_debug=True)
    bname = bnames[0]
    res = fitter.fit(curves_o.get(bname), curves_mdl.get(bname))
    # tshift, tsigma = res['dt']
    # dm, dmsigma = res['dm']

    print(""" Results: """)
    print("Band: {0} time shift  = {1:.2f}+/-{2:.4f}".format(bname, res.tshift, res.tsigma))

    curves_o.set_tshift(res.tshift)
    # curves.set_tshift(0.)

    plot_curves(curves_mdl, curves_o)
    # print(" dm_abs      = {0:.4f}+/-{1:.4f}".format(dm, dmsigma))


if __name__ == '__main__':
    main()
