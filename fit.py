#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import logging
import os
import sys
from collections import OrderedDict

import numpy as np

import pystella.util.callback as cb
from pystella.fit.fit_mcmc import FitLcMcmc
from pystella.fit.fit_mpfit import FitMPFit
from pystella.rf import band
from pystella.rf import light_curve_func as lcf
from pystella.rf.band import band_is_exist
from pystella.util.path_misc import get_model_names
from pystella.util.string_misc import str2interval

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
    print(" Observational data could be loaded with plugin, ex: -c lcobs:filedata:tshift:mshift")

    parser.add_argument('-b', '--band',
                        required=False,
                        dest="bnames",
                        help="-b <bands>: string, default: U-B-V-R-I, for example g-i-r-UVW2")
    parser.add_argument('-c', '--call',
                        required=True,
                        type=str,
                        dest='call',
                        help='Call observational data')
    parser.add_argument('-g', '--engine',
                        required=False,
                        type=str,
                        default='mpfit',
                        dest='engine',
                        help='The fitter algorithm-engine: [{}]. Default: mpfit'.format(', '.join(engines())))
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
                        dest="input",
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
    parser.add_argument('-dt', '--dtshift',
                        required=False,
                        type=str,
                        default=None,
                        dest="dtshift",
                        help="The range of tshift in model LC. Default: None (any time). Format: {0}".format('2:50'))
    parser.add_argument('-q', '--quiet',
                        action='store_const',
                        const=True,
                        dest="is_quiet",
                        help="Just show result, no plots, no info")
    return parser


def engines(nm=None):
    switcher = {
        'mpfit': FitMPFit(),
        'mcmc': FitLcMcmc(),
    }
    if nm is not None:
        return switcher.get(nm)
    return list(switcher.keys())


def plot_curves(curves_o, res_models, res_sorted, **kwargs):
    from pystella.rf import light_curve_plot as lcp
    from matplotlib import pyplot as plt

    font_size = kwargs.get('font_size', 10)

    xlim = None
    ylim = None
    num = len(res_sorted)
    nrow = int(num / 2.1) + 1
    ncol = 2 if num > 1 else 1
    fig = plt.figure(figsize=(12, nrow * 4))
    plt.matplotlib.rcParams.update({'font.size': font_size})

    i = 0
    for k, v in res_sorted.items():
        i += 1
        if i > num:
            break
        ax = fig.add_subplot(nrow, ncol, i)

        curves = res_models[k]
        lcp.curves_plot(curves, ax=ax, figsize=(12, 8), linewidth=1, is_legend=False)

        lt = {lc.Band.Name: 'o' for lc in curves_o}
        lcp.curves_plot(curves_o, ax, lt=lt, markersize=3, is_legend=False)

        ax.text(0.99, 0.94, k, horizontalalignment='right', transform=ax.transAxes)
        ax.text(0.98, 0.85, "dt={:.2f}".format(v.tshift), horizontalalignment='right', transform=ax.transAxes)
        ax.text(0.01, 0.05, "$\chi^2: {:.2f}$".format(v.measure), horizontalalignment='left', transform=ax.transAxes, bbox=dict(facecolor='green', alpha=0.3))
        # ax.text(0.9, 0.9, "{:.2f}".format(v.measure), horizontalalignment='right', transform=ax.transAxes, bbox=dict(facecolor='green', alpha=0.3))

        # fix axes
        if i % ncol == 0:
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position("right")
            # ax.set_ylabel('')
            # ax.set_yticklabels([])
            # ax2.set_ylabel('Magnitude')
        else:
            ax.yaxis.tick_left()
            ax.yaxis.set_label_position("left")
            # ax.set_ylabel('Magnitude')
        ax.yaxis.set_ticks_position('both')
    fig.subplots_adjust(wspace=0, hspace=0)
    # return fig
    plt.show()


def main():
    model_ext = '.ph'
    t_diff = 1.01
    name = None
    # z = 0.
    # distance = 10  # pc
    bnames = ['U', 'B', 'V', 'R', "I"]
    Nbest = 6
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
    if args.bnames:
        bnames = []
        for bname in str(args.bnames).split('-'):
            if not band.band_is_exist(bname):
                print('No such band: ' + bname)
                parser.print_help()
                sys.exit(2)
            bnames.append(bname)

    # if args.distance:
    # distance = args.distance

    if args.call:
        callback = cb.lc_wrapper(str(args.call), method='load')
    else:
        print('No obs data. Use key: -c: ')
        parser.print_help()
        sys.exit(2)

    # Get observations
    curves_o = callback.load({'is_debug': not args.is_quiet})

    if len(bnames) == 0:
        bnames = [bn for bn in curves_o.BandNames if band_is_exist(bn)]

    # Get models
    times = (0, None)

    if args.times:
        times = list(map(float, args.times.split(':')))

    # The fit engine
    fitter = engines(args.engine)
    fitter.is_debug = not args.is_quiet  # fitter = FitMPFit(is_debug=not args.is_quiet)

    # The filter results by tshift
    if args.dtshift:
        dtshift = str2interval(args.dtshift, llim=float("-inf"), rlim=float('inf'))
    else:
        dtshift = (float("-inf"), float("inf"))

    # tshift = 0.
    res_models = {}
    res_chi = {}
    if len(names) == 1:
        name = names[0]
        if not args.is_quiet:
            if times is not None:
                print("Fitting for model %s %s for %s moments" % (path, name, args.times))
            else:
                print("Fitting for model %s %s " % (path, name))
        curves = lcf.curves_compute(name, path, bnames, z=args.redshift, distance=args.distance,
                                    t_beg=times[0], t_end=times[1], t_diff=t_diff)

        res = fitter.fit_curves(curves_o, curves)
        print("{}: time shift  = {:.2f}+/-{:.4f} Measure: {:.4f}".format(name, res.tshift, res.tsigma, res.measure))
        best_tshift = res.tshift
        res_models[name] = curves
        res_chi[name] = res
    elif len(names) > 1:
        i = 0
        for name in names:
            i += 1
            if args.is_quiet:
                sys.stdout.write((u"\u001b[1000D Fitting for model %30s  [%d/%d]" % (name, i, len(names))))
                sys.stdout.flush()
            else:
                if times is not None:
                    print("Fitting for model %s %s for %s moments" % (path, name, args.times))
                else:
                    print("Fitting for model %s %s " % (path, name))
            curves = lcf.curves_compute(name, path, bnames, z=args.redshift, distance=args.distance,
                                        t_beg=times[0], t_end=times[1], t_diff=t_diff)
            res = fitter.fit_curves(curves_o, curves)
            res_models[name] = curves
            res_chi[name] = res

        # select with dtshift
        res_chi_sel = {}
        for k, v in res_chi.items():
            if dtshift[0] < v.tshift < dtshift[1]:
                res_chi_sel[k] = v
        res_chi = res_chi_sel
        # sort with measure
        res_sorted = OrderedDict(sorted(res_chi.items(), key=lambda kv: kv[1].measure))
        print("\n Results (tshift in range:{:.2f} -- {:.2f}".format(dtshift[0], dtshift[1]))
        print("{:40s} ||{:18s}|| {:10}".format('Model', 'dt+-t_err', 'Measure'))
        for k, v in res_sorted.items():
            # if dtshift[0] < v.tshift < dtshift[1]:
            print("{:40s} || {:.2f}+/-{:.4f} || {:.4f}".format(k, v.tshift, v.tsigma, v.measure))

        best_mdl = list(res_sorted)[0]
        res = list(res_sorted.values())[0]
        # best_mdl, res = res_sorted.popitem(last=False)
        print("Best fit model:")
        print("{}: time shift  = {:.2f}+/-{:.4f} Measure: {:.4f}".format(best_mdl, res.tshift, res.tsigma, res.measure))
        best_tshift = res.tshift
    else:
        print("No any data about models. Path: {}".format(path))
        parser.print_help()
        sys.exit(2)

    if not args.is_quiet:
        curves_o.set_tshift(best_tshift)
        # show only Nbest modeles
        while len(res_sorted) > Nbest:
            res_sorted.popitem()
        # plot the best models
        plot_curves(curves_o, res_models, res_sorted)
        # print(" dm_abs      = {0:.4f}+/-{1:.4f}".format(dm, dmsigma))


if __name__ == '__main__':
    main()
