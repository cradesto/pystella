#!/usr/bin/env python3
# #!/usr/bin/python3

import getopt
import numpy as np
import sys

import matplotlib
# matplotlib.use("Agg")
matplotlib.rcParams['backend'] = "Qt4Agg"
import matplotlib.pyplot as plt
from matplotlib import gridspec

from pystella.rf import spectrum, rad_func
from pystella.util.phys_var import phys


def plot_bb(Trad, Tcol, W1, W2, wl_lim=(100., 1e5, 100)):

    # freq
    # nf, start, end = 100, 10., 1e5
    wl = np.exp(np.linspace(np.log(wl_lim[1]), np.log(wl_lim[0]), wl_lim[2]))
    freq = rad_func.val_to_hz(wl, inp="A")

    sp_1 = spectrum.SpectrumDilutePlanck(freq, Trad, W1)
    sp_2 = spectrum.SpectrumDilutePlanck(freq, Tcol, W2)

    # setup figure
    plt.matplotlib.rcParams.update({'font.size': 14})
    fig = plt.figure(num=None, figsize=(7, 11), dpi=100, facecolor='w', edgecolor='k')

    gs1 = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs1[0, 0])
    gs1.update(wspace=0.3, hspace=0.3, left=0.15, right=0.95)

    # plot Trad
    x = sp_1.Wl * phys.cm_to_angs
    y = sp_1.FluxWl
    ax.semilogy(x, y, color='red', ls=":", linewidth=2.5, label='Diluted 1: Tcol=%6.1f W=%4.2f' % (sp_1.T, sp_1.W))

    # plot Tcol & W
    xx = sp_2.Wl * phys.cm_to_angs
    yy = sp_2.FluxWl
    ax.semilogy(xx, yy, color='blue', ls="--", linewidth=2.5, label='Diluted 2: Tcol=%6.1f W=%4.2f' % (sp_2.T, sp_2.W))

    ymax = np.max([y, yy])
    ylim = [ymax*1e-15, ymax*1e1]
    ax.set_ylim(ylim)

    ax.set_ylabel(r'$F_\lambda, \, [erg\, s^{-1} cm^2]$')
    ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set_xlabel(r'$\lambda, \, [\AA]$')
    ax.legend(prop={'size': 9}, loc=4, borderaxespad=0.)
    plt.grid()
    plt.show()


def usage():
    print("Usage:")
    print("  plot_bb.py [params]")
    print("  -c <Tcol1>: color temperature for W*B(tcol)")
    print("  -r <Tcol2>: color temperature for W*B(Trad)")
    print("  -w <dilution1>: dilution for W*B(tcol)")
    print("  -g <dilution2>: dilution for W*B(tcol)")
    print("  -n <nfreq>: the number of frequencies")
    print("  -l <lambda1:lambda2>: wave length range")
    print("  -h  print usage")


def main():
    Trad = 6e3
    Tcol = Trad
    W1 = 1.
    W2 = 1.
    nfreq = 100
    wl_lim = None

    try:
        opts, args = getopt.getopt(sys.argv[1:], "c:g:l:n:r:w:")
    except getopt.GetoptError as err:
        print(str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(opts) == 0:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-n':
            nfreq = int(arg)
            continue

    for opt, arg in opts:
        if opt == '-c':
            Tcol = float(arg)
            continue
        if opt == '-r':
            Trad = float(arg)
            continue
        if opt == '-w':
            W1 = float(arg)
            continue
        if opt == '-g':
            W2 = float(arg)
            continue
        if opt == '-l':
            wl_lim = tuple(map(float, arg.split(':'))) + (nfreq,)
            continue
        elif opt == '-h':
            usage()
            sys.exit(2)

    if wl_lim is None:
        wl_lim = [1e2, 5e4, nfreq]

    plot_bb(Trad, Tcol, W1, W2, wl_lim=wl_lim)


if __name__ == '__main__':
    main()
