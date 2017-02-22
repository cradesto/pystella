#!/usr/bin/python3
# -*- coding: utf-8 -*-
import getopt

import sys

from pystella.util.phys_var import phys
import matplotlib.pyplot as plt
import pystella.rf.band as band

__author__ = 'bakl'


def plot_griuz():
    plt.title('griuz filter response')
    bands = dict(g='g', i='r+', r='r', u='o', z='*')
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl * phys.cm_to_angs, b.resp_wl, v, label=k)

    plt.legend()
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def plot_UBVRI():
    plt.title('UBVRI filter response')
    bands = dict(U='b', B='c', V='g', R='r', I='p')
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl * phys.cm_to_angs, b.resp_wl, v, label=k)
    plt.legend()
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def plot_JHK():
    plt.title('JHK filter response')
    bands = dict(J='b', H='c', K='r')
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl * phys.cm_to_angs, b.resp_wl, v, label=k)
    plt.legend()
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def plot_SWIFT():
    plt.title('SWIFT filter response')
    bands = dict(UVM2='m', UVW1='r', UVW2='b', SwiftU='k', SwiftB='c', SwiftV='g')
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl * phys.cm_to_angs, b.resp_wl, v, label=k)
    plt.legend()
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def plot_PS1():
    plt.title('The Pan-STARRS1 Photometric  filter responses')
    bands = dict(PS1g='m', PS1i='r', PS1r='b', PS1z='y', PS1y='g', PS1w='p')

    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl * phys.cm_to_angs, b.resp_wl, v, label=k)
    plt.legend()
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def plot_HSC():
    plt.title('The Hyper Suprime-Cam(HSC) Photometric  filter responses')

    bands = dict(HSCg='m', HSCr='b', HSCi='r', HSCz='y', HSCY='g')
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl * phys.cm_to_angs, b.resp_wl, v, label=k, linewidth=2)
    plt.legend(loc=4)
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def plot_HST():
    plt.title('The Hubble Space Telescope  filter responses')

    bands = dict(F105W="blue", F125W="g", F435W="skyblue", F140W="orange", F160W="r", F606W="cyan", F814W="magenta")
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl * phys.cm_to_angs, b.resp_wl, v, label=k, linewidth=2)
    plt.legend(loc=4)
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def plot_Kepler():
    plt.title('The Kepler Space Telescope  filter responses')

    bands = dict(Kepler="magenta")
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl * phys.cm_to_angs, b.resp_wl, v, label=k, linewidth=2)
    plt.legend(loc=4)
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def plot_bands(bands, color_dic=None):
    if color_dic is None:
        color_dic = band.bands_colors()

    for bname in bands:
        b = band.band_by_name(bname)
        plt.plot(b.wl * phys.cm_to_angs, b.resp_wl, color_dic[bname], label=bname, linewidth=2)
    plt.legend(loc=4)
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def usage():
    print("Usage:")
    print("  plot_bands.py [params]")
    print("  -b <bands>: string, default: U-B-V-R-I, for example U-B-V")
    print("  -h  print usage")
    print("   --- ")
    band.print_bands()


def main():
    band.Band.load_settings()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hb:")
    except getopt.GetoptError as err:
        print(str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    bands = []
    for opt, arg in opts:
        if opt == '-b':
            bands = str(arg).split('-')
            continue
        elif opt == '-h':
            usage()
            sys.exit(2)

    if len(bands) > 0:
        plot_bands(bands)
    else:
        plot_UBVRI()
        plot_JHK()
        plot_griuz()
        plot_SWIFT()
        plot_PS1()
        plot_HSC()
        plot_HST()
        plot_Kepler()


if __name__ == '__main__':
    main()
