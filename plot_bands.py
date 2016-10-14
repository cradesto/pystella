#!/usr/bin/python
# -*- coding: utf-8 -*-

from pystella.util.phys_var import phys
import matplotlib.pyplot as plt
import pystella.rf.band as band

__author__ = 'bakl'


def plot_griuz():
    plt.title('griuz filter response')
    bands = dict(g='g', i='r+', r='r', u='o', z='*')
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl*phys.cm_to_angs,  b.resp, v, label=k)

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
        plt.plot(b.wl*phys.cm_to_angs,  b.resp, v, label=k)
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
        plt.plot(b.wl*phys.cm_to_angs,  b.resp, v, label=k)
    plt.legend()
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def plot_SWIFT():
    plt.title('SWIFT filter response')
    bands = dict(UVM2='m', UVW1='r', UVW2='b', UVOTU='k', UVOTB='c', UVOTV='g')
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl*phys.cm_to_angs,  b.resp, v, label=k)
    plt.legend()
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def plot_PS1():
    plt.title('The Pan-STARRS1 Photometric  filter responses')
    bands = dict(PS1g='m', PS1i='r', PS1r='b', PS1z='y', y='g', w='p')

    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl*phys.cm_to_angs,  b.resp, v, label=k)
    plt.legend()
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def plot_HSC():
    plt.title('The Hyper Suprime-Cam(HSC) Photometric  filter responses')

    bands = dict(HSCg='m', HSCr='b', HSCi='r', HSCz='y', HSCy='g')
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl*phys.cm_to_angs,  b.resp, v, label=k, linewidth=2)
    plt.legend(loc=4)
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def plot_HST():
    plt.title('The Hubble Space Telescope  filter responses')

    bands = dict(F105W="blue",  F125W="g", F435W="skyblue",  F140W="orange", F160W="r", F606W="cyan", F814W="magenta")
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl*phys.cm_to_angs,  b.resp, v, label=k, linewidth=2)
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
        plt.plot(b.wl*phys.cm_to_angs,  b.resp, v, label=k, linewidth=2)
    plt.legend(loc=4)
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def main():
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
