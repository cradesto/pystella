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
    bands = dict(UVM2='m', UVW1='r', UVW2='b', U_UVOT='k', B_UVOT='c', V_UVOT='g')
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
    bands = dict(gps1='m', ips1='r', rps1='b', zps1='y', y='g', w='p')
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl*phys.cm_to_angs,  b.resp, v, label=k)
    plt.legend()
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

if __name__ == '__main__':
    main()