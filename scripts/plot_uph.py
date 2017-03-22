#!/usr/bin/python3
# -*- coding: utf-8 -*-

import getopt
import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# plot both results


def read_pystella_swd_uph(fname, path=None):
    if path is not None:
        fname = os.path.join(path, fname)
    if not os.path.isfile(fname):
        logger.error(' No uph-data for %s' % fname)
        raise ValueError(' No uph-data for %s' % fname)
    # return None
    logger.info(' Load uph-data from  %s' % fname)
    col_names = "time zone R V "
    dt = np.dtype({'names': col_names.split(), 'formats': np.repeat('f8', len(col_names))})
    data = np.loadtxt(fname, comments='#', dtype=dt)
    return data


def read_uph(fname, path=None):
    if path is not None:
        fname = os.path.join(path, fname)
    if not os.path.isfile(fname):
        logger.error(' No uph-data for %s' % fname)
        raise ValueError(' No uph-data for %s' % fname)
    # return None
    logger.info(' Load uph-data from  %s' % fname)
    col_names = "time R V"
    dt = np.dtype({'names': col_names.split(), 'formats': np.repeat('f8', len(col_names))})
    data = np.loadtxt(fname, comments='#', dtype=dt)
    return data


def read_ph_swd(fname, path=None):
    if path is not None:
        fname = os.path.join(path, fname)
    if not os.path.isfile(fname):
        logger.error(' No ph-swd-data for %s' % fname)
        raise ValueError(' No ph-swd-data for %s' % fname)
    # return None
    logger.info(' Load ph-swd-data from  %s' % fname)
    col_names = "time R V M T"
    dt = np.dtype({'names': col_names.split(), 'formats': np.repeat('f8', len(col_names))})
    data = np.loadtxt(fname, comments='#', dtype=dt)
    return data


def usage():
    print("Usage:")
    print("  plot_uph.py  uph-file ph-swd-pystella   ph-swd-file ")


def main():
    lw = 2
    fname_pystella = None  # 's15s7b2v1z532E1.swd.ph'
    fname_swd = None  # 's15s7b2v1z532E1.swd.ph'
    fname_uph = None  # 'uph_s15s7b2v1z532E1.txt'

    try:
        opts, args = getopt.getopt(sys.argv[1:], "h")
    except getopt.GetoptError as err:
        print(str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(args) > 2:
        fname_uph = args[0]
        fname_pystella = args[1]
        fname_swd = args[2]
    elif len(args) > 1:
        fname_uph = args[0]
        fname_swd = args[1]
    elif len(args) > 0:
        fname_uph = args[0]
    elif len(opts) == 0:
        usage()
        sys.exit(2)

    # setup plot
    plt.matplotlib.rcParams.update({'font.size': 12})
    fig = plt.figure(figsize=(8, 8))

    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('Time [day]')
    ax.set_ylabel(r'Velocity [$\times 10^3$ km/s]')
    #    ax.set_title(fname_uph)

    # read and plot ph-swd-data
    if fname_swd is not None:
        dswd = read_ph_swd(fname_swd)
        x = dswd['time']
        y = dswd['V']
        ax.plot(x, y, label=fname_swd, color='blue', ls="-", linewidth=lw)

    # read and plot uph-data
    if fname_uph is not None:
        dswd = read_uph(fname_uph)
        x = dswd['time']
        y = dswd['V']
        ax.plot(x, y, label=fname_uph, color='orange', ls="--", linewidth=lw)

    # read and plot uph-data
    if fname_pystella is not None:
        dswd = read_pystella_swd_uph(fname_pystella)
        x = dswd['time']
        y = dswd['V']
        ax.plot(x, y, label=fname_pystella, color='red', ls=":", linewidth=lw)

    ax.legend()
    plt.show()


if __name__ == '__main__':
    main()
