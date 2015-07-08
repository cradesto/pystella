import os
from os.path import dirname
from model.stella import Stella
from rf import band
from util.arr_dict import dict_save
import numpy as np

import matplotlib.pyplot as plt

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))



def plot_bands(dict_mags, bands):
    plt.title(''.join(bands)+' filter response')
    lblbands = dict(U='U', B='B', V='V', R='R', I='I')
    bandshift = dict(U=6.9, B=3.7, V=0, R=-2.4, I=-4.7, UVM2=11.3, UVW1=10, UVW2=13.6, g=2.5, r=-1.2, i=-3.7)
    colors = dict(U="magenta", B="black", V="red", R="magenta", I="black", UVM2="green", UVW1="red", UVW2="blue", g="blue", r="cyan", i="green")
    dist = 24e6  # distance to SN 2013ab
    dm = 5*np.log10(dist)-5  # distance module
    x = dict_mags['time']
    lims = (-12, -23)
    lims = (39, 9)
    for n in bands:
        y = dict_mags[n]
        y += dm + bandshift[n]
        plt.plot(x, y, label=n, color=colors[n])
        # lims = (max(y), min(y))

    plt.gca().invert_yaxis()
    plt.ylim(lims)
    plt.legend()
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def compute_mag(bands):
    path = os.path.join(ROOT_DIRECTORY, 'tests', 'data', 'stella')
    model = Stella("cat_R1000_M15_Ni007_E15", path=path)
    model.show_info()

    serial_spec = model.read_serial_spectrum(t_diff=1.05)

    dict_mags = dict((k, None) for k in bands)
    dict_mags['time'] = serial_spec.times
    for n in bands:
        b = band.band_by_name(n)
        mags = serial_spec.flux_to_mags(b)
        if mags is not None:
            dict_mags[n] = mags

    return dict_mags


def main():
    bands = ['U', 'B', 'V', 'R', "I", 'UVM2', "UVW1", "UVW2", 'g', "r", "i"]
    dict_mags = compute_mag(bands)

    dict_save(dict_mags, 'test_mags_' + ''.join(bands) + '.txt')
    plot_bands(dict_mags, bands)


if __name__ == '__main__':
    main()