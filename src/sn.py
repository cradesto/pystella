import os
from os.path import dirname
from model.stella import Stella
from rf import band
from util.arr_dict import dict_save

import matplotlib.pyplot as plt

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))



def plot_bands(dict_mags, bands):
    plt.title(''.join(bands)+' filter response')
    lblbands = dict(U='U', B='B', V='V', R='R', I='I')
    x = dict_mags['time']
    lims = (-12, -23)
    for n in bands:
        y = dict_mags[n]
        plt.plot(x, y, label=n)
        lims = (max(y), min(y))

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