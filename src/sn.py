import os
from os.path import dirname
from io.stella import Stella
from rf import band
from util.arr_dict import dict_save

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))


def compute_mag():
    path = os.path.join(ROOT_DIRECTORY, 'tests', 'data', 'stella')
    model = Stella("cat_R1000_M15_Ni007_E15", path=path)
    model.show_info()

    serial_spec = model.read_serial_spectrum()

    bands = ['U', 'B', 'V', 'R', "I"]
    dict_mags = dict((k, None) for k in bands)
    dict_mags['time'] = serial_spec.times
    for n in bands:
        b = band.band_by_name(n)
        mags = serial_spec.flux_to_mag(b)
        if mags is not None:
            dict_mags[n] = mags

    dict_save(dict_mags, 'test_mags_'+''.join(bands)+'.txt')


def main():
    compute_mag()

if __name__ == '__main__':
    main()