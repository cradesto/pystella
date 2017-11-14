import numpy as np
import unittest
import h5py
import unittest
from os.path import join, dirname, abspath

import matplotlib.pyplot as plt

import pystella.rf.light_curve_func as lcf
from pystella.model.stella import Stella

__author__ = 'bakl'


class TestHDF5(unittest.TestCase):
    def test_hdf5_simple(self):
        f = h5py.File("tmp.hdf5", "w")
        dset = f.create_dataset("mydataset", (100,), dtype='i')
        grp = f.create_group("subgroup")
        dset2 = grp.create_dataset("another_dataset", (50,), dtype='f')
        dset3 = f.create_dataset('subgroup2/dataset_three', (10,), dtype='i')
        dset.attrs['temperature'] = 99.5

        f.close()

    def test_curves_save(self):
        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        bands = ('UVW1', 'UVW2', 'UVM2')
        # bands = ('UVW1', 'UVW2', 'UVM2', 'U', 'B', 'R', 'I')
        mdl = Stella(name, path=path)
        curves = mdl.curves(bands)

        with h5py.File("tmp_curves.hdf5", "w") as f:
            grp = f.create_group("curves")
            for lc in curves:
                ds = grp.create_dataset(lc.Band.Name, data=lc.toarray())
                ds.attrs['z'] = lc.attrs('z')
                ds.attrs['d'] = lc.attrs('d')

            self.assertEqual(curves.Length, len(grp.items()))
