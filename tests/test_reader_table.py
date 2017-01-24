# coding=utf-8
import numpy as np
import unittest

from pystella.rf import band
from pystella.rf.lc import LightCurve
from pystella.util.reader_table import read_table_header_float, table2curves

__author__ = 'bakl'


def lc_create(b, m=-19, dt=0.):
    n = 10
    time = np.linspace(0. + dt, 200. + dt, n)
    mags = m * np.ones(n)
    return LightCurve(b, time, mags)


class TestReaderTable(unittest.TestCase):
    def test_read_table_header_float(self):
        fname = 'data/stella/cat_R500_M15_Ni006_E12.gri'
        data = read_table_header_float(fname)
        cols = len(data.dtype.names)
        self.assertTrue(cols == 15,
                        msg="The number of colums in the data should be 15, but it's : %d." % cols)

    def test_read_table_header_float_skiprows(self):
        fname = 'data/stella/rednova_R3.2_M6_Ni0_E0.25.tt'
        data = read_table_header_float(fname, skip=87)
        cols = len(data.dtype.names)
        self.assertTrue(cols == 14,
                        msg="The number of colums in [%s] should be 14, but it's : %d." % (fname, cols))

    def test_table2curves_no_bands(self):
        band.Band.load_settings()
        fname = 'data/stella/rednova_R3.2_M6_Ni0_E0.25.tt'
        data = read_table_header_float(fname, skip=87)
        data.dtype.names = [col.replace('M', '') for col in data.dtype.names]
        curves = table2curves('test', data)
        for bname in curves.BandNames:
            self.assertTrue(bname in data.dtype.names,
                            msg="No band %s in [%s] after table2curves." % (bname, ''.join(data.dtype.names)))
