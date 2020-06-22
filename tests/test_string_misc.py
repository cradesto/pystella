# coding=utf-8
import numpy as np
import unittest

import pystella as ps

__author__ = 'bakl'


class TestStringMisc(unittest.TestCase):
    def test_pd2np(self):
        import pandas as pd

        d = {'dt': [22, 25, 27, 35],
             'dterrp': [2, 3, 4, 5],
             'dterrm': [2.2, 2.5, 2.7, 3.5],
             'dm': [-1, -1.25, -2.27, -2.35],
             'dmerr': [-0.1, -0.25, -0.27, -0.35]
             }
        t_init = tuple(sorted(d.keys()))
        df = pd.DataFrame(d)

        print(df)
        arr = ps.util.pd2np(df)
        names = arr.dtype.names
        t_np = tuple(sorted(names))
        print('numpy array: ', names)
        print(arr)
        self.assertEqual(len(d), len(names))
        self.assertTupleEqual(t_init, t_np)

    def test_np2tex(self):
        import pandas as pd

        d = {'dt': [22, 25, 27, 35],
             'dterrp': [2, 3, 4, 5],
             'dterrm': [2.2, 2.5, 2.7, 3.5],
             'dm': [-1, -1.25, -2.27, -2.35],
             'dmerr': [-0.1, -0.25, -0.27, -0.35]
             }
        df = pd.DataFrame(d)

        arr = ps.util.pd2np(df)
        ss = ps.util.np2tex(arr)
        for s in ss:
            print(s)
