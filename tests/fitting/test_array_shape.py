import unittest

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

from pystella.util.math import shrink, portion_index


class TestArrayShape(unittest.TestCase):
    # @unittest.skip("just for plot")
    @staticmethod
    def gen_array(start=1, end=10, nums=(200, 100, 200)):
        y = []
        for n in nums:
            t = np.linspace(start, end, n)
            start, end = end, start
            y = np.append(y, t)
        return y

    def test_spline(self):
        # y = np.linspace(1, 10, 200)
        # y = np.append(y, np.linspace(10, 2, 100))
        # y = np.append(y, np.linspace(2, 5, 200))
        y = self.gen_array()
        x = np.arange(len(y))

        new_length = 25
        xx = np.linspace(x.min(), x.max(), new_length)
        yy = interp1d(x, y, kind='cubic')(xx)
        # plot
        plt.plot(x, y)
        plt.plot(xx, yy, marker='o')
        plt.show()

    def test_point_rm(self):
        # y = np.linspace(1, 10, 200)
        # y = np.append(y, np.linspace(10, 2, 100))
        # y = np.append(y, np.linspace(2, 5, 200))
        y = self.gen_array()
        x = np.arange(len(y))

        idxs = self.point_rm(y, mode='g')
        xx = x[idxs]
        yy = y[idxs]
        # plot
        plt.plot(x, y, ls='-.')
        plt.plot(xx, yy, marker='o')
        plt.show()

    def test_shrunk(self):
        # y = np.linspace(1, 10, 200)
        # y = np.append(y, np.linspace(10, 2, 100))
        # y = np.append(y, np.linspace(2, 5, 200))
        y = self.gen_array(nums=(200, 100, 200))
        idxs = shrink(y)
        # self.assertEqual(end-start, len(idxs), msg="(end-start)={} len(idxs)={}".format(end-start, len(idxs)))

        yy = y[idxs]
        # plot
        plt.plot(y, ls='-.')
        plt.plot(idxs, yy, marker='o')
        plt.show()

    def test_portion_index(self):
        # y = np.linspace(1, 10, 200)
        # y = np.append(y, np.linspace(10, 2, 100))
        # y = np.append(y, np.linspace(2, 5, 200))
        y = self.gen_array(nums=(200, 100, 200))

        # trivial
        start, end = 0, None
        idxs = portion_index(y, lambda i, x: i, start=start, end=end)
        self.assertEqual(len(y), len(idxs), msg="len(y)={} len(idxs)={}".format(len(y), len(idxs)))

        # left boundary
        start, end = 117, None
        idxs = portion_index(y, lambda i, x: i, start=start, end=end)
        self.assertEqual(len(y), len(idxs), msg="len(y)={} len(idxs)={}".format(len(y), len(idxs)))

        # right boundary
        start, end = 0, 99
        idxs = portion_index(y, lambda i, x: i, start=start, end=end)
        self.assertEqual(len(y), len(idxs), msg="len(y)={} len(idxs)={}".format(len(y), len(idxs)))

        # right negative boundary
        start, end = 0, -99
        idxs = portion_index(y, lambda i, x: None, start=start, end=end)
        self.assertEqual(abs(end), len(idxs), msg="len(y)={} len(idxs)={}".format(abs(end), len(idxs)))

        #  False boundary
        start, end = 44, 99
        idxs = portion_index(y, lambda i, x: None, start=start, end=end)
        self.assertEqual(len(y) - (end - start), len(idxs),
                         msg="(end-start)={} len(idxs)={}".format(end - start, len(idxs)))

        y = y[idxs]

        #  False boundary
        start, end = 74, 99
        idxs = portion_index(y, lambda i, x: None, start=start, end=end)
        self.assertEqual(len(y) - (end - start), len(idxs),
                         msg="(end-start)={} len(idxs)={}".format(end - start, len(idxs)))

        yy = y[idxs]
        # plot
        plt.plot(y, ls='-.')
        plt.plot(idxs, yy, marker='o')
        plt.show()

    @staticmethod
    def point_rm(a, diff=1.1, mode='l'):
        """
        l - linear progression
        g - geom progression
        """
        res = []
        prev = a[0]
        for i, x in enumerate(a):
            if prev != 0:
                if mode == 'g':
                    if x >= prev and abs(x / prev) < diff:
                        continue
                    if x < prev and abs(prev / x) < diff:
                        continue
                else:
                    if abs(x - prev) < diff * prev:
                        continue
            res.append(i)
            prev = x
        return res
