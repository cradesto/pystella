##
#  Interface to read stella h5-file
#  Usage example: >>> python3 h5stella.py  stella-model.h5
#

import os
import logging

import numpy as np
import sys

logging.basicConfig()
logger = logging.getLogger(__name__)
# logger.setLevel(logging.INFO)
# logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

try:
    import h5py
except ImportError:
    logging.debug('h5py failed to import', exc_info=True)
    pass

__author__ = 'bakl'

"""
Tools for reading HDF5 Stella's output 
"""


class H5Stella(object):
    def __init__(self, fname):
        """Load Stella model data from h5.  Required parameters:  filename."""
        self._fname = fname

    @property
    def Name(self):
        return os.path.splitext(os.path.basename(self.Fname))[0]

    @property
    def Fname(self):
        return self._fname

    @property
    def Freqmn(self):
        with h5py.File(self.Fname, "r") as h5f:
            logging.debug('Freqmn from {}'.format(self.Fname))
            freqmn = np.array(h5f.get('/parameters/freqmn'))
        return freqmn

    @property
    def Fh(self):
        logging.debug('Fh from {}'.format(self.Fname))
        with h5py.File(self.Fname, "r") as h5f:
            res = H5Fh('Fh').fill(h5f)
        return res

    @property
    def Fj(self):
        logging.debug('Fj from {}'.format(self.Fname))
        with h5py.File(self.Fname, "r") as h5f:
            res = H5Fj('Fj').fill(h5f)
        return res

    @property
    def Tau(self):
        logging.debug('Tau from {}'.format(self.Fname))
        with h5py.File(self.Fname, "r") as h5f:
            res = H5Tau('Tau').fill(h5f)
        return res

    @property
    def Iray(self):
        logging.debug('Iray from {}'.format(self.Fname))
        with h5py.File(self.Fname, "r") as h5f:
            res = H5Iray('Iray').fill(h5f)
        return res

    @property
    def Hyd(self):
        logging.debug('Hyd from {}'.format(self.Fname))
        with h5py.File(self.Fname, "r") as h5f:
            res = H5Hyd('Hyd').fill(h5f)
        return res

    @property
    def Yabun(self):
        logging.debug('Yabun from {}'.format(self.Fname))
        with h5py.File(self.Fname, "r") as h5f:
            yabun = np.array(h5f.get('/presn/Yabun'))  
        return yabun
    
    @property
    def Xisotopes(self):
        logging.debug('Xisotopes from {}'.format(self.Fname))
        with h5py.File(self.Fname, "r") as h5f:
            yabun = np.array(h5f.get('/presn/Xisotopes'))  
        return yabun    
    
    @property
    def M(self):
        logging.debug('Mass from {}'.format(self.Fname))
        with h5py.File(self.Fname, "r") as h5f:
            yabun = np.array(h5f.get('/presn/M'))  
        return yabun

    def ds_write(self, path, ds, attrs=None):
        """
        Write to h5-file dataset.
        :param path: The path inside h5-file
        :param ds: data to write
        :param attrs: Additional information. Default: None
        :return: void
        """
        if path == '':
            raise ValueError('path')
        if ds is None:
            raise ValueError('ds')
        print('Writing to {} dset [{}] with shape [{}]'.
              format(self.Fname, path, ':'.join(map(str, ds.shape))))
        with h5py.File(self.Fname, 'a') as h5f:
            if path in h5f:
                del h5f[path]
            dset = h5f.create_dataset(path, data=ds)
            if attrs is not None:
                for k, v in attrs.items():
                    dset.attrs[k] = v

            
class H5TimeElement(object):
    def __init__(self, name, path):
        self._name = name
        self._path = path
        self._attrs = None  # Attributes
        self._s = None  # step
        self._t = None  # time
        self._v = None  # value

    @property
    def Path(self):
        return self._path

    @property
    def Name(self):
        return self._name

    @property
    def Attrs(self):
        return self._attrs

    @property
    def S(self):
        """
        Step
        :return:
        """
        return self._s

    @property
    def Time(self):
        """
        Time
        :return:
        """
        return self._t

    @property
    def Val(self):
        """Value"""
        return self._v

    @property
    def Shape(self):
        """The shape of data"""
        return self._v.shape

    @property
    def Ntime(self):
        return self.Shape[0]

    @property
    def Nzon(self):
        return self.Shape[1]

    def fill(self, h5f):
        self._s = np.array(h5f.get(self.Path + '/step'))  # step
        self._t = np.array(h5f.get(self.Path + '/time'))  # time
        self._v = np.array(h5f.get(self.Path + '/value'))  # value
        attrs = h5f.get(self.Path).attrs
        if attrs is not None:
            self._attrs = {}
            for k, v in attrs.items():
                self._attrs[k] = v
        return self

    def IdxByTime(self, time):
        for idx, t in enumerate(self.Time):
            if t >= time:
                return idx
        return -1

    def __getitem__(self, item):
        return self.ValueByTime(item)

    def IdxValueByTime(self, time):
        idx = self.IdxByTime(time)
        if idx >= 0:
            return idx, self.Val[idx]
        return None, None

    def ValueByTime(self, nt):
        return self.Val[nt, :, :]

    def ValueByZone(self, nz):
        return self.Val[:, :, nz]

    def Info(self):
        print("Ntime={} self.Nzon={}".format(self.Ntime, self.Nzon))


class H5FreqTimeElement(H5TimeElement):
    def __init__(self, name, path):
        super(H5FreqTimeElement, self).__init__(name, path)

    @property
    def Nfreq(self):
        return self.Shape[1]

    @property
    def Nzon(self):
        return self.Shape[2]

    def Info(self):
        print("Ntime={} self.Nzon={} Nfreq={} ".format(self.Ntime, self.Nzon, self.Nfreq))


class H5Fh(H5FreqTimeElement):
    def __init__(self, name):
        self._name = name
        super(H5Fh, self).__init__(name, path='/timing/Fh')


class H5Fj(H5FreqTimeElement):
    def __init__(self, name):
        self._name = name
        super(H5Fj, self).__init__(name, path='/timing/Fj')


class H5Tau(H5FreqTimeElement):
    def __init__(self, name):
        self._name = name
        super(H5Tau, self).__init__(name, path='/timing/Tau')


class H5Iray(H5FreqTimeElement):
    def __init__(self, name):
        self._name = name
        super(H5Iray, self).__init__(name, path='/timing/Iray')

    @property
    def Nrays(self):
        return self.Shape[-1]


class H5Hyd(H5TimeElement):
    def __init__(self, name):
        self._name = name
        super(H5Hyd, self).__init__(name, path='/timing/Hyd')

    @property
    def Nvars(self):
        return self.Shape[-1]

    @property
    def Columns(self):
        return self.Attrs['columns'].decode()

    def Var(self, ncol):
        return self.Val[:, :, ncol]

    def __getattr__(self, name):
        """
        Extract variable as property with name from the "column" attribute.
        :param name: name of  variable
        :return: array of variable
        """
        if self.Attrs is None:
            raise ValueError('There are no any Attributes.')
        s = 'columns'
        if s not in self.Attrs:
            raise ValueError('There is no "{}" in  Attrs.'.format(s))
        columns = self.Attrs['columns'].decode()
        # return columns
        #
        # columns = columns
        if name not in columns:
            raise ValueError('There is no key: "{}" in  columns: {}.'.format(name, s))

        cols = columns.split()
        for k, c in enumerate(cols):
            if c == name:
                return self.Val[:, :, k]
        return None


def main():
    import matplotlib.pyplot as plt

    if len(sys.argv) < 2:
        print('Missed  h5-file')
        exit()

    fname = sys.argv[1]  # 'dumm030307mh5.h5'
    h5 = H5Stella(fname)

    # read frequencies vars
    Fj = h5.Fj
    Fj.Info()

    # plot Hyd
    hyd = h5.Hyd
    hyd.Info()

    nums = [0, 2, 5, 10, 20, 50]  # time steps
    columns = hyd.Attrs['columns'].split()
    plt.figure(figsize=(12, 12))
    for i, var in enumerate(columns):
        for num in nums:
            plt.subplot(2, 2, i + 1)
            x = hyd.Val[num, :, 0]
            y = hyd.Val[num, :, i]
            plt.plot(x, y, '-', label='{} t={:.3f}'.format(str(var, 'utf-8'), hyd.Time[num]))
            plt.xscale('log')
            plt.yscale('log')
        plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
