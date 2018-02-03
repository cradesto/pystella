import os
import logging

import h5py
import numpy as np
import sys

__author__ = 'bakl'

"""
Tools for reading HDF5 Stella's output 
"""

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)


class H5Stella(object):
    def __init__(self, fname):
        """Load Stella model data from h5.  Required parameters:  filename."""
        self._fname = fname

    @property
    def Name(self):
        return os.path.basename(self.Fname)[0]

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

    def ds_write(self, path, ds, attrs=None):
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
    def S(self):
        return self._s

    @property
    def T(self):
        return self._t

    @property
    def V(self):
        return self._v

    @property
    def Shape(self):
        return self._v.shape

    @property
    def Ntime(self):
        return self.Shape[0]

    @property
    def Nfreq(self):
        return self.Shape[1]

    @property
    def Nzon(self):
        return self.Shape[2]

    def fill(self, h5f):
        self._s = np.array(h5f.get(self.Path + '/step'))  # step
        self._t = np.array(h5f.get(self.Path + '/time'))  # time
        self._v = np.array(h5f.get(self.Path + '/value'))  # value
        return self

    def IdxByTime(self, time):
        for idx, t in enumerate(self.T):
            if t >= time:
                return idx
        return -1

    def IdxValueByTime(self, time):
        idx = self.IdxByTime(time)
        if idx >= 0:
            return idx, self.V[idx]
        return None, None

    def ValueByTime(self, nt):
        return self.V[nt, :, :]

    def ValueByZone(self, nz):
        return self.V[:, :, nz]

    def Info(self):
        print("Ntime={} self.Nzon={} Nfreq={} ".format(self.Ntime, self.Nzon, self.Nfreq))


class H5Fh(H5TimeElement):
    def __init__(self, name):
        self._name = name
        super(H5Fh, self).__init__(name, path='/timing/Fh')


class H5Fj(H5TimeElement):
    def __init__(self, name):
        self._name = name
        super(H5Fj, self).__init__(name, path='/timing/Fj')
