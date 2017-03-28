import scipy as sp
import scipy.interpolate
import numpy as np


def log_interp1d(x, y, kind='linear'):
    """
    The logarithmic interpolation.  See comment here: http://stackoverflow.com/a/29359275
    :param x: 
    :param y: 
    :param kind: (‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic, ‘cubic’) 
    See: https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.interpolate.interp1d.html#scipy.interpolate.interp1d
    :return: 
    """
    logx = np.log10(x)
    logy = np.log10(y)
    lin_interp = sp.interpolate.interp1d(logx, logy, kind=kind)
    return lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
