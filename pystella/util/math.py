import scipy.interpolate
import numpy as np


def log_interp1d(x, y, kind='linear'):
    """
    The logarithmic interpolation.  See comment here: http://stackoverflow.com/a/29359275
    :param x: 
    :param y: 
    :param kind: (‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic, ‘cubic’) 
    See: https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.interpolate.interp1d.html
    :return:
    """
    logx = np.log10(x)
    logy = np.log10(y)
    lin_interp = scipy.interpolate.interp1d(logx, logy, kind=kind)
    return lambda zz: np.power(10.0, lin_interp(np.log10(zz)))


def is_number(s):
    try:
        float(s)
        return True
    except:
        return False


def portion_value(a, func, start=0, end=None):
    """
    l - linear progression
    g - geom progression
    """
    if start < 0:
        raise ValueError("start < 0")
    if end is None:
        end = len(a)
    if start > end:
        raise ValueError("start > end")
    if end > len(a):
        raise ValueError("end > len(a)")

    aa = func(a[start:end])
    res = np.concatenate((a[:start], aa, a[end + 1:]))
    return res


def portion_index(a, where, start=0, end=None, isByEl=True):
    """
    :param a:  array
    :type where: condition(i, a)
    :param start: the left boundary to work with condition
    :param end: the right boundary to work with condition
    :param isByEl:
    :return: the element indexes satisfying the condition
    """
    if start < 0:
        raise ValueError("start < 0")
    if end is None or end > len(a):
        end = len(a)
    if end < 0:
        end = len(a) + end
    if start > end:
        raise ValueError("start > end")

    if isByEl:
        idxs = []
        for i in np.arange(start, end, dtype=int):
            if not where(i, a[i]) is None:
                idxs.append(i)
    else:
        idxs = where(a[start:end])
        idxs = np.add(idxs, start)  # add start index

    if len(idxs) > 0:
        res = np.concatenate((np.arange(start, dtype=int), idxs, np.arange(end, len(a), dtype=int)))
    else:
        res = np.concatenate((np.arange(start, dtype=int), np.arange(end, len(a), dtype=int)))

    # if not np.all(np.diff(idxs) > 0):
    #     raise ArithmeticError("not np.all(np.diff(idxs) > 0): {}".format(idxs[np.diff(idxs) < 0.]))
    return res


def shrink(a, diff=1.1, mode='g'):
    """
    l - linear progression
    g - geom progression
    """
    res = [0]
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
