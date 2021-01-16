import os
import numpy as np

from .. import read_obs_table_header, SetLightCurve, LightCurve, Lum2MagBol
from ..rf import band

__author__ = 'bakl'


def load(name, path='.', tcol='time', cols=('L_ubvri', 'L_bol'), is_lum=False):
    """
    Load the bolometric light curves from lbol-file.
    :param name: model (file without extension)
    :param path: the model directory
    :param tcol: column with time, default: 'time'
    :param cols: columns with Lum, default: ('L_ubvri', 'L_bol')
    :param is_lum: if is_lum=True return log10(Lum) else abs.magnitudes. Default: False
    :return: SetLightCurve
    """
    fname = os.path.join(path, f'{name}.lbol')
    # print(f"Loading {fname}")
    d, header = read_obs_table_header(fname)
    # print(header)
    time = np.array(d[tcol])
    cut = time > 0
    t = time[cut]
    res = SetLightCurve(f'L_bol {name}')
    for col in cols:  # header.items():
        b = band.BandUni(name=f"{col}")
        lum = np.array(d[col])[cut]
        if not is_lum:  # convert to magnitudes
            lum = Lum2MagBol(10. ** lum)
        lc = LightCurve(b, t, lum)
        res.add(lc)
    return res
