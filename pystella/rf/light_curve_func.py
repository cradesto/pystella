import csv
import os

import numpy as np
import pystella.rf.rad_func as rf
from pystella.model.stella import Stella
from pystella.rf import extinction
from pystella.rf.lc import LightCurve, SetLightCurve
from pystella.rf.reddening import ReddeningLaw

__author__ = 'bakl'


def compute_mag(name, path, bands, ext=None, z=0., distance=10., magnification=1., t_diff=1.05, is_show_info=True,
                is_save=False):
    """
        Compute magnitude in bands for the 'name' model.
    :param name: the name of a model and data files
    :param path: the directory with data-files
    :param bands: photometric bands
    :param ext: extinction
    :param z: redshift, default 0
    :param distance: distance to star in parsec, default 10 pc
    :param magnification: gravitational lensing magnification
    :param t_diff:  depression between time points
    :param is_show_info: flag to write some information, default True
    :param is_save: flag to save result in file, default False
    :return: dictionary with keys = bands, value = star's magnitudes
    """
    model = Stella(name, path=path)
    if is_show_info:
        print('')
        model.info()

    if not model.is_ph:
        model.info()
        print("Error: No data for: " + str(model))
        return None

    serial_spec = model.get_ph(t_diff=t_diff)
    mags = serial_spec.mags_bands(bands, z=z, d=rf.pc_to_cm(distance), magnification=magnification)

    if mags is not None:
        fname = os.path.join(path, name + '.ubv')
        if is_save:
            mags_save(mags, bands, fname)
            print("Magnitudes have been saved to " + fname)

    if is_show_info:
        # print the time of maximum LC
        tmin = 2.0
        t = mags['time']
        for n in bands:
            t_min = t[t > tmin][mags[n][t > tmin].argmin()]
            print("t_max(%s) = %f" % (n, t_min))

    if ext is not None:  # add extinction
        for n in bands:
            mags[n] = mags[n] + ext[n]

    return mags


def curves_save_mags(curves, fname, sep='\t'):
    """
       Save curves to CSV-format. It required for correct operation the common time for all LC.
    :param curves:
    :param fname:
    :param sep:
    :return:
    """
    if curves.Length == 0:
        print("Nothing to save: curves.Length=%d" % curves.Length)
        return False

    if not curves.IsCommonTime:
        print("The curves has different time arrays.")
        print("  LC  |  len(time)  ")
        for lc in curves:
            print("{:40s}   {:d}".format(lc.Name, len(lc.Time)))
        return False

    with open(fname, 'w') as f:
        writer = csv.writer(f, delimiter=sep, quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['{:^8s}'.format(x) for x in ['time'] + curves.BandNames])
        for i, (row) in enumerate(zip(curves.TimeCommon, *[curves.get(b).Mag for b in curves.BandNames])):
            # row = row[-1:] + row[:-1]  # make time first column
            writer.writerow(['{:8.3f}'.format(x) for x in row])
            # writer.writerow(['{:3.4e}'.format(x) for x in row])

    return True


def curves_save(curves, fname, sep='\t'):
    """
       Save curves to CSV-format. It required for correct operation the common time for all LC.
    :param curves:
    :param fname:
    :param sep:
    :return:
    """
    if curves.Length == 0:
        print("Nothing to save: curves.Length=%d" % curves.Length)
        return False

    if not curves.IsCommonTime:
        print("The curves has different time arrays.")
        print("  LC  |  len(time)  ")
        for lc in curves:
            print("{:40s}   {:d}".format(lc.Name, len(lc.Time)))
        return False

    arr = curves2nparray(curves)
    fmt_header = "%10s " * len(arr.dtype.names)
    header = fmt_header % arr.dtype.names
    fmt = "%10.5f  " * len(arr.dtype.names)
    np.savetxt(fname, arr, delimiter=sep, header=header, comments='', fmt=fmt)
    return True


def curves_read(fname, is_out=False):
    """
       Read curves from file with header-format. It required the correct format of header, ex:
       time          B       errB          I          V       errV          U       errU
    :param fname: file with data
    :param is_out: show additional information during parsing
    :return: curves
    """
    from pystella.util.reader_table import read_obs_table_header, table2curves
    if not os.path.exists(fname):
        raise("No file to read: {}".format(fname))

    from pystella.rf import band
    band.Band.load_settings()
    tbl, cols_data = read_obs_table_header(fname, is_out=is_out)
    curves = table2curves(os.path.basename(fname), tbl)
    return curves


def curves2nparray(curves):
    """
       Convert curves to numpy array.
    :param curves:
    :return:
    """
    if curves.Length == 0:
        print("Nothing to convert: curves.Length=%d" % curves.Length)
        return None

    if not curves.IsCommonTime:
        print("The curves has different time arrays.")
        print("  LC  |  len(time)  ")
        for lc in curves:
            print("{:40s}   {:d}".format(lc.Name, len(lc.Time)))
        return None

    dtype = {'names': ['time']}
    data = [curves.TimeCommon]
    for lc in curves:
        # res = np.hstack((res, lc.Mag.reshape(lc.Length, 1)))
        # res = np.column_stack((res, lc.Mag))
        dtype['names'].append(lc.Band.Name)
        data.append(lc.Mag)
        if lc.IsErr:
            # res = np.hstack((res, lc.Err.reshape(lc.Length, 1)))
            # res = np.column_stack(res, lc.Err)
            dtype['names'].append('err'+lc.Band.Name)
            data.append(lc.Err)

    res = np.column_stack(data)
    # res = np.dstack(data)
    dtype['formats'] = [np.float64] * len(dtype['names'])
    res.dtype = dtype
    return res.ravel()


def curves_compute(name, path, bands, z=0., distance=10., magnification=1.,
                   **kwargs):
    """
        Compute magnitude in bands for the 'name' model.
    :param name: the name of a model and data files
    :param path: the directory with data-files
    :param bands: photometric bands
    :param z: redshift, default 0
    :param distance: distance to star in parsec, default 10 pc
    :param magnification: gravitational lensing magnification
    :return: dictionary with keys = bands, value = star's magnitudes
    """
    t_beg = kwargs.get("t_beg", 0.)
    t_end = kwargs.get("t_end", float('inf'))
    is_show_info = kwargs.get("is_show_info", False)
    t_diff = kwargs.get("t_diff", 1.05)

    if len(bands) == 0:
        raise ValueError("You have not set any bands for model: " + str(name))

    model = Stella(name, path=path)
    if not model.is_ph:
        model.info()
        raise ValueError("Error: No spectral data for: " + str(model))

    if is_show_info:
        print('')
        model.info()

    serial_spec = model.get_ph(t_diff=t_diff, t_beg=t_beg, t_end=t_end)
    curves = serial_spec.flux_to_curves(bands, z=z, d=rf.pc_to_cm(distance), magnification=magnification)

    if is_show_info:
        # print the time of maximum LC
        tmin = 2.0
        t = curves.TimeCommon
        for n in bands:
            t_min = t[t > tmin][curves[n][t > tmin].argmin()]
            print("t_max(%s) = %f" % (n, t_min))

    return curves


def flux_reddening(freq, flux, ebv, law=ReddeningLaw.MW):
    """
    Apply extinction curves to flux(freq) values
    :param freq:  [Hz]
    :param flux:  [ergs s^-1 cm^-2 Hz^-1]
    :param ebv: E(B-V)
    :param law: type of extinction curves (MW, LMC, SMC)
    :return: reddening flux
    """
    from pystella.util.phys_var import phys

    flux_wl = rf.Fnu2Fwl(freq, flux)
    # flux_wl = flux * freq ** 2 / phys.c
    wl = np.array(phys.c / freq) * phys.cm_to_angs
    res = flux_reddening_wl(wl, flux_wl, ebv, law)
    res = rf.Fwl2Fnu(freq, res)
    # res = flux_reddening_wl(wl, flux_wl, ebv, law) * phys.c / freq**2
    return res


def flux_reddening_wl(wl, flux_wl, ebv, law=ReddeningLaw.MW):
    """
    Apply extinction curves to flux(lambda) values
    :param wl:  [A]
    :param flux_wl:  [ergs s^-1 cm^-2 A^-1]
    :param ebv: E(B-V)
    :param law: type of extinction curves (MW, LMC, SMC)
    :return: reddening flux
    """
    from pystella.rf.reddening import ReddeningLaw
    A_lambda = ReddeningLaw.Almd(wl, ebv, law)
    res = flux_wl * 10 ** (-0.4 * A_lambda)
    return res


def series_spec_reddening(series, ebv, law=ReddeningLaw.MW):
    from pystella.rf.spectrum import Spectrum, SeriesSpectrum
    res = SeriesSpectrum(series.Name)
    for k, t in enumerate(series.Time):
        s = series.get_spec(k)
        freq = s.Freq
        flux_red = flux_reddening(freq, s.Flux, ebv, law=law)
        ss = Spectrum(s.Name, freq, flux=flux_red)
        res.add(t, ss)
    return res


# todo Implement direct reddening from spectra
# http://webast.ast.obs-mip.fr/hyperz/hyperz_manual1/node10.html
def curves_reddening(curves, ebv, z=None, law=extinction.law_default, is_info=True):
    if ebv < 0.:
        raise ValueError("ebv should be > 0")
    if ebv == 0.:
        return curves

    is_instance_lc = isinstance(curves, LightCurve)
    if is_instance_lc:
        lc = curves
        curves = SetLightCurve(lc.Name)
        curves.add(lc)
        is_instance_lc = True

    bands = curves.BandNames

    if z is not None and z > 0.1:
        ext = extinction.reddening_z(ebv=ebv, bands=bands, z=z, law=law, is_info=is_info)
    else:
        ext = extinction.reddening(ebv=ebv, bands=bands, law=law, is_info=is_info)

    res = SetLightCurve(curves.Name)
    for lc in curves:
        lc_red = lc.copy()
        lc_red.M += ext[lc.Band.Name]
        res.add(lc_red)
    # for n in bands:
    #     curves[n].mshift = curves[n].mshift + ext[n]

    if is_instance_lc:
        return res.get(bands[0])
    return res


def mags_save(dictionary, bands, fname):
    with open(fname, 'wb') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['{:^8s}'.format(x) for x in ['time'] + bands])
        for i, (row) in enumerate(zip(*[dictionary[k] for k in 'time'.split() + bands])):
            # row = row[-1:] + row[:-1]  # make time first column
            writer.writerow(['{:8.3f}'.format(x) for x in row])
            # writer.writerow(['{:3.4e}'.format(x) for x in row])
