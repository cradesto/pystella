import csv
import os

import pystella.rf.rad_func as rf
from pystella.model.stella import Stella
from pystella.rf import extinction
from pystella.rf.lc import LightCurve, SetLightCurve

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
        model.show_info()

    if not model.is_ph_data:
        model.show_info()
        print("Error: No data for: " + str(model))
        return None

    # serial_spec = model.read_serial_spectrum(t_diff=0.)
    serial_spec = model.read_series_spectrum(t_diff=t_diff)
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


def curves_save(curves, fname, sep='\t'):
    """
       Save curves to CSV-format. It required for correct operation the common time for all LC.
    :param curves:
    :param fname:
    :return:
    """
    if curves.Length > 0:
        with open(fname, 'w') as f:
            writer = csv.writer(f, delimiter=sep, quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['{:^8s}'.format(x) for x in ['time'] + curves.BandNames])
            for i, (row) in enumerate(zip(curves.TimeDef, *[curves.get(b).Mag for b in curves.BandNames])):
                # row = row[-1:] + row[:-1]  # make time first column
                writer.writerow(['{:8.3f}'.format(x) for x in row])
                # writer.writerow(['{:3.4e}'.format(x) for x in row])
    else:
        print("Nothing to save: curves.Length=%d" % curves.Length)


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
    :param t_end:
    :param t_beg:
    :param is_show_info: flag to write some information, default True
    :return: dictionary with keys = bands, value = star's magnitudes
    """
    t_beg = kwargs.get("t_beg", 0.)
    t_end = kwargs.get("t_end", None)
    is_show_info = kwargs.get("is_show_info", False)
    t_diff = kwargs.get("t_diff", 1.05)

    if len(bands) == 0:
        raise ValueError("You have not set any bands for model: " + str(name))

    model = Stella(name, path=path)
    if not model.is_ph_data:
        model.show_info()
        raise ValueError("Error: No spectral data for: " + str(model))

    if is_show_info:
        print('')
        model.show_info()

    # serial_spec = model.read_serial_spectrum(t_diff=0.)
    serial_spec = model.read_series_spectrum(t_diff=t_diff, t_beg=t_beg, t_end=t_end)
    curves = serial_spec.flux_to_curves(bands, z=z, d=rf.pc_to_cm(distance), magnification=magnification)
    # curves = SetLightCurve(name)
    # for n in bands:
    #     b = band.band_by_name(n)
    #     lc = serial_spec.flux_to_curve(b, z=z, dl=rf.pc_to_cm(distance), magnification=magnification)
    #     # time = serial_spec.times * (1. + z)
    #     # lc = LightCurve(b, time, mags)
    #     curves.add(lc)

    if is_show_info:
        # print the time of maximum LC
        tmin = 2.0
        t = curves.TimeDef
        for n in bands:
            t_min = t[t > tmin][curves[n][t > tmin].argmin()]
            print("t_max(%s) = %f" % (n, t_min))

    return curves


def curves_reddening(curves, ebv, z=None, law=extinction.law_default, is_info=True):
    if ebv < 0.:
        raise ValueError("ebv should be > 0")
    if ebv == 0.:
        return curves

    is_instance_lc = isinstance(curves, LightCurve)
    if is_instance_lc:
        lc = curves
        # bands = tuple(curves.Band.Name
        curves = SetLightCurve(lc.Name)
        curves.add(lc)
        is_instance_lc = True

    bands = curves.BandNames

    if z is not None and z > 0.1:
        ext = extinction.reddening_law_z(ebv=ebv, bands=bands, z=z, law=law, is_info=is_info)
    else:
        ext = extinction.reddening_law(ebv=ebv, bands=bands, law=law, is_info=is_info)

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

