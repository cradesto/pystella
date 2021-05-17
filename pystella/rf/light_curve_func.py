import csv
import os

import numpy as np

from .rad_func import Fnu2Fwl, Fwl2Fnu
from pystella.util.phys_var import phys
from pystella.model.stella import Stella
from . import extinction
from .lc import LightCurve, SetLightCurve
from .reddening import ReddeningLaw
from .reddening import LawFitz

__author__ = 'bakl'

# __all__ = ['curves_read', 'curves_read_mix']


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
    mags = serial_spec.mags_bands(bands, z=z, d=phys.pc2cm(distance), magnification=magnification)

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


def curves_save(curves, fname, sep='\t', is_mix=False):
    """
    Save curves to CSV-format. For band-column format it's required the  common time for all LC.
    :param curves: saved curves
    :param fname: output file
    :param sep: string or character separating columns.
    :param is_mix: format flag. If True, save rows format
    :return: True if Success, otherwise False
    """
    if curves.Length == 0:
        print("Nothing to save: curves.Length=%d" % curves.Length)
        return False

    if curves.IsCommonTime and not is_mix:
        arr = curves2nparray(curves)
        fmt_header = "%10s  " * len(arr.dtype.names)
        header = fmt_header % arr.dtype.names
        fmt = "%10.3e  " # time
        fmt += "%10.4f  " * (len(arr.dtype.names)-1)  # mags
    else:
        print("The curves has different time arrays.")
        print("  LC  |  len(time)  ")
        for lc in curves:
            print("{:40s}   {:d}".format(lc.Name, len(lc.Time)))
        arr = curves2nparraymix(curves)
        fmt_header = "%8s " * len(arr.dtype.names)
        header = fmt_header % arr.dtype.names
        fmt = "%10.5f %10s %10.5f"
        if len(arr.dtype.names) > 3:
            fmt += ' %10.5f'

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


def curves_read_mix(fname, dtype=None,
                    skiprows=1, comments='#', is_full=False, is_out=False):
    """
    Reader data-file with mix bname data, like:
    >>% head photometry.txt
    jd filter mag mage
    2457059.6228778586 V 17.493766309999998 0.0592200135089
    2457059.6244578934 V 17.539956019999998
    0.0542402986717 2457059.6261980557 g 17.782871193345898
    0.0454000142503 2457059.6287036575 g 17.7782469177482 0.0395424488201

    :param is_out:
    :param fname: file with data
    :param dtype: You should define the colums: time, mag, err and filter,
                  example:  (('time', '<f4'), ('filter', 'S1'), ('mag', '<f4'), ('err', '<f4'))
    :param skiprows: skip rows, default: 1 for header
    :param comments: skip comments, default: #
    :param is_full: return also table data, default: False
    :param is_out: print first line of data
    :return: curves
    """
    from pystella.rf import band
    band.Band.load_settings()

    if dtype is None:
        dtype = [('time', np.float), ('filter', 'S8'), ('mag', np.float), ('err', np.float)]

    lc_data = np.loadtxt(fname, skiprows=skiprows, dtype=dtype, comments=comments)  # jd filter mag mage
    if is_out:
        print(lc_data.dtype)
        print(lc_data[0])
        print(lc_data[1])

    b_tot = lc_data['filter']
    bnames = np.unique(b_tot)

    curves = SetLightCurve()
    for bname in bnames:
        bname_s = bname
        try:
            bname_s = bname_s.decode()
        except AttributeError:
            pass
        if band.is_exist(bname_s):
            # filter of the current band
            is_good = np.array(list(map(lambda x: x == bname, b_tot)), dtype=bool)
            t = lc_data['time'][is_good]
            m = lc_data['mag'][is_good]
            e = lc_data['err'][is_good]
            # add light curve
            b = band.band_by_name(bname_s)
            lc = LightCurve(b, t, m, e)
            curves.add(lc)
        else:
            print('Could read the light curve. There is no band: {}. '
                  'You may try to add it to dir data/bands'.format(bname))
    if is_full:
        return curves, lc_data

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


def curves2nparraymix(curves):
    """
    Convert curves to numpy array with mix format.
    >>% head photometry.txt
    jd filter mag mage
    2457059.6228778586 V 17.493766309999998 0.0592200135089
    :param curves:
    :return:
    """
    if curves.Length == 0:
        print("Nothing to convert: curves.Length=%d" % curves.Length)
        return None

    is_err = True
    for lc in curves:
        if not lc.IsErr:
            is_err = False
            break
    dtype = [('time', '<f4'), ('b', 'U8'), ('mag', '<f4')]
    if is_err:
        dtype.append(('err', '<f4'))

    rows = np.sum([lc.Length for lc in curves])
    data = np.zeros(rows, dtype=dtype)
    i = 0
    for lc in curves:
        b = lc.Band.Name
        data['time'][i:i+lc.Length] = lc.Time
        data['b'][i:i+lc.Length] = [b] * lc.Length
        data['mag'][i:i+lc.Length] = lc.Mag
        if is_err:
            data['err'][i:i + lc.Length] = lc.MagErr
        i += lc.Length

    # res = np.column_stack(data)
    # res = np.dstack(data)
    # dtype['formats'] = [np.float64] * len(dtype['names'])
    # res.dtype = dtype
    return data
    # return res.ravel()


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
    curves = serial_spec.flux_to_curves(bands, z=z, d=phys.pc2cm(distance), magnification=magnification)

    if is_show_info:
        # print the time of maximum LC
        tmin = 2.0
        t = curves.TimeCommon
        for n in bands:
            t_min = t[t > tmin][curves[n][t > tmin].argmin()]
            print("t_max(%s) = %f" % (n, t_min))

    return curves


def flux_reddening(freq, flux, ebv, Rv=None, law=LawFitz, mode=ReddeningLaw.MW):
    """
    Apply extinction curves to flux(freq) values
    :param freq:  [Hz]
    :param flux:  [ergs s^-1 cm^-2 Hz^-1]
    :param ebv: E(B-V)
    :param Rv: R_V
    :param law: the variant of extinction curves
    :param mode: type of extinction curves (MW, LMC, SMC)
    :return: reddening flux
    """
    from pystella.util.phys_var import phys

    flux_wl = Fnu2Fwl(freq, flux)
    # flux_wl = flux * freq ** 2 / phys.c
    wl = np.array(phys.c / freq) * phys.cm_to_angs
    res = flux_reddening_wl(wl, flux_wl, ebv, Rv=Rv, law=law, mode=mode)
    res = Fwl2Fnu(freq, res)
    # res = flux_reddening_wl(wl, flux_wl, ebv, law) * phys.c / freq**2
    return res


def flux_reddening_wl(wl, flux_wl, ebv, Rv=None, law=LawFitz, mode=ReddeningLaw.MW):
    """
    Apply extinction curves to flux(lambda) values
    :param wl:  [A]
    :param flux_wl:  [ergs s^-1 cm^-2 A^-1]
    :param ebv: E(B-V)
    :param Rv: R_V
    :param law: the variant of extinction curves
    :param mode: type of extinction curves (MW, LMC, SMC)
    :return: reddening flux
    """
    if Rv is None:
        Rv = law.Rv[mode]
    A_lambda = law.Almd(wl, ebv, Rv=Rv)
    res = flux_wl * 10 ** (-0.4 * A_lambda)
    return res


def series_spec_reddening(series, ebv, Rv=None, law=LawFitz, mode=ReddeningLaw.MW):
    from pystella.rf.spectrum import Spectrum, SeriesSpectrum
    res = SeriesSpectrum(series.Name)
    for k, t in enumerate(series.Time):
        s = series.get_spec(k)
        freq = s.Freq
        flux_red = flux_reddening(freq, s.Flux, ebv, Rv=Rv, law=law, mode=mode)
        ss = Spectrum(s.Name, freq, flux=flux_red)
        res.add(t, ss)
    return res


# todo Implement direct reddening from spectra
# http://webast.ast.obs-mip.fr/hyperz/hyperz_manual1/node10.html
def curves_reddening(curves, ebv, z=None, law=extinction.law_default, is_info=True):
    # is_unred = ebv < 0.
    # if ebv < 0.:
    #     ebv = abs(ebv)
    #     # raise ValueError("ebv should be > 0")
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
        # if is_unred:
        #     lc_red.M -= ext[lc.Band.Name]
        # else:
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
