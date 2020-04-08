# origin:  https://github.com/iancze/Pysplotter/blob/master/photometry.py

import numpy as np
from ..util.phys_var import phys as p, phys


def distance_modulus(distance):
    """Given a distance in parsecs, returns the value :math:`m - M`, a characteristic value used to convert
    from apparent magnitude :math:`m` to absolute magnitude, :math:`M`.
    Uses the formula from Carroll and Ostlie, where :math:`d` is in parsecs
    .. math::   m - M = 5 \log_{10}(d) - 5 = 5 \log_{10}(d/10)"""
    return 5 * np.log(distance / 10.0) / np.log(10.)


def distance_from_modulus(mag, abs_mag):
    """Given the distance modulus, return the distance to the source, in parsecs.
    Uses Carroll and Ostlie's formula,
    .. math:: d = 10^{(m - M + 5)/5}"""
    return 10.0 ** (mag - abs_mag + 5) / 5


def val_to_hz(val, inp="Hz"):
    """Takes value in a wavelength in Angstroms and converts it to Hz"""
    if inp == "Hz":
        return val
    elif inp == "A":
        return p.c / (val * 1e-8)
    elif inp == "cm":
        return p.c / val
    else:
        raise ValueError("inp must be 'Hz', 'A','cm'")


def val_to_wl(val, inp="Hz"):
    """Takes value and converts it to wavelength in cm
    :param inp:
    :param val:
    """
    if inp == "Hz":
        return p.c / val
    elif inp == "A":
        return val * 1e-8
    elif inp == "cm":
        return val
    else:
        raise ValueError("inp must be 'Hz', 'A','cm'")


def au_to_cm(dist):
    """Takes in a measurement in AU and returns cm"""
    return dist * p.AU


def planck_partition(nu, temperature):
    """
    Planck partition function
    :param nu: [Hz]
    :param temperature: [K]
    :return: partition
    """
    x = np.exp(-p.h * nu / (p.k * temperature))
    return x / (1. - x)
    # return 1. / (np.exp(p.h*nu / (p.k*temperature)) - 1)


def planck(x, temperature, inp="Hz", out="freq"):
    """General Planck blackbody function.
    Specify `inp` for whether one desires input as wavelength [Angstroms(A), cm] or frequency [Hz].
    Specify `out` for whether one desires output as :math:`B_\lambda` or :math:`B_\nu`."""
    if out == "freq":
        nu = val_to_hz(x, inp)
        return 2 * p.h / p.c ** 2 * nu ** 3 * planck_partition(nu, temperature)
    elif out == 'wave':
        wl = val_to_wl(x, inp)
        return 2 * p.h * p.c ** 2 / wl ** 5 * planck_partition(p.c / wl, temperature)
    else:
        raise ValueError("out must be 'freq' or 'wave' ")


def bb_luminosity_bolometric(temperature, radius):
    """Assume a spherical blackbody at constant temperature.
    Determine the integrated luminosity in ergs/s. All input units cgs"""
    return 4 * np.pi * p.sigma_SB * radius ** 2 * temperature ** 4


def bb_flux_dist(nu, temperature, r, radius):
    """Use to model the flux from an astrophysical source as from a blackbody.
    Takes in frequency in Hz, Temperature in Kelvin, radius and distance in cm.
    Assumes spherical source and F = pi B (R/r)**2 from Rybicki.
    Returns flux"""
    return np.pi * planck(nu, temperature, inp="Hz") * radius ** 2 / r ** 2


def bb_fit_flux_dist(nu, r):
    """Takes in known values and returns a function useful for fitting temperature and radius"""
    return lambda temperature, radius: bb_flux_dist(nu, temperature, r, radius)


def fit_planck(x, inp="Hz", out="freq"):
    return lambda temperature: planck(x, temperature, inp, out)


def compute_x_bb():
    """
    Compute const x_bb for nu weighted with Planck.
    H nu / k T = x_bb --- const
    :return: const
    """
    from scipy.integrate import quad
    f = lambda x: x ** 4 * np.exp(-x) / (1. - np.exp(-x))
    a, err1 = quad(f, 0., np.inf)
    b, err2 = quad(lambda x: x ** 3 * np.exp(-x) / (1. - np.exp(-x)), 0, np.inf)
    return a / b


def temp_wien(wl):
    """
    Find temperature of Spectrum as Wien's law
    :return: Wien's temperature
    """
    b = 0.28977721  # [cm]
    Twien = b / wl
    return Twien


def Lum2MagBol(l):
    """Convert bolometric luminosity to abs. bol. magnitude"""
    mag = -2.5 * np.log10(l / phys.L_sun) + phys.Mag_sun
    return mag


def MagBol2Lum(mag):
    """Convert abs. bol. magnitude to bolometric luminosity"""
    lum = 10. ** (0.4 * (phys.Mag_sun - mag)) * phys.L_sun
    return lum


def Flux2MagAB(f):
    """Convert monochromatic flux [ erg sec^-1 cm^-2 Hz^-1] to AB magnitudes"""
    # ab = 3631e-23  # erg
    mag = -2.5 * np.log10(f) + phys.ZP_AB
    return mag


def Fnu2Fwl(nu, flux):
    """Convert monochromatic flux [ erg sec^-1 cm^-2 Hz^-1] to [ erg sec^-1 cm^-2 A^-1]"""
    return flux * nu ** 2 / phys.c


def Fwl2Fnu(nu, flux_wl):
    """Convert monochromatic flux [ erg sec^-1 cm^-2 A^-1] to [ erg sec^-1 cm^-2 Hz^-1]"""
    return flux_wl * phys.c / nu ** 2


def MagAB2Flux(ab):
    """Convert AB magnitudes to monochromatic flux [ erg sec^-1 cm^-2 Hz^-1]
    see http://www.jstor.org/stable/10.1086/429382
    """
    zp_def = 3631  # in Jy
    f = zp_def * phys.jy_to_erg * 10 ** (-0.4 * ab)
    return f


def kcorrection(series, z, bn_rest, bn_obs, is_verbose=False):
    """
    Compute the K-correction for the set of spectrum.
    See  https://ned.ipac.caltech.edu/level5/Sept02/Hogg/frames.html
    and https://www.osti.gov/biblio/823198

    :param series:  Set of Spectrum
    :param z:  redshift
    :param bn_rest:  rest band
    :param bn_obs:  obs band
    :param is_verbose: if True print info
    :return: iterator for tuple of (t, k)
    """
    from ..rf.band import band_is_exist, band_by_name

    if not band_is_exist(bn_rest):
        raise ValueError('No such rest band: {}'.format(bn_rest))
    if not band_is_exist(bn_obs):
        raise ValueError('No such obs band: {}'.format(bn_obs))

    br = band_by_name(bn_rest)
    bo = band_by_name(bn_obs)

    if is_verbose:
        print("Run k-correction with z={} band_rest={}  band_obs={} ".format(z, bn_rest, bn_obs))

    for i, (t, spec) in enumerate(series):
        k = kcorrection_spec(spec, z, br, bo)
        if is_verbose:
            print(" t={}  kcorr={}".format(t, k))
        yield t, k


def kcorrection_spec_2band(spec, z, br, bo, is_photons=True, is_vega=False):
    """
    Compute the K-correction for the Spectrum
    See  https://ned.ipac.caltech.edu/level5/Sept02/Hogg/frames.html
    and https://www.osti.gov/biblio/823198

    :param spec:  the Spectrum
    :param z:  redshift
    :param br:  rest band
    :param bo:  obs band
    :param is_photons:
    :param is_vega:
    :return:  kcorrection
    """
    from scipy import integrate
    from ..util.math import log_interp1d

    log_interp = log_interp1d(spec.Wl, spec.FluxWl)
    # obs
    flo = log_interp(bo.wl)  # band grid of frequencies
    bo_resp_wl = bo.resp_wl * bo.wl

    # rest
    br_wl = br.wl / (1. + z)
    flr = log_interp(br_wl)
    br_resp_wl = br.resp_wl * br.wl

    if is_photons:
        bo_resp_wl = bo_resp_wl * bo.wl
        br_resp_wl = br_resp_wl * br_wl

    a_obs = integrate.simps(flo * bo_resp_wl, bo.wl)
    a_rest = integrate.simps(flr * br_resp_wl, br.wl)

    if is_vega:
        from pystella import Spectrum
        sp_vg = Spectrum.vega()
        vg_interp = log_interp1d(sp_vg.Wl, sp_vg.FluxWl)
        vg_flo = vg_interp(bo.wl)  # band grid of frequencies
        vg_flr = vg_interp(br.wl)  # band grid of frequencies

        b_obs = integrate.simps(bo_resp_wl * vg_flo / bo.wl ** 2, bo.wl)
        b_rest = integrate.simps(br_resp_wl * vg_flr / br_wl ** 2, br_wl)
    else:
        b_obs = integrate.simps(bo_resp_wl, bo.wl)
        b_rest = integrate.simps(br_resp_wl, br.wl)

    # k = b_rest/b_obs * a_rest/a_obs  # todo check
    k = b_rest / b_obs * a_obs / a_rest  # todo check
    # print(z, b_obs, b_rest)
    k = -2.5 * np.log10(k / (1. + z))
    return k


def kcorrection_spec(spec, z, b, is_photons=True, is_vega=False):
    """
    Compute the K-correction for the Spectrum
    See  https://ned.ipac.caltech.edu/level5/Sept02/Hogg/frames.html
    and https://www.osti.gov/biblio/823198

    :param spec:  the Spectrum
    :param z:  redshift
    :param b:  obs and rest band
    :param is_photons:
    :param is_vega:
    :return:  kcorrection
    """
    from scipy import integrate
    from ..util.math import log_interp1d

    log_interp = log_interp1d(spec.Wl, spec.FluxWl)
    # obs
    flo = log_interp(b.wl)  # band grid of frequencies
    bo_resp_wl = b.resp_wl * b.wl

    # rest
    br_wl = b.wl / (1. + z)
    flr = log_interp(br_wl)
    br_resp_wl = b.resp_wl * b.wl

    if is_photons:
        bo_resp_wl = bo_resp_wl * b.wl
        br_resp_wl = br_resp_wl * br_wl

    a_obs = integrate.simps(flo * bo_resp_wl, b.wl)
    a_rest = integrate.simps(flr * br_resp_wl, b.wl)

    # k = b_rest/b_obs * a_rest/a_obs  # todo check
    k = a_obs / a_rest  # todo check
    # print(z, b_obs, b_rest)
    k = 2.5 * np.log10(1. + z) + 2.5 * np.log10(k)
    return k
