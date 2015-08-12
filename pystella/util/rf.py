# origin:  https://github.com/iancze/Pysplotter/blob/master/photometry.py

import numpy as np
from pystella.util.phys_var import phys as p


def distance_modulus(distance):
    """Given a distance in parsecs, returns the value :math:`m - M`, a characteristic value used to convert
    from apparent magnitude :math:`m` to absolute magnitude, :math:`M`.
    Uses the formula from Carroll and Ostlie, where :math:`d` is in parsecs
    .. math::   m - M = 5 \log_{10}(d) - 5 = 5 \log_{10}(d/10)"""
    return 5 * np.log(distance / 10.0) / np.log(10.)


def distance_from_modulus(m, M):
    """Given the distance modulus, return the distance to the source, in parsecs.
    Uses Carroll and Ostlie's formula,
    .. math:: d = 10^{(m - M + 5)/5}"""
    return 10.0 ** (m - M + 5) / 5


def lamb_to_hz(lamb):
    """Takes in a wavelength in Angstroms and converts it to Hz"""
    return p.c / (lamb * 1e-8)


def hz_to_lamb(Hz):
    """Takes in a frequency in Hz and converts it to wavelength in Angstroms
    :param Hz:
    """
    return 1e8 * p.c / Hz


def pc_to_cm(parsec):
    """Takes in a measurement in parsecs and returns cm"""
    return parsec * p.pc


def au_to_cm(dist):
    """Takes in a measurement in AU and returns cm"""
    return dist * p.AU


def astro_blackbody(nu, T, r, R):
    """Use to model the flux from an astrophysical source as from a blackbody.
    Takes in frequency in Hz, Temperature in Kelvin, radius and distance in cm.
    Assumes spherical sounce and F = pi B (R/r)**2 from Rybicki.
    Returns flux in Jy"""
    Fnu = 2 * np.pi * p.h * nu ** 3 * R ** 2 / (r ** 2 * p.c ** 2 * (np.exp(p.h * nu / (p.k * T)) - 1))
    Jy = 1e23 * Fnu
    return Jy


def astro_blackbody_bolometric_luminosity(T, R):
    """Assume a spherical blackbody at constant temperature.
    Determine the integrated luminosity in ergs/s. All input units cgs"""
    return 4 * np.pi * p.sigma_SB * R ** 2 * T ** 4


def fit_astro_blackbody(nu, r):
    """Takes in known values and returns a function useful for fitting temperature and radius"""
    func = lambda T, R: astro_blackbody(nu, T, r, R)
    return func


def planck_partition(nu, temperature):
    """
    Planck partition function
    :param nu: [Hz]
    :param temperature: [K]
    :return: partition
    """
    return 1. / (np.exp(p.h * nu / (p.k * temperature)) - 1)


def planck(x, T, inp="waveA", out="freq"):
    """General Planck blackbody function.
    Specify `inp` for whether one desires input as wavelength (Angstroms) or frequency (Hz).
    Specify `out` for whether one desires output as :math:`B_\lambda` or :math:`B_\nu`."""
    # TODO: better error checking here on inp and out.
    if out == "freq":
        if inp == "waveA":
            nu = p.c / (x * 1e-8)
        elif inp == "wave":
            nu = p.c / x
        else:
            nu = x
        return ((2 * p.h * nu ** 3) / p.c ** 2) * planck_partition(nu, T)
    else:
        if inp == "waveA":
            lamb = x * 1e-8
        elif inp == "wave":
            lamb = x
        else:
            lamb = p.c / x
        return 2 * p.h * p.c**2 / lamb**5 * planck_partition(p.c/lamb, T)


def fit_planck(x, inp="waveA", out="freq"):
    return lambda T: planck(x, T, inp, out)

