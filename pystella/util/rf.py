# origin:  https://github.com/iancze/Pysplotter/blob/master/photometry.py

import numpy as np
from pystella.util.phys_var import phys as p


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


def pc_to_cm(parsec):
    """Takes in a measurement in parsecs and returns cm"""
    return parsec * p.pc


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
    x = np.exp(-p.h*nu / (p.k*temperature))
    return x/(1.-x)
    # return 1. / (np.exp(p.h*nu / (p.k*temperature)) - 1)


def planck(x, temperature, inp="Hz", out="freq"):
    """General Planck blackbody function.
    Specify `inp` for whether one desires input as wavelength [Angstroms(A), cm] or frequency [Hz].
    Specify `out` for whether one desires output as :math:`B_\lambda` or :math:`B_\nu`."""
    if out == "freq":
        nu = val_to_hz(x, inp)
        return 2*p.h/p.c**2 * nu**3 * planck_partition(nu, temperature)
    elif out == 'wave':
        wl = val_to_wl(x, inp)
        return 2 * p.h * p.c**2 / wl**5 * planck_partition(p.c/wl, temperature)
    else:
        raise ValueError("out must be 'freq' or 'wave' ")


def bb_luminosity_bolometric(temperature, radius):
    """Assume a spherical blackbody at constant temperature.
    Determine the integrated luminosity in ergs/s. All input units cgs"""
    return 4*np.pi * p.sigma_SB * radius**2 * temperature**4


def bb_flux_dist(nu, temperature, r, radius):
    """Use to model the flux from an astrophysical source as from a blackbody.
    Takes in frequency in Hz, Temperature in Kelvin, radius and distance in cm.
    Assumes spherical source and F = pi B (R/r)**2 from Rybicki.
    Returns flux"""
    return np.pi * planck(nu, temperature, inp="Hz") * radius**2 / r**2


def bb_fit_flux_dist(nu, r):
    """Takes in known values and returns a function useful for fitting temperature and radius"""
    func = lambda temperature, radius: bb_flux_dist(nu, temperature, r, radius)
    return func


def fit_planck(x, inp="Hz", out="freq"):
    return lambda temperature: planck(x, temperature, inp, out)
