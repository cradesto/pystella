import numpy as np
from scipy.integrate import quad

__author__ = 'bakl'


#  CGS

class phys:
    h = 6.626068e-27  # erg s
    c = 2.9979245800e10  # cm/s
    k = 1.3806504e-16  # erg K^-1
    sigma_SB = 5.6704e-5  # erg cm^-2 s^-1 K^-4, Stefan-Boltzman Constant
    H0 = 68  # Hubble constant [km/c/Mpc]

    # conversions
    angs_to_cm = 1.e-8
    cm_to_angs = 1. / angs_to_cm
    ZP_AB = -48.6  # zero point AB magnitude for nu
    ZP_AB_lmb = -21.10  # zero point AB magnitude  for lambda

    # units
    AU = 1.4959787066e13  # cm
    pc = 206265 * AU
    R_sun = 6.957e10  # cm
    M_sun = 1.99e33  # g
    L_sun = 3.9e33  # ergs

    d2s = 24.*60.*60.  # convert days to seconds


def cosmology_D_by_z(z):
    Omega_m = 0.31
    Omega_e = 0.69
    c = 2.998e5
    H0 = 67.7
    D = (1. + z) * c / H0 * \
        quad(lambda zz: 1 / np.sqrt(Omega_m * (1. + zz) ** 3 + Omega_e), 0, z)[0]
    return D
