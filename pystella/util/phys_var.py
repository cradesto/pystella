__author__ = 'bakl'


#  CGS
class phys:
    h = 6.626068e-27  # erg s
    c = 2.9979245800e10  # cm/s
    k = 1.3806504e-16  # erg K^-1
    sigma_SB = 5.6704e-5  # erg cm^-2 s^-1 K^-4, Stefan-Boltzman Constant
    H0 = 68  # Hubble constant [km/c/Mpc]
    G = 6.6743e-8  # Newton's gravitational constant cm3 g-1 s-2

    echarg = 4.8032042000e-10
    avogar = 6.0221419900e+23

    # conversions
    angs_to_cm = 1.e-8
    cm_to_angs = 1. / angs_to_cm
    ZP_AB = -48.6  # zero point AB magnitude for nu
    ZP_AB_lmb = -21.10  # zero point AB magnitude  for lambda
    jy_to_erg = 1.e-23  # 1 Jy = 10^-23 erg sec^-1 cm^-2 Hz^-1
    jy_to_photon = 1.51e3  # 1 Jy = 1.51e7 photons sec^-1 m^-2 (dlambda/lambda)^-1

    # units
    AU = 1.4959787066e13  # cm
    pc = 206265 * AU
    R_sun = 6.957e10  # cm
    M_sun = 1.99e33  # g
    L_sun = 3.8270e33  # ergs  # see https://sites.google.com/site/mamajeksstarnotes/bc-scale
    Mag_sun = 4.62  # https://ui.adsabs.harvard.edu/abs/1938ApJ....88..429K/abstract
    # Mag_sun = 4.7554

    FOE = 1.e51  # ergs

    d2s = 24. * 60. * 60.  # convert days to seconds
    ev2erg = 1.6021764630e-12  # convert eV to erg

    @staticmethod
    def pc2cm(parsec):
        """Takes in a measurement in parsecs and returns cm"""
        return parsec * phys.pc

    @staticmethod
    def cosmology_D_by_z(*args, **kwargs): # clone
        return cosmology_D_by_z(*args, **kwargs)

    @staticmethod
    def dist2MD(d):  # clone
        return dist2MD(d)


def dist2MD(d):
    import math
    return 5*math.log10(d) - 5.


def cosmology_D_by_z(z, H0=67.7, Omega_m=0.31, Omega_e=0.69):
    """Compute the photometric distance for Lambda-CDM model of cosmology

    Returns
    -------
    D : float
        Distance [Mpc]
    """
    from scipy.integrate import quad
    import numpy as np

    c = 2.998e5
    D = (1. + z) * c / H0 * \
        quad(lambda zz: 1 / np.sqrt(Omega_m * (1. + zz) ** 3 + Omega_e), 0, z)[0]
    return D
