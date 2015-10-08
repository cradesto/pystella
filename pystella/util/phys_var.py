__author__ = 'bakl'


#  CGS

class phys:
    h = 6.626068e-27  # erg s
    c = 2.9979245800e10  # cm/s
    k = 1.3806504e-16  # erg K^-1
    sigma_SB = 5.6704e-5  # erg cm^-2 s^-1 K^-4, Stefan-Boltzman Constant

    # conversions
    angs_to_cm = 1.e-8
    cm_to_angs = 1. / angs_to_cm
    ZP_AB = -48.6  # zero point AB magnitude for nu
    ZP_AB_lmb = -21.10  # zero point AB magnitude  for lambda

    # units
    AU = 1.4959787066e13  # cm
    pc = 206265 * AU
    M_sun = 1.99e33  # g
    L_sun = 3.9e33  # ergs

