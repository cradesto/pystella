import numpy as np
import unittest
import matplotlib.pyplot as plt

from scipy import interpolate

from pystella.model.popov import Popov
from plugin import sn1999em
from pystella.fit import mpfit

__author__ = 'bakl'


def myfit(mags, lc, is_verbose=True, xtol=1e-10, ftol=1e-10, gtol=1e-10):
    if is_verbose:
        quiet = 0
    else:
        quiet = 1
    t0 = 0.

    mags_interp = interpolate.splrep(lc.Time, mags, s=0)

    def leastsq(p, fjac):
        lc.tshift = p[0]
        m = interpolate.splev(lc.Time, mags_interp)
        return 0, lc.Mag - m

    parinfo = [{'value': t0}]
    result = mpfit.mpfit(leastsq, parinfo=parinfo, quiet=quiet, maxiter=200,
                         ftol=ftol, gtol=gtol, xtol=xtol)
    if result.status == 5:
        print 'Maximum number of iterations exceeded in mangle_spectrum'

    tshift = result.params[0]
    if is_verbose:
        print "The final tshift are: %f" % tshift
    return tshift


def popov_fit(lc, is_verbose=True, xtol=1e-10, ftol=1e-10, gtol=1e-10):
    if is_verbose:
        quiet = 0
    else:
        quiet = 1

    R0 = 100.
    M0 = 10.
    Mni0 = 0.1
    E0 = 1.
    time = lc.Time
    def leastsq(p, fjac):
        ppv = Popov('test', R=p[0], M=p[1], Mni=p[2], E=p[3]*1e51)
        m = ppv.MagBol(time)
        return 0, lc.Mag - m

    parinfo = [{'value': R0, 'limited': [1, 1], 'limits': [10., 150e0]},
               {'value': M0, 'limited': [1, 1], 'limits': [1., 150.]},
               {'value': Mni0, 'limited': [1, 1], 'limits': [0.1, 50.]},
               {'value': E0, 'limited': [1, 1], 'limits': [0.1, 50.]}]
    result = mpfit.mpfit(leastsq, parinfo=parinfo, quiet=quiet, maxiter=200,
                         ftol=ftol, gtol=gtol, xtol=xtol)
    if result.status == 5:
        print 'Maximum number of iterations exceeded in mangle_spectrum'

    tshift = result.params[0]
    ppv = Popov('test', R=result.params[0], M=result.params[1], Mni=result.params[2], E=result.params[3]*1e51)

    if is_verbose:
        print "The final params are: R=%f M=%f Mni=%f E=%f " % (result.params[0],result.params[1],result.params[2],result.params[3])
    return ppv


class TestFit(unittest.TestCase):
    def test_fit_time_popov_SN1999em(self):
        jd_shift = 20.
        dm = -29.38  # D = 7.5e6 pc
        # dm = -30.4  # D = 12.e6 pc
        curves = sn1999em.read_curves()
        lc = curves.get('V')
        lc.mshift = dm

        time = lc.Time - lc.tmin
        # time = np.exp(np.linspace(np.log(start), np.log(end), n))
        popov = Popov('test', R=450., M=15., Mni=0.04, E=0.7e51)
        mags = popov.MagBol(time)


        # fit
        tshift = myfit(mags, lc)

        # plot
        ax = popov.plot_Lbol(time)
        x = lc.Time + tshift
        # x = lc.Time + jd_shift + res
        y = lc.Mag
        ax.plot(x, y, label='%s SN 1999em' % lc.Band.Name,
                ls=".", color='red', markersize=8, marker="o")
        plt.show()

    def test_fit_popov_SN1999em(self):
        ##  todo fit M, E, Mni
        dm = -29.38  # D = 7.5e6 pc
        # dm = -30.4  # D = 12.e6 pc
        curves = sn1999em.read_curves()
        lc = curves.get('V')
        lc.mshift = dm
        lc.tshift = -lc.tmin

        # fit
        ppv = popov_fit(lc)

        # plot
        ax = ppv.plot_Lbol(lc.Time)
        x = lc.Time
        # x = lc.Time + jd_shift + res
        y = lc.Mag
        ax.plot(x, y, label='%s SN 1999em' % lc.Band.Name,
                ls=".", color='red', markersize=8, marker="o")
        plt.show()