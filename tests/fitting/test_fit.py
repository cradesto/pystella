import unittest
from os.path import dirname, abspath, join

import matplotlib.pyplot as plt
import numpy as np
# from altair.vega.v2 import transform
from scipy import interpolate

from plugin import sn1999em, sn87a, rednova
# from pystella import FitMCMC
# from pystella.fit import mpfit
# from pystella.fit.fit_mpfit import FitMPFit
# from pystella.util.phys_var import phys

import pystella as ps
from pystella.model.popov import Popov

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
    result = ps.mpfit.mpfit(leastsq, parinfo=parinfo, quiet=quiet, maxiter=200,
                            ftol=ftol, gtol=gtol, xtol=xtol)
    if result.status == 5:
        print('Maximum number of iterations exceeded in mangle_spectrum')

    tshift = result.params[0]
    if is_verbose:
        print("The final tshift are: %f" % tshift)
    return tshift


def popov_fit(lc, R0, M0, Mni0=None, E0=None, dt0=None, is_verbose=True, xtol=1e-10, ftol=1e-10, gtol=1e-10):
    if is_verbose:
        quiet = 0
    else:
        quiet = 1

    time = lc.Time

    def leastsq(p, fjac):
        mdl = Popov('test', R=p[0], M=p[1], Mni=p[2], E=p[3])
        l_dt = p[4]
        t = time + l_dt
        m = mdl.MagBol(t)
        res = (lc.Mag - m) ** 2 / m
        w = np.exp(-(max(abs(lc.Mag)) - abs(lc.Mag)) * 2)  # weight
        w = 1.
        # w = w / max(w)
        if lc.IsErr:
            res = res * w / lc.Err
        return 0, res

    parinfo = [{'value': R0, 'limited': [1, 1], 'limits': [10., 1500e0]},
               {'value': M0, 'limited': [1, 1], 'limits': [1., 150.]}
               # {'value': Mni0, 'limited': [1, 1], 'limits': [0.0000001, 0.00001]},
               # , {'value': E0,   'limited': [1, 1], 'limits': [0.01, 5.]}
               # , {'value': dt0,  'limited': [1, 1], 'limits': [-200., 250.]}
               ]
    if Mni0 is not None:
        parinfo.append({'value': Mni0, 'limited': [1, 1], 'limits': [0., 15.]})
    else:
        parinfo.append({'value': 0.01, 'fixed': 1})
    if E0 is not None:
        parinfo.append({'value': E0, 'limited': [1, 1], 'limits': [0.01, 5.]})
    else:
        parinfo.append({'value': 1., 'fixed': 1})
    if dt0 is not None:
        parinfo.append({'value': dt0, 'limited': [1, 1], 'limits': [-200., 250.]})
    else:
        parinfo.append({'value': 0., 'fixed': 1})

    result = ps.mpfit.mpfit(leastsq, parinfo=parinfo, quiet=quiet, maxiter=200,
                            ftol=ftol, gtol=gtol, xtol=xtol)
    if result.status == 5:
        print('Maximum number of iterations exceeded in mangle_spectrum')

    ppv = Popov('test', R=result.params[0], M=result.params[1], Mni=result.params[2], E=result.params[3])
    tshift = result.params[4]

    if is_verbose:
        print("The final params are: R=%f M=%f Mni=%f E=%f; dt = %f " % (
            result.params[0], result.params[1], result.params[2], result.params[3], tshift))
    return ppv, tshift


class TestFit(unittest.TestCase):
    # @unittest.skip("just for plot")
    def test_fit_time_popov_SN1999em(self):
        jd_shift = 20.
        dm = -29.38  # D = 7.5e6 pc
        # dm = -30.4  # D = 12.e6 pc
        curves = sn1999em.read_curves()
        lc = curves.get('V')
        lc.mshift = dm

        time = lc.Time - lc.Tmin
        # time = np.exp(np.linspace(np.log(start), np.log(end), n))
        popov = Popov('test', R=450., M=15., Mni=0.04, E=0.7)
        mags = popov.MagBol(time)

        # fit
        tshift = myfit(mags, lc)

        # plot
        ax = popov.plot_Lbol(time)
        x = lc.Time + tshift
        # x = lc.Time + jd_shift + res
        y = lc.Mag
        ax.plot(x, y, label='%s SN 1999em' % lc.Band.Name,
                ls=":", color='red', markersize=8, marker="o")
        plt.show()

    # @unittest.skip("just for plot")
    def test_fit_popov_SN1999em(self):
        D = 11.5e6  # pc
        dm = -5. * np.log10(D) + 5
        # dm = -30.4  # D = 12.e6 pc
        curves = sn1999em.read_curves()
        lc = curves.get('R')
        lc.mshift = dm
        # lc.tshift = -lc.tmin

        # fit
        ppv, tshift = popov_fit(lc, R0=1200., M0=10., Mni0=0.01, E0=1., dt0=0.)
        # print
        txt = '{0:10} {1:.4e} R_sun\n'.format('R:', ppv.R0 / ps.phys.R_sun) + \
              '{0:10} {1:.4e} M_sun\n'.format('Mtot:', ppv.Mtot / ps.phys.M_sun) + \
              '{0:10} {1:.4e} M_sun\n'.format('Mni:', ppv.Mni / ps.phys.M_sun) + \
              '{0:10} {1} ergs\n'.format('Etot:', ppv.Etot) + \
              '{0:10} {1} d'.format('tshift:', tshift)
        print(txt)
        # plot model
        # time = lc.Time - tshift
        ax = ppv.plot_Lbol(lc.Time)
        # plot obs
        lc.tshift = lc.tshift + tshift
        x = lc.Time
        # x = lc.Time + jd_shift + res
        y = lc.Mag
        ax.plot(x, y, label='%s SN 1999em' % lc.Band.Name,
                ls="", color='red', markersize=8, marker="o")
        plt.show()

    # @unittest.skip("just for plot")
    def test_fit_popov_rednova(self):
        D = 7.5e5  # pc
        dm = 5. * np.log10(D) - 5
        curves = rednova.read_curves_kurtenkov()
        lc = curves.get('R')
        lc.mshift = -dm
        lc.tshift = -lc.Tmin

        # fit
        ppv, tshift = popov_fit(lc, R0=10., M0=3., E0=1.e-2, dt0=-10)
        # ppv, tshift = popov_fit(lc, R0=10., M0=3., Mni0=0.00001, E0=1.e-2, dt0=0.)

        # plot
        ax = ppv.plot_Lbol(lc.Time)
        x = lc.Time + tshift
        # x = lc.Time + jd_shift + res
        y = lc.Mag
        ax.plot(x, y, label='%s M31LRN ' % lc.Band.Name,
                ls="", color='red', markersize=8, marker="o")
        plt.show()

    # @unittest.skip("just for plot")
    def test_fit_popov_SN1987A(self):
        #  todo fit M, E, Mni
        D = 5e4  # pc
        dm = -5. * np.log10(D) + 5  # D = 7.5e6 pc
        # dm = -30.4  # D = 12.e6 pc
        curves = sn87a.read_curves()
        lc = curves.get('R')
        lc.mshift = dm
        lc.tshift = -lc.Tmin

        # fit
        ppv, tshift = popov_fit(lc, R0=1000., M0=10., Mni0=0.01, E0=1., dt0=0.)

        # plot
        ax = ppv.plot_Lbol(lc.Time)
        x = lc.Time + tshift
        # x = lc.Time + jd_shift + res
        y = lc.Mag
        ax.plot(x, y, label='%s SN 1987A' % lc.Band.Name,
                ls="", color='red', markersize=8, marker="o")
        plt.show()

    # @unittest.skip("just for plot")
    def test_fit_lc_Stella_SN1999em(self):
        from pystella.rf import light_curve_func as lcf
        from pystella.rf import light_curve_plot as lcp
        # matplotlib.rcParams['backend'] = "TkAgg"
        # matplotlib.rcParams['backend'] = "Qt4Agg"
        # from matplotlib import pyplot as plt
        # Get observations
        D = 11.5e6  # pc
        dm = -5. * np.log10(D) + 5
        # dm = -30.4  # D = 12.e6 pc
        curves_obs = sn1999em.read_curves()
        curves_obs.set_mshift(dm)

        # Get model
        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(dirname(abspath(__file__))), 'data', 'stella')

        curves_mdl = lcf.curves_compute(name, path, curves_obs.BandNames)

        # fit
        # fitter = FitLcMcmc()
        fitter = ps.FitMPFit()
        fitter.is_debug = True
        lc_o = curves_obs.get('V')
        res = fitter.fit_lc(lc_o, curves_mdl.get(lc_o.Band.Name), dt0=0.)

        # print
        print("tshift:   {:.2f} +- {:.4f} ".format(res.tshift, res.tsigma))
        # txt = '{0:10} {1:.4e} \n'.format('tshift:', res.tshift) + \
        #       '{0:10} {1:.4e} \n'.format('tsigma:', res.tsigma)
        # print(txt)
        # plot model
        curves_obs.set_tshift(res.tshift)
        # curves_mdl.set_tshift(0.)
        ax = lcp.curves_plot(curves_mdl)

        lt = {lc.Band.Name: 'o' for lc in curves_obs}
        lcp.curves_plot(curves_obs, ax, lt=lt, xlim=(-10, 300), is_line=False)
        plt.show()

    # @unittest.skip("just for plot")
    def test_fit_mpfit_curves_Stella_SN1999em(self):
        from pystella.rf import light_curve_func as lcf
        from pystella.rf import light_curve_plot as lcp
        # matplotlib.rcParams['backend'] = "TkAgg"
        # matplotlib.rcParams['backend'] = "Qt4Agg"
        # from matplotlib import pyplot as plt
        # Get observations
        D = 11.5e6  # pc
        dm = -5. * np.log10(D) + 5
        # dm = -30.4  # D = 12.e6 pc
        curves_obs = sn1999em.read_curves()
        curves_obs.set_mshift(dm)

        # Get model
        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(dirname(abspath(__file__))), 'data', 'stella')

        curves_mdl = lcf.curves_compute(name, path, curves_obs.BandNames)

        # fit
        # fitter = FitLcMcmc()
        fitter = ps.FitMPFit()
        fitter.is_debug = True
        res = fitter.fit_curves(curves_obs, curves_mdl)

        # print
        # txt = '{0:10} {1:.4e} \n'.format('tshift:', res.tshift) + \
        #       '{0:10} {1:.4e} \n'.format('tsigma:', res.tsigma) + \
        #       '{0}\n'.format(res.comm)
        txt = '{:10s} {:.4f}+-{:.4f} \n'.format('tshift:', res.tshift, res.tsigma) + \
              '{:10s} {:.4f}+-{:.4f}\n'.format('msigma:', res.mshift, res.msigma) + \
              '{0}\n'.format(res.comm)

        print(txt)
        # plot model
        curves_obs.set_tshift(res.tshift)
        # curves_mdl.set_tshift(0.)
        ax = lcp.curves_plot(curves_mdl)

        lt = {lc.Band.Name: 'o' for lc in curves_obs}
        lcp.curves_plot(curves_obs, ax, lt=lt, xlim=(-10, 300), is_line=False)
        ax.text(0.1, 0.1, txt, transform=ax.transAxes)
        plt.show()

    # @unittest.skip("just for plot")
    def test_FitMPF_best_curves_SN1999em(self):
        from pystella.rf import light_curve_func as lcf
        from pystella.rf import light_curve_plot as lcp
        # Get observations
        D = 11.5e6  # pc
        ebv_sn = 0.
        # ebv_sn = 0.1
        # dm = -5. * np.log10(D) + 5
        # dm = -30.4  # D = 12.e6 pc
        curves_obs = sn1999em.read_curves()

        # Get model
        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(dirname(abspath(__file__))), 'data', 'stella')

        curves_mdl = ps.Stella(name, path=path).curves(curves_obs.BandNames, distance=D, ebv=ebv_sn)

        # fit
        fitter = ps.FitMPFit()
        fitter.is_info = True
        fitter.is_debug = True
        fitter.is_quiet = True
        fit_result, res, dum = fitter.best_curves(curves_mdl, curves_obs, dt0=0., )
        # fit_result, res, dum = fitter.best_curves(curves_mdl, curves_obs, dt0=0., dm0=0.)
        dt, dtsig = res['dt'], res['dtsig']
        dm, dmsig = res['dm'], res['dmsig']
        # plot model

        curves_mdl.set_tshift(dt)
        ax = lcp.curves_plot(curves_mdl)

        lt = {lc.Band.Name: 'o' for lc in curves_obs}
        lcp.curves_plot(curves_obs, ax, lt=lt, xlim=(-10, 300), is_line=False)
        # print
        txt = '{:10} {:.2f} +- {:.4f}'.format('dt:', res['dt'], res['dtsig'])
        txt += '\n{:10} {:.2f} +- {:.4f}'.format('dm:', res['dm'], res['dmsig'])
        print(txt)
        title_font = {'size': '12', 'color': 'black', 'weight': 'normal', 'verticalalignment': 'bottom'}
        ax.text(0.03, 0.95, txt, transform=ax.transAxes,
                bbox={'facecolor': 'none', 'alpha': 0.2, 'edgecolor': 'none'}, **title_font)
        plt.show()

    # @unittest.skip("just for plot")
    def test_FitMPFit_best_curves_gp_SN1999em(self):
        from pystella.rf import light_curve_func as lcf
        from pystella.rf import light_curve_plot as lcp
        # matplotlib.rcParams['backend'] = "TkAgg"
        # matplotlib.rcParams['backend'] = "Qt4Agg"
        # from matplotlib import pyplot as plt
        # Get observations
        D = 11.5e6  # pc
        ebv_sn = 0.1
        dm = -5. * np.log10(D) + 5
        # dm = -30.4  # D = 12.e6 pc
        curves_obs = sn1999em.read_curves()
        curves_obs.set_mshift(dm)

        # Get model
        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(dirname(abspath(__file__))), 'data', 'stella')

        curves_mdl = lcf.curves_compute(name, path, curves_obs.BandNames)

        # fit
        # fitter = FitLcMcmc()
        fitter = ps.FitMPFit()
        fitter.is_info = True
        fitter.is_debug = True
        res = fitter.best_curves_gp(curves_mdl, curves_obs, dt0=0., dm0=0.)
        # print
        txt = '{0:10} {1:.4e} \n'.format('tshift:', res['dt']) + \
              '{0:10} {1:.4e} \n'.format('tsigma:', res['dtsig'])
        print(txt)
        # plot model
        curves_obs.set_tshift(res['dt'])
        # curves_mdl.set_tshift(0.)
        ax = lcp.curves_plot(curves_mdl)

        lt = {lc.Band.Name: 'o' for lc in curves_obs}
        lcp.curves_plot(curves_obs, ax, lt=lt, xlim=(-10, 300), is_line=False)
        plt.show()

    # @unittest.skip("just for plot")
    def test_FitMCMC_best_curves_dtdm_SN1999em(self):
        from pystella.rf import light_curve_func as lcf
        from pystella.rf import light_curve_plot as lcp
        # Get observations
        D = 11.5e6  # pc
        # MD0 = -5. * np.log10(D) + 5
        ebv_sn = 0.1
        # dm = -30.4  # D = 12.e6 pc
        curves_obs = sn1999em.read_curves()
        # curves_obs.set_mshift(MD0)

        # Get model
        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(dirname(abspath(__file__))), 'data', 'stella')

        curves_mdl = ps.Stella(name, path=path).curves(curves_obs.BandNames, distance=D, ebv=ebv_sn)

        # fit
        is_debug = True  # True
        # fitter = FitLcMcmc()
        fitter = ps.FitMCMC()
        fitter.is_info = True
        if is_debug:
            fitter.is_debug = is_debug
            fitter.nwalkers = 100  # number of MCMC walkers
            fitter.nburn = 20  # "burn-in" period to let chains stabilize
            fitter.nsteps = 200  # number of MCMC steps to take
        threads = 3

        # fitter.is_debug = False
        # fit_result, res, samples \
        fit_result, res, (th, ep, em), sampler = fitter.best_curves(curves_mdl, curves_obs, dt0=0., dm0=0.,
                                                                    threads=threads, is_sampler=True)
        nb = curves_obs.Length
        dt, dm, sigs = fitter.theta2arr(th, nb)
        ep_dt, ep_dm, ep_sigs = fitter.theta2arr(ep, nb)
        em_dt, em_dm, em_sigs = fitter.theta2arr(em, nb)

        samples = sampler.flatchain[fitter.nburn:, :]
        fig = fitter.plot_corner(samples, bnames=curves_obs.BandNames, bins=55, alpha=0.5, verbose=fitter.is_info)

        # print
        txt = 'chi2= {:.1f}  BIC= {:.1f} AIC= {:.1f} measure= {:.3f}\n'.\
            format(res['chi2'], res['bic'], res['aic'], res['measure'])
        txt += '\n dt= {:.1f} ^{{{:.1f}}}_{{{:.1f}}} '.format(dt, ep_dt, em_dt)
        txt += '\n dm= {:.2f} ^{{{:.2f}}}_{{{:.2f}}} '.format(dm, ep_dm, em_dm)
        for i, bname in enumerate(curves_obs.BandNames):
            txt += '\n sigma_{{{:s}}} = {:.2f}^{{{:.2f}}}_{{{:.2f}}} '. \
                format(bname, sigs[i], ep_sigs[i], em_sigs[i])
        print(txt)

        curves_mdl.set_tshift(dt)
        curves_mdl.set_mshift(dm)

        ax = None
        errs = {}
        for i, bn in enumerate(curves_mdl.BandNames):
            errs[bn] = [sigs[i]] * curves_mdl[bn].Length  # The error are the same for all points
        curves_clone = curves_mdl.clone(err=errs)
        ax = ps.curves_plot(curves_clone, ax=ax, is_legend=False, is_fill=True, alpha=0.2)
        # Obs
        lt = {lc.Band.Name: 'o' for lc in curves_obs}
        lcp.curves_plot(curves_obs, ax, lt=lt, xlim=(-10, 300), is_line=False)

        ax.text(0.05, 0.05, txt, transform=ax.transAxes)
        plt.show()

    # @unittest.skip("just for plot")
    def test_FitMCMC_best_curves_dt_SN1999em(self):
        from pystella.rf import light_curve_func as lcf
        from pystella.rf import light_curve_plot as lcp
        # Get observations
        D = 11.5e6  # pc
        dm = -5. * np.log10(D) + 5
        # dm = -30.4  # D = 12.e6 pc
        curves_obs = sn1999em.read_curves()
        curves_obs.set_mshift(dm)

        # Get model
        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(dirname(abspath(__file__))), 'data', 'stella')

        curves_mdl = lcf.curves_compute(name, path, distance=10, bands=curves_obs.BandNames)

        # fit
        is_debug = True  # True
        # fitter = FitLcMcmc()
        fitter = ps.FitMCMC()
        fitter.is_info = True
        if is_debug:
            fitter.is_debug = is_debug
            fitter.nwalkers = 100  # number of MCMC walkers
            fitter.nburn = 20  # "burn-in" period to let chains stabilize
            fitter.nsteps = 200  # number of MCMC steps to take
        threads = 1

        # fitter.is_debug = False
        fit_result, res, (th, e1, e2), samples = fitter.best_curves(curves_mdl, curves_obs, dt0=0., threads=threads,
                                                                    is_sampler=True)
        fig = fitter.plot_corner(samples, labels=('dt',), bnames=curves_obs.BandNames)

        # print
        # '{:10s} {:.4f} ^{:.4f}_{:.4f}\n'.format('lnf:', res['lnf'], res['lnfsig2'], res['lnfsig1']) + \
        txt = '{:10s} {:.4f} ^{:.4f}_{:.4f} \n'.format('tshift:', th.dt, e1.dt, e2.dt) + \
              '{:10s} chi2= {:.1f}  BIC= {:.1f} AIC= {:.1f} dof= {} accept= {:.3f}\n'. \
                  format('stat:', res['chi2'], res['bic'], res['aic'], res['dof'], res['acceptance_fraction'])
        print(txt)
        # plot model
        curves_obs.set_tshift(-th.dt)
        ax = lcp.curves_plot(curves_mdl)

        lt = {lc.Band.Name: 'o' for lc in curves_obs}
        lcp.curves_plot(curves_obs, ax, lt=lt, xlim=(-10, 300), is_line=False)
        ax.text(0.1, 0.1, txt, transform=ax.transAxes)
        plt.show()

    # @unittest.skip("just for plot")
    def test_FitMCMC_best_curves_dt_sigmas_SN1999em(self):
        from pystella.rf import light_curve_func as lcf
        from pystella.rf import light_curve_plot as lcp
        # Get observations
        D = 11.5e6  # pc
        dm = -5. * np.log10(D) + 5
        # dm = -30.4  # D = 12.e6 pc
        curves_obs = sn1999em.read_curves()
        curves_obs.set_mshift(dm)

        # Get model
        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(dirname(abspath(__file__))), 'data', 'stella')

        curves_mdl = lcf.curves_compute(name, path, distance=10, bands=curves_obs.BandNames)

        # fit
        is_debug = True  # True
        # fitter = FitLcMcmc()
        fitter = ps.FitMCMC()
        fitter.is_info = True
        if is_debug:
            fitter.is_debug = is_debug
            fitter.nwalkers = 100  # number of MCMC walkers
            fitter.nburn = 20  # "burn-in" period to let chains stabilize
            fitter.nsteps = 200  # number of MCMC steps to take
        threads = 1

        # fitter.is_debug = False
        bnames = curves_obs.BandNames
        res, samples = fitter.best_curves_sigmas(curves_mdl, curves_obs, bnames=bnames,
                                                 dt0=0., threads=threads, is_sampler=True)
        fig = fitter.plot_corner(samples, bnames=bnames)

        # print
        txt = ''
        txt += r'{:10s} {:.4f} ^{{{:.4f}}}_{{{:.4f}}} \n'.format('tshift:', res['dt'], res['dtsig2'], res['dtsig1'])
        for i, bname in enumerate(bnames):
            txt += r'sig+{:7s}: {:.4f}^{{{:.4f}}}_{{{:.4f}}}\n'. \
                format(bname, res['lnf'][i], res['lnfsig2'][i], res['lnfsig1'][i])

        txt += '{:10s} chi2= {:.1f}  BIC= {:.1f} AIC= {:.1f} dof= {} accept= {:.3f}\n'. \
            format('stat:', res['chi2'], res['bic'], res['aic'], res['dof'], res['acceptance_fraction'])
        print(txt)

        # plot model
        curves_obs.set_tshift(res['dt'])
        ax = lcp.curves_plot(curves_mdl)

        lt = {lc.Band.Name: 'o' for lc in curves_obs}
        lcp.curves_plot(curves_obs, ax, lt=lt, xlim=(-10, 300), is_line=False)
        ax.text(0.1, 0.02, txt, transform=ax.transAxes)
        plt.show()

    # @unittest.skip("just for plot")
    def test_FitMCMC_best_curves_dtdm_sigmas_SN1999em(self):
        from pystella.rf import light_curve_func as lcf
        from pystella.rf import light_curve_plot as lcp
        # Get observations
        D = 11.5e6  # pc
        dm = -5. * np.log10(D) + 5
        # dm = -30.4  # D = 12.e6 pc
        curves_obs = sn1999em.read_curves()
        curves_obs.set_mshift(dm)

        # Get model
        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(dirname(abspath(__file__))), 'data', 'stella')

        curves_mdl = lcf.curves_compute(name, path, distance=10, bands=curves_obs.BandNames)

        # fit
        is_debug = True  # True
        # fitter = FitLcMcmc()
        fitter = ps.FitMCMC()
        fitter.is_info = True
        if is_debug:
            fitter.is_debug = is_debug
            fitter.nwalkers = 100  # number of MCMC walkers
            fitter.nburn = 20  # "burn-in" period to let chains stabilize
            fitter.nsteps = 200  # number of MCMC steps to take
        threads = 1

        # fitter.is_debug = False
        bnames = curves_obs.BandNames
        res, samples = fitter.best_curves_sigmas(curves_mdl, curves_obs, bnames=bnames,
                                                 dt0=0., dm0=0., threads=threads, is_sampler=True)
        fig = fitter.plot_corner(samples, bnames=bnames)

        # print
        txt = ''
        txt += '{:10s} {:.4f} ^{:.4f}_{:.4f} \n'.format('tshift:', res['dt'], res['dtsig2'], res['dtsig1'])
        txt += '{:10s} {:.4f} ^{:.4f}_{:.4f}\n'.format('msigma:', res['dm'], res['dmsig2'], res['dmsig1'])
        for i, bname in enumerate(bnames):
            txt += 'sig+{:7s}: {:.4f} ^{:.4f}_{:.4f}\n'. \
                format(bname, res['lnf'][i], res['lnfsig2'][i], res['lnfsig1'][i])

        txt += '{:10s} chi2= {:.1f}  BIC= {:.1f} AIC= {:.1f} dof= {} accept= {:.3f}\n'. \
            format('stat:', res['chi2'], res['bic'], res['aic'], res['dof'], res['acceptance_fraction'])
        print(txt)
        # plot model
        curves_obs.set_tshift(res['dt'])
        curves_obs.set_mshift(dm + res['dm'])
        # curves_mdl.set_tshift(0.)
        ax = lcp.curves_plot(curves_mdl)

        lt = {lc.Band.Name: 'o' for lc in curves_obs}
        lcp.curves_plot(curves_obs, ax, lt=lt, xlim=(-10, 300), is_line=False)
        ax.text(0.1, 0.02, txt, transform=ax.transAxes)
        plt.show()

    # @unittest.skip("just for plot")
    def test_FitMCMC_fake_curves(self):
        # Get observations
        # D = 11.5e6  # pc
        D = 10.5e6  # pc
        dm = -5. * np.log10(D) + 5
        curves_o = sn1999em.read_curves()
        curves_o.set_mshift(dm)

        # Model
        fname = 'cat_R500_M15_Ni006_E12.tt.ubv'
        path = join(dirname(dirname(abspath(__file__))), 'data', 'stella')
        # curves_m = pd.read_csv(join(path, fname), header=0, delim_whitespace=True)
        curves_m = ps.curves_read(join(path, fname))

        # Filter
        bnames = curves_m.BandNames
        for bname in bnames:
            if bname not in curves_o.BandNames:
                curves_m.pop(bname)

        # curves_m = curves_m[curves_m.TimeCommon > 0]

        # fit
        is_debug = True
        fitter = ps.FitMCMC()
        fitter.is_info = True  # True False
        fitter.is_debug = is_debug
        if is_debug:
            fitter.nwalkers = 100  # number of MCMC walkers
            fitter.nburn = 20  # "burn-in" period to let chains stabilize
            fitter.nsteps = 200  # number of MCMC steps to take
        else:
            fitter.nwalkers = 200  # number of MCMC walkers
            fitter.nburn = 100  # "burn-in" period to let chains stabilize
            fitter.nsteps = 500  # number of MCMC steps to take

        print('is_debug ', is_debug)

        res, samples = fitter.best_curves(curves_m, curves_o, dt0=0., dm0=0., is_sampler=True)

        curves_o.set_tshift(0.)
        curves_o.set_mshift(dm)

        curves_m.set_tshift(-res['dt'])
        curves_m.set_mshift(-res['dm'])

        # Plot
        fig, (ax, axhist) = plt.subplots(2, 1, figsize=(9, 9))
        fig = fitter.plot_corner(samples, axhist=axhist)

        ps.curves_plot(curves_m, ax=ax)
        ps.curves_plot(curves_o, ax=ax, is_line=False, markersize=2)

        # print
        txt = '{:10s} {:.4f} ^{:.4f}_{:.4f} \n'.format('tshift:', res['dt'], res['dtsig2'], res['dtsig1']) + \
              '{:10s} {:.4f} ^{:.4f}_{:.4f}\n'.format('msigma:', res['dm'], res['dmsig2'], res['dmsig1']) + \
              '{:10s} {:.4f} ^{:.4f}_{:.4f}\n'.format('lnf:', res['lnf'], res['lnfsig2'], res['lnfsig1']) + \
              '{:10s} chi2= {:.4f} dof= {} accept= {:.3f}\n'. \
                  format('stat:', res['chi2'], res['dof'], res['acceptance_fraction'])
        print(txt)
        ax.text(0.1, 0.02, txt, transform=ax.transAxes)
        xlim = ax.get_xlim()
        xlim = -30, xlim[1]
        ax.set_xlim(xlim)
        fig.tight_layout()
        plt.show()
