#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import os
from os.path import dirname
from scipy import interpolate

from pystella.model.stella import Stella
from pystella.rf.ts import TimeSeries, SetTimeSeries

__author__ = 'bakl'

ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))


class SetVelocityCurve(SetTimeSeries):
    """Set of the Velocity Curves"""
    def __init__(self, name=''):
        """Creates a Set of Light Curves."""
        super().__init__(name)
        self._loop = 0

    @classmethod
    def Merge(cls, vels1, vels2):
        if vels1 is None:
            return vels2
        if vels2 is None:
            return vels1

        res = SetVelocityCurve("{}+{}".format(vels1.Name, vels2.Name))
        for vel1 in vels1:
            vel2 = vels2.get(vel1.Name)
            if vel2 is None:
                res.add(vel1)
            else:
                vel = VelocityCurve.Merge(vel1, vel2)
                res.add(vel)
        for vel in vels2:
            if not res.IsName(vel.Name):
                res.add(vel)

        return res


class VelocityCurve(TimeSeries):
    def __init__(self, name, time, vel, errs=None, tshift=0., vshift=0.):
        """Creates a Velocity Time Series instance.  Required parameters:  name, time, vel."""
        super().__init__(name, time, vel, errs, tshift=tshift)

        self._vshift = vshift

    @property
    def Vel(self):
        return self.V * self.vshift

    @property
    def vshift(self):
        return self._vshift

    @vshift.setter
    def vshift(self, shift):
        self._vshift = shift

    def copy(self, name=None, f=None):
        errs = None
        if name is None:
            name = self.Name

        if f is not None:
            is_good = np.where(f(self))
            # is_good = np.where((self.Time >= tlim[0]) & (self.Time <= tlim[1]))
            t = self.T[is_good]
            v = self.V[is_good]
            if self.IsErr:
                errs = self.Err[is_good]
        else:
            t = self.T
            v = self.V
            if self.IsErr:
                errs = self.Err

        new = VelocityCurve(name, t, v, errs)
        new.tshift = self.tshift
        new.vshift = self.vshift
        return new


def plot_vels_sn87a(ax, z=0):
    print("Plot the velocities of Sn 87A ")
    d = os.path.expanduser('~/Sn/Release/svn_kepler/stella/branches/lucy/run/res/sncurve/sn1987a')

    jd_shift = 2446850  # moment of explosion SN 1987A, Hamuy 1988, doi:10.1086/114613

    # Blanco's data from plot
    fs = {'Halpha': os.path.join(d, 'Halpha_blanco.csv'), 'Hbeta': os.path.join(d, 'Hbeta_blanco.csv'),
          'Hgamma': os.path.join(d, 'Hgamma_blanco.csv'), 'NaID': os.path.join(d, 'NaID_blanco.csv'),
          'FeII5018': os.path.join(d, 'FeII5018_blanco.csv'), 'FeII5169': os.path.join(d, 'FeII5169_blanco.csv')
          }
    elcolors = {'Halpha': "black", 'Hbeta': "cyan", 'Hgamma': "orange", 'NaID': "red", 'FeII5018': "orange",
                'FeII5169': "magenta"}
    elmarkers = {'Halpha': u's', 'Hbeta': u'x', 'Hgamma': u'd', 'NaID': u'+', 'FeII5018': u'D', 'FeII5169': u'o'}

    for el, fname in fs.items():
        data = np.loadtxt(fname, comments='#')
        x = data[:, 0] - jd_shift
        x *= 1. + z  # redshift
        y = data[:, 1]
        ax.plot(x, y, label='%s, SN 87A' % el, ls=".", color=elcolors[el], markersize=6, marker=elmarkers[el])


def plot_vels_models(ax, models_dic, xlim=None, ylim=None):
    is_x_lim = xlim is None
    is_y_lim = ylim is None

    # t_points = [0.2, 1, 2, 3, 4, 5, 10, 20, 40, 80, 150]

    lw = 1.
    mi = 0
    x_max = []
    y_mid = []
    for mname, mdic in models_dic.items():
        mi += 1
        x = mdic['time']
        y = mdic['vel'] / 1e8
        ax.plot(x, y, label='Vel  %s' % mname, color='blue', ls="-", linewidth=lw)
        if is_x_lim:
            x_max.append(np.max(x))
        if is_y_lim:
            y_mid.append(np.max(y))

    if is_x_lim:
        xlim = [-10, np.max(x_max) + 10.]
    ax.set_xlim(xlim)

    if is_y_lim:
        ylim = [1e-1, np.max(y_mid) + 5]
        # ylim = [np.min(y_mid) + 7., np.min(y_mid) - 2.]
    ax.set_ylim(ylim)

    ax.set_ylabel('Velocity')
    ax.set_xlabel('Time [days]')


def plot_vel(ax, vel, xlim=None, ylim=None):
    is_x_lim = xlim is None
    is_y_lim = ylim is None

    lw = 1.
    x_max = []
    y_mid = []
    x = vel['time']
    y = vel['vel'] / 1e8
    ax.plot(x, y, label='Velocity', color='blue', ls="-", linewidth=lw)
    if is_x_lim:
        x_max.append(np.max(x))
    if is_y_lim:
        y_mid.append(np.max(y))

    if is_x_lim:
        xlim = [-10, np.max(x_max) + 10.]
    ax.set_xlim(xlim)

    if is_y_lim:
        ylim = [1e-1, np.max(y_mid) + 1]
        # ylim = [np.min(y_mid) + 7., np.min(y_mid) - 2.]
    ax.set_ylim(ylim)

    ax.set_ylabel('Velocity')
    ax.set_xlabel('Time [days]')
    ax.grid()


def compute_vel_swd(name, path):
    model = Stella(name, path=path)
    # check data
    if not model.is_swd_data:
        raise ValueError("There are no swd-file for %s in the directory: %s " % (name, path))

    swd = model.get_swd().load()
    data = swd.params_ph(cols=['V'])

    res = np.array(np.zeros(len(data['V'])),
                   dtype=np.dtype({'names': ['time', 'vel'], 'formats': [np.float] * 2}))
    res['time'] = data['time']
    res['vel'] = data['V']

    return res


def compute_vel_res_tt(name, path, z=0., t_beg=1., t_end=None, t_diff=1.05):
    model = Stella(name, path=path)
    # check data
    if not model.is_res_data:
        raise ValueError("There are no res-file for %s in the directory: %s " % (name, path))
    if not model.is_tt_data:
        raise ValueError(("There are no tt-file for %s in the directory: %s " % (name, path)))

    if t_end is None:
        t_end = float('inf')

    res = model.get_res()
    tt = model.get_tt().read()
    tt = tt[tt['time'] >= t_beg]  # time cut  days

    radiuses = list()
    vels = list()
    times = list()
    Rph_spline = interpolate.splrep(tt['time'], tt['Rph'], s=0)
    for nt in range(len(tt['time'])):
        t = tt['time'][nt]
        if t > t_end:
            break
        if t < t_beg or np.abs(t / t_beg < t_diff):
            continue
        t_beg = t
        radius = interpolate.splev(t, Rph_spline)
        if np.isnan(radius):
            radius = np.interp(t, tt['time'], tt['Rph'], 0, 0)  # One-dimensional linear interpolation.
        block = res.read_at_time(time=t)
        if block is None:
            break

        if True:
            vel = np.interp(radius, block['R14']*1e14, block['V8'], 0, 0)  # One-dimensional linear interpolation.
            vels.append(vel * 1e8)
        else:
            idx = np.abs(block['R14'] - radius / 1e14).argmin()
            vels.append(block['V8'][idx] * 1e8)

        radiuses.append(radius)
        times.append(t * (1. + z))  # redshifted time

    # show results
    res = np.array(np.zeros(len(vels)),
                   dtype=np.dtype({'names': ['time', 'vel', 'r'],
                                   'formats': [np.float] * 3}))
    res['time'] = times
    res['vel'] = vels
    res['r'] = radiuses
    return res
