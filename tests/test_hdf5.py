import unittest
from os.path import join, dirname, abspath
import numpy as np
import h5py
import pystella as ps

import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.cm as cm
except ImportError:
    mpl_logger.warn("You have not got the module matplotlib. Please, install it.")
    plt = None


__author__ = 'bakl'


class TestHDF5(unittest.TestCase):
    def test_hdf5_simple(self):
        f = h5py.File("tmp.hdf5", "w")
        dset = f.create_dataset("mydataset", (100,), dtype='i')
        grp = f.create_group("subgroup")
        dset2 = grp.create_dataset("another_dataset", (50,), dtype='f')
        dset3 = f.create_dataset('subgroup2/dataset_three', (10,), dtype='i')
        dset.attrs['temperature'] = 99.5

        f.close()

    def test_curves_save(self):
        name = 'cat_R500_M15_Ni006_E12'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        bands = ('UVW1', 'UVW2', 'UVM2')
        # bands = ('UVW1', 'UVW2', 'UVM2', 'U', 'B', 'R', 'I')
        mdl = ps.Stella(name, path=path)
        curves = mdl.curves(bands)

        with h5py.File("tmp_curves.hdf5", "w") as f:
            grp = f.create_group("curves")
            for lc in curves:
                ds = grp.create_dataset(lc.Band.Name, data=lc.toarray())
                ds.attrs['z'] = lc.attrs('z')
                ds.attrs['d'] = lc.attrs('d')

            self.assertEqual(curves.Length, len(grp.items()))

    def test_rho_plot(self):
        name = 'lvl14E5R50M26Ni2m2b1m3Z01h5.h5'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        snh5 = ps.H5Stella(join(path, name))  # ('lvl14E5R50M26Ni2m2b1m3Z01h5.h5')
        hyd = snh5.Hyd

        times = hyd.T
        x = hyd.R
        y = hyd.Rho

        colorparams = times
        colormap = cm.viridis
        # colormap = cm.jet
        normalize = mcolors.LogNorm(vmin=np.min(colorparams), vmax=np.max(colorparams))

        fig, ax = plt.subplots(figsize=(12, 8))
        for i, t in enumerate(times):
            color = colormap(normalize(t))
            ax.plot(x[i], y[i], color=color)

        # Colorbar setup
        s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
        s_map.set_array(colorparams)

        # Use this to emphasize the discrete color values
        cbar = fig.colorbar(s_map, spacing='proportional', format='%3.1f')  # format='%2i' for integer

        cbar.set_label(r'Times', fontsize=20)

        ax.set_xlabel(r'Radius', fontsize=20)
        ax.set_ylabel('Density', fontsize=20)

        ax.set_xscale('log')
        ax.set_yscale('log')

        ylims = ax.get_ylim()
        ax.set_ylim(top=1.1 * ylims[1])

        plt.show()

    def test_tau_plot(self):
        name = 'lvl14E5R50M26Ni2m2b1m3Z01h5.h5'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        snh5 = ps.H5Stella(join(path, name))  # ('lvl14E5R50M26Ni2m2b1m3Z01h5.h5')
        hyd = snh5.Hyd
        freq = snh5.Freqmn
        nfreq = len(freq)

        tau = snh5.Tau

        times = hyd.T
        it = 15
        zones = np.arange(0, tau.Nzon, 1)
        print('Plot the flux t[{}]= {} '.format(it, times[it]))

        # colorparams = times
        colormap = plt.get_cmap('gray')
        # colormap = cm.viridis

        waves = ps.phys.c / freq * ps.phys.cm_to_angs
        r, theta = np.meshgrid(zones, waves)
        values = np.log10(tau[it])

        # -- Plot... ------------------------------------------------
        levels = np.linspace(np.min(values), np.max(values), 10)
        fig, ax = plt.subplots(1, 1, figsize=(12, 8), constrained_layout=True)
        # ax1.contourf(theta, r,  values, cmap=colormap)
        CS = ax.contour(theta, r,  values)  #, levels=levels)  # cmap=colormap)
        ax.clabel(CS, inline=1, fontsize=10)
        pc = ax.pcolor(theta, r,  values, cmap=colormap)
        ax.set_xscale('log')

        xlim = (1e1, 5e4)
        ax.set_xlim(xlim)
        ax.set_xlabel('Wave [A]')
        ax.set_ylabel('Zone')

        # make a colorbar for the contour lines
        cbar = plt.colorbar(pc, shrink=0.8, extend='both')
        #
        # s_map = cm.ScalarMappable(cmap=colormap)
        # s_map.set_array(levels)
        # cbar = fig.colorbar(s_map,  format='%3.2e')
        # # cbar = fig.colorbar(s_map, spacing='proportional', format='%3.2e')

        cbar.set_label(r'log10(Tau)')  # , fontsize=20)

        ax.text(0.12, 0.98, 't[{:d}]= {:.2f} d'.format(it, times[it]),
                horizontalalignment='right', transform=ax.transAxes)

        plt.show()

    def test_flux_plot(self):
        from scipy.interpolate import griddata

        name = 'lvl14E5R50M26Ni2m2b1m3Z01h5.h5'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        snh5 = ps.H5Stella(join(path, name))  # ('lvl14E5R50M26Ni2m2b1m3Z01h5.h5')
        hyd = snh5.Hyd
        freq = snh5.Freqmn
        nfreq = len(freq)

        fh = snh5.Tau

        times = hyd.T
        it = 20
        r = np.log10(hyd.R)
        zones = np.arange(0, fh.Nzon, 1)
        # y = hyd.Rho
        print('Plot the flux t[{}]= {} '.format(it, times[it]))

        colorparams = times
        colormap = cm.viridis
        # colormap = cm.jet
        normalize = mcolors.LogNorm(vmin=np.min(colorparams), vmax=np.max(colorparams))

        # -- Generate Data -----------------------------------------
        # Using linspace so that the endpoint of 360 is included...
        # azimuths = np.radians(np.linspace(0, 360, nfreq))
        # zeniths = zones

        r, theta = np.meshgrid(zones, freq)
        values = np.log10(fh[it])
        area = 1
        # values = np.random.random((azimuths.size, zeniths.size))
        points = values

        # -- Plot... ------------------------------------------------
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        ax.contourf(theta, r, values, cmap=colormap)
        # Theta, R = np.meshgrid(theta_edges, r_edges)
        # ax.pcolormesh(r, theta, values)
        # ax.scatter(theta, r, c=values, s=2, edgecolor='')

        # grid_r, grid_theta = np.meshgrid(r, theta)
        # data = griddata(points, values, (grid_r, grid_theta), method='cubic', fill_value=0)
        # ax.pcolormesh(theta, r, data.T)

        # c = ax.scatter(theta, r, c=values, s=area, cmap=colormap, alpha=0.75)
        # ax.pcolormesh(theta, r, values, cmap=colormap)
        # ax1.pcolormesh(df.index, r, df.values.T)
        # Colorbar setup
        s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
        s_map.set_array(colorparams)

        # Use this to emphasize the discrete color values
        cbar = fig.colorbar(s_map, spacing='proportional', format='%3.1f')  # format='%2i' for integer
        cbar.set_label(r'Times', fontsize=20)

        ax.text(0.09, 0.98, 't[{:d}]= {:.2f} d'.format(it, times[it]),
                horizontalalignment='right', transform=ax.transAxes)
        plt.show()
