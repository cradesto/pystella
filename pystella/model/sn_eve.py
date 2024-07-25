import logging
# import pprint

import numpy as np
import os

try:
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    is_matplotlib = True
except:
    is_matplotlib = False

from pystella.util.phys_var import phys

logger = logging.getLogger(__name__)
# logger.setLevel(logging.INFO)

__author__ = 'bakl'

eve_elements = ("Ni56", "H", "He", "C", "N", "O", "Ne", "Na", "Mg", "Al"
                , "Si", "S", "Ar", "Ca", "Fe", "Ni")

eve_colors = dict(Ni56="red", H="blue", He="cyan", C="darkorange", N="coral",
                  O="violet", Ne="green", Na="sandybrown",
                  Mg="skyblue", Si="olive", Al="lime",
                  S="indigo", Ar="brown", Ca="purple",
                  Fe='maroon', Ni='magenta',
                  Fe52='blue', Cr48='cyan',
                  Z='black', )
eve_lntypes = dict((k, '--') for k, v in eve_colors.items())  # no y-shift
eve_lntypes['H'] = '-'
eve_lntypes['He'] = '-'
eve_lntypes['O'] = '-'
eve_lntypes['C'] = '-'
eve_lntypes['Ni56'] = '-'
eve_lntypes['Z'] = '-'  # metals

eve_el_m = {'H': 1.008, 'He': 4.003, 'C': 12.011, 'N': 14.007, 'O': 15.999,
            'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305,
            'Al': 26.982, 'Si': 28.086, 'P': 30.974, 'S': 32.066,
            'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
            'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996,
            'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693
            }
eve_el_m['Ni56'] = eve_el_m['Ni']


class PreSN(object):
    """
    A class that holds data of presupernova
    """
    sRho = 'Rho'
    sM = 'M'
    sMcore = 'm_core'
    sT = 'T'
    sR = 'R'
    sV = 'V'
    presn_hydro = (sM, sR, sT, sRho, sV)
    stl_elements = ("H", "He", "C", "N", "O", "Ne", "Na", "Mg", "Al"
                    , "Si", "S", "Ar", "Ca", "Fe", "Ni", "Ni56")
    stl_elements_iso = ("H", "He", "C", "N", "O", "Ne", "Na", "Mg", "Al"
                        , "Si", "S", "Ar", "Ca", "Fe", "Ni", "Ni56", 'Fe52', 'Cr48')

    def __init__(self, name, nzon, elements=stl_elements):
        """Creates a PreSN model instance.  Required parameters:  name, nzon"""
        self._name = name
        self._nzon = nzon
        self._elements = elements
        self._data_hyd = np.empty(nzon, dtype=PreSN.dtype_hyd())
        self._data_chem = np.empty(nzon, dtype=self.dtype_chem())
        self._params = {}
        self._loads = []

    @staticmethod
    def dtype_hyd():
        dt = np.dtype({'names': PreSN.presn_hydro, 'formats': np.repeat('f8', len(PreSN.presn_hydro))})
        return dt

    def dtype_chem(self):
        dt = np.dtype({'names': self.Elements, 'formats': np.repeat('f8', self.Nelements)})
        return dt

    def show_info(self):
        print("-" * 20)
        print(" Name: %s    nzon: %d" % (self.Name, self.nzon))
        print(" m_tot: {:5.3f} r_cen: {:12.6e}".format(self.m_tot, self.r_cen))

    @property
    def Name(self):
        return self._name

    @property
    def Elements(self):
        return self._elements

    @property
    def Nelements(self):
        return len(self.Elements)

    @property
    def nzon(self):
        """ time start"""
        return self._nzon

    @property
    def time_start(self):
        """ time start"""
        return self.par('time_start', 0.)

    @property
    def r_cen(self):
        """Center radius"""
        p = 'r_cen'
        if self.is_set(PreSN.sR):
            d = self.hyd(PreSN.sR)[0] * 0.99  # / 2.  # todo check Rcen
        else:
            d = 0.
        return self.par(p, d)

    @property
    def m_core(self):
        """Core mass"""
        p = PreSN.sMcore
        if self.is_set(PreSN.sM):
            d = self.hyd(PreSN.sM)[0]
        else:
            d = 0.
        return self.par(p, d)

    @property
    def m_tot(self):
        """Total mass"""
        p = 'm_tot'
        if self.is_set(PreSN.sM):
            d = self.hyd(PreSN.sM)[-1]
        else:
            d = 0.
        return self.par(p, d)

    @property
    def rho_cen(self):
        """ Center density"""
        p = 'rho_cen'
        if self.is_set(PreSN.sRho):
            d = self.hyd(PreSN.sRho)[0]
        else:
            d = 0.
        return self.par(p, d)

    @property
    def lg_rho(self):
        """logarithmic density"""
        return np.log10(self.rho)

    @property
    def rho(self):
        """density"""
        return self.hyd(PreSN.sRho)

    @property
    def lg_r(self):
        """logarithmic radius"""
        return np.log10(self.r)

    @property
    def r(self):
        """logarithmic radius"""
        return self.hyd(PreSN.sR)

    @property
    def m(self):
        """Mass"""
        return self.hyd(PreSN.sM)

    @property
    def lgT(self):
        """Log T"""
        return np.log10(self.T)

    @property
    def T(self):
        """Temperature"""
        return self.hyd(PreSN.sT)

    @property
    def V(self):
        """Velocity"""
        return self.hyd(PreSN.sV)

    @property
    def lgV(self):
        """Log Velocity"""
        return np.log10(self.V)

    def hyd(self, v):
        """Hydro data"""
        if v not in self._loads:
            raise ValueError("There is no information about the parameter [%s]. You should set it." % v)
        return self._data_hyd[v]

    @property
    def chem(self):
        """Full hydro data"""
        return self._data_chem

    @property
    def params_keys(self):
        return self._params.keys()

    def par(self, name, d=None):
        return self._params.get(name, d)

    def is_par(self, key):
        return key in self._params

    def set_par(self, name, v):
        self._params[name] = v

    def copy_par(self, src, keys=None):
        if keys is None:
            keys = src.params_keys
        for k in keys:
            try:
                self.set_par(k, getattr(src, k))
            except AttributeError:
                self.set_par(k, src.par(k))

    def lg_el(self, el):
        return np.log10(self.el(el))

    def mass_tot_el(self, el=None, is_diff=False):
        """
        Compute the total mass of element el. Return dict of elements with total mass ff el = None
        :param el: the name of element. Default: None
        :param is_diff: if True use  np.sum(self.el(e)*np.diff(self.m))
        :return: the total mass of the element el
        """

        def m_el(e):
            return np.trapz(self.el(e), self.m)

        def m_el_diff(e):
            dmass = np.diff(self.m)
            dmass = np.insert(dmass, -1, dmass[-1])
            return np.sum(self.el(e) * dmass)

        fm = m_el
        if is_diff:
            fm = m_el_diff

        elements = self.Elements
        if el is not None:
            if isinstance(el, str):
                return fm(el)
                # return m_el_diff(el)
            else:
                elements = el

        mass = {}
        for el in elements:
            mass[el] = fm(el)

        return mass

    def mass_tot_rho(self):
        """
        Compute the total mass via Radius and Density
        """

        dm = np.zeros(self.nzon)
        dm[0] = 4. * np.pi / 3. * (self.r[0] ** 3 - self.r_cen ** 3) * self.rho[0]
        for i in range(1, self.nzon):
            dm[i] = 4. / 3. * np.pi * (self.r[i] ** 3 - self.r[i - 1] ** 3) * self.rho[i]
        # print(f'  M_tot(Density) =  {np.sum(dm)/phys.M_sun:.3f}')
        return np.sum(dm)

    def abund(self, k=None):
        """
        Abundances in k-zone.  k in [1, Nzon]
        :param k: zone. If None, return 2d array for all zones
        :return: array
        """
        if k is None:
            abun = np.zeros((self.nzon, len(self.Elements)))
            for i, ename in enumerate(self.Elements):
                abun[:, i] = self.el(ename)
            return abun
        else:
            abun = [self.el(e)[k - 1] for e in self.Elements]
            return abun

    def chem_norm(self, k=None, norm=None):
        if k is None:
            for j in range(self.nzon):
                self.chem_norm(k=j + 1, norm=norm)
            return

        if norm is None:
            norm = np.sum(self.abund(k))
        for e in self.Elements:
            self._data_chem[e][k-1] = self.el(e)[k-1] / norm

    def el(self, el):
        """
        Get abundances for the element
        :param el: the Element name
        :return: array
        """
        if el not in self.Elements:
            raise ValueError("There is no  element [%s] in elements" % el)
        if el not in self._loads:
            raise ValueError("There is no information about the element [%s]. You should set it." % el)
        return self._data_chem[el]

    def xyz(self, k=-1, xy=('H', 'He'), is_norm=False):
        """
        Compute XYZ for chemical abundances
        :param k: zone, default: -1, last zone. If k = None, return the array for all zones
        :param xy: array Not-metal elements, default: ('H', 'He')
        :param is_norm: normalize to 1, default: False
        :return: XYZ value(s) for zone(s)
        """
        if any([el not in self.Elements for el in xy]):
            raise ValueError("There is no  elements of xy [{}] in Elements [{}]".format(xy, self.Elements))

        norm = 1.
        metals = [el for el in self.Elements if el not in xy]
        if k is None:
            if is_norm:
                norm = np.sum([self.el(ze) for ze in self.Elements], axis=0)
            ed = {el: self.el(el) / norm for el in xy}
            y = np.sum([self.el(ze) for ze in metals], axis=0)
            ed['Z'] = y / norm
        else:
            if is_norm:
                norm = np.sum([self.el(ze)[k] for ze in self.Elements], axis=0)
            ed = {el: self.el(el)[k] / norm for el in xy}
            y = np.sum([self.el(ze)[k] for ze in metals])
            ed['Z'] = y / norm
        return ed

    def write_hyd(self, fname):
        """
        Code readheger.trf:
          BM1=cutmass; -- core Mass
          r(0)=Rcen;
          dum=0.;
          write(12,'(1p,e12.3,i6,2e13.5)') timeStart, NzoneHyd, BM1, Rcen;
          do km=1,NzoneHyd;
            write(12,'(1x,i4,1p,e12.4,e18.10,3e15.7,2e12.4)')
               km, dum, rHyd(km), rhoHyd(km), TpHyd(km), uHyd(km), aMr(km), dum;
          end do;
        :return:
        """
        dum = np.zeros(self.nzon)
        logger.info(' Write hyd-data to %s' % fname)
        zones = range(1, self._nzon + 1)
        with open(fname, 'w') as f:
            f.write('{:12.3e} {:6d} {:13.2e} {:13.5e} {:13.5e}\n'
                    .format(self.time_start, self.nzon, self.m_core / phys.M_sun, self.r_cen, self.rho_cen))
            # a = '#No. Mr  dM  R  dR  Rho PRE T   V'.split()
            # f.write('           '.join(a)+'\n')
            # for _ in zip(zones, self.m/phys.M_sun, dum, self.r, dum, self.rho, dum, self.T, self.V):
            #     f.write(' %4d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e \n' % _)
            # 'evehyd.trf: idum,dum,Radius(j),RHOeve(j),TMPR(j),VELOC(j), dum,dum; '
            a = '#No. M  R  Rho T   V   M  dum '.split()
            f.write('  ' + (' '*18).join(a) + '\n')
            # for _ in zip(zones, self.m / phys.M_sun, self.r, np.log10(self.rho), np.log10(self.T), self.V, self.m / phys.M_sun, dum):
            for _ in zip(zones, self.m / phys.M_sun, self.r, self.rho, self.T, self.V, self.m / phys.M_sun, dum):
                f.write(' %4d %24.16e %24.16e %24.16e %15.7e %15.7e  %24.16e  %8.1e\n' % _)
                # f.write(' %4d %15.5e %15.5e %15.5e %15.5e %15.5e  %15.5e  %8.1e\n' % _)
        return os.path.isfile(fname)

    def plot_chem(self, x='m', elements=eve_elements, ax=None, xlim=None, ylim=None, **kwargs):
        """
        Plot the chemical composition.

        ls = kwargs.get('ls', eve_lntypes), if ls is str then ls is the same for all elements
        colors = kwargs.get('colors', eve_colors), if colors is str then colors is the same for all elements
        loc = kwargs.get('leg_loc', 'best')
        leg_ncol = kwargs.get('leg_ncol', 4)
        lw = kwargs.get('lw', 2), if lw is number then lw is the same for all elements
        marker = kwargs.get('marker', None)
        markersize = kwargs.get('markersize', 4)
        alpha = kwargs.get('alpha', 1)
        figsize = kwargs.get('figsize', (8, 8))
        fontsize = kwargs.get('fontsize', 14)
        is_legend = kwargs.get('is_legend', True)
        """
        if not is_matplotlib:
            return
        # elements = kwargs.get('elements', eve_elements)
        # lntypes = kwargs.get('lntypes', eve_lntypes)
        lntypes = kwargs.get('ls', eve_lntypes)
        if isinstance(lntypes, str):
            lntypes = {el: lntypes for el in elements}
        colors = kwargs.get('colors', eve_colors)
        if isinstance(colors, str):
            colors = {el: colors for el in elements}
        lw = kwargs.get('lw', 2)
        if isinstance(lw, (int, float)):
            lw = {el: lw for el in elements}
        loc = kwargs.get('leg_loc', 'best')
        leg_ncol = kwargs.get('leg_ncol', 4)
        marker = kwargs.get('marker', None)
        markersize = kwargs.get('markersize', 4)
        alpha = kwargs.get('alpha', 1)
        figsize = kwargs.get('figsize', (8, 8))
        fontsize = kwargs.get('fontsize', 14)
        is_legend = kwargs.get('is_legend', True)

        if isinstance(lntypes, str):
            tmp = lntypes
            lntypes = {e: tmp for e in elements}

        is_new_plot = ax is None
        # setup figure
        if is_new_plot:
            plt.matplotlib.rcParams.update({'font.size': fontsize})
            fig = plt.figure(num=None, figsize=figsize, dpi=100, facecolor='w', edgecolor='k')

            gs1 = gridspec.GridSpec(1, 1)
            # gs1.update(wspace=0.1, hspace=0.1, top=0.97, left=0.12, right=0.98)
            gs1.update(wspace=0.1, hspace=0.1, top=0.97, left=0.12, right=0.87)
            ax = fig.add_subplot(gs1[0, 0])

        is_x_lim = xlim is not None
        is_y_lim = ylim is not None

        if x.lower() == 'rsun':
            x = self.r / phys.R_sun
            ax.set_xlabel(r'R [$\mathrm{R}_\odot$]')
        elif x.lower() == 'lgr':
            x = self.r
            ax.set_xscale('log')
            ax.set_xlabel(r'R [cm]')
        elif x.lower() == 'm':
            x = self.m / phys.M_sun
            ax.set_xlabel(r'M [$\mathrm{M}_\odot$]')
        elif x.lower() == 'v':
            x = self.V / 1e5  # to km/s
            ax.set_xlabel(r'V [$km\, s^{-1}$]')
        elif x.lower() == 'z':  # zones
            x = np.arange(0, stop=self.nzon, dtype=np.int) + 1
            ax.set_xlabel(r'Zone')
        else:
            x = self.r
            ax.set_xlabel(r'R [cm]')

        y_min = []
        y_max = []
        for el in elements:
            if self.is_set(el):
                # y = self.lg_el(el)
                y = self.el(el)
                # x = y[np.nonzero(y)]
                # y = y[np.nonzero(y)]
                # y[y<=0] == 1e-15
                ax.plot(x, y, label='{0}'.format(el), color=colors[el], ls=lntypes[el], linewidth=lw[el]
                        , marker=marker, markersize=markersize, alpha=alpha)
                # ax.semilogy(x, y, label='{0}'.format(el), color=colors[el], ls=lntypes[el], linewidth=lw
                #             , marker=marker, markersize=markersize)

                if not is_y_lim:
                    y_min.append(np.min(y))
                    y_max.append(np.max(y))

        if not is_y_lim and len(y_min) > 0:
            ylim = [np.min(y_min), np.max(y_min)]

        if not is_x_lim:
            xlim = np.min(x), np.max(x)

        if is_x_lim or not is_new_plot:
            ax.set_xlim(xlim)

        if is_y_lim or not is_new_plot:
            ax.set_ylim(ylim)
        ax.set_yscale('log')
        if is_new_plot:
            ax.set_ylabel(r'$X_i$')

        if is_legend:
            ax.legend(prop={'size': 9}, loc=loc, ncol=leg_ncol, fancybox=False, frameon=False,
                      markerscale=0, handlelength=3)
            # ax.legend(prop={'size': 9}, loc=3, ncol=4, fancybox=True, shadow=True)
            # plt.grid()
            # plt.show()
        return ax

    def write_abn(self, fname, is_header=False, is_dum=False):
        """
        Write data to file in abn format.
        See code readheger.trf:
         _do km=1,NzoneHyd;
            write(13,'(i4,1p,19e10.3)')km,dum,dum,dum,
        -- No.  Mr    X(H He C N O Ne Na Mg Al Si S Ar Ca Fe Co Ni 56Ni)
              bh(km), bhe(km),bc(km),
              bn(km),bo(km),bne(km),bna(km),bmg(km),bal(km),bsi(km),
              bs(km),bar(km),bca(km),bfe(km),
              bni58(km),bni56(km); -- with Ni58 separated
         _od;
        :return:
        """
        dum = 0.
        logger.info(' Write abn-data to %s' % fname)
        with open(fname, 'w') as f:
            # f.write('%d\n' % self.nzon_abn)
            if is_header:
                if is_dum:
                    s = '%4s  %10s %10s %10s' % ('# zn', ' ', ' ', ' ')
                else:
                    s = '%4s' % '# zn'
                for ename in self.Elements:
                    s += ' %10s' % ename
                f.write('%s\n' % s)
            for i in range(self.nzon):
                if is_dum:
                    s = '%4d  %10.3e %10.3e %10.3e' % (i + 1, dum, dum, dum)
                else:
                    s = '%4d' % (i + 1)
                for ename in self.Elements:
                    s += ' %10.3e' % self.el(ename)[i]
                f.write('%s\n' % s)
        return os.path.isfile(fname)

    def plot_rho(self, x='m', ax=None, xlim=None, ylim=None, **kwargs):
        if not is_matplotlib:
            return
        lw = kwargs.get('lw', 2)
        ls = kwargs.get('ls', '-')
        label = kwargs.get('label', '')
        color = kwargs.get('color', 'black')
        xnorm = kwargs.get('xnorm', 1)
        marker = kwargs.get('marker', None)
        markersize = kwargs.get('markersize', 4)
        alpha = kwargs.get('alpha', 1)

        is_new_plot = ax is None
        # setup figure
        if is_new_plot:
            plt.matplotlib.rcParams.update({'font.size': 14})
            fig = plt.figure(num=None, figsize=(9, 5), dpi=100, facecolor='w', edgecolor='k')

            gs1 = gridspec.GridSpec(1, 1)
            gs1.update(wspace=0.1, hspace=0.1, top=None, left=0.13, right=0.98)
            ax = fig.add_subplot(gs1[0, 0])
            ax.set_ylabel(r'$\rho, [g/cm^3]$ ')

            # if x == 'r':
            #     ax.set_xlabel(r'R [cm]')
            # elif x == 'm':
            #     ax.set_xlabel(r'M [$M_\odot$]')
            # elif x == 'v':
            #     ax.set_xlabel(r'V [$km\, s^{-1}$]')
            # elif x == 'z':
            #     ax.set_xlabel(r'Zone')
            # else:
            #     ax.set_xscale('log')
            #     ax.set_xlabel(r'R [cm]')

        is_x_lim = xlim is not None
        is_y_lim = ylim is not None

        if x == 'm':
            xi = self.m / phys.M_sun * xnorm
            ax.set_xlabel(r'M [$\mathrm{M}_\odot$]')
        elif x == 'v':
            xi = self.V * xnorm
            ax.set_xlabel(r'V [$km\, s^{-1}$]')
        elif x == 'z':
            xi = np.arange(0, self.nzon, dtype=np.int) + 1
            ax.set_xlabel(r'Zone')
        else:
            xi = self.r * xnorm
            ax.set_xlabel(r'R [cm]')

        y = self.rho
        ax.semilogy(xi, y, color=color, ls=ls, linewidth=lw,
                    marker=marker, markersize=markersize, label=label, alpha=alpha)

        if is_new_plot:
            if not is_x_lim and len(xi) > 0:
                xlim = [np.min(xi), np.max(xi)]
                ax.set_xlim(xlim)
            if not is_y_lim and len(y) > 0:
                ylim = [np.min(y), np.max(y)]
                ax.set_ylim(ylim)

        return ax

    def plot_structure(self, elements=eve_elements, xlimR=None, xlimM=None, ylimRho=None, ylimChem=None,
                       title=None, figsize=(12, 8)):
        def set_xlim(ax, lim):
            if lim is not None:
                # ax.set_xlim(lim[0] * 0.5, lim[1] * 2.)
                ax.set_xlim(lim)

        def set_ylim(ax, lim):
            if lim is not None:
                # ax.set_ylim(lim[0]*0.1, lim[1]*10.)
                ax.set_ylim(lim)

        def lims(ain, aout, lim):
            res = np.interp(lim, ain, aout)
            return res

        # if xlimR is not None and xlimM is None:
        #     xlimM = lims(self.r, self.m / phys.M_sun, xlimR)
        #     print("xlimM = {} for xlimR={}".format(xlimM, xlimR))
        # elif xlimM is not None and xlimR is None:
        #     xlimR = lims(self.m / phys.M_sun, self.r, xlimM)
        #     print("xlimR = {} for xlimM={}".format(xlimR, xlimM))

        # Set up the axes with gridspec
        fig = plt.figure(figsize=figsize)
        # fig.subplots_adjust(hspace=0.4, wspace=0.4)
        grid = plt.GridSpec(2, 3, hspace=0.2, wspace=0.4)
        axR = fig.add_subplot(grid[0, 0:2])
        axM = fig.add_subplot(grid[1, 0:2])
        axRhoR = fig.add_subplot(grid[0, 2])
        axRhoM = fig.add_subplot(grid[1, 2])

        self.plot_chem(ax=axR, x='lgR', elements=elements)
        axR.set_xlabel('R, cm')
        axR.set_ylabel(r'$X_i$')
        axR.set_xscale('log')
        axR.legend(frameon=False, ncol=4)
        set_xlim(axR, xlimR)
        set_ylim(axR, ylimChem)

        self.plot_chem(ax=axM, x='m', elements=elements)
        # axM.set_xlabel(r'$M [M_\odot$')
        axM.set_ylabel(r'$X_i$')
        set_xlim(axM, xlimM)
        set_ylim(axM, ylimChem)

        self.plot_rho(ax=axRhoR, x='lgR')
        axRhoR.set_xlabel('R, cm')
        axRhoR.set_xscale('log')
        axRhoR.set_ylabel(r'$\rho, g/cm^3$')
        set_xlim(axRhoR, xlimR)
        set_ylim(axRhoR, ylimRho)

        self.plot_rho(ax=axRhoM, x='m')
        # axRhoM.set_xlabel(r'$M, M_\odot$')
        axRhoM.set_ylabel(r'$\rho, g/cm^3$')
        set_xlim(axRhoM, xlimM)
        set_ylim(axRhoM, ylimRho)

        if title is not None:
            axR.text(0.5, 1.07, title, transform=axR.transAxes, fontsize=14)
        return fig

    def is_set(self, name):
        return name in self._loads

    def set_hyd(self, name, vec, is_exp=False):
        if len(vec) != self.nzon:
            raise ValueError("The length of vector [%d] should be %d" % (len(vec), self.nzon))
        if name not in self._loads:
            self._loads.append(name)

        if is_exp:
            self._data_hyd[name] = 10. ** vec
        else:
            self._data_hyd[name] = vec

    def set_chem(self, name, vec, is_exp=False):
        if len(vec) != self.nzon:
            raise ValueError("The length of vector [%d] should be %d" % (len(vec), self.nzon))
        if name not in self._loads:
            self._loads.append(name)
        if is_exp:
            self._data_chem[name] = 10. ** vec
        else:
            self._data_chem[name] = vec

    def zone_reduce(self, by=sM, diff=1.01, start=0, end=None, mode='g'):
        """

        :param by: 'Rho' 'M' 'T' 'R' 'V', default: 'M'
        :param diff: geom progression, default: 1.01
        :param start:
        :param end:
        :param mode:
        :return:
        """
        from pystella.util.math import shrink, portion_index

        x = self.hyd(by)

        def where(a):
            return shrink(a, diff=diff, mode=mode)

        idxs = portion_index(x, where, start=start, end=end, isByEl=False)

        newPreSN = PreSN(self.Name, len(idxs), elements=self.Elements)
        # hyd reshape
        for v in PreSN.presn_hydro:
            old = self.hyd(v)
            new = old[idxs]
            newPreSN.set_hyd(v, new)

        # abn reshape
        for el in self.Elements:
            old = self.el(el)
            new = old[idxs]
            newPreSN.set_chem(el, new)

        # copy parameters
        # copy parameters
        newPreSN.copy_par(self)  # keys=['time_start', 'm_tot',  'm_core', 'r_cen'])

        return newPreSN

    def set_composition(self, zones, sample=None, is_add=True, is_normalize=True):
        """
        Set abundances with solar composition
        :return:
        """
        if sample is None:
            sample = sample_sol()

        # abn reshape
        for el, Xi in sample.items():
            y = self.el(el)
            for k in zones:
                if is_add:
                    y[k - 1] += Xi
                else:
                    y[k - 1] = Xi
            self.set_chem(el, y)

        if is_normalize:
            for k in zones:
                self.chem_norm(k)

    def bad_zone_reduce(self, diff=1.05, start=0, end=None, mode='g'):
        from pystella.util.math import shrink
        x = self.m
        if end is None:
            end = len(x)
        idxs = np.arange(len(x))

        xx = x[start:end]
        idx = shrink(xx, diff=diff, mode=mode)
        idxs = np.concatenate((np.arange(start), idx, np.arange(end, len(x))))

        if start > 0:
            idxs = idxs[:start - 1]
        else:
            idxs = []

        newPreSN = PreSN(self.Name, len(idxs), elements=self.Elements)
        # hyd reshape
        for v in PreSN.presn_hydro:
            old = self.hyd(v)
            new = old[idxs]
            newPreSN.set_hyd(v, new)

        # abn reshape
        for el in self.Elements:
            old = self.el(el)
            new = old[idxs]
            newPreSN.set_chem(el, new)

        return newPreSN

    def cut(self, name=None, start=0, end=None, elements=None, pars=None):
        """
        Cut zones in the envelope between nstart:nend
        @param name: the name of new PreSN. Take from parent, if it's None.
        @param start: zone number of the left edge. Default: 0 (first zone)
        @param end: zone number of the right edge. Default: None,  (equal last zone)
        @param elements: the elements to be left on hold.. Take from parent, if it's None.
        @return: new PreSN
        """
        if name is None:
            name = self.Name
        if end is None:
            end = self.nzon
        if elements is None:
            elements = self.Elements

        nznew = end - start
        newPreSN = PreSN(name, nznew, elements=elements)

        for v in PreSN.presn_hydro:
            old = self.hyd(v)
            new = old[start:end]
            newPreSN.set_hyd(v, new)

        # abn reshape
        for el in elements:
            old = self.el(el)
            new = old[start:end]
            newPreSN.set_chem(el, new)

        # copy parameters
        newPreSN.copy_par(self)
        # for p in ['m_core', 'r_cen']:
        #     v = getattr(newPreSN, p)
        #     newPreSN.set_par(p, v)
        # # newPreSN.copy_par(self)  # keys=['time_start', 'm_tot',  'm_core', 'r_cen'])
        # newPreSN.copy_par(self, keys=['time_start', 'm_tot'])
        return newPreSN

    def reshape(self, nz: int, name=None, start: int = 0, end: int = None, axis: str = sM, xmode: str = 'resize', kind='np'):
        """
        Reshape parameters of envelope from nstart to nend to nz-zones
        :param nz: new zones
        :param name: the name of new PreSN. Take from parent, if it's <=0.
        :param start: zone number to start reshaping. Default: 0 (first zone)
        :param end: zone number to end reshaping. Default: None,  (equal last zone)
        :param axis: [M OR R OR V] - reshape along mass or radius or velocity coordinate. Default: M
        :param xmode: [lin OR rlog OR resize] - linear OR reversed log10 OR add/remove points. Default: resize
        :param kind: [np OR interp1d(..kind)], kind is  ('np=np.interp', 'linear', 'nearest', 'zero', 'slinear', 'quadratic, 'cubic'). Default: np
        :return: new preSN with reshaping zones
        """
        from scipy.interpolate import interp1d
        # from scipy.interpolate import splev, splrep
        from scipy.interpolate import UnivariateSpline
        from scipy.ndimage import gaussian_filter1d

        def rlogspace(s, e, n):
            r = np.exp(np.linspace(np.log(s), np.log(e), n))
            r = (e - r + s)
            return r[::-1]

        def add_point(x, mode='lin'):  # 'lin' 'log'
            """
            Find max interval in x and insert the new point in the middle (lin or geom) of it
            :param x: array
            :param mode:
            :return:
            """
            dif = np.diff(x) / x[1:]
            idx = np.argmax(dif)
            if mode == 'lin':
                p = (x[idx] + x[idx + 1]) / 2.
            elif mode == 'geom':
                p = np.sqrt(x[idx] * x[idx + 1])
            else:
                raise ValueError('Mode should be "lin" lor "geom"')
            logger.info('         To interval {}[{:.6e} - {:.6e}] added {} '.format(idx, x[idx], x[idx + 1], p))
            xn = np.insert(x, idx + 1, p)
            return xn

        def remove_point(x):  # 'lin' 'log'
            """
            Find min delta and remove the right point
            """
            dif = np.diff(x) / x[1:]
            idx = np.argmin(dif)
            xn = np.delete(x, idx + 1)
            if idx+2 == len(x):
                logger.debug('         Remove right point {} {:.8e} [dx= {:.8e}] '.format(idx+1, x[idx + 1], dif[idx]))
            else:
                logger.debug('         Remove point {} {:.8e} [dx= {:.8e}] next: {:.8e}'.format(idx+1, x[idx + 1], dif[idx], x[idx + 2]))  
            return xn

        def resize_points(x, n: int, mode: str = 'lin'):
            """
            Add or remove points in the array x
            :param x:  the array is not changed
            :param n: number points to add or remove
            :param mode:  should be "lin" or "geom". Default: lin
            :return: the resized array
            """
            n_old = len(x)
            xn = np.copy(x)
            if n == n_old:  # nothing to do, return the copy
                return xn

            if n > n_old:
                f = lambda xxx: add_point(xxx, mode=mode)
            else:
                f = lambda xxx: remove_point(xxx)

            for i in range(abs(n - n_old)):
                xn = f(xn)
            return xn

        def x_reshaped(x, n):
            if xmode == 'lin':
                res = np.linspace(x[0], x[-1], n)
            elif xmode == 'rlog':
                res = rlogspace(x[0], x[-1], n)
            elif xmode == 'resize':
                res = resize_points(x, n)
            else:
                raise ValueError('Such xmode "{}" is not supported.'.format(xmode))
            return res

        def interp(xn, x, v, s: int, e: int, kind: str, is_log: bool = False):
            res = []
            if s > 0:
                res = v[:s]  # save points before start
            xi = x[s:e]
            yi = v[s:e]
            if is_log:
                yi = np.log10(yi)

            if kind == 'np':
                yy = np.interp(xn, xi, yi)
            elif kind == 'spline':
                spl = UnivariateSpline(xi, yi)
                yy = spl(xn)
            elif kind == 'gauss':
                yii = gaussian_filter1d(yi, 3)
                yy = np.interp(xn, xi, yii)
            else:
                interp_linear = interp1d(xi, yi, kind=kind)
                yy = interp_linear(xn)

            if is_log:
                yy = 10. ** yy
            res = np.append(res, yy)
            return res

        if nz <= 0:
            nz = self.nzon

        nznew = start + nz
        if name is None:
            name = self.Name

        newPreSN = PreSN(name, nznew, elements=self.Elements)

        if end is None:
            end = self.nzon

        logger.info(f'axis= {axis} nz= {nz} nznew= {nznew} start= {start} end= {end}')

        # hyd reshape
        if axis == PreSN.sM:
            xx = self.m
        elif axis == PreSN.sR:
            xx = self.r
        elif axis == PreSN.sV:
            xx = self.V
        else:
            raise ValueError('Such axis "{}" is not supported.'.format(axis))

        xx = xx / max(abs(xx))  # norm
        xxx = x_reshaped(xx, nz)
        if np.any(np.diff(xxx) < 0.):
            for i, dx in enumerate(np.diff(xxx)):
                if dx <= 0:
                    logger.error("ERROR reshaped: {} xxx= {}  dx= {} ".format(i, xxx[i], dx))
            raise ValueError('The interval beetween some of {} elements is < 0.'.format(len(xxx)))

        # from pprint import pprint

        for vv in PreSN.presn_hydro:
            old = self.hyd(vv)
            new = interp(xxx, xx, old, s=start, e=end, kind=kind, is_log=False)
            # if vv == PreSN.sRho:
            #    rho_new = interp(xxx, xx, old, s=start, e=end, kind='next') #, is_log=True)
            # else:
            #    new = interp(xxx, xx, old, s=start, e=end, kind=kind)
            newPreSN.set_hyd(vv, new)
            # print(f'{vv} before: old[{len(xx)}-1]= {old[len(xx)-2]:12.7e} new[{len(xxx)}-1]= {new[len(xxx)-2]:12.7e}')
            logger.info(f'{vv} before: old[0]= {old[0]:12.7e} new[0]= {new[0]:12.7e}')
            logger.info(f'{vv} before: old[{len(xx)}]= {old[len(xx) - 1]:12.7e} new[{len(xxx)}]= {new[len(xxx) - 1]:12.7e}')
            # print(f'\n{vv} before: {len(xx)}')
            # pprint(list(zip(range(1, len(xx)+1), xx, old)))
            # print(f'{vv} after:  {len(xxx)}')
            # pprint(list(zip(range(1, len(xxx)+1), xxx, new)))

        # Density Normalization: m_tot(NEW) should be equal m_tot(OLD)
        m_rho = newPreSN.mass_tot_rho() + newPreSN.m_core
        rho = newPreSN.rho * newPreSN.m_tot / m_rho
        newPreSN.set_hyd(PreSN.sRho, rho)

        # abn reshape
        for el in self.Elements:
            old = self.el(el)
            new = interp(xxx, xx, old, s=start, e=end, kind='np')
            # new = interp(xxx, xx, old, s=start, e=end, kind=kind, is_log=True)
            newPreSN.set_chem(el, new)

        # copy parameters
        newPreSN.copy_par(self)  # keys=['time_start', 'm_tot',  'm_core', 'r_cen'])

        return newPreSN

    def clone(self):
        presn = PreSN(self.Name, self.nzon, elements=self.Elements)
        presn.copy_par(self)

        # hydro
        for k in self.presn_hydro:
            presn.set_hyd(k, self.hyd(k))

        # chem
        for ename in self.Elements:
            presn.set_chem(ename, self.el(ename))

        return presn

    def boxcar(self, box_dm: float = 0.5, n: int = 4, el_included=None, m_uplim: float=np.inf, is_info: bool = False):
        """
        The function runs a boxcar average to emulate the mixing of chemical composition.
        :param box_dm: float. The boxcar width. Default value is 0.5 Msun.
        :param n: int. The number of repeats. Default value is 4
        :param el_included: the tuple of included elements. If None = all elements are included. Default: None
        :param m_uplim: float. The mass limit from above. Default value is np.inf
        :param is_info: bool. Prints some debug information. Default value is False
        """
        clone = self.clone()
        abund = clone.abund()
        #     abun = np.zeros((clone.nzon, len(clone.Elements)))
        if el_included is None:
            el_included = clone.Elements
        if m_uplim is None:
            m_uplim = np.inf
        m = clone.m / phys.M_sun
        dmass = np.diff(m)
        dmass = np.insert(dmass, -1, dmass[-1])

        # todo Check left boundary condition fo Fe, Si
        for l in range(n):  # the iteration number
            logger.debug(f'Attempt # {l}')
            for k in range(clone.nzon):
                if m[k] > m_uplim:
                    logger.warn(f'Exit from zone cycle m[k] > m_uplim: k= {k} m= {m[k]:.4f} m_uplim= {m_uplim:.4f}')
                    break #  stop because mass is over the limit
                kk = k + 1
                dm = dmass[k]
                while dm < box_dm and kk <= clone.nzon:
                    kk += 1
                    dm = np.sum(dmass[k:kk])

                logger.debug(f'{k}: kk= {kk} dm= {dm:.4f} m= {m[k]:.4f}')
                if dm > 1e-6:
                    for i, ename in enumerate(clone.Elements):
                        if ename in el_included:
                            dm_e = np.dot(abund[k:kk, i], dmass[k:kk])
                            abund[k, i] = dm_e / dm
            #             abun[k,i] = x[k]
        for i, ename in enumerate(clone.Elements):
            clone.set_chem(ename, abund[:, i])        
            logger.debug(f"{ename}: {clone.el(ename)}")

        return clone

    def smooth(self, window_length: int = 15, polyorder: int = 2, mode: str = 'interp'):
        from scipy.signal import savgol_filter
        """
        The function runs a Savitzky-Golay filter to smooth density distribution.
        :param window_length: int. Default value is 5
        :param polyorder: int. Default value is 2
        # :param el_included: the tuple of included elements. If None = all elements are included. Default: None
        :param is_info: bool. Prints some debug information. Default value is False
        """
        clone = self.clone()
        rho_s = savgol_filter(np.log(clone.rho), window_length, polyorder, mode=mode)
        rho_s = np.exp(rho_s)
        r = clone.r
        # np.set_printoptions(precision=2)
        logger.debug('Rho')
        logger.debug(list(zip(clone.rho, rho_s)))
        m_s = np.zeros(clone.nzon)
        # RHO(1) = DM(1) / (Ry(1) ** 3 - Rce ** 3);
        m_s[0] = clone.m_core + 4. / 3. * np.pi * clone.rho_cen * (r[0] ** 3 - clone.r_cen ** 3)

        for k in range(1, clone.nzon):
            dm = 4. / 3. * np.pi * rho_s[k] * (r[k] ** 3 - r[k - 1] ** 3)
            m_s[k] = m_s[k-1] + dm
        
        logger.debug('Mass: ')
            # pprint.pprint(list(zip(clone.m, m_s)))
        logger.debug(list(zip(clone.m, m_s)))
        clone.set_hyd('Rho', rho_s)
        clone.set_hyd('M', m_s)
        # print(clone.rho)
        return clone


# ==============================================
def load_rho(fname, path: str = None):
    import linecache
    if path is not None:
        fname = os.path.join(path, fname)
    if not os.path.isfile(fname):
        logger.error(' No rho-data for %s' % fname)
        raise ValueError(' No rho-data for %s' % fname)
        # return None
    logger.info(' Load rho-data from  %s' % fname)

    # try to get column names from header
    header = linecache.getline(fname, 1).split()
    logger.info('header(len= {}): {} '.format(len(header), ' '.join(header)))

    if len(header) < 22 or len(header) >= 25:  # default header
        header = "zone mass lgR lgTp lgRho u Ni56 H He C N O Ne Na  Mg  Al  Si  S  Ar  Ca  Fe  Ni".split()

    usecols, col_names = [], []
    for i, e in enumerate(header):
        if not e in col_names:
            col_names.append(e)
            usecols.append(i)

    logger.info(' col_names=  {}'.format(' '.join(col_names)))
    # logger.info(' usecols=  {}'.format(usecols))
    dt = np.dtype({'names': col_names, 'formats': np.repeat('f8', len(col_names))})

    data = np.loadtxt(fname, comments='#', skiprows=2, dtype=dt, usecols=usecols)

    nz = len(data['lgR'])
    ###
    name = os.path.basename(os.path.splitext(fname)[0])
    col_map = {PreSN.sR: 'lgR', PreSN.sM: 'mass', PreSN.sT: 'lgTp', PreSN.sRho: 'lgRho', PreSN.sV: 'u'}
    presn = PreSN(name, nz)
    presn.set_hyd(PreSN.sV, np.zeros(nz))
    for k, v in col_map.items():
        try:
            presn.set_hyd(k, data[v], is_exp=v.startswith('lg'))
        except ValueError as ex:
                logger.warn('No col: {} (set {} = 0.) in data from: {}. '.format(v, k, fname))
                presn.set_hyd(k, np.zeros(nz))


    # CGS
    presn.set_hyd('M', presn.m * phys.M_sun)

    for ename in presn.Elements:
        presn.set_chem(ename, data[ename], is_exp=True)

    return presn


def load_hyd_abn(name, path='.', abn_elements=PreSN.stl_elements, skiprows=0, comments='#',
                 is_rho=False, is_dm=True, is_dum=False):
    """
    Load progenitor from hyd- + abn- files.

    is_dm: if True, the column 2 is used as dM. Default: False, he column 2 is used as M.
    if is_dum:
        col_names = ("zone dum1 dum2 dum3 " + ' '.join(abn_elements)).split()
    else:
        col_names = ("zone " + ' '.join(abn_elements)).split()

    Code readheger.trf:
      BM1=cutmass; -- core Mass
      r(0)=Rcen;
      dum=0.;
      write(12,'(1p,e12.3,i6,2e13.5)') timeStart, NzoneHyd, BM1, Rcen;
      do km=1,NzoneHyd;
        write(12,'(1x,i4,1p,e12.4,e18.10,3e15.7,2e12.4)')
           km, dum, rHyd(km), rhoHyd(km), TpHyd(km), uHyd(km), aMr(km), dum;
      enddo;
    :return: PreSN
    """
    # abn_elements = 'H He C N O Ne Na  Mg  Al  Si  S  Ar  Ca  Fe  Ni Ni56'.split()

    # hydro
    ext_hyd = '.hyd'
    hyd_file = os.path.join(path, name + ext_hyd)
    if not os.path.isfile(hyd_file):
        logger.error(' No file for %s' % hyd_file)
        return None

    logger.info(' Load hyd-data from  %s' % hyd_file)

    def set_params(pre, a):
        if len(a) > 0:
            if len(a) == 5:
                time_start, nzon, m_core, r_cen, rho_cen = a
                pre.set_par('time_start', time_start)
                pre.set_par('m_core', m_core * phys.M_sun)
                pre.set_par('r_cen', r_cen)
                pre.set_par('rho_cen', rho_cen)
            elif len(a) == 4:
                time_start, nzon, m_core, r_cen = a
                pre.set_par('time_start', time_start)
                pre.set_par('m_core', m_core * phys.M_sun)
                pre.set_par('r_cen', r_cen)
            elif len(a) == 2:
                time_start, nzon = a
                pre.set_par('time_start', time_start)
        return pre

    # read table data
    if is_dm:
        col_names = "zone dm R Rho T V M".split()
    else:
        col_names = "zone M R Rho T V M2".split()

    a = []
    # Load header
    with open(hyd_file, 'r') as f:
        header_line = f.readline()
    if len(header_line) > 0:
        a = [float(x) for x in header_line.split()]

    # Load data
    dt = np.dtype({'names': col_names,
                   'formats': ['i4'] + list(np.repeat('f8', len(col_names) - 1))})

    data_hyd = np.loadtxt(hyd_file, comments='#', skiprows=1, dtype=dt, usecols=np.arange(len(col_names)))

    nz = len(data_hyd['R'])

    presn = PreSN(name, nz, elements=abn_elements)
    set_params(presn, a)

    col_map = {PreSN.sR, PreSN.sT, PreSN.sRho, PreSN.sV}
    for v in col_map:
        presn.set_hyd(v, data_hyd[v], is_exp=v.startswith('lg'))

    # Set header data
    set_params(presn, a)

    # Set Mass
    if is_rho:
        r = presn.r
        rho = presn.rho
        r = np.insert(r, 0, presn.r_cen)
        #       rho = np.insert(rho, 0, presn.rho_cen)
        dm = np.zeros(nz)
        for i in range(nz):
            dm[i] = (r[i + 1] ** 3 - r[i] ** 3) * rho[i] * 4. / 3. * np.pi
        #            dm[i] = (r[i + 1] ** 3 - r[i] ** 3) * rho[i + 1] * 4. * np.pi / 3.
        m = np.cumsum(dm)
        m += presn.m_core
    else:
        m = data_hyd[PreSN.sM] * phys.M_sun

    presn.set_hyd(PreSN.sM, m)

    # Set chemical composition
    ext_abn = '.abn'
    abn_file = os.path.join(path, name + ext_abn)
    if not os.path.isfile(abn_file):
        logger.error(' No file for %s' % abn_file)
        return None

    logger.info(' Load abn-data from  %s' % abn_file)
    col_names = ("zone " + ' '.join(abn_elements)).split()
    if is_dum:
        col_names = ("zone dum1 dum2 dum3 " + ' '.join(abn_elements)).split()

    # dt = np.dtype({'names': col_names, 'formats': np.repeat('f8', len(col_names))})
    dt = np.dtype({'names': col_names,
                   'formats': ['i4'] + list(np.repeat('f8', len(col_names) - 1))})
    # logger.info(dt)
    data_chem = np.loadtxt(abn_file, comments=comments, skiprows=skiprows, dtype=dt)

    for ename in abn_elements:
        presn.set_chem(ename, data_chem[ename])

    return presn


def sample_sol():
    el = dict(H=7.0600E-01, He=2.7500E-01, C=3.0700E-03, N=1.1100E-03, O=9.6100E-03, Ne=1.7500E-03, Na=3.3400E-05,
              Mg=6.6000E-04, Al=5.8100E-05, Si=7.1100E-04, S=4.1800E-04, Ar=9.2800E-05, Ca=6.2000E-05,
              Fe=1.3700E-03, Ni=7.3400e-05)
    norm = sum(el.values())
    el = {e: v / norm for e, v in el.items()}
    return el
