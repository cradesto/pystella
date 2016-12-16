import os
import warnings
import numpy as np
import logging

import matplotlib.pyplot as plt
from matplotlib import gridspec

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

__author__ = 'bakl'

eve_elements = map(str.strip, "Ni56 H He C N O Ne Na  Mg  Al  Si  S  Ar  Ca  Fe  Ni".split())
eve_colors = dict(Ni56="red", H="blue", He="cyan", C="darkorange", N="coral",
                  O="violet", Ne="green", Na="sandybrown",
                  Mg="skyblue", Si="olive", Al="lime",
                  S="indigo", Ar="brown", Ca="purple",
                  Fe='maroon', Ni='magenta')
eve_lntypes = dict((k, '--') for k, v in eve_colors.items())  # no y-shift
eve_lntypes['H'] = '-'
eve_lntypes['He'] = '-'
eve_lntypes['O'] = '-'
eve_lntypes['C'] = '-'
eve_lntypes['Ni56'] = '-'


# eve_lntypes['Si'] = '-'
# eve_lntypes['Fe'] = '-'


class StellaEve:
    def __init__(self, name, path='./'):
        """Creates a StellaEve model instance.  Required parameters:  name."""
        self.name = name
        self.path = path  # path to model files
        self._data = None

    def show_info(self):
        ext = ('eve', 'rho', 'xni')
        for e in ext:
            fname = os.path.join(self.path, self.name + '.' + e)
            if os.path.isfile(fname):
                print "Exist %s-file: %s" % (e, fname)
            else:
                print "No %s-file: %s" % (e, fname)

    @property
    def lg_rho(self):
        """logarithmic density"""
        return self._data['lgRho']

    @property
    def rho(self):
        """density"""
        return 10. ** self.lg_rho

    @property
    def lg_r(self):
        """logarithmic radius"""
        return self._data['lgR']

    @property
    def r(self):
        """logarithmic radius"""
        return 10. ** self.lg_r

    @property
    def mass(self):
        """Mass"""
        return self._data['mass']

    @property
    def lgT(self):
        """Log T"""
        return self._data['lgTp']

    @property
    def T(self):
        """Temperature"""
        return 10. ** self.lgT

    @property
    def data(self):
        """Full data"""
        return self._data

    @property
    def is_load(self):
        """Check if data has been loaded."""
        warnings.warn("deprecated", DeprecationWarning)
        return self.is_load_rho

    @property
    def is_load_rho(self):
        """Check if data has been loaded."""
        return self._data is not None

    @property
    def rho_file(self):
        return os.path.join(self.path, self.name + '.rho')

    @property
    def is_rho_data(self):
        ext = ['rho']
        return any(map(os.path.isfile, [os.path.join(self.path, self.name + '.' + e) for e in ext]))

    def lg_el(self, el):
        if el not in eve_elements:
            raise ValueError("There is no  element [%s] in eve" % el)
        if not self.is_load:
            raise StandardError("Eve-Data has not been loaded. Check and load from %s" % self.rho_file)

        return self.data[el]

    def el(self, el):
        return 10. ** self.lg_el(el)

    def load(self):  # for backward compatibility
        warnings.warn("deprecated", DeprecationWarning)
        return self.load_rho()

    def load_rho(self):
        if not self.is_rho_data:
            logger.error(' No rho-data for %s' % self.rho_file)
            return None
        logger.info(' Load rho-data from  %s' % self.rho_file)
        col_names = "zone mass lgR lgTp lgRho lgDm Ni56 H He C N O Ne Na  Mg  Al  Si  S  Ar  Ca  Fe  Ni"
        dt = np.dtype({'names': map(str.strip, col_names.split()), 'formats': np.repeat('f8', len(col_names))})
        # dt = np.dtype({'names': map(str.strip, col_names.split()), 'formats': ['i4']+np.repeat('f8', len(col_names))})

        data = np.loadtxt(self.rho_file, comments='#', skiprows=2, dtype=dt)

        self._data = data
        return self
        # return self._data

    def plot_chem(self, x='m', ax=None, elements=None, xlim=None, ylim=None,
                  leg_loc=3, leg_ncol=4, lw=2, lntypes=None, is_save=False):
        if elements is None:
            elements = eve_elements
        if lntypes is None:
            lntypes = eve_lntypes

        is_new_plot = ax is None
        # setup figure
        if is_new_plot:
            plt.matplotlib.rcParams.update({'font.size': 14})
            fig = plt.figure(num=None, figsize=(12, 12), dpi=100, facecolor='w', edgecolor='k')

            gs1 = gridspec.GridSpec(1, 1)
            gs1.update(wspace=0.1, hspace=0.1, top=None, left=0.1, right=0.98)
            ax = fig.add_subplot(gs1[0, 0])

        is_x_lim = xlim is not None
        is_y_lim = ylim is not None

        ib = 0
        x_min = []
        x_max = []
        y_min = []
        y_max = []
        if x == 'r':
            x = self.r
            ax.set_xlabel(r'R [cm]')
        if x == 'lgr':
            x = self.r
            ax.set_xscale('log')
            ax.set_xlabel(r'R [cm]')
        else:
            x = self.mass
            ax.set_xlabel(r'M [$M_\odot$]')

        for el in elements:
            ib += 1
            y = self.el(el)
            ax.plot(x, y, label='%s' % el, color=eve_colors[el], ls=lntypes[el], linewidth=lw)

            if not is_x_lim:
                x_min.append(np.min(x))
                x_max.append(np.max(x))
            if not is_y_lim:
                y_min.append(np.min(y))
                y_max.append(np.max(y))

        if not is_x_lim:
            xlim = [np.min(x_min), np.max(x_max)]
        if not is_y_lim:
            ylim = [np.min(y_min), np.max(y_min)]

        ax.set_xlim(xlim)

        ax.set_yscale('log')
        ax.set_ylim(ylim)
        ax.set_ylabel(r'$log10(X_i)$')

        if is_new_plot:
            ax.legend(prop={'size': 9}, loc=leg_loc, ncol=leg_ncol, fancybox=False, frameon=True)
            # ax.legend(prop={'size': 9}, loc=3, ncol=4, fancybox=True, shadow=True)
            # plt.grid()
            # plt.show()

        if is_save:
            fsave = os.path.join(os.path.expanduser('~/'), 'chem_%s.pdf' % self.name)
            logger.info(" Save plot to %s " % fsave)
            ax.get_figure().savefig(fsave, bbox_inches='tight', format='pdf')

        return ax


class StellaHydAbn:
    abn_elements = 'H He C N O Ne Na  Mg  Al  Si  S  Ar  Ca  Fe  Ni Ni56'.split()

    def __init__(self, name, path='./'):
        """Creates a StellaEve model instance.  Required parameters:  name."""
        self.name = name
        self.path = path  # path to model files
        self._data_hyd = None
        self._data_chem = None

    def show_info(self):
        ext = ('hyd', 'abn')
        for e in ext:
            fname = os.path.join(self.path, self.name + '.' + e)
            if os.path.isfile(fname):
                print "Exist %s-file: %s" % (e, fname)
            else:
                print "No %s-file: %s" % (e, fname)

    @property
    def zone_hyd(self):
        """logarithmic density"""
        return self._data_hyd['zone']

    @property
    def zone_abn(self):
        """logarithmic density"""
        return self._data_chem['zone']

    @property
    def lg_rho(self):
        """logarithmic density"""
        return np.log10(self.rho)

    @property
    def rho(self):
        """density"""
        return self._data_hyd['Rho']

    @property
    def lg_r(self):
        """logarithmic radius"""
        return np.log10(self.r)

    @property
    def r(self):
        """logarithmic radius"""
        return self._data_hyd['R']

    @property
    def mass(self):
        """Mass"""
        return self._data_hyd['Mr']

    @property
    def lgT(self):
        """Log T"""
        return np.log10(self.T)

    @property
    def T(self):
        """Temperature"""
        return self._data_hyd['Tp']

    @property
    def V(self):
        """Velocity"""
        return self._data_hyd['V']

    @property
    def lgV(self):
        """Log Velocity"""
        return np.log10(self.V)

    @property
    def hyd(self):
        """Full hydro data"""
        return self._data_hyd

    @property
    def chem(self):
        """Full hydro data"""
        return self._data_chem

    @property
    def is_load_hyd(self):
        """Check if data has been loaded."""
        return self._data_hyd is not None

    @property
    def is_load_chem(self):
        """Check if data has been loaded."""
        return self._data_chem is not None

    @property
    def hyd_file(self):
        return os.path.join(self.path, self.name + '.hyd')

    @property
    def abn_file(self):
        return os.path.join(self.path, self.name + '.abn')

    @property
    def is_file_hyd(self):
        return os.path.isfile(os.path.join(self.path, self.name + '.hyd'))

    @property
    def is_file_abn(self):
        return os.path.isfile(os.path.join(self.path, self.name + '.abn'))

    @property
    def nzon(self):
        if self.is_load_hyd and self.is_load_chem:
            if self.nzon_abn != self.nzon_hyd:
                raise (ValueError("Zones numbers should be the same in hyd-[%d] and abn-[%d] files."
                                  % (self.nzon_hyd, self.nzon_abn)))
            return self.nzon_abn
        elif self.is_load_hyd:
            return self.nzon_hyd
        elif self.is_load_chem:
            return self.nzon_abn
        else:
            return None

    @property
    def nzon_hyd(self):
        if not self.is_load_hyd:
            raise StandardError("Eve-hyd has not been loaded. Check and load from %s" % self.hyd_file)
        return len(self.r)

    @property
    def nzon_abn(self):
        if not self.is_load_chem:
            raise StandardError("Eve-chem has not been loaded. Check and load from %s" % self.abn_file)
        return len(self.el('H'))

    def lg_el(self, el):
        return np.log10(self.el(el))

    def el(self, el):
        if el not in StellaHydAbn.abn_elements:
            raise ValueError("There is no  element [%s] in abn_elements" % el)
        if not self.is_load_chem:
            raise StandardError("Eve-Data has not been loaded. Check and load from %s" % self.abn_file)

        return self._data_chem[el]

    def load_hyd(self):
        """
        Code readheger.trf:
          BM1=cutmass; -- core Mass
          r(0)=Rcen;
          dum=0.;
          write(12,'(1p,e12.3,i6,2e13.5)') timeStart, NzoneHyd, BM1, Rcen;
          do km=1,NzoneHyd;
            write(12,'(1x,i4,1p,e12.4,e18.10,3e15.7,2e12.4)')
               km, dum, rHyd(km), rhoHyd(km), TpHyd(km), uHyd(km), aMr(km), dum;
          enddo;
        :return:
        """
        if not self.is_file_hyd:
            logger.error(' No hyd-data for %s' % self.hyd_file)
            return None
        logger.info(' Load hyd-data from  %s' % self.hyd_file)
        dt = self.dtype_hyd()
        self._data_hyd = np.loadtxt(self.hyd_file, comments='#', skiprows=1, dtype=dt)
        return self

    @staticmethod
    def dtype_hyd():
        col_names = "zone dum1 R Rho Tp V Mr dum2"
        dt = np.dtype(dict(names=map(str.strip, col_names.split()),
                           formats=['i4'] + np.repeat('f8', len(col_names) - 1).tolist()))
        return dt

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
          enddo;
        :return:
        """
        if not self.is_load_hyd:
            raise ValueError(' No loaded hyd-data for %s')

        dum = np.zeros(self.nzon_hyd)
        logger.info(' Write hyd-data to %s' % fname)
        with open(fname, 'a') as f:
            f.write('%d\n' % self.nzon_hyd)
            for _ in zip(self.zone_hyd, dum, self.r, self.rho, self.T, self.V, self.mass, dum):
                f.write(' %4d %12.4e %18.10e %15.7e %15.7e %15.7e %12.4e %12.4e\n' % _)
        return os.path.isfile(fname)

    def load_abn(self):
        """
        Code readheger.trf:
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
        if not self.is_file_abn:
            logger.error(' No abn-data for %s' % self.abn_file)
            return None
        logger.info(' Load abn-data from  %s' % self.abn_file)
        dt = self.dtype_abn()
        self._data_chem = np.loadtxt(self.abn_file, comments='#', dtype=dt)
        return self

    @staticmethod
    def dtype_abn():
        col_names = "zone dum1 dum2 dum3 " + ' '.join(StellaHydAbn.abn_elements)
        dt = np.dtype({'names': map(str.strip, col_names.split()), 'formats': np.repeat('f8', len(col_names))})
        return dt

    def write_abn(self, fname):
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
        if not self.is_load_chem:
            raise ValueError(' No loaded abn-data for %s')

        dum = 0.
        logger.info(' Write abn-data to %s' % fname)
        with open(fname, 'w') as f:
            # f.write('%d\n' % self.nzon_abn)
            for i in range(self.nzon_abn):
                s = '%4d  %10.3e %10.3e %10.3e' % (i+1, dum, dum, dum)
                for ename in StellaHydAbn.abn_elements:
                    s += ' %10.3e' % self.el(ename)[i]
                f.write('%s\n' % s)
        return os.path.isfile(fname)

    def plot_chem(self, x='m', ax=None, elements=None, xlim=None, ylim=None,
                  leg_loc=3, leg_ncol=4, lw=2, lntypes=None, is_save=False):
        if elements is None:
            elements = eve_elements
        if lntypes is None:
            lntypes = eve_lntypes

        is_new_plot = ax is None
        # setup figure
        if is_new_plot:
            plt.matplotlib.rcParams.update({'font.size': 14})
            fig = plt.figure(num=None, figsize=(12, 12), dpi=100, facecolor='w', edgecolor='k')

            gs1 = gridspec.GridSpec(1, 1)
            gs1.update(wspace=0.1, hspace=0.1, top=None, left=0.1, right=0.98)
            ax = fig.add_subplot(gs1[0, 0])

        is_x_lim = xlim is not None
        is_y_lim = ylim is not None

        ib = 0
        x_min = []
        x_max = []
        y_min = []
        y_max = []
        if x == 'r':
            x = self.r
            ax.set_xlabel(r'R [cm]')
        if x == 'lgr':
            x = self.r
            ax.set_xscale('log')
            ax.set_xlabel(r'R [cm]')
        else:
            x = self.mass
            ax.set_xlabel(r'M [$M_\odot$]')

        for el in elements:
            ib += 1
            y = self.el(el)
            ax.plot(x, y, label='%s' % el, color=eve_colors[el], ls=lntypes[el], linewidth=lw)

            if not is_x_lim:
                x_min.append(np.min(x))
                x_max.append(np.max(x))
            if not is_y_lim:
                y_min.append(np.min(y))
                y_max.append(np.max(y))

        if not is_x_lim:
            xlim = [np.min(x_min), np.max(x_max)]
        if not is_y_lim:
            ylim = [np.min(y_min), np.max(y_min)]

        ax.set_xlim(xlim)

        ax.set_yscale('log')
        ax.set_ylim(ylim)
        ax.set_ylabel(r'$log10(X_i)$')

        if is_new_plot:
            ax.legend(prop={'size': 9}, loc=leg_loc, ncol=leg_ncol, fancybox=False, frameon=True)
            # ax.legend(prop={'size': 9}, loc=3, ncol=4, fancybox=True, shadow=True)
            # plt.grid()
            # plt.show()

        if is_save:
            fsave = os.path.join(os.path.expanduser('~/'), 'chem_%s.pdf' % self.name)
            logger.info(" Save plot to %s " % fsave)
            ax.get_figure().savefig(fsave, bbox_inches='tight', format='pdf')

        return ax
