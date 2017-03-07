import os
import numpy as np
import logging
import matplotlib.pyplot as plt
from matplotlib import gridspec

from pystella.model.sn_eve import PreSN
from pystella.util.phys_var import phys

__author__ = 'bakl'

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# 681     15
# 1.0d0 1.0d0 4.0d0 12.0d0 16.0d0 20.0d0 24.0d0 28.0d0 32.0d0 36.0d0 40.0d0 44.0d0 48.0d0 52.0d0 56.0d0
# 0.0d0 1.0d0 2.0d0 6.0d0 8.0d0 10.0d0 12.0d0 14.0d0 16.0d0 18.0d0 20.0d0 22.0d0 24.0d0 26.0d0 28.0d0

snec_elements = "NN H He C O Ne Mg Si S Ar Ca Ti Cr Fe Ni".split()
snec_el_colors = dict(NN="yellow", H="blue", He="cyan", C="darkorange",
                      O="violet", Ne="green", Mg="skyblue", Si="olive",
                      S="indigo", Ar="brown", Ca="purple", Ti="hotpink",
                      Cr="m", Fe='maroon', Ni='magenta')

snec_el_lntypes = dict((k, '--') for k, v in snec_el_colors.items())  # no y-shift
snec_el_lntypes['H'] = '-'
snec_el_lntypes['He'] = '-'
snec_el_lntypes['O'] = '-'
snec_el_lntypes['C'] = '-'
snec_el_lntypes['Ni56'] = '-'

snec_profile_cols = "i M R T Rho V".split()


class Problem:
    def __init__(self, name):
        """Creates a Problem instance.  It's initial conditions for SNEC. Required parameters:  name."""
        self._name = name
        self._chem_file = None
        self._chem = None
        self._profile_file = None
        self._profile = None

    @property
    def name(self):
        return self._name

    @property
    def chem_file(self):
        return self._chem_file

    @property
    def r(self):
        """radius"""
        return self._chem['R']

    @property
    def nzon(self):
        """Number of zones"""
        return len(self.r)

    @property
    def m(self):
        """Mass"""
        return self._chem['M']

    @property
    def is_chem_load(self):
        """Check if data has been loaded."""
        return self._chem is not None

    @property
    def chem(self):
        """Full data"""
        return self._chem

    def el(self, el):
        if el not in snec_elements:
            raise ValueError("There is no such element [%s]." % el)
        if not self.is_chem_load:
            raise Exception("SNEC chem-data has not been loaded. Check and load from %s" % self._chem_file)
        return self._chem[el]

    def load_chem(self, fname):
        if not os.path.isfile(fname):
            logger.error(' No snec chem-data for %s' % fname)
            return None
        self._chem_file = fname
        logger.info('Load chemical data from  %s' % self.chem_file)

        names = ['M', 'R'] + snec_elements
        print("Names: %s" % ' '.join(names))
        dtype = np.dtype({'names': names, 'formats': [np.float64] * len(names)})
        self._chem = np.loadtxt(fname, skiprows=3, dtype=dtype, comments='#')
        return self

    # Profile structure
    @property
    def profile_file(self):
        return self._profile_file

    @property
    def profile(self):
        """Full data"""
        return self._profile

    @property
    def is_profile_load(self):
        """Check if data has been loaded."""
        return self._profile is not None

    @property
    def pmass(self):
        """Mass"""
        return self.hyd('M')

    @property
    def pradius(self):
        """Radius"""
        return self.hyd('R')

    @property
    def ptemp(self):
        """Temperature"""
        return self.hyd('T')

    @property
    def prho(self):
        """Density"""
        return self.hyd('Rho')

    @property
    def pvel(self):
        """Velocity"""
        return self.hyd('V')

    def hyd(self, v):
        """Hydro data"""
        if v not in snec_profile_cols:
            raise ValueError("There is no information about the parameter [%s]. You should set it." % v)
        return self._profile[v]

    def load_profile(self, fname):
        if not os.path.isfile(fname):
            logger.error(' No snec chem-data for %s' % fname)
            return None
        self._profile_file = fname
        logger.info('Load profile data from  %s' % self.profile_file)

        use_cols = list(range(0, len(snec_profile_cols)))
        dtype = np.dtype({'names': snec_profile_cols, 'formats': [np.float64] * len(snec_profile_cols)})
        self._profile = np.loadtxt(self.profile_file, skiprows=1, dtype=dtype, usecols=use_cols)
        return self

    # Plotting
    def plot_chem(self, x='m', ax=None, xlim=None, ylim=None, **kwargs):
        elements = kwargs.get('elements', snec_elements)
        lntypes = kwargs.get('lntypes', snec_el_lntypes)
        colors = kwargs.get('colors', snec_el_colors)
        loc = kwargs.get('leg_loc', 3)
        font_size = kwargs.get('font_size', 14)
        leg_ncol = kwargs.get('leg_ncol', 4)
        lw = kwargs.get('lw', 2)
        is_save = kwargs.get('is_save', False)

        is_new_plot = ax is None
        # setup figure
        if is_new_plot:
            plt.matplotlib.rcParams.update({'font.size': font_size})
            fig = plt.figure(num=None, figsize=(12, 12), dpi=100, facecolor='w', edgecolor='k')

            gs1 = gridspec.GridSpec(1, 1)
            gs1.update(wspace=0.1, hspace=0.1, top=None, left=0.1, right=0.98)
            ax = fig.add_subplot(gs1[0, 0])

            if x == 'r':
                ax.set_xlabel(r'R [cm]')
            elif x == 'm':
                ax.set_xlabel(r'M [$M_\odot$]')
            else:
                ax.set_xlabel(r'R [cm]')
                ax.set_xscale('log')

        is_x_lim = xlim is not None
        is_y_lim = ylim is not None

        if x == 'r':
            x = self.r
        elif x == 'm':
            x = self.m / phys.M_sun
        else:
            x = self.r

        y_min = []
        y_max = []
        for el in elements:
            y = self.el(el)
            # y = np.log10(self.el(el))
            ax.semilogy(x, y, label='%s' % el, color=colors[el], ls=lntypes[el], linewidth=lw)
            # ax.plot(x, y, label='%s' % el, color=colors[el], marker='o', ls=':', markersize=3)

            if not is_y_lim:
                y_min.append(np.min(y))
                y_max.append(np.max(y))

        if not is_y_lim:
            ylim = [np.min(y_min), np.max(y_min)]

        if not is_x_lim:
            xlim = np.min(x), np.max(x)

        ax.set_xlim(xlim)

        # ax.set_yscale('log')
        # ax.set_ylim(ylim)
        # ax.set_ylabel(r'$log10(X_i)$')
        ax.set_ylim(ylim)
        ax.set_ylabel(r'$X_i$')

        if is_new_plot:
            ax.legend(prop={'size': 9}, loc=loc, ncol=leg_ncol, fancybox=False, frameon=True)

        if is_save:
            fsave = os.path.join(os.path.expanduser('~/'), 'chem_%s.pdf' % self._name)
            logger.info(" Save plot to %s " % fsave)
            ax.get_figure().savefig(fsave, bbox_inches='tight', format='pdf')

        return ax


def to_presn(p):
    if not p.is_profile_load:
        raise ValueError("There are no data in SNEC problem. "
                         "Probably, You should run: load_profile and load_chem.")
    presn = PreSN(p.name, p.nzon)
    col_map = {'R', 'M', 'T', 'Rho', 'V'}
    for v in col_map:
        presn.set_hyd(v, p.hyd(v))

    for e in presn.presn_elements:
        if e in snec_elements:
            presn.set_chem(e, p.el(e))
        else:
            presn.set_chem(e, np.zeros(presn.nzon))

    # todo check with Viktoriya: in SNEC Ni used Ni as Ni56
    presn.set_chem('Ni56', presn.el('Ni'))

    return presn
