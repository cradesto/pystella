import os
import numpy as np
import logging

import matplotlib.pyplot as plt
from matplotlib import gridspec

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

__author__ = 'bakl'

eve_elements = map(str.strip, "Ni56 H He C N O Ne Na  Mg  Al  Si  S  Ar  Ca  Fe  Ni".split())
eve_colors = dict(Ni56="black", H="blue", He="magenta", C="coral", N="red",
                  O="cyan", Ne="green", Na="sandybrown",
                  Mg="skyblue", Si="violet", Al="lime",
                  S="indigo", Ar="olive", Ca="purple",
                  Fe='maroon', Ni='steelblue')
eve_lntypes = dict((k, '-') for k, v in eve_colors.items())  # no y-shift
eve_lntypes['O'] = '--'
eve_lntypes['Ni'] = '--'
eve_lntypes['Si'] = '--'
eve_lntypes['Fe'] = '--'


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
    def mass(self):
        """Mass"""
        return self._data['mass']

    @property
    def data(self):
        """Full data"""
        return self._data

    @property
    def rho_file(self):
        return os.path.join(self.path, self.name + '.rho')

    @property
    def is_rho_data(self):
        ext = ['rho']
        return any(map(os.path.isfile, [os.path.join(self.path, self.name + '.' + e) for e in ext]))

    def rho_load(self):
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

    def plot_chem(self, elements=None, xlim=None, ylim=None, is_save=False):
        if elements is None:
            elements = eve_elements

        # setup figure
        plt.matplotlib.rcParams.update({'font.size': 14})
        fig = plt.figure(num=None, figsize=(7, 7), dpi=100, facecolor='w', edgecolor='k')

        gs1 = gridspec.GridSpec(1, 1)
        gs1.update(wspace=0.1, hspace=0.1, top=None, left=0.1, right=0.98)
        ax = fig.add_subplot(gs1[0, 0])

        is_x_lim = xlim is not None
        is_y_lim = ylim is not None

        lw = 2.
        ib = 0
        x_max = []
        y_mid = []
        x = self.mass
        for el in elements:
            ib += 1
            y = self._data[el]
            ax.plot(x, y, label='%s' % el, color=eve_colors[el], ls=eve_lntypes[el], linewidth=lw)

            if not is_x_lim:
                x_max.append(np.max(x))
            if not is_y_lim:
                y_mid.append(np.min(y))

        if not is_x_lim:
            xlim = [1.5, np.max(x_max)]
        if not is_y_lim:
            ylim = [np.min(y_mid), 0.]

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_ylabel(r'$log10(X_i)$')
        ax.set_xlabel(r'M [$M_\odot$]')

        ax.legend(prop={'size': 9}, loc=3, ncol=4, fancybox=False, frameon=True)
        # ax.legend(prop={'size': 9}, loc=3, ncol=4, fancybox=True, shadow=True)
        plt.grid()
        plt.show()

        if is_save:
            d = os.path.expanduser('~/')
            fsave = os.path.join(d, 'chem_%s.pdf' % self.name)
            logger.info(" Save plot to %s " % fsave)
            fig.savefig(fsave, bbox_inches='tight', format='pdf')

        return fig
