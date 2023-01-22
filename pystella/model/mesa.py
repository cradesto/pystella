import os
import numpy as np
import logging


import pystella as ps
# from pystella.model.sn_eve import PreSN
# from pystella.util.phys_var import phys

logger = logging.getLogger(__name__)

try:
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    is_matplotlib = True
except ImportError:
    logging.debug('matplotlib failed to import', exc_info=True)
    is_matplotlib = False
    pass
# logger.setLevel(logging.INFO)
# logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)


__author__ = 'bakl'

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# 681     15
# 1.0d0 1.0d0 4.0d0 12.0d0 16.0d0 20.0d0 24.0d0 28.0d0 32.0d0 36.0d0 40.0d0 44.0d0 48.0d0 52.0d0 56.0d0
# 0.0d0 1.0d0 2.0d0  6.0d0  8.0d0 10.0d0 12.0d0 14.0d0 16.0d0 18.0d0 20.0d0 22.0d0 24.0d0 26.0d0 28.0d0

# mesa_elements = "NN H He C O Ne Mg Si S Ar Ca Ti Cr Fe Ni".split()
mesa_elements = ("h1", "he3", 'he4', "c12", "n14", "o16", "ne20", "mg24", "si28", 's32'
                , 'ar36', "ca40", 'ti44','cr48','cr56','fe52','fe54','fe56', "ni56")

mesa_el_colors = dict(h1="blue", he4="cyan", c12="darkorange", n14="yellow", 
                      o16="violet", ne20="green", mg24="skyblue", si28="olive",
                      s32="indigo", ar36="brown", ca40="purple", ti44="hotpink",
                      cr48="m", fe56='maroon', ni56='red')
for k in mesa_elements:
    if k not in mesa_el_colors:
        mesa_el_colors[k] = 'grey'

mesa_el_lntypes = dict((k, '--') for k, v in mesa_el_colors.items())  # no y-shift
mesa_el_lntypes['h1'] = '-'
mesa_el_lntypes['he4'] = '-'
mesa_el_lntypes['o16'] = '-'
mesa_el_lntypes['c12'] = '-'
mesa_el_lntypes['ni56'] = '-'


class Mesa:

    def __init__(self, name, elements=mesa_elements):
        """Creates a Problem instance.  It's initial conditions for MESA. Required parameters:  name."""
        self._profile_file = None
        self._profile = None
        self._header = None
        self._name = name
        self._elements = elements

    @property
    def Name(self):
        return self._name

    @property
    def Nzone(self):
        return self.profile.shape[0]

    # @property
    # def elements(self):
    #     return self._elements

    @property
    def profile(self):
        return self._profile

    @property
    def header(self):
        return self._header

    @property
    def profile_file(self):
        return self._profile_file

    @property
    def m_init(self):
        return float(self.header['initial_mass'])

    @property
    def m_core(self):
        return float(self.header['neutron_rich_core_mass'])

    @property
    def m_tot(self):
         return float(self.header['star_mass'])

    @property
    def Elements(self):
        """Elements"""
        els = []
        keys = self.profile.dtype.names
        for el in mesa_elements:
            if el in keys:
                els.append(el)
        return els

    @property
    def is_load(self):
        return self._header is not None and self._header is not None

    def load(self, fname):
        if not os.path.isfile(fname):
            logger.error(f' No mesa-data file for {fname}')
            return None
        self._profile_file = fname
        logger.info('Load profile data from  %s' % self.profile_file)

        def read_head(fn):
            keys = fn.readline().split()
            values = fn.readline().split()
            return dict(zip(keys,values))

        # dic_header = {}
        # Read header
        with open(fname, 'r') as fn: 
            line = fn.readline()
        #     print(line)
            ncol = len(line.split())
            print('Ncol= {}'.format(ncol))
            
            # header1
            header1 = read_head(fn)
            line, line = fn.readline(), fn.readline() # empty line & nums
            cols = fn.readline().split()
            
        dtype = dtype={'names': cols, 'formats': [int] + [float]*(len(cols)-1)}    
        data = np.loadtxt(fname, skiprows=6, dtype=dtype)    
        # data = np.flip(data)

        self._profile = np.flip(data)
        self._header = header1
        return self

    def el(self, el):
        if not self.is_load:
            raise Exception("The data has not been loaded. Check and load from {}".format(self._profile_file))
        if el not in self.profile.dtype.names:
            raise ValueError("There is no such element [%s]." % el)
        return self.profile[el]

    def to_presn_stella_el(self, is_mcut=True, chem_lim=1e-15):
        if not self.is_load:
            raise Exception("The data has not been loaded. Check and load from {}".format(self._profile_file))
            
        # d = self.profile
        if len(self.profile) <= 0:
            raise ValueError("There are no data i. ")
        nzonDvd = self.profile.shape[0]
        print(f'nzonDvd= {nzonDvd}')
        # Hyd
        #     print(f'ADD m_core= {m_core}')
        # Set cut inner core
        mass_mdl = m_cut = self.profile['mass'].flatten()
        cut_core = mass_mdl > 0.
        # Set cutting core
        if is_mcut:
            cut_core = mass_mdl > self.m_core
            m_cut = mass_mdl[cut_core]

        # Create PreSN
        presn = ps.PreSN(self.Name, len(m_cut))
        map_col = {'logRho':('Rho',1.), 'radius': ('R', ps.phys.R_sun), 'mass':('M',ps.phys.M_sun), 'temperature':('T', 1.), 'velocity':('V', -1.)}
        for cfrom, a in map_col.items():
            cto, unit = a
            # print(f'{cfrom} => {cto}')
            val = self.profile[cfrom].flatten()
            if 'log' in cfrom:
                val = 10.**val
            val *= unit
            if cto == 'M':
                # print(f'ADD m_core= {m_core}')
                # val += m_core*unit
                # print(m_core, val)
                pass
            # cut 
            val = val[cut_core]
            presn.set_hyd(cto, val)

        # Abn
    #     map_el = 'h1' 'he3' 'he4' 'c12' 'n14' 'o16' 'ne20' 'mg24' 'si28' 's32' 'ar36' 'ca40' 'ti44' 'cr48' 'cr56' 'fe52' 'fe54' 'fe56' 'ni56'
        map_el = {'H':['h1'], 'He':['he3','he4'], 'C':['c12'], 'N':['n14'], 'O':['o16']
                , 'Ne':['ne20'], 'Mg':['mg24'], 'Si':['si28'], 'S':['s32'], 'Ca':['ca40']
                #  , 'Fe':['fe56']
                , 'Fe':['ti44','cr48','cr56','fe52','fe54','fe56']
                , 'Ni56':['ni56']
                }
        for e in presn.Elements:
            val = np.zeros(presn.nzon)
            if e in map_el:
                xi = np.zeros(nzonDvd)
                s = ''
                for i, edvd in enumerate(map_el[e]):                                                                                                                                         
                    # print('  add edvd= ', edvd)
                    xi += self.profile[edvd]
                    s += f' {edvd}'
                # cut
                xi[xi<chem_lim] = 0.
                val = xi[cut_core]
                print(f'Set {e} from {s}')
            presn.set_chem(e, val)

        print("{} nzon={}  m_init= {:.2f}  m_tot= {:.2f}  m_core= {:.2f}".
              format(self.Name, presn.nzon, self.m_init, self.m_tot, self.m_core))

        return presn

    def to_presn(self, is_mcut=False, chem_lim=1e-15):
        if not self.is_load:
            raise Exception("The data has not been loaded. Check and load from {}".format(self._profile_file))
            
        # d = self.profile
        if len(self.profile) <= 0:
            raise ValueError("There are no data i. ")
        nzonDvd = self.profile.shape[0]
        print(f'nzonDvd= {nzonDvd}')
        # Hyd
        #     print(f'ADD m_core= {m_core}')
        # Set cut inner core
        mass_mdl = m_cut = self.profile['mass'].flatten()
        cut_core = mass_mdl > 0.
        # Set cutting core
        if is_mcut:
            cut_core = mass_mdl > self.m_core
            m_cut = mass_mdl[cut_core]

        # Create PreSN
        presn = ps.PreSN(self.Name, len(m_cut), elements=self.Elements)
        map_col = {'logRho':('Rho',1.), 'radius': ('R', ps.phys.R_sun), 'mass':('M',ps.phys.M_sun), 'temperature':('T', 1.), 'velocity':('V', -1.)}
        for cfrom, a in map_col.items():
            cto, unit = a
            # print(f'{cfrom} => {cto}')
            val = self.profile[cfrom].flatten()
            if 'log' in cfrom:
                val = 10.**val
            val *= unit
            if cto == 'M':
                # print(f'ADD m_core= {m_core}')
                # val += m_core*unit
                # print(m_core, val)
                pass
            # cut 
            val = val[cut_core]
            presn.set_hyd(cto, val)

        # Abn
    #     map_el = 'h1' 'he3' 'he4' 'c12' 'n14' 'o16' 'ne20' 'mg24' 'si28' 's32' 'ar36' 'ca40' 'ti44' 'cr48' 'cr56' 'fe52' 'fe54' 'fe56' 'ni56'
        # map_el = {'H':['h1'], 'He':['he3','he4'], 'C':['c12'], 'N':['n14'], 'O':['o16']
        #         , 'Ne':['ne20'], 'Mg':['mg24'], 'Si':['si28'], 'S':['s32'], 'Ca':['ca40']
        #         #  , 'Fe':['fe56']
        #         , 'Fe':['ti44','cr48','cr56','fe52','fe54','fe56']
        #         , 'Ni56':['ni56']
        #         }
        for e in self.Elements:
            xi = self.profile[e]
            # cut
            xi[xi<chem_lim] = 0.
            val = xi[cut_core]
            print(f'Set {e}')
            presn.set_chem(e, val)

        print("{} nzon={}  m_init= {:.2f}  m_tot= {:.2f}  m_core= {:.2f}".
              format(self.Name, presn.nzon, self.m_init, self.m_tot, self.m_core))

        return presn        

    # Plotting
    @staticmethod
    def plot_chem(presn, x='m', ax=None, xlim=None, ylim=None, **kwargs):
        elements = kwargs.get('elements', mesa_elements)
        lntypes = kwargs.get('lntypes', mesa_el_lntypes)
        if isinstance(lntypes, str):
            lntypes = {el: lntypes for el in elements}
        colors = kwargs.get('colors', mesa_el_colors)
        loc = kwargs.get('leg_loc', 3)
        font_size = kwargs.get('font_size', 14)
        leg_ncol = kwargs.get('leg_ncol', 4)
        lw = kwargs.get('lw', 2)
        is_save = kwargs.get('is_save', False)
        alpha = kwargs.get('alpha', 1.)
        figsize = kwargs.get('figsize', (8, 8))


        is_new_plot = ax is None
        # setup figure
        if is_new_plot:
            plt.matplotlib.rcParams.update({'font.size': font_size})
            fig = plt.figure(num=None, figsize=figsize, dpi=100, facecolor='w', edgecolor='k')

            gs1 = gridspec.GridSpec(1, 1)
            gs1.update(wspace=0.1, hspace=0.1, top=0.97, left=0.12, right=0.87)
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
            x = presn.r
        elif x == 'm':
            x = presn.m / ps.phys.M_sun
        else:
            x = presn.r

        y_min = []
        y_max = []
        for el in elements:
            y = presn.el(el)
            # y = np.log10(self.el(el))
            ax.semilogy(x, y, label='%s' % el, color=colors[el], ls=lntypes[el], linewidth=lw, alpha=alpha)
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
            fsave = os.path.join(os.path.expanduser('~/'), 'chem_%s.pdf' % presn._name)
            logger.info(" Save plot to %s " % fsave)
            ax.get_figure().savefig(fsave, bbox_inches='tight', format='pdf')

        return ax



