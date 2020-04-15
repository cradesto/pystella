import os
import logging
from itertools import islice
import numpy as np

from pystella.util.phys_var import phys
from pystella.rf import rad_func as rf

__author__ = 'bakl'

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

LEGEND_MASK_None = 0x00
LEGEND_MASK_Rho = 0x01
LEGEND_MASK_Vars = 0x10


class StellaTauDetail:
    """
    Reader for tau-files
    """

    def __init__(self, name, path='./'):
        """Creates a StellaTauDetail instance.  Required parameters:  name."""
        fname, ext = os.path.splitext(name)
        if ext is "tau":
            name = fname
        self._name = name
        self._path = path  # path to model files
        self._times = None
        self._nzon = None
        self._freq = None
        self._block_positions = None

    def __str__(self):
        return "%s, path: %s" % (self._name, self._path)

    def __repr__(self):
        return "%s, path: %s" % (self._name, self._path)
        # return "%s" % self.name

    @property
    def Nzon(self):
        return self._nzon

    @property
    def Freq(self):
        return self._freq

    @property
    def Wl(self):
        return rf.val_to_wl(self.Freq)
        # return np.array(map(lambda x: phys.c/x, self.Freq))
        # return phys.c / self.Freq

    @property
    def Wl2angs(self):
        return self.Wl * phys.cm_to_angs

    @property
    def Times(self):
        return self._times

    @property
    def Ntimes(self):
        return len(self.Times)

    @property
    def NFreq(self):
        return len(self.Freq)

    @property
    def Positions(self):
        return self._block_positions

    @property
    def Name(self):
        return self._name

    @property
    def FName(self):
        return os.path.expanduser(os.path.join(self._path, self._name + ".tau"))

    @staticmethod
    def parse_start_end_blocks(fname):
        """1
        Find the line indexes of the block data for all time moments if n is None.
        :param fname:
        :return: array( (s1,e1),...), array(times)
        """
        # PROPER   T = 0.00000000E+00
        # NZON        R14.  V8.    T5.          13.801        13.848
        start = 3
        indexes, times = [], []
        count = 0
        with open(fname, "r") as f:
            for i, line in enumerate(f):
                if line.lstrip().startswith("PROPER"):
                    count += 1
                    if count > 1:
                        indexes.append((start, i))
                        start = i+3
                    # times
                    t = float(line.split()[2])
                    times.append(t)
            else:
                indexes.append((start, i+1))  # last block

        return np.array(indexes), np.array(times)

    @staticmethod
    def parse_header(fname, nline):
        with open(fname) as f:
            try:
                line = next(islice(f, nline, nline+1))
            except StopIteration:
                print('Not header line in file {}'.format(fname))
        # NZON        R14.  V8.    T5.          13.801        13.848
        arr = line.split()
        freq = np.array(list(map(float, arr[6:])))
        freq = 10.**freq
        return freq

    @staticmethod
    def parse_block(fname, start, end, nfreq):
        colstr = " NZON    R14  V8    T5  " + ' '.join(list(map(str, range(1, nfreq+1))))
        cols = colstr.split()
        dt = np.dtype({'names': cols, 'formats': [np.float] * len(cols)})

        with open(fname, "rb") as f:
            from itertools import islice
            b = np.genfromtxt(islice(f, start-1, end), names=cols, dtype=dt)

        return b

    def load(self):
        """
        Load data from tau-file
        :return: self
        """
        se, times = self.parse_start_end_blocks(self.FName)
        if not len(se) == len(times):
            raise "Error in load: len(se)=len(times), " \
                  "but should be len(se)= {}, but len(times)= {}".format(len(se), len(times))
        if len(se) == 0:
            raise 'No data in file: {}'.format(self.FName)

        self._block_positions = se
        s, e = self._block_positions[0]
        self._nzon = e - s + 1
        self._times = times
        #  Freq
        self._freq = self.parse_header(self.FName, s - 2)

        logger.debug("Load data from  %s " % self.FName)
        return self

    def time_nearest(self, time):
        idx = (np.abs(self.Times - time)).argmin()
        return idx, self.Times[idx]

    def block_nearest(self, time):
        idx = self.time_nearest(time)
        block = self[idx]
        return BlockTau(self.Times[idx], block)

    def __getitem__(self, idx):
        b, e = self.Positions[idx]
        block = self.parse_block(self.FName, b, e, self.NFreq)
        return BlockTau(self.Times[idx], block)

    def evolution(self, var, nz):
        """Return the time evolution VAR in zone NZ"""
        if nz < 1 or nz > self.Nzon:
            raise ValueError('NZ should be in range: 0 - {}'.format(self.Nzon))

        if var not in vars(BlockTau):
            raise ValueError('var should be property of BlockTau, like Zon, M, R, Vel, T, Trad, Rho, Lum, Cappa, M')

        # x = []
        # y = []
        for idx, time in enumerate(self.Times):
            b = self.block_nearest(time)
            v = getattr(b, var)[nz]
            yield time, v
        # return False
        #     x.append(time)
        #     y.append(v)
        # return x, y

    def taus(self):
        """Compute tau for each moment of time
        Return: 2d-array[i,k], where i - time index, k - zone index.
        """
        taus = np.zeros((self.Ntimes, self.Nzon))
        for i, time in enumerate(self.Times):
            s = self[i]
            # tau = s.Tau
            taus[i, :] = s.Tau[:]
            # for k in range(self.Nzon-1, 1, -1):
            #     tau += s.Cappa[k] * s.Rho[k] * (s.R[k] - s.R[k-1])
            #     taus[i, k] = tau
        return taus

    def params_ph(self, cols=('R', 'T', 'V'), tau_ph=2./3.):
        is_str = isinstance(cols, str)
        if is_str:
            cols = [cols]
        res = {k: np.zeros(self.Ntimes) for k in ['time', 'zone'] + list(cols)}
        taus = self.taus()
        for i, time in enumerate(self.Times):
            s = self[i]
            kph = 1  # todo check the default value
            for k in range(self.Nzon-1, 1, -1):
                if taus[i][k] >= tau_ph:
                    kph = k
                    break
            res['time'][i] = time
            res['zone'][i] = kph
            for v in cols:
                res[v][i] = getattr(s, v)[kph]
        return res

    def vel_ph(self, tau_ph=2./3., z=0.):
        v_dat = self.params_ph(tau_ph=tau_ph, cols=['V'])
        t = v_dat['time'] * (1. + z)  # redshifted time
        v = v_dat['V']
        z = v_dat['zone']
        return t, v, z


class BlockTau:
    """
        Block tau-data for given time
        line: " NZON    R14  V8    T5   + freqs"
    """

    def __init__(self, time, block):
        self._time = time
        self._block = block

    @property
    def Time(self):
        return self._time

    @property
    def Nzon(self):
        return len(self.Zon)

    @property
    def Zon(self):
        return self._block['NZON']

    @property
    def R(self):
        return  self._block['R14'] * 1e14 # [cm]

    def Rsun(self):
        return self.R / phys.R_sun  # [Rsun]

    @property
    def V(self):
        return self._block['V8'] * 1e8  # [cm/s]

    @property
    def V8(self):
        return self._block['V8']  # [1000 km/s]

    @property
    def T(self):
        return self._block['T5'] * 1e5  # [K]

    @property
    def Tau(self):
        return self._block[4:, :]

    @property
    def Size(self):
        return self._block.size


# ==================================================================
def isfloat(value):
    try:
        float(value)
        return True
    except:
        return False


def plot_tau(ax, b, **kwargs):
    xlim = kwargs.get('xlim', None)
    ylim = kwargs.get('ylim', None)
    legmask = kwargs.get('legmask', 0x11)  # mask 0x11 - both, 0x01-rho, 0x10 - pars
    is_frameon = kwargs.get('is_frameon', False)
    islim = kwargs.get('islim', True)
    is_xlabel = kwargs.get('is_xlabel', True)
    is_yllabel = kwargs.get('is_yllabel', True)
    is_yrlabel = kwargs.get('is_yrlabel', True)
    is_day = kwargs.get('is_day', True)
    is_grid = kwargs.get('is_grid', False)
    rnorm = kwargs.get('rnorm', 'm')
    text_posy = kwargs.get('text_posy', 1.01)

    vnorm = kwargs.get('vnorm', 1.e8)
    lumnorm = kwargs.get('lumnorm', 1.e40)

    lw = 1.

    if rnorm == 'sun':
        rnorm = phys.R_sun
        x, xlabel = b.R / rnorm, r'Ejecta Radius, [$\mathtt{R}_\odot$]'
    elif rnorm == 'lgr':
        x, xlabel = b.R, r'Ejecta Radius, [cm]'
    elif rnorm == 'm':
        x, xlabel = b.M, r'Ejecta Mass [$\mathtt{M}_\odot$]'
    elif isfloat(rnorm):
        x, xlabel = b.R / float(rnorm), r'Ejecta Radius, [$\times 10^{%d}$ cm]' % int(np.log10(float(rnorm)))
    else:
        x, xlabel = b.M, r'Ejecta Mass [$\mathtt{M}_\odot$]'

    y = np.log10(b.Rho)
    if rnorm == 'lgr':
        ax.semilogx(x, y, label='Rho', color='black', ls="-", linewidth=lw)
    else:
        ax.plot(x, y, label='Rho', color='black', ls="-", linewidth=lw)

    if is_day:
        ax.text(.02, text_posy, r'$%5.2f^d$' % b.Time, horizontalalignment='left',
                transform=ax.transAxes)
    # ax.text(.5, 1.01, '%5.2f days' % b.Time, horizontalalignment='center', transform=ax.transAxes)

    if islim:
        if xlim is None:
            xlim = [min(x), max(x) * 1.2]
        if ylim is None:
            ylim = [np.min(y), np.max(y) + 2]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    if is_xlabel:
        ax.set_xlabel(xlabel)
    else:
        pass
        # ax.set_xticklabels([])

    # Right axe
    ax2 = ax.twinx()
    ax2.set_ylim((0., 9.9))

    if is_yllabel:
        ax.set_ylabel(r'$\log_{10}(\rho)$')
    else:
        ax.set_yticklabels([])

    if is_yrlabel:
        ax2.set_ylabel('T, Vel, Lum, Tau')
    else:
        ax2.set_yticklabels([])

    # y2 = b.Tau
    y2 = np.ma.masked_where(b.Tau <= 0., b.Tau)
    ax2.plot(x, y2, 'r-', label='tau')

    y2 = np.log10(b.T)
    ax2.plot(x, y2, 'g-', label=r'$\log(T)$')

    y2 = np.ma.log10(b.Lum) - np.log10(lumnorm)
    y22 = np.ma.log10(-b.Lum) - np.log10(lumnorm)
    ax2.plot(x, y2, color='orange', ls="-", label='Lum{0:d}'.format(int(np.log10(lumnorm))))
    ax2.plot(x, y22, color='brown', ls="--", label='-Lum{0:d}'.format(int(np.log10(lumnorm))))

    y2 = b.V / vnorm
    ax2.plot(x, y2, 'b-', label='V{0:d}'.format(int(np.log10(vnorm))))

    if legmask & LEGEND_MASK_Rho:
        ax.legend(loc=1, prop={'size': 8}, frameon=is_frameon)
    if legmask & LEGEND_MASK_Vars:
        ax2.legend(loc=1, prop={'size': 8}, ncol=3, frameon=is_frameon)

    if is_grid:
        ax2.grid(linestyle=':')
