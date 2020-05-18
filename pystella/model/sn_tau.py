import os
import logging
from itertools import islice
import numpy as np

from pystella.util.phys_var import phys

# from pystella.rf import rad_func as rf

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
    col_zon = 'zon'

    def __init__(self, name, path='./'):
        """Creates a StellaTauDetail instance.  Required parameters:  name."""
        fname, ext = os.path.splitext(name)
        if ext is "tau":
            name = fname
        self._name = name
        self._path = path  # path to model files
        self._times = None
        self._nzon = None
        self._block_positions = None
        self._blocks = None

    def __str__(self):
        return "%s, path: %s" % (self._name, self._path)

    def __repr__(self):
        return "%s, path: %s" % (self._name, self._path)
        # return "%s" % self.name

    @property
    def Nzon(self):
        return self._nzon

    @property
    def Times(self):
        return self._times

    @property
    def Ntimes(self):
        return len(self.Times)

    @property
    def TimesObs(self):
        res = []
        for block in self:
            res.append(block.TimeObs)
        return np.array(res)

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
    def parse_start_end_blocks(fname, tnorm=24.*3600.):
        """
        Find the line indexes of the block data for all time moments if n is None.
        :param fname: filename
        :param tnorm: normalisation of time. Default: 86400, days
        :return: array( (s1,e1),...), array(times)
        """
        # PROPER   T = 0.00000000E+00
        # NZON        R14.  V8.    T5.          13.801        13.848
        start = 1
        indexes, times = [], []
        count = 0
        with open(fname, "r") as f:
            for i, line in enumerate(f):
                if line.lstrip().startswith("PROPER"):
                    count += 1
                    if count > 1:
                        indexes.append((start, i - 1))
                        start = i + 1
                    # times
                    t = float(line.split()[2])
                    times.append(t / tnorm)  # in days
            else:
                indexes.append((start, i))  # last block

        # check
        check = [e - b for b, e in indexes]
        if not len(np.unique(check)) == 1:
            print("      Start    End     NZon-1")
            for i, (b, e) in enumerate(indexes, 1):
                print("{:3d}  {:5d}  - {:5d} = {:5d} ".format(i, b, e, e - b))
            print("Diff zone number: ", np.unique(check))
            raise ValueError("There are different number of zone in file {}".format(fname))

        return np.array(indexes), np.array(times)

    @staticmethod
    def parse_header(fname, nline):
        with open(fname, "r") as f:
            try:
                line = next(islice(f, int(nline), int(nline + 1)))
            except StopIteration:
                print('Not header line in file {}'.format(fname))
        # NZON        R14.  V8.    T5.          13.801        13.848
        arr = line.split()
        freq = np.array(list(map(float, arr[6:])))
        freq = 10. ** freq
        return freq

    @staticmethod
    def parse_block(fname, start, nzon):
        from itertools import islice

        #  Freq
        freq = StellaTauDetail.parse_header(fname, start)

        b = []
        with open(fname, "r") as f:
            # b = np.genfromtxt(islice(f, start-1, start-1+nzon), names=cols, dtype=dt)
            for i, line in enumerate(islice(f, int(start + 1), int(start + 1 + nzon))):
                items = [float(v) for v in line.split()]
                b.append(items)

        # print(start, start-1+nzon, nzon)
        # b = np.loadtxt(fname, dtype=dt, skiprows=start-1, max_rows=nzon)

        return np.array(b), freq

    def load(self, is_info=False):
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
        self._blocks = [None] * len(times)
        s, e = self._block_positions[0]
        self._nzon = e - s
        self._times = times

        logger.debug("Load data from  %s " % self.FName)
        if is_info:
            logger.info('Loaded Tau from {} \n    '
                        'nzone= {} ntimes= {}'.format(self.FName, self.Nzon, self.Ntimes))

        return self

    def time_nearest(self, time):
        idx = (np.abs(self.Times - time)).argmin()
        return idx, self.Times[idx]

    def block_nearest(self, time):
        idx, t = self.time_nearest(time)
        block = self[idx]
        return block

    def __getitem__(self, idx):
        if idx in self._blocks:
            return self._blocks[idx]

        b, e = self.Positions[idx]
        block, freq = self.parse_block(self.FName, b, self.Nzon)
        # block = self.parse_block(self.FName, b, e, self.NFreq)
        block = BlockTau(self.Times[idx], freq, block)
        self._blocks[idx] = block
        return block

    def __iter__(self):
        for i, t in enumerate(self.Times):
            yield self[i]

    def params_ph(self, pars, moments, tau_ph=2. / 3.):
        """
        Compute  photosphere as  Func(nu). Maybe: R, V, V8, T
        :param pars: photosphere parameters. Maybe: R, V, V8, T. Example: [logR:V8:T]
        :param moments: time moments
        :param tau_ph:  the optical depth at the photosphere. Default: 2/3
        :return: yield the Param(t, nu)
        """

        res = {k: [] for k in [self.col_zon] + list(pars)}
        for j, time in enumerate(moments):
            # b = self.block_nearest(time)
            idx = (np.abs(self.TimesObs - time)).argmin()
            b = self[idx]
            idxs = b.ph_indexes(tau_ph=tau_ph)
            if b.NFreq != len(idxs):
                raise ValueError("Error in photo. indexes: len(wl)= {}  len(idxs)= {}".format(b.NFreq, len(idxs)))
            # Add the zone of photosphere
            res[self.col_zon].append((b.Time, b.Freq, idxs))

            for i, p in enumerate(pars, 1):
                arr = getattr(b, p)
                y = [arr[idx] for idx in idxs]
                res[p].append((b.Time, b.Freq, y))
        return res

    @staticmethod
    def data_save(fwrite_prefix, tau_data, pars, fmt_nu='6.3f', fmt_def='6.3f'):
        pars = [StellaTauDetail.col_zon] + pars
        formats = {
            StellaTauDetail.col_zon: '6d',
            'T': '6.1f',
            'V': '6.3e',
            'V8': '6.3f',
            'R': '6.3e',
        }

        for i, p in enumerate(pars):
            is_log = p.startswith('log')
            p_data = p.replace('log', '') if is_log else p
            fwrite = "{}{}.dat".format(fwrite_prefix, p_data)
            with open(fwrite, 'w') as f:
                for j, (t, freq, y) in enumerate(tau_data[p_data]):
                    n = int(5 - np.log10(max(1e-03, abs(t))))  # label format
                    freq = np.log10(freq)
                    s = "Time Nfreq: {:.{}f}   {}".format(t, n, len(freq))
                    print(s, file=f)
                    s = "  ".join(["{:{}}".format(nu, fmt_nu) for nu in freq])
                    print(s, file=f)
                    fmt = formats.get(p_data, fmt_def)
                    if is_log:
                        y = np.log10(y)
                    s = "  ".join(["{:{}}".format(v, fmt) for v in y])
                    print(s, file=f)
            print(' Saved {} to {}'.format(p, fwrite))


class BlockTau:
    """
        Block tau-data for given time
        line: " NZON    R14  V8    T5   + freqs"
    """
    cols = {'Zon': 0, 'R14': 1, 'V8': 2, 'T5': 3}

    def __init__(self, time, freq, block):
        self._time = time
        self._freq = freq
        self._block = block

    @property
    def Time(self):
        return self._time

    @property
    def TimeObs(self):
        t_prop = self.Time
        tretard = self.R[self.Nzon-1] / phys.c / (24.*3600.)  # to days
        t_ob = t_prop - tretard
        return t_ob

    @property
    def Nzon(self):
        return len(self.Zon)

    @property
    def Freq(self):
        return self._freq

    @property
    def NFreq(self):
        return len(self.Freq)

    @property
    def Wl(self):
        return phys.c / self.Freq

    @property
    def Wl2angs(self):
        return self.Wl * phys.cm_to_angs

    @property
    def Zon(self):
        return self.ColByName('Zon')

    @property
    def R(self):
        return self.ColByName('R14') * 1e14  # [cm]

    @property
    def Rsun(self):
        return self.R / phys.R_sun  # [Rsun]

    @property
    def V(self):
        return self.V8 * 1e8  # [cm/s]

    @property
    def V8(self):
        return self.ColByName('V8')  # [1000 km/s]

    @property
    def T(self):
        return self.ColByName('T5') * 1e5  # [K]

    @property
    def Tau(self):
        return self._block[:, len(self.cols):]

    @property
    def Size(self):
        return self._block.size

    def Col(self, idx):
        return self._block[:, idx]

    def ColByName(self, name):
        return self._block[:, self.cols[name]]

    def ph_indexes(self, tau_ph=2. / 3.):
        tau = self.Tau
        nzon, nfreq = tau.shape
        idxs = np.zeros(nfreq, dtype=int)
        for k in range(nfreq):
            array = tau[:, k]
            idxs[k] = (np.abs(array - tau_ph)).argmin()
        return idxs


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
