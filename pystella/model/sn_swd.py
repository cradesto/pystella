import numpy as np
import os
import logging

from pystella.util.phys_var import phys

__author__ = 'bakl'

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class StellaShockWaveDetail:
    """
    Reader for swd-files
    """

    def __init__(self, name, path='./'):
        """Creates a StellaShockWaveDetail instance.  Required parameters:  name."""
        fname, ext = os.path.splitext(name)
        if ext is "swd":
            name = fname
        self._name = name
        self._path = path  # path to model files
        self._times = None
        self._nzon = None
        self._data = None

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
    def Data(self):
        return self._data

    @property
    def Name(self):
        return self._name

    def load(self):
        """
        Load datd from swd-file
        :return:
        """
        fname = os.path.expanduser(os.path.join(self._path, self._name + ".swd"))
        # tout, Km,log10(AMPR),log10(UR*Ry(Km)),Uy(Km)*1.D+6/(UTIME*CRAP),
        # log10(max(UTP*Ty(Km),1.d0)),log10(max(UTP*TpRAD,1.d0)),
        # PLLOG,PLOG,QVLOG,log10(max(eng,1.d-50)),Flum*1.d-40,WRKX(Km);
        colstr = "tday km lgM lgR14 V8 lgT lgTrad lgDm6 lgP7  lgQv lgEng Flum40 cap"
        cols = [s.strip() for s in colstr.split()]
        dt = np.dtype({'names': cols, 'formats': [np.float] * len(cols)})
        data = np.loadtxt(fname, dtype=dt)

        self._nzon = int(np.max(data['km']))
        times = np.unique(data['tday'])
        if data['tday'][0] != 0.:
            times = np.delete(times, np.where(times == 0.))
        self._times = times

        # self.ntimes = len(self.times)   # len(data['tday']) // self.nzon
        self._data = data
        logger.debug("Read data from  %s " % fname)
        return self

    def time_nearest(self, time):
        idx = (np.abs(self.Times - time)).argmin()
        return idx, self.Times[idx]

    def block_nearest(self, time):
        idx = (np.abs(self._data['tday'] - time)).argmin()
        b = idx
        e = b + self._nzon
        block = self._data[:][b:e]
        return BlockSwd(self._data['tday'][idx], block)

    def __getitem__(self, idx):
        b = idx*self._nzon
        e = b + self._nzon
        block = self._data[:][b:e]
        return BlockSwd(self._data['tday'][idx], block)


class BlockSwd:
    """
        Block swd-data for given time
        line: "tday km lgM lgR14 V8 lgT lgTrad lgDm6 lgP7  lgQv lgEng Flum40 cap"
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
        return self._block['km']

    @property
    def R(self):
        return 10. ** self._block['lgR14']  # [cm]

    def Rsun(self):
        return self.R / phys.R_sun  # [Rsun]

    @property
    def Vel(self):
        return self._block['V8'] * 1e8  # [cm/s]

    @property
    def Vel8(self):
        return self._block['V8']  # [1000 km/s]

    @property
    def T(self):
        return 10. ** self._block['lgT']  # [K]

    @property
    def Trad(self):
        return 10. ** self._block['lgTrad']  # [K]

    @property
    def Rho(self):
        return 10. ** (self._block['lgDm6'] - 6.)  # [g*cm-3]

    @property
    def Lum(self):
        return self._block['Flum40'] * 1e40  # [erg/s]  ??

    @property
    def Cappa(self):
        return self._block['cap']  # [cm^2/g]   ??

    @property
    def Mtot(self):
        return 10 ** self._block['lgM'][0]  # [Msun]   ??

    @property
    def M(self):
        return self.Mtot - 10 ** self._block['lgM']  # [Msun]   ??

    @property
    def Tau(self):
        tau = np.zeros(self.Nzon)
        for i in range(self.Nzon - 2, 0, -1):
            tau[i] = tau[i + 1] + self.Cappa[i] * self.Rho[i] * (self.R[i + 1] - self.R[i])
        tau[self.Nzon - 1] = tau[self.Nzon - 2] / 2.
        return tau


# ==================================================================
def isfloat(value):
    try:
        float(value)
        return True
    except:
        return False


def plot_swd(ax, b, **kwargs):
    xlim = kwargs.get('xlim', None)
    ylim = kwargs.get('ylim', None)
    is_legend = kwargs.get('is_legend', True)
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
        x, xlabel = b.M,  r'Ejecta Mass [$\mathtt{M}_\odot$]'
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
        ax.text(.05, text_posy, '%5.2f days' % b.Time, horizontalalignment='left',
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

    y2 = b.Vel / vnorm
    ax2.plot(x, y2, 'b-', label='V{0:d}'.format(int(np.log10(vnorm))))

    if is_legend:
        ax.legend(loc=2, prop={'size': 8})
        ax2.legend(loc=1, prop={'size': 8}, ncol=2)

    if is_grid:
        ax.grid(linestyle=':')
