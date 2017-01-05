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
        fname = os.path.abspath(os.path.join(self._path, self._name + ".swd"))
        # tout, Km,log10(AMPR),log10(UR*Ry(Km)),Uy(Km)*1.D+6/(UTIME*CRAP),
        # log10(max(UTP*Ty(Km),1.d0)),log10(max(UTP*TpRAD,1.d0)),
        # PLLOG,PLOG,QVLOG,log10(max(eng,1.d-50)),Flum*1.d-40,WRKX(Km);
        colstr = "tday km lgM lgR14 V8 lgT lgTrad lgDm6 lgP7  lgQv lgEng Flum40 cap"
        cols = map(str.strip, colstr.split())
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
        return idx, self.Times[time]

    def block_nearest(self, time):
        idx = (np.abs(self._data['tday'] - time)).argmin()
        b = idx
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
def plot_swd(ax, b, **kwargs):
    lw = 1.
    x, xlabel = b.M, 'Ejecta Mass [Msun]'

    if 'rnorm' in kwargs:
        rnorm = kwargs['rnorm']

        if rnorm == 'sun':
            rnorm = phys.R_sun
            x, xlabel = b.R / rnorm, 'Ejecta Radius, [Rsun]'
        elif rnorm == 'm':
            x, xlabel = b.M, 'Ejecta Mass [Msun]'
        else:
            x, xlabel = b.R / rnorm, 'Ejecta Radius, [x10^%d cm]' % int(np.log10(rnorm))

    y = np.log10(b.Rho)
    ax.plot(x, y, label='Rho', color='black', ls="-", linewidth=lw)
    ax.text(.5, 1.01, '%5.2f days' % b.Time, horizontalalignment='center', transform=ax.transAxes)

    if 'xlim' in kwargs:
        xlim = kwargs['xlim']
    else:
        xlim = [min(x), max(x) * 1.1]

    if 'ylim' in kwargs:
        ylim = kwargs['ylim']
    else:
        ylim = [np.min(y), np.max(y) + 2]

    if 'is_xlabel' in kwargs:
        if kwargs['is_xlabel']:
            ax.set_xlabel(xlabel)
    else:
        ax.set_xlabel(xlabel)
    if 'is_ylabel' in kwargs:
        if kwargs['is_ylabel']:
            ax.set_ylabel('log10(Rho)')
    else:
        ax.set_ylabel('log10(Rho)')

    # Right axe
    ax2 = ax.twinx()
    ax2.set_ylim((0., 10.))

    y2 = b.Tau
    ax2.plot(x, y2, 'g-', label='tau')

    y2 = np.log10(b.T)
    ax2.plot(x, y2, 'r-', label='T')

    lumnorm = 1.e40
    if 'lumnorm' in kwargs:
        lumnorm = kwargs['lumnorm']
    # y2 = np.ma.log10(b.Lum)
    y2 = np.ma.log10(b.Lum) - np.log10(lumnorm)
    ax2.plot(x, y2, color='orange', ls="-", label='Lum%d' % int(np.log10(lumnorm)))

    vnorm = 1.e8
    if 'vnorm' in kwargs:
        vnorm = kwargs['vnorm']
    y2 = b.Vel / vnorm
    ax2.plot(x, y2, 'b-', label='V%d' % int(np.log10(vnorm)))

    # for tl in ax2.get_yticklabels():
    #     tl.set_color('r')
    is_legend = True
    if 'is_legend' in kwargs:
        is_legend = kwargs['is_legend']
    if is_legend:
        ax.legend(loc=2, prop={'size': 9})
        ax2.legend(prop={'size': 9})
    ax2.set_ylabel('T, Vel, Lum, Tau')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # ax.grid()
