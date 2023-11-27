import numpy as np
import os
import logging

from pystella.util.phys_var import phys

__author__ = 'bakl'

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

LEGEND_MASK_None = 0x00
LEGEND_MASK_Rho = 0x01
LEGEND_MASK_Vars = 0x10


class StellaShockWaveDetail:
    """
    Reader for swd-files
    """

    def __init__(self, name, path='./'):
        """Creates a StellaShockWaveDetail instance.  Required parameters:  name."""
        fname, ext = os.path.splitext(name)
        if ext == "swd":
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
        dt = np.dtype({'names': cols, 'formats': [float] * len(cols)})
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
        b = idx * self._nzon
        e = b + self._nzon
        block = self._data[:][b:e]
        return BlockSwd(self._data['tday'][idx], block)

    def evolution(self, var, nz):
        """Return the time evolution VAR in zone NZ"""
        if nz < 1 or nz > self.Nzon:
            raise ValueError('NZ should be in range: 0 - {}'.format(self.Nzon))

        if var not in vars(BlockSwd):
            raise ValueError('var should be property of BlockSwd, like Zon, M, R, Vel, T, Trad, Rho, Lum, Cappa, M')

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

    def params_ph(self, cols=('R', 'M', 'T', 'V', 'Rho'), tau_ph=2. / 3.):
        is_str = isinstance(cols, str)
        if is_str:
            cols = [cols]
        res = {k: np.zeros(self.Ntimes) for k in ['time', 'zone'] + list(cols)}
        taus = self.taus()
        for i, time in enumerate(self.Times):
            s = self[i]
            kph = 1  # todo check the default value
            for k in range(self.Nzon - 1, 1, -1):
                if taus[i][k] >= tau_ph:
                    kph = k
                    break
            res['time'][i] = time
            res['zone'][i] = kph
            for v in cols:
                res[v][i] = getattr(s, v)[kph]
            # res['R'][i] = s.R[kph]
            # res['M'][i] = s.M[kph]
            # res['T'][i] = s.T[kph]
            # res['V'][i] = s.V[kph] / 1e8  # to 1000 km/s
        return res

    def vel_ph(self, tau_ph=2. / 3., z=0.):
        v_dat = self.params_ph(tau_ph=tau_ph, cols=['V'])
        t = v_dat['time'] * (1. + z)  # redshifted time
        v = v_dat['V']
        z = v_dat['zone']
        return t, v, z


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
    def V(self):
        return self._block['V8'] * 1e8  # [cm/s]

    @property
    def V8(self):
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


def plot_swd(axs, b, **kwargs):
    name = kwargs.get('name', '')
    xlim = kwargs.get('xlim', None)
    ylim_rho = kwargs.get('ylim_rho', None)
    ylim_par = kwargs.get('ylim_par', (0., 9.9))
    legmask = kwargs.get('legmask', 0x11)  # mask 0x11 - both, 0x01-rho, 0x10 - pars
    is_frameon = kwargs.get('is_frameon', False)
    islim = kwargs.get('islim', True)
    is_xlabel = kwargs.get('is_xlabel', True)
    is_yllabel = kwargs.get('is_yllabel', True)
    is_yrlabel = kwargs.get('is_yrlabel', True)
    is_day = kwargs.get('is_day', True)
    is_grid = kwargs.get('is_grid', False)
    text_posy = kwargs.get('text_posy', 1.01)

    axeX = kwargs.get('axeX', 'm')
    tnorm = kwargs.get('tnorm', None)
    # tnorm = kwargs.get('tnorm', 1e3)
    vnorm = kwargs.get('vnorm', 1e8)
    lumnorm = kwargs.get('lumnorm', 1e40)

    ls = kwargs.get('ls', '-')
    lw = kwargs.get('lw', 1.)

    if axeX == 'sun':
        axeX = phys.R_sun
        x, xlabel = b.R / axeX, r'Ejecta Radius, [$\mathtt{R}_\odot$]'
    elif axeX == 'lgr':
        x, xlabel = b.R, r'Ejecta Radius, [cm]'
    elif axeX == 'z':
        x, xlabel = b.Zon, r'Zone'
    elif axeX == 'm':
        x, xlabel = b.M, r'Ejecta Mass [$\mathtt{M}_\odot$]'
    elif isfloat(axeX):
        x, xlabel = b.R / float(axeX), r'Ejecta Radius, [$\times 10^{%d}$ cm]' % int(np.log10(float(axeX)))
    else:
        x, xlabel = b.M, r'Ejecta Mass [$\mathtt{M}_\odot$]'

    y = np.log10(b.Rho)
    axrho, axpar = axs
    is_axpar_none = axpar is None
    if axeX == 'lgr':
        axrho.semilogx(x, y, label=rf'$\rho$ {name}', color='black', ls=ls, linewidth=lw)
    else:
        axrho.plot(x, y, label=rf'$\rho$ {name}', color='black', ls=ls, linewidth=lw)

    if is_day:
        axrho.text(.02, text_posy, r'$%5.2f^d$' % b.Time, horizontalalignment='left',
                   transform=axrho.transAxes)

    if islim:
        if xlim is None:
            xlim = [min(x), max(x) * 1.2]
        if ylim_rho is None:
            ylim_rho = [np.min(y), np.max(y) + 2]
        axrho.set_xlim(xlim)
        axrho.set_ylim(ylim_rho)

    if is_xlabel:
        axrho.set_xlabel(xlabel)
    else:
        pass
        # ax.set_xticklabels([])

    # Right axe
    if is_axpar_none:
        axpar = axrho.twinx()
        axpar.set_ylim(ylim_par)
        if is_yrlabel:
            axpar.set_ylabel(r'T, Vel, Lum, $\tau$')
        else:
            axpar.set_yticklabels([])

    if is_yllabel:
        axrho.set_ylabel(r'$\log_{10}(\rho)$')
    else:
        axrho.set_yticklabels([])

    # y2 = b.Tau
    y2 = np.ma.masked_where(b.Tau <= 0., b.Tau)
    axpar.plot(x, y2, color='red', ls=ls, label=r'$\tau$')

    if tnorm is None:
        y2 = np.log10(b.T)
        axpar.plot(x, y2, color='green', ls=ls, label=r'$\log(T)$')
    else:
        y2 = b.T / tnorm
        axpar.plot(x, y2, color='green', ls=ls, label='T{0:d}'.format(int(np.log10(tnorm))))

    y2 = np.ma.log10(b.Lum) - np.log10(lumnorm)
    y22 = np.ma.log10(-b.Lum) - np.log10(lumnorm)
    axpar.plot(x, y2, color='orange', ls=ls, label=r'$Lum_{{{0:d}}}$'.format(int(np.log10(lumnorm))))
    axpar.plot(x, y22, color='brown', ls="--", label=r'$-Lum_{{{0:d}}}$'.format(int(np.log10(lumnorm))))

    y2 = b.V / vnorm
    axpar.plot(x, y2, color='blue', ls=ls, label=r'$V_{{{0:d}}}$'.format(int(np.log10(vnorm))))

    if legmask & LEGEND_MASK_Rho:
        axrho.legend(loc=1, prop={'size': 8}, frameon=is_frameon, handlelength=2.5)
    if (legmask & LEGEND_MASK_Vars) and is_axpar_none:
        axpar.legend(loc=1, prop={'size': 8}, ncol=1, frameon=is_frameon, handlelength=1)

    if is_grid:
        axpar.grid(linestyle=':')
    return axrho, axpar
