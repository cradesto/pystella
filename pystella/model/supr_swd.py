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


class SupremnaShockWaveDetail:
    """
    Reader for swd-files
    """

    def __init__(self, name, path='./'):
        """Creates a SupremnaShockWaveDetail instance.  Required parameters:  name."""
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
        # tout, Km,log10(AMPR),log10(Ry)),Uy(Km)*1.D+6/(u%Utime*CRAP),
        # log10(max(UTP*Tye(Km),1.d0)),log10(max(UTP*Tyi(Km),1.d0)),
        # PLLOG,PeLOG,PiLOG,QVLOG,CRLOG,log10(max(ABS(eng)*1.d12,1.d-50)),Flum*1.d-40,
        # Clum*FL1unsat*1.d-40,Clum*Fsat*1.d-40,WRKX(Km),XEswd;
        colstr = "tyear km lgM lgRpc V8 lgTe lgTi lgPl lgPe lgPi  lgQv lgCR lgEng Flum40 FL1unsat40 Fsat40 rhoke XEswd"
        cols = [s.strip() for s in colstr.split()]
        dt = np.dtype({'names': cols, 'formats': [float] * len(cols)})
        data = np.loadtxt(fname, dtype=dt)

        self._nzon = int(np.max(data['km']))
        times = np.unique(data['tyear'])
        if data['tyear'][0] != 0.:
            times = np.delete(times, np.where(times == 0.))
        self._times = times

        # self.ntimes = len(self.times)   # len(data['tyear']) // self.nzon
        self._data = data
        logger.debug("Read data from  %s " % fname)
        return self

    def time_nearest(self, time):
        idx = (np.abs(self.Times - time)).argmin()
        return idx, self.Times[idx]

    def block_nearest(self, time):
        idx = (np.abs(self._data['tyear'] - time)).argmin()
        b = idx
        e = b + self._nzon
        block = self._data[:][b:e]
        return BlockSwd(self._data['tyear'][idx], block)

    def __getitem__(self, idx):
        b = idx * self._nzon
        e = b + self._nzon
        block = self._data[:][b:e]
        return BlockSwd(self._data['tyear'][idx], block)

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

class BlockSwd:
    """
        Block swd-data for given time
        line: "tyear km lgM lgRpc V8 lgT lgTrad lgDm6 lgP7  lgQv lgEng Flum40 cap"
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
        return 10. ** self._block['lgRpc']  # [cm]

    def Rsun(self):
        return self.R / phys.R_sun  # [Rsun]

    @property
    def V(self):
        return self._block['V8'] * 1e8  # [cm/s]

    @property
    def V8(self):
        return self._block['V8']  # [1000 km/s]

    @property
    def Te(self):
        return 10. ** self._block['lgTe']  # [K]

    @property
    def Ti(self):
        return 10. ** self._block['lgTi']  # [K]

    @property
    def Rho(self):
        return 10. ** self._block['lgPl']  # [g*cm-3]

    @property
    def Lum(self):
        return self._block['Flum40'] * 1e40  # [erg/s]  ??

    @property
    def lgPl(self):
        return self._block['lgPl']  
    
    @property
    def lgPe(self):
        return self._block['lgPe']  
    
    @property
    def lgPi(self):
        return self._block['lgPi']  
    
    @property
    def lgEng(self):
        return self._block['lgEng']  
    
    @property
    def Mtot(self):
        return 10 ** self._block['lgM'][0]  # [Msun]   ??

    @property
    def M(self):
        return self.Mtot - 10 ** self._block['lgM']  # [Msun]   ??


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
    ylim_par = kwargs.get('ylim_par', (0., 15.9))
    legmask = kwargs.get('legmask', 0x11)  # mask 0x11 - both, 0x01-rho, 0x10 - pars
    is_frameon = kwargs.get('is_frameon', False)
    islim = kwargs.get('islim', True)
    is_xlabel = kwargs.get('is_xlabel', True)
    is_yllabel = kwargs.get('is_yllabel', True)
    is_yrlabel = kwargs.get('is_yrlabel', True)
    is_day = kwargs.get('is_day', True)
    is_grid = kwargs.get('is_grid', False)
    text_posy = kwargs.get('text_posy', 1.01)

    rnorm = kwargs.get('rnorm', 'm')
    tnorm = kwargs.get('tnorm', None)
    # tnorm = kwargs.get('tnorm', 1e3)
    vnorm = kwargs.get('vnorm', 1e8)
    lumnorm = kwargs.get('lumnorm', 1e40)

    ls = kwargs.get('ls', '-')
    lw = kwargs.get('lw', 1.)

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
    axrho, axpar = axs
    is_axpar_none = axpar is None
    if rnorm == 'lgr':
        axrho.semilogx(x, y, label=rf'$\rho$ {name}', color='black', ls=ls, linewidth=lw)
    else:
        axrho.plot(x, y, label=rf'$\rho$ {name}', color='black', ls=ls, linewidth=lw)

    if is_day:
        axrho.text(.02, text_posy, r'$%5.2f^y$' % b.Time, horizontalalignment='left',
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
            axpar.set_ylabel(r'$T_e, T_i, Vel, Lum$')
        else:
            axpar.set_yticklabels([])

    if is_yllabel:
        axrho.set_ylabel(r'$\log_{10}(\rho)$')
    else:
        axrho.set_yticklabels([])

    if tnorm is None:
        y2 = np.log10(b.Te)
        axpar.plot(x, y2, color='green', ls=ls, label=r'$\log(Te)$')
        y2 = np.log10(b.Ti)
        axpar.plot(x, y2, color='red', ls=ls, label=r'$\log(Ti)$')
    else:
        y2 = b.Te / tnorm
        axpar.plot(x, y2, color='green', ls=ls, label='Te{0:d}'.format(int(np.log10(tnorm))))
        y2 = b.Ti / tnorm
        axpar.plot(x, y2, color='red', ls=ls, label='Ti{0:d}'.format(int(np.log10(tnorm))))

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
