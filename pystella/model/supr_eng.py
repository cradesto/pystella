import numpy as np
import os
import logging

from pystella.util.phys_var import phys

__author__ = 'bakl'

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

LEGEND_MASK_None = 0x00
LEGEND_MASK_Eng = 0x01


class SupremnaEngDetail:
    """
    Reader for eng-files with ENGemiss engxedot engeion
    """
    def __init__(self, name, path='./'):
        """Creates a SupremnaEngDetail instance.  Required parameters:  name."""
        fname, ext = os.path.splitext(name)
        if ext == "eng":
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
        Load datd from eng-file
        :return:
        """
        fname = os.path.expanduser(os.path.join(self._path, self._name + ".eng"))
        colstr = "tday km engemiss engxedot engeion"   #  log10() ergs
        cols = [s.strip() for s in colstr.split()]
        dt = np.dtype({'names': cols, 'formats': [float] * len(cols)})
        data = np.loadtxt(fname, dtype=dt)

        self._nzon = int(np.max(data['km']))
        times = np.unique(data['tday'])
        if data['tday'][0] != 0.:
            times = np.delete(times, np.where(times == 0.))
        self._times = times

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
        return BlockEng(self._data['tday'][idx], block)

    def __getitem__(self, idx):
        b = idx * self._nzon
        e = b + self._nzon
        block = self._data[:][b:e]
        return BlockEng(self._data['tday'][idx], block)

    def evolution(self, var, nz):
        """Return the time evolution VAR in zone NZ"""
        if nz < 1 or nz > self.Nzon:
            raise ValueError('NZ should be in range: 0 - {}'.format(self.Nzon))

        if var not in vars(BlockEng):
            raise ValueError('var should be property of BlockEng, like engemiss engxedot engeion')

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

    # def params_eng(self, cols=('engemiss', 'engxedot', 'engeion')):
    #     '''
    #         --     em(j) = emissivity in erg/(s cm3 eV sr)
    #         --     cool = 4*pi*energy integration of em(j) (= total cooling)
    #         emiss = ENGrad = cool/(Pl*URho)/UEng;
    #     '''        
    #     is_str = isinstance(cols, str)
    #     if is_str:
    #         cols = [cols]
    #     res = {k: np.zeros(self.Ntimes) for k in ['time', 'zone'] + list(cols)}
    #     for i, time in enumerate(self.Times):
    #         s = self[i]
    #         kph = 1  # todo check the default value
    #         for k in range(self.Nzon - 1, 1, -1):
    #         res['time'][i] = time
    #         res['zone'][i] = kph
    #         for v in cols:
    #             res[v][i] = getattr(s, v)[kph]
    #         # res['R'][i] = s.R[kph]
    #         # res['M'][i] = s.M[kph]
    #         # res['T'][i] = s.T[kph]
    #         # res['V'][i] = s.V[kph] / 1e8  # to 1000 km/s
    #     return res


class BlockEng:
    """
        Block eng-data for given time
        line: "tday km engemiss engxedot engeion"
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
    def engemiss(self):
        return 10. ** self._block['engemiss']  # [ergs]

    @property
    def engxedot(self):
        return 10. ** self._block['engxedot']  # [ergs]

    @property
    def engeion(self):
        return 10. ** self._block['engeion']  # [ergs]


# ==================================================================
def isfloat(value):
    try:
        float(value)
        return True
    except:
        return False


def plot_eng(axeng, b, **kwargs):
    name = kwargs.get('name', '')
    xlim = kwargs.get('xlim', None)
    ylim_eng = kwargs.get('ylim_eng', None)
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

    axeX = kwargs.get('axeX', 'z')
    engnorm = kwargs.get('engnorm', 1.)

    ls = kwargs.get('ls', '-'); #('-','--','-.'))
    lw = kwargs.get('lw', 1.)
    colors = ('black', 'red', 'blue')
    x, xlabel = b.Zon, r'Zone'    

    for i, eng in enumerate(('engemiss', 'engxedot', 'engeion')):
            y = getattr(b, eng)
        # y = np.log10(b.engemiss)
            axeng.plot(x, y, label=rf'{eng} {name}', ls=ls, linewidth=lw, color=colors[i])

    if is_day:
        axeng.text(.02, text_posy, r'$%5.2f^y$' % b.Time, horizontalalignment='left',
                   transform=axeng.transAxes)

    if islim:
        if xlim is None:
            xlim = [min(x), max(x) * 1.2]
        if ylim_eng is None:
            ylim_eng = [min(y), max(y)*1.2]
        axeng.set_xlim(xlim)
        axeng.set_ylim(ylim_eng)

    if is_xlabel:
        axeng.set_xlabel(xlabel)

    if is_yllabel:
        axeng.set_ylabel(r'$\log_{10}(Eng)$')
    else:
        axeng.set_yticklabels([])

    if legmask & LEGEND_MASK_Eng:
        axeng.legend(loc=1, prop={'size': 8}, frameon=is_frameon, handlelength=2.5)

    return axeng
