##
#  Callbacks
##

import numpy as np

from pystella.model.popov import Popov
from pystella.rf.rad_func import Lum2MagBol
from pystella.util.phys_var import phys


def plot(ax, dic=None):
    arg = []
    if dic is not None and 'args' in dic:
        arg = dic['args']

    r_init = 450.  # M_sun
    m_tot = 15.  # M_sun
    m_ni = 0.07  # M_sun
    e_tot = 0.7  # FOE
    if len(arg) > 0:
        r_init = float(arg.pop(0))
    if len(arg) > 0:
        m_tot = float(arg.pop(0))
    if len(arg) > 0:
        e_tot = float(arg.pop(0))
    if len(arg) > 0:
        m_ni = float(arg.pop(0))

    n = 100
    start, end = map(lambda x: max(x, 0.1), ax.get_xlim())
    time = np.exp(np.linspace(np.log(start), np.log(end), n))
    ppv = Popov('plugin', R=r_init, M=m_tot, Mni=m_ni, E=e_tot)

    mags = ppv.MagBol(time)
    mags_ni = Lum2MagBol(ppv.e_rate_ni(time))

    lbl = 'R M E Mni: %4.1f %2.1f %2.1f %3.2f' % (r_init, m_tot, e_tot, m_ni)
    print("Plot Popov model:  %s " % lbl)
    times = {'Diffusion time [d]': ppv.t_d,
             'Expansion time [d]': ppv.t_e,
             'T surf > Tion [d]': ppv.t_i,
             'Max bol time [d]': ppv.t_max,
             'Plateau duration time [d]': ppv.t_p}
    for k, v in times.items():
        print("  %25s: %8.2f " % (k, v / phys.d2s))

    ax.plot(time, mags,  color='blue', ls='-.', label=lbl, lw=2.5)
    ax.plot(time, mags_ni, color='red', ls='-.', label='Ni56 & Co56', lw=2.)
