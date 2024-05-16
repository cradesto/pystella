#!/usr/bin/env python3

import logging
import os
from itertools import cycle
import h5py

import numpy as np
import pystella as ps

import sys
try:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
except ImportError:
    plt = None
    mcolors, cm = None, None


__author__ = 'bakl'

# ROOT_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))
logging.basicConfig()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.ERROR)


def get_parser(time="-1"):
    import argparse
    from argparse import RawTextHelpFormatter

    parser = argparse.ArgumentParser(description='Process the date from hdf-file.',
                                     formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('-i', '--input', action='append', nargs=1,
                        metavar='Model name', help='Key -i can be used multiple times.')
    
    parser.add_argument('-p', '--path',
                    required=False,
                    type=str,
                    default='./',
                    dest="path",
                    help="Model directory")
    parser.add_argument('--plot_rho',
                    action='store_const',
                    const=True,
                    dest="is_plot_rho",
                    help="Plot densites for all times.")
    parser.add_argument('--save_hyd',
                    action='store_const',
                    const=True,
                    dest="is_save_hyd",
                    help="Save hyd&abn at the current moment.")

    parser.add_argument('-t', '--time',
                        required=False,
                        type=str,
                        default=time,
                        dest="time",
                        help="""Selected time points for shock wave snapshots.
If value is integer (like 10), it used as an index for points in time.
It can be negative, like -2.
If value is float (like 10.), it used as time [days].
if value is string, like None, all times used for plot.
Default: {} = last saved moment".format(time))    
""")
    return parser

def hdf2presn(fnameh5, idx_time):
    if idx_time is None:
        raise ValueError("Bad value for time. Try key -t -1 ")
    
    abn_elements = ps.PreSN.stl_elements
    # Load 
    snh5 = ps.H5Stella(fnameh5)         # 
    print(snh5.Name)
    hyd = snh5.Hyd
    print(hyd.Attrs)
    times = hyd.Time
    t_idx = times[idx_time] * 86400 # tday to seconds

    presn = ps.PreSN(snh5.Name, hyd.Nzon, elements=abn_elements)
    # time_start, nzon, m_core, r_cen, rho_cen = a
    presn.set_par('time_start', t_idx)
    # presn.set_par('m_core', m_core * phys.M_sun)
    # presn.set_par('r_cen', r_cen)
    # presn.set_par('rho_cen', rho_cen)
    
    # Hyd
    print("hyd.Columns: ", hyd.Columns)
    print("snh5.M : ", snh5.M )
    presn.set_hyd(ps.PreSN.sM, snh5.M * ps.phys.M_sun)
    presn.set_hyd(ps.PreSN.sR, hyd.R[idx_time])
    presn.set_hyd(ps.PreSN.sT, hyd.T[idx_time])
    presn.set_hyd(ps.PreSN.sRho, hyd.Rho[idx_time])
    presn.set_hyd(ps.PreSN.sV, hyd.V[idx_time]*1e8)

    # # Set Mass
    # if is_rho:
    #     r = presn.r
    #     rho = presn.rho
    #     r = np.insert(r, 0, presn.r_cen)
    #     #       rho = np.insert(rho, 0, presn.rho_cen)
    #     dm = np.zeros(nz)
    #     for i in range(nz):
    #         dm[i] = (r[i + 1] ** 3 - r[i] ** 3) * rho[i] * 4. / 3. * np.pi
    #     #            dm[i] = (r[i + 1] ** 3 - r[i] ** 3) * rho[i + 1] * 4. * np.pi / 3.
    #     m = np.cumsum(dm)
    #     m += presn.m_core
    # else:
    #     m = data_hyd[PreSN.sM] * phys.M_sun

    # presn.set_hyd(PreSN.sM, m)

    # Set chemical composition 
    if idx_time == 0:
        logger.info(' Load abn-data from  /presn/Yabun')
        data_chem = snh5.Yabun
        x_isotopes = snh5.Xisotopes
        for i, ename in enumerate(abn_elements):
            if ename == 'Ni56':
                print(i, ename, x_isotopes[:,0])
                presn.set_chem(ename, x_isotopes[:,0])
            else:
                print(i, ename)
                presn.set_chem(ename, data_chem[:,i]*ps.AZ[ename])
    else:
        logger.info(' Load abn-data from  /timing/AbunIso')
        abun = snh5.Abun
        elements = abun.Columns.split()
        # print("Colunms: ", abun.Columns)
        for i, ename in enumerate(elements):
            el = abun.el(ename)
            print(i, ename, el.shape)
            presn.set_chem(ename.strip(), el[idx_time,:]*ps.AZ[ename])
    # print("dbg stop")
    # exit(1)
    return presn


def presn_save_hyd(fname, eve):
    f = fname + '.abn'
    if eve.write_abn(f, is_header=True):
        print(" abn has been saved to {}".format(f))
    else:
        print("Error with abn saving to {}".format(f))
    # f = fname + '.eve.hyd'
    f = fname + '.hyd'
    if eve.write_hyd(f):
        print(" hyd has been saved to {}".format(f))
    else:
        print("Error with hyd saving to {}".format(f))


def plot_density(hyd, idx_time=None):
    times = hyd.Time
    x = hyd.R
    y = hyd.Rho

    fig, ax = plt.subplots(figsize=(12, 8))
    colorparams = times
    colormap = cm.viridis
    # colormap = cm.jet
    normalize = mcolors.LogNorm(vmin=np.min(colorparams), vmax=np.max(colorparams))

    if idx_time is not None:
        t = times[idx_time]
        color = colormap(normalize(t))
        ax.plot(x[idx_time], y[idx_time], color=color, label=f't= {t}')
        ax.legend()
    else:
        print("Times: ", [(i, t) for i, t in enumerate(times)])
        for i, t in enumerate(times):
            color = colormap(normalize(t))
            ax.plot(x[i], y[i], color=color)

        # Colorbar setup
        s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
        s_map.set_array(colorparams)

        # Use this to emphasize the discrete color values
        cbar = fig.colorbar(s_map, ax=ax, spacing='proportional', format='%3.1f')  # format='%2i' for integer

        cbar.set_label(r'Times', fontsize=20)

    ax.set_xlabel(r'Radius', fontsize=20)
    ax.set_ylabel('Density', fontsize=20)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ylims = ax.get_ylim()
    ax.set_ylim(top=1.1 * ylims[1])

    plt.show()


def main():
    parser = get_parser()
    args, unknownargs = parser.parse_known_args()    


    if args.path:
        pathDef = os.path.expanduser(args.path)
    else:
        pathDef = os.getcwd()

    # Set model names
    names = []
    if args.input:
        for nm in args.input:
            names.append(nm[0])  # remove extension
    else:
        if len(unknownargs) > 0:
            names.append(unknownargs[0])

    if len(names) == 0:
        logger.error(" No data. Use key '-i' ")
        parser.print_help()
        sys.exit(2)

    #
    path, fname = os.path.split(ps.first(names))
    if len(path) == 0:
        path = pathDef
    fname_full = os.path.join(path, fname)            
    # name = fname.replace('.h5', '')  # remove extension

    snh5 = ps.H5Stella(fname_full)         # 
    hyd = snh5.Hyd
    times = np.asarray(hyd.Time)

    # Time 
    idx_time = None
    if args.time.replace('-','').isdigit():
        idx_time = int(args.time)
        time = times[idx_time]
    elif args.time.lower() == 'none':
        idx_time = None
        time = -1
    else:
        try:
            time = float(args.time)
            idx_time = (np.abs(times - time)).argmin()
            time = times[idx_time]
        except:
            idx_time = -1
            time = -1

    # if time < 0:
    #     idx_time = -1

    print("Loading h5-data from {} at t= {:12.4f} [idx={}]".format(fname_full, time, idx_time))        
    
    if args.is_plot_rho:
        plot_density(hyd, idx_time=idx_time)

    if args.is_save_hyd:
        print("Save h5-data to hyd&abn format at t= {:12.4f} [idx={}]".format(time, idx_time))        
        presn = hdf2presn(fname_full, idx_time=idx_time)
        # ax = presn.plot_chem(ylim=(1e-12, 1.))
        # plt.show()
        fname = '~/{}_t{:.0f}'.format(presn.Name, time)
        fname = os.path.expanduser(fname)
        presn_save_hyd(fname,presn)

if __name__ == '__main__':
    main()