#!/usr/bin/env python3

import logging
# import numpy as np

import pystella.model.sn_eve as sneve
import pystella.rf.light_curve_plot as lcp
from pystella import phys

try:
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines

    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)
except ImportError as ex:
    import os
    import sys

    exc_type, exc_obj, exc_tb = sys.exc_info()
    fn = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    logging.info(exc_type, fn, exc_tb.tb_lineno, ex)
    logging.info('  Probably, you should install module: {}'.format('matplotlib'))
    print()
    #    print(ex)
    plt = None
    mlines = None

# matplotlib.use("Agg")
# import matplotlib
# matplotlib.rcParams['backend'] = "TkAgg"
# matplotlib.rcParams['backend'] = "Qt5Agg"
# matplotlib.rcParams['backend'] = "Qt4Agg"

__author__ = 'bakl'

logging.basicConfig()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

markers = {u'x': u'x', u'o': u'circle', u'v': u'triangle_down', u'd': u'thin_diamond',
           u'+': u'plus', u'*': u'star', u'<': u'triangle_left'}
markers_style = list(markers.keys())
lines_style = lcp.linestyles


def get_parser():
    import argparse

    parser = argparse.ArgumentParser(description='Process PreSN configuration.')

    parser.add_argument('-b', '--box',
                        required=False,
                        type=str,
                        default=None,
                        dest='box',
                        help='Make boxcar average, for example: '
                             'Delta_mass:Number:[True, if info], -b 0.5:4 . '
                             'Use key -e _ELEM [-e _Ni56]to exclude elements')

    parser.add_argument('-r', '--rho', nargs="?",
                        required=False,
                        const=True,
                        dest="rho",
                        metavar="<r OR m>",
                        help="Plot Rho-figure")

    parser.add_argument('--is_dum', nargs="?",
                        required=False,
                        const=True,
                        dest="is_dum",
                        help="Set is_dum = TRUE to parse abn-file with dum columns")

    parser.add_argument('-x',
                        required=False,
                        dest="x",
                        default='m',
                        metavar="<m OR r OR lgR OR rsun OR m OR z>",
                        help="Setup abscissa: radius or lg(R) OR mass OR zone")

    parser.add_argument('-s', '--save',
                        required=False,
                        type=bool,
                        default=False,
                        dest="is_save_plot",
                        help="save plot to pdf-file, default: False")

    # parser.add_argument('-c', '--chem',
    #                     required=False,
    #                     type=bool,
    #                     default=True,
    #                     dest="is_chem",
    #                     help="Show chemical composition, default: True")
    parser.add_argument('--structure', dest='is_structure', action='store_true',
                        help="Show the chemical composition and rho with R/M coordinates.")
    parser.add_argument('--chem', dest='is_chem', action='store_true', help="Show chemical composition [default].")
    parser.add_argument('--no-chem', dest='is_chem', action='store_false', help="Not show chemical composition")
    parser.set_defaults(is_chem=True)

    parser.add_argument('-i', '--input', action='append', nargs=1,
                        metavar='model name', help='Key -i can be used multiple times')

    parser.add_argument('-p', '--path',
                        required=False,
                        type=str,
                        default='./',
                        dest="path",
                        help="Model directory")

    parser.add_argument('-e', '--elements',
                        required=False,
                        type=str,
                        default='H:He:C:O:Si:Fe:Ni:Ni56',
                        dest="elements",
                        help="Elements directory. \n   Available: {0}".format(':'.join(sneve.eve_elements)))
    parser.add_argument('--reshape',
                        required=False,
                        type=str,
                        default=None,
                        dest="reshape",
                        help="Reshape parameters of envelope from nstart to nend to nz-zones."
                             "\n Format: --reshape NZON:AXIS:XMODE:START:END:KIND. You may use * to set default value."
                             "\n NZON: value of zones between START and END. "
                             "If < 0 Nzon is the same as Nzon of the initial model "
                             "\n AXIS: [M* OR R OR V] - reshape along mass or radius or velocity coordinate."
                             "\n XMODE: [lin OR rlog* OR resize] - linear OR reversed log10 OR add/remove points. "
                             "\n START: zone number to start reshaping. Default: 0 (first zone)"
                             "\n END: zone number to end reshaping. Default: None,  (equal last zone)"
                             "\n KIND: [np OR interp1d(..kind)], kind is  ('np=np.interp', 'linear', 'nearest', "
                             "'zero', 'slinear', 'quadratic, 'cubic', "
                             "'spline' = UnivariateSpline, 'gauss' = gaussian_filter1d). Default: np "
                        )
    # parser.add_argument('-w', '--write',
    #                     action='store_const',
    #                     const=True,
    #                     dest="is_write",
    #                     help="To write the data to hyd-, abn-files")
    parser.add_argument('-w', '--write',
                        required=False,
                        type=str,
                        default=False,
                        dest="write_to",
                        help="To write the data to hyd-, abn-files")
    return parser


def print_masses(presn):
    m_el_tot = 0.
    for ii, el in enumerate(presn.Elements):
        m = presn.mass_tot_el(el) / phys.M_sun
        m_el_tot += m
        print(f'  {el:3}:  {m:.3e}')
    print(f'  M_full(Elements) =  {m_el_tot:.3f}')
    print(f'  M_total =  {presn.m_tot / phys.M_sun:.3f}')
    # via density
    print(f'  M_tot(Density) =  {presn.mass_tot_rho() / phys.M_sun:.3f}')


def main():
    import os
    import sys
    from itertools import cycle

    def get(arr, i, default):
        if i < len(arr):
            if a[i] != '*':
                return a[i]
        return default

    parser = get_parser()
    args, unknownargs = parser.parse_known_args()
    eve_prev = None
    markersize = 6
    fig = None

    if args.path:
        pathDef = os.path.expanduser(args.path)
    else:
        pathDef = os.getcwd()

    # if args.elements:
    if '_' in args.elements:
        elements = list(sneve.eve_elements)
        excluded = args.elements.split(':')
        for e in excluded:
            if not e.startswith('_'):
                logger.error('For excluded mode all elements should be starts from _. Even element: ' + e)
                sys.exit(2)
            e = e[1:]
            if e not in sneve.eve_elements:
                logger.error('No such element: ' + e)
                sys.exit(2)
            elements.remove(e)
    else:
        elements = args.elements.split(':')
        for e in elements:
            if e not in sneve.eve_elements:
                logger.error('No such element: ' + e)
                sys.exit(2)

    # Set model names
    names = []
    if args.input:
        for nm in args.input:
            names.append(nm[0])  # remove extension
    else:
        if len(unknownargs) > 0:
            names.append(unknownargs[0])

    if len(names) == 0:
        # logger.error(" No data. Use key '-i' ")
        parser.print_help()
        sys.exit(2)

    if len(names) > 1:  # special case
        markers_cycler = cycle(markers_style)
        lines_cycler = cycle(lines_style)
    else:
        markers_cycler = cycle([None])
        lines_cycler = cycle(['-'])

    ax = None
    ax2 = None
    handles_nm = []
    for nm in names:
        # print("Run eve-model %s" % nm)
        path, fullname = os.path.split(nm)
        if len(path) == 0:
            path = pathDef
        # print("Run eve-model %s in %s" % (name, path))

        if fullname.endswith('hyd') or fullname.endswith('abn'):
            name = fullname.replace('.hyd', '')  # remove extension
            name = name.replace('.abn', '')  # remove extension
            try:
                # With header
                eve = sneve.load_hyd_abn(name=name, path=path, is_dm=False, is_dum=args.is_dum)
            except ValueError:
                # No header
                eve = sneve.load_hyd_abn(name=name, path=path, is_dm=False, is_dum=args.is_dum, skiprows=0)
        else:
            name = fullname.replace('.rho', '')  # remove extension
            rho_file = os.path.join(path, name + '.rho')
            eve = sneve.load_rho(rho_file)

        if args.reshape is not None:
            a = args.reshape.split(':')
            nz, axis, xmode = get(a, 0, eve.nzon), get(a, 1, 'M'), get(a, 2, 'resize')  # rlog
            start, end = get(a, 3, 0), get(a, 4, None)
            kind = get(a, 5, 'np')
            start = int(start)
            if end.upper() in 'NONE':
                end = None
            if end is not None:
                end = int(end)
            nz = int(nz)
            print(f'Resize: before Nzon={eve.nzon}')
            print(f'Resize parameters: nznew= {nz}  axis={axis}  xmode={xmode}  '
                  f'start= {start}  end= {end} kind= {kind}')
            print("The element masses: before Resize")
            print_masses(eve)
            eve = eve.reshape(nz=nz, axis=axis, xmode=xmode, start=start, end=end, kind=kind)
            eve.chem_norm()
            # eve = eve_resize
            print(f'Resize: after Nzon={eve.nzon}')
            print("The element masses: after Resize")
            print_masses(eve)

        # Boxcar
        if args.box is not None:
            is_info = False
            s = args.box.split(':')
            dm, n = float(s[0]), int(s[1])
            if len(s) == 3:
                is_info = bool(s[2])
            print(f'Running boxcar average: dm= {dm} Msun  Repeats= {n}')
            print("The element masses: Before boxcar")
            print_masses(eve)
            eve_box = eve.boxcar(box_dm=dm, n=n, el_included=elements, is_info=is_info)
            print("The element masses: After boxcar")
            print_masses(eve_box)
            eve, eve_prev = eve_box, eve

        if args.write_to:
            fname = os.path.expanduser(args.write_to)
            # fname = os.path.join(path, name)
            # f = fname + '.eve.abn'
            fname = fname.replace('.rho', '')
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

            continue

        marker = next(markers_cycler)
        ls = next(lines_cycler)

        if args.is_structure:
            fig = eve.plot_structure(elements=elements, title=name, ylimChem=(1e-8, 1.))
        else:
            if args.is_chem:
                # print "Plot eve-model %s" % name
                ax = eve.plot_chem(elements=elements, ax=ax, x=args.x, ylim=(1e-8, 1.), marker=marker,
                                   markersize=markersize, leg_loc='lower center')
                if eve_prev is not None:
                    eve_prev.plot_chem(elements=elements, ax=ax, x=args.x, ylim=(1e-8, 1.), marker=marker,
                                       markersize=max(1, markersize - 2), alpha=0.5, leg_loc='lower center')
                    # ax.set_title('{}: before boxcar'.format(eve_prev.Name))

            if args.rho:
                if args.is_chem:
                    if ax2 is None:
                        ax2 = ax.twinx()
                        ax2.set_ylabel(r'$\rho, [g/cm^3]$ ')
                else:
                    ax2 = ax
                ax2 = eve.plot_rho(x=args.x, ax=ax2, ls=ls, marker=marker)
                if eve_prev is not None:
                    eve_prev.plot_rho(x=args.x, ax=ax2, ls=ls, markersize=max(1, markersize - 2), alpha=0.5)
            else:
                ls = 'None'

            handle = mlines.Line2D([], [], color='black', marker=marker,
                                   markersize=markersize, label=name, linestyle=ls)
            handles_nm.append(handle)
            if len(names) > 1:
                if ax2 is None:
                    ax2 = ax.twinx()
                ax2.legend(handles=handles_nm, loc=4, fancybox=False, frameon=False)

    if not args.write_to:
        if args.is_save_plot:
            if args.rho:
                fsave = os.path.join(os.path.expanduser('~/'), 'rho_%s.pdf' % names[0])
            else:
                fsave = os.path.join(os.path.expanduser('~/'), 'chem_%s.pdf' % names[0])
            logger.info(" Save plot to %s " % fsave)
            if fig is None:
                fig = ax.get_figure()
            fig.savefig(fsave, bbox_inches='tight')
        else:
            plt.show()


if __name__ == '__main__':
    main()
