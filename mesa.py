#!/usr/bin/env python3

import os
import sys
import logging

import pystella as ps

try:
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines

    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.ERROR)
except ImportError as ex:
    import os
    import sys

    exc_type, exc_obj, exc_tb = sys.exc_info()
    fn = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    logging.info(exc_type, fn, exc_tb.tb_lineno, ex)
    logging.info('  Probably, you should install module: {}'.format('matplotlib'))
    #    print(ex)
    plt = None
    mlines = None

logger = logging.getLogger(__name__)
level = logging.INFO

markers = {u'x': u'x', u'o': u'circle', u'v': u'triangle_down', u'd': u'thin_diamond',
           u'+': u'plus', u'*': u'star', u'<': u'triangle_left'}
markers_style = list(markers.keys())
lines_style = ps.rf.light_curve_plot.linestyles

def get_parser():
    import argparse
    from argparse import RawTextHelpFormatter

    parser = argparse.ArgumentParser(
              description='Convert MESA to STELLA presupernova model.', 
              formatter_class=RawTextHelpFormatter)

    parser.add_argument('-r', '--rho', nargs="?",
                        required=False,
                        const=True,
                        dest="rho",
                        metavar="<r OR m>",
                        help="Plot Rho-figure")

    parser.add_argument('-x',
                        required=False,
                        dest="x",
                        default='m',
                        metavar="<m OR r OR lgR OR rsun OR m OR z>",
                        help="Setup abscissa: radius or lg(R) OR mass OR zone")

    parser.add_argument('-i', '--input', action='append', nargs=1,
                        metavar='model name', help='Key -i can be used multiple times')

    parser.add_argument('-p', '--path',
                        required=False,
                        type=str,
                        default='./',
                        dest="path",
                        help="Model directory")
    parser.add_argument('-es', '--elements_stella',
                        required=False,
                        type=str,
                        default='H:He:C:O:Si:Fe:Ni56',
                        dest="elements_stella",
                        help="Elements in STELLA. \n   Available: {0}"
                             .format(':'.join(ps.eve.eve_elements)))
    parser.add_argument('-em', '--elements_mesa',
                        required=False,
                        type=str,
                        default='h1:he4:c12:o16:si28:fe56:ni56',
                        dest="elements_mesa",
                        help="Elements in MESA. \n   Available: {0}"
                             .format(':'.join(ps.model.mesa.mesa_elements)))
    parser.add_argument('--log',
                        required=False,
                        type=str,
                        default='INFO',
                        dest="log",
                        help="Set logging level: "
                             "CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET"
                        )
    parser.add_argument('--is_mcut',
                        action='store_const',
                        const=True,
                        dest="is_mcut",
                        help="Cut m_core")
    parser.add_argument('--verb',
                        action='store_const',
                        const=True,
                        dest="is_verb",
                        help="To enable verbose output")
    parser.add_argument('--stella',
                        action='store_const',
                        const=True,
                        dest="is_stella",
                        help="Convert MESA to Stella elements")
    parser.add_argument('-s', '--save',
                        required=False,
                        type=str,
                        default=False,
                        dest="save_plot",
                        help="save plot to pdf-file, default: False")
    parser.add_argument('-w', '--write',
                        required=False,
                        type=str,
                        default=False,
                        dest="write_to",
                        help="To write the data to hyd-, abn-files")
    return parser


# def usage():
#     print("Convert MESA to STELLA presupernova model.")
#     print("Usage:")
#     print("  mesa  -i  <datafile>")
#     print("  -p <model directory>, default: ./")
#     print("  -h  print usage")
#     print("   --- ")


def main():
    from itertools import cycle

    def is_level(lvl):
        levels = ['CRITICAL', 'FATAL','ERROR','WARN','WARNING','INFO','DEBUG','NOTSET']
        return lvl.upper() in levels

    def get_elements(elements_in, elements_default):
        if '_' in elements_in:
            elements = list(elements_default)
            excluded = elements_in.split(':')
            for e in excluded:
                if not e.startswith('_'):
                    logger.error('For excluded mode all elements should be starts from _. Even element: ' + e)
                    sys.exit(2)
                e = e[1:]
                if e not in elements_default:
                    logger.error('No such element: ' + e)
                    sys.exit(2)
                elements.remove(e)
        else:
            elements = elements_in.split(':')
            for e in elements:
                if e not in elements_default:
                    logger.error('No such element: ' + e)
                    sys.exit(2)
        return elements

    file_profile = None
    markersize = 6

    parser = get_parser()
    args, unknownargs = parser.parse_known_args()

    if args.path:
        pathDef = os.path.expanduser(args.path)
    else:
        pathDef = os.getcwd()

    level = logging.INFO
    if args.log is not None and is_level(args.log):
        level = logging.getLevelName(args.log.upper())
    else:
        logger.error(f"ERROR: Bad value for log: {level}. See help.")
        sys.exit(2)

    if args.is_stella:
        elements = get_elements(args.elements_stella, ps.eve.eve_elements)
        lntypes = ps.model.sn_eve.eve_lntypes
        colors = ps.model.sn_eve.eve_colors
    else:
        elements = get_elements(args.elements_mesa, ps.model.mesa.mesa_elements)
        lntypes = ps.model.mesa.mesa_el_lntypes 
        colors = ps.model.mesa.mesa_el_colors

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
        sys.exit()
    
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
        # logger.info("Run mesa-model %s" % nm)
        path, fullname = os.path.split(nm)
        if len(path) == 0:
            path = pathDef
        name = fullname.replace('.data', '')  # remove extension
        # print(f"Load mesa model {name} in {path}")
        file_profile = os.path.join(path, fullname)
        if not os.path.isfile(file_profile):
            raise FileNotFoundError(f"No file {file_profile}")

        prb_mesa = ps.Mesa(name).load(file_profile)
        if args.is_stella:
            presn = prb_mesa.to_presn_stella_el(is_mcut=args.is_mcut)
        else:
            presn = prb_mesa.to_presn(is_mcut=args.is_mcut)

        if args.is_verb:
            m_tot = 0.
            m_ni56 = 0.
            print('{:22s}: '.format(presn.Name))   
            for n, m in presn.mass_tot_el().items():
                m_tot += m
                print(f'{n} {m/ps.phys.M_sun:.3e}')
            print('  m_tot= {:6.3f}    m_tot(El)= {:6.3f} '.format(presn.m_tot/ps.phys.M_sun, m_tot/ps.phys.M_sun))   
        
        if args.write_to:
            fname = presn.Name + '.hyd'
            # fname = os.path.join(os.path.dirname(file_profile), presn.name + '.hyd')
            res = presn.write_hyd(fname)
            if res:
                print(f"Save in {fname}")
            else:
                print("Error while saving in %s" % fname)

            fname = presn.Name + '.abn'
            # fname = os.path.join(os.path.dirname(file_profile), presn.name + '.abn')
            res = presn.write_abn(fname, is_header=True)
            if res:
                print(f"Save in {fname}")
            else:
                print("Error while saving in %s" % fname)

            continue

        # Plot
        # print "Plot eve-model %s" % name
        marker = next(markers_cycler)
        ls = next(lines_cycler)

        ax = presn.plot_chem(ax=ax, elements=elements, x=args.x, ylim=(1e-8, 1.)
                            , marker=marker, colors=colors, ls=lntypes
                            , markersize=markersize, leg_loc='lower center')
        # else:
        #     ax = ps.Mesa.plot_chem(presn, elements=elements, ax=ax, x=args.x, ylim=(1e-8, 1.), marker=marker,
        #                         markersize=markersize, leg_loc='lower center')

        if args.rho:
            if ax2 is None:
                ax2 = ax.twinx()
                ax2.set_ylabel(r'$\rho, [g/cm^3]$ ')
            ax2 = presn.plot_rho(x=args.x, ax=ax2, ls=ls, marker=marker)
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
        if args.save_plot:
            fsave = os.path.expanduser(args.save_plot)
            if not fsave.endswith('.pdf'):
                fsave += '.pdf'
            logger.info(" Save plot to %s " % fsave)
            fig = ax.get_figure()
            fig.savefig(fsave, bbox_inches='tight')
            print(f'The plot has been saved to {fsave}')
        else:
            plt.show()

if __name__ == '__main__':
    main()
