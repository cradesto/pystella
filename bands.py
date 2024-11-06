#!/usr/bin/env python3
# #!/usr/bin/python3

import numpy as np

from pystella.util.phys_var import phys
import pystella.rf.band as band

try:
    import matplotlib.pyplot as plt
except ImportError as ex:
    import os
    import sys
    # import traceback

    exc_type, exc_obj, exc_tb = sys.exc_info()
    fn = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fn, exc_tb.tb_lineno, ex)
    print('  Probably, you should install module: {}'.format('matplotlib'))
    #    print(ex)
    plt = None
    pass

__author__ = 'bakl'


def plot_griuz():
    plt.title('griuz filter response')
    bands = dict(g='g', i='r+', r='r', u='o', z='*')
    plot_bands(bands.keys(), bands)


def plot_UBVRI():
    plt.title('UBVRI filter response')
    bands = dict(U='b', B='c', V='g', R='r', I='p')
    plot_bands(bands.keys(), bands)


def plot_JHK():
    plt.title('JHK filter response')
    bands = dict(J='b', H='c', K='r')
    plot_bands(bands.keys(), bands)


def plot_SWIFT():
    plt.title('SWIFT filter response')
    bands = dict(UVM2='m', UVW1='r', UVW2='b', SwiftU='k', SwiftB='c', SwiftV='g')
    plot_bands(bands.keys(), bands)


def plot_PS1():
    plt.title('The Pan-STARRS1 Photometric  filter responses')
    bands = dict(PS1g='m', PS1i='r', PS1r='b', PS1z='y', PS1y='g', PS1w='p')
    plot_bands(bands.keys(), bands)


def plot_HSC():
    plt.title('The Hyper Suprime-Cam(HSC) Photometric  filter responses')

    bands = dict(HSCg='m', HSCr='b', HSCi='r', HSCz='y', HSCY='g')
    plot_bands(bands.keys(), bands)


def plot_HST():
    plt.title('The Hubble Space Telescope  filter responses')

    bands = dict(F105W="blue", F125W="g", F435W="skyblue", F140W="orange", F160W="r", F606W="cyan", F814W="magenta")
    plot_bands(bands.keys(), bands)


def plot_Misc():
    plt.title('The filter responses of the different sources')
    bands = dict(AtlasC="cyan", AtlasO="orange")
    plot_bands(bands.keys(), bands)

    # for k, v in bands.items():
    #     b = band.band_by_name(k)
    #     plt.plot(b.wl * phys.cm_to_angs, b.resp_wl, v, label=k, linewidth=2)
    # plt.legend(loc=4)
    # plt.ylabel('Amplitude Response')
    # plt.xlabel('Wave [A]')
    # plt.grid(linestyle=':')
    # plt.show()


def plot_Kepler():
    plt.title('The Kepler Space Telescope  filter responses')

    bands = dict(Kepler="magenta")
    plot_bands(bands.keys(), bands)
    # for k, v in bands.items():
    #     b = band.band_by_name(k)
    #     plt.plot(b.wl * phys.cm_to_angs, b.resp_wl, v, label=k, linewidth=2)
    # plt.legend(loc=4)
    # plt.ylabel('Amplitude Response')
    # plt.xlabel('Wave [A]')
    # plt.grid(linestyle=':')
    # plt.show()


def plot_bands(bands, color_dic=None, is_norm=False):   
    for bname in bands:
        b = band.band_by_name(bname)
        if color_dic is None:
            c = band.colors(bname)
        else:
            c  = color_dic[bname]
        resp = b.resp_wl 
        if is_norm:
            resp = (resp-np.min(resp))/(np.max(resp)-np.min(resp))
        plt.plot(b.wl * phys.cm_to_angs, resp, c, label=bname, linewidth=2)
        # np.savetxt(bname+'.dat', np.c_[b.wl * phys.cm_to_angs, resp], fmt='%12.0f %12.8f')

    plt.legend(loc=4)
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid(linestyle=':')
    plt.show()

def get_parser(times='1:4:15:65', bnames='U:B:V:R', tau_ph=2. / 3):
    import argparse
    from argparse import RawTextHelpFormatter

    description = 'Plot and Configuration for the bands. \n'
    description += band.print_bands(is_print=False)
    parser = argparse.ArgumentParser(description=description,formatter_class=RawTextHelpFormatter)
    parser.add_argument('-b', '--band',
                        required=False,
                        type=str,
                        dest="bnames",
                        help="<bnames>: string. Ex: U:B ")
    parser.add_argument('--add',
                        # nargs='?',
                        required=False,
                        type=str,
                        dest="add",
                        help="<nm:leftWv:rightWv:[color]:[lntype]> "
                            ", where nm - the name of new band, leftWv and rightWv are left and right wave range in AA. "
                            "Ex: f1150:1150:1900 . Default color is black, ln is solid line. "
                            "See colors https://matplotlib.org/stable/gallery/color/named_colors.html")
    
    parser.add_argument('--norm', nargs="?",
                            required=False,
                            const=True,
                            dest="is_norm",
                            metavar="normalize bands",
                            help="Normalize bands to max-min values")
    return parser

# def usage():
#     print("Usage:")
#     print("  bands.py [params]")
#     print("  -b <bands>: string, default: U-B-V-R-I, for example U-B-V")
#     print("  --add nm:leftWv:rightWv:color:lntype  string. "
#           "nm - the name of new band. leftWv and rightWv are left and right wave range in AA. "
#           "for example f1150:1150:1900. Default color is black, ln is solid line. "
#           "See colors https://matplotlib.org/stable/gallery/color/named_colors.html")
#     print("  -h  print usage")
#     print("   --- ")
#     band.print_bands()


def main():
    # import getopt
    import sys

    band.Band.load_settings()

    parser = get_parser()
    args, unknownargs = parser.parse_known_args()

    # if len(names) == 0:
    #     # logger.error(" No data. Use key '-i' ")
    #     parser.print_help()
    #     sys.exit(2)

    # try:
    #     opts, args = getopt.getopt(sys.argv[1:], "hb:")
    # except getopt.GetoptError as err:
    #     print(str(err))  # will print something like "option -a not recognized"
    #     usage()
    #     sys.exit(2)

    bands = []
    if args.bnames:
        bands = str(args.bnames).split(':')

    if args.add:
        c, ln = None, None
        opt = str(args.add).split(':')
        if len(opt) < 3:
            print('To add new band you should set up it correctly. See -h') 
            parser.print_help()
        else:
            name, lWv, rWv = opt[0:3]
            lWv, rWv = float(lWv), float(rWv)
            if len(opt) > 3:
                c = opt[3]
            if len(opt) > 4:
                ln = opt[4]
            band.band_add_new_bol(name, lWv, rWv, c=c, ln=ln, length=300)    
            band.print_bands()

            plot_bands([name])
        sys.exit(2)
    
    # for opt, arg in opts:
    #     if opt == '-b':
    #         bands = str(arg).split('-')
    #         continue
    #     elif opt == '-h':
    #         parser.print_help()
    #         sys.exit(2)

    band.print_bands()
    print("-"*50)

    if len(bands) > 0:
        try:
            plot_bands(bands, is_norm=args.is_norm)
        except AttributeError:
            parser.print_help()
            sys.exit(2)
    else:
        plot_UBVRI()
        plot_JHK()
        plot_griuz()
        plot_SWIFT()
        plot_PS1()
        plot_HSC()
        plot_HST()
        plot_Kepler()
        plot_Misc()


if __name__ == '__main__':
    main()
