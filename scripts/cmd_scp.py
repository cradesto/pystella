#!/usr/bin/env python3
import os
import sys
import argparse


# read arguments
def get_parser():
    pfrom = "gray:/home/bakl/Sn/Release/seb_git/run/87a/strad/"
    parser = argparse.ArgumentParser(description=' Get files.')
    parser.add_argument('-f', '--from',
                        required=False,
                        dest="path_from",
                        default=pfrom,
                        help="-f <path_from>. Default: "+pfrom)

    parser.add_argument('-t', '--to',
                        required=False,
                        dest="path_to",
                        default="../",
                        help="-f <path_to>. Default: ../")
    parser.add_argument('-e', '--extension',
                        required=False,
                        dest="ext",
                        default=".ph",
                        help="-e <extension with point>. Default: .ph")
    return parser


def main():
    parser = get_parser()
    args, unknownargs = parser.parse_known_args()

    # Set model names
    names = []

    if len(unknownargs) > 0:
        path, name = os.path.split(unknownargs[0])
        path = os.path.expanduser(path)
        name = name.replace(args.ext, '')
        names.append(name)

    if len(names) == 0:
        parser.print_help()
        sys.exit(2)
    
    # print "This is the name of the script: ", sys.argv[0]
    print("Number of files: ", len(names))

    for i, fname in enumerate(names):
        cmd = "scp {}.{{swd,tt,ph}} {}".format(os.path.join(args.path_from, fname), args.path_to)
        print(cmd)
        os.system(cmd)
        fname = os.path.join(args.path_to, fname) + args.ext
        if os.path.isfile(fname):
            cmd = "ln -s {}".format(fname)
            print(cmd)
            os.system(cmd)
        else:
            print("no file {}".format(fname))
            #        print( subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout.read())

                
if __name__ == '__main__':
    main()
