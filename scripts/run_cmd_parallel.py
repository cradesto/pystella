#!/usr/bin/python3
# -*- coding: utf-8 -*-


import argparse
import os
import subprocess
import sys


def get_parser():
    parser = argparse.ArgumentParser(description='Run <command> in parallel.\n'
                                     'Ex: ./%(prog)s ./xronfictd.exe -n 11',
                                     usage='%(prog)s <command> [options]')
    
    parser.add_argument('-n', '--node',
                        required=False,
                        type=int,
                        default=1,
                        dest="nodes",
                        help="-n <nodes>: number of processes ")
    parser.add_argument('-m', '--max',
                        required=False,
                        type=int,
                        default=6,
                        dest="nmax",
                        help="--max <nmax>: the maximum number of concurrent processors. It could be less <nodes>. ")
    return parser


def main():
    n = 1
    max_processes = 6
    command = None
    #    command = "./xronfictd.exe"

    parser = get_parser()
    args, unknownargs = parser.parse_known_args()
    if len(unknownargs) > 0:
        command = 'nice -n 20 ' + unknownargs[0]

    if args.nodes:
        n = int(args.nodes)

    if args.nmax:
        max_processes = int(args.nmax)

    if command is None:
        parser.print_help()
        sys.exit(2)
                        
    processes = set()

    print("mother pid {0}".format(os.getpid()))
    for i in range(n):
        a = [command, '-i', str(i+1)]
        print("run # {0}".format(' '.join(a)))
        processes.add(subprocess.Popen(' '.join(a), shell=True))
        if len(processes) >= max_processes:
            os.wait()
            processes.difference_update([
                p for p in processes if p.poll() is not None])
            print("mother still alive with pid {0}".format(os.getpid()))


if __name__ == '__main__':
    main()
