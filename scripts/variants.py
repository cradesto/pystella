#!/usr/bin/python3
# -*- coding: utf-8 -*-
import argparse
import sys
import numpy as np


class Tickets(object):
    def __init__(self, n):
        self._n = n
        self._tickets = []

    @property
    def n(self):
        return self._n

    @property
    def length(self):
        return len(self._tickets)

    def get(self, i):
        return self._tickets[i]

    def add(self, new):
        is_same = False
        for i, t in enumerate(self._tickets):
            diff = t - new
            if len(diff[diff==0]) == len(new):
                is_same = True
                break
        if not is_same:
            self._tickets.append(new)
        return is_same

    def remove(self, i):
        self._tickets.pop(i)

    def diff(self, new):
        weight = 0
        n = 0
        for i, t in enumerate(self._tickets):
            w = Tickets.ticket_diff(t, new)
            if weight < w:
                weight = w
                n = i

        return n, weight

    def max_diff(self):
        weight = 0
        n_i, n_k = 0, 0
        for i, t1 in enumerate(self._tickets):
            for k, t2 in enumerate(self._tickets):
                if i != k:
                    w = Tickets.ticket_diff(t1, t2)
                    if w > weight:
                        weight = w
                        n_i, n_k = i, k
        return weight, n_i, n_k

    def print(self):
        print("Print the {:d}. ".format(self.n))
        # print(self._tickets)
        print('\n'.join([''.join(['{:2d}'.format(item) for item in row])
                         for row in self._tickets]))
        w, i, k = self.max_diff()
        print("Max weight: {}  for pair: {:d}-{:d}".format(w, i, k))

    @staticmethod
    def ticket_diff(t1, t2):
        diff = t1 - t2
        return len(diff[diff==0])


def variants(n_vars, n_tasks):
    """Returns variant."""
    i = 1
    v = np.ones(n_tasks, dtype=int)
    while i < n_vars**n_tasks:
        yield v
        for k in range(n_tasks):
            if v[k] == n_vars:
                v[k] = 1
            else:
                v[k] += 1
                break
        i = i + 1


def get_parser():
    parser = argparse.ArgumentParser(description='Calculate the variants of tasks.')

    parser.add_argument('-n', '--n_tasks',
                        required=False,
                        type=int,
                        default=6,
                        dest="n_tasks",
                        help="-n <tasks>: number of tasks ")
    parser.add_argument('-t', '--n_tickets',
                        required=False,
                        type=int,
                        default=36,
                        dest="n_tickets",
                        help="-t <tickets>: number of tickets ")
    parser.add_argument('-v', '--n_vars',
                        required=False,
                        type=int,
                        default=6,
                        dest="n_vars",
                        help="-v <variants>: number of variants ")
    # parser.add_argument('-h', '--help',
    #                     action='store_const',
    #                     const=True,
    #                     dest="is_help",
    #                     help="Show this help")
    return parser


def run(n_tickets, n_vars, n_tasks):

    tickets = Tickets(n_tickets)
    vars_possible = variants(n_vars, n_tasks)
    weight = n_tasks
    count = 0
    for var in vars_possible:
        count += 1
        if count % 1000 == 0:
            print("{:d} var: {} , weight={}".format(count, var, w))
        if tickets.length < n_tickets:
            tickets.add(var.copy())
        else:
            i_worse, w = tickets.diff(var)
            if w <= weight:
                weight = w
                tickets.remove(i_worse)
                tickets.add(var.copy())
    return tickets


def main():
    parser = get_parser()
    args, unknownargs = parser.parse_known_args()

    n_tasks = args.n_tasks
    n_vars = args.n_vars
    n_tickets = args.n_tickets
    #
    # if args.:
    #     parser.print_help()
    #     sys.exit(2)

    print("Search the best {:d} tickets in {:d} variants of {:d} tasks. ".format(n_tickets, n_vars, n_tasks))

    tickets = run(n_tickets, n_vars, n_tasks)
    tickets.print()


if __name__ == '__main__':
    main()
