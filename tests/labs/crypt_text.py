#!/usr/bin/python
# -*- coding: utf-8 -*-

import getopt
import sys, os
from os.path import basename, isfile, join, dirname


def decrypt_text(crypt_text, dic_key):
    if len(dic_key) == 0:
        return crypt_text

    mapping = {k: '<%s>' % i for i, k in zip(range(len(dic_key)), dic_key.keys())}

    # заменить все ключи на теги: <1>
    for k, v in mapping.iteritems():
        crypt_text = crypt_text.replace(k, v)

    # заменить все вхождения ключей в текст на их значение
    for k, v in dic_key.iteritems():
        r = mapping[k]
        crypt_text = crypt_text.replace(r, v)

    return crypt_text


def usage():
    print "Usage:"
    print "  %s [params]" % basename(__file__)
    print "  -k <keys>: string, default: a-b_b_c_c-d"
    print "  -i <text file name>.  Example: puzzle_25.txt"


def main():
    fname = 'puzzle_25.txt'
    path = './'
    mapping = {}

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:k:")
    except getopt.GetoptError as err:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    if len(opts) == 0:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-k':
            set_pairs = str(arg).split('_')
            for pair in set_pairs:
                k, v = pair.split('-')
                k, v = map(lambda x: x.decode('utf-8'), (k, v))
                mapping[k] = v
            break
        if opt == '-i':
            fname = os.path.expanduser(str(arg))
            if not isfile(fname):
                print "No such file: " + fname
                sys.exit(2)
            continue

    with open(fname, 'r') as f:
        read_text = f.read()
        f.closed

    # key = {'a': 'b', 'b': 'c'}
    res = decrypt_text(read_text.decode('utf-8').lower(), mapping)
    print res


if __name__ == '__main__':
    main()
    # main(name="cat_R1000_M15_Ni007_E15")
