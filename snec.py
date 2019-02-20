#!/usr/bin/env python3
# #!/usr/bin/python3

import getopt
import os
import sys

from pystella.model import snec


def usage():
    print("Convert SNEC to STELLA presupernova model.")
    print("Usage:")
    print("  snec.py  -c <chemfile> -e <profile>")
    print("  -c <chemfile> ")
    print("  -e <profile> ")
    print("  -p <model directory>, default: ./")
    print("  -h  print usage")
    print("   --- ")


def main():
    file_profile = None
    file_chem = None
    path = os.getcwd()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hc:p:e:")
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    if len(opts) == 0:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-e':
            file_profile = str.strip(arg)
            continue
        if opt == '-c':
            file_chem = str.strip(arg)
            continue
        if opt == '-p':
            path = os.path.expanduser(str(arg))
            if not (os.path.isdir(path) and os.path.exists(path)):
                print("No such directory: " + path)
                sys.exit(2)
            continue
        elif opt == '-h':
            usage()
            sys.exit(2)

    if not os.path.isfile(file_profile):
        file_profile = os.path.join(path, file_profile)
    if not os.path.isfile(file_chem):
        file_chem = os.path.join(path, file_chem)

    name = os.path.splitext(os.path.basename(file_profile))[0]

    prb_snec = snec.Problem(name).load_chem(file_chem)
    prb_snec.load_profile(file_profile)

    presn = snec.to_presn(prb_snec)

    fname = presn.Name + '.hyd'
    # fname = os.path.join(os.path.dirname(file_profile), presn.name + '.hyd')
    res = presn.write_hyd(fname)
    if res:
        print("Save in %s" % fname)
    else:
        print("Error while saving in %s" % fname)

    fname = presn.Name + '.abn'
    # fname = os.path.join(os.path.dirname(file_profile), presn.name + '.abn')
    res = presn.write_abn(fname)
    if res:
        print("Save in %s" % fname)
    else:
        print("Error while saving in %s" % fname)


if __name__ == '__main__':
    main()
