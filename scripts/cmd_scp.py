#!/usr/bin/env python3
import os
import sys

# print "This is the name of the script: ", sys.argv[0]
print("Number of files: ", len(sys.argv) - 1)

if len(sys.argv) < 2:
    print("Usage:")
    print(" >{}   file[s] ".format(sys.argv[0]))
    sys.exit()

path_from = "gray:/home/bakl/Sn/Release/seb_git/run/87a/strad/"
path_to = "../"

for i, fname in enumerate(sys.argv):
    if i > 0:
        cmd = "scp {}.{{swd,tt,ph}} {}".format(os.path.join(path_from, fname), path_to)
        print(cmd)
        os.system(cmd)
        cmd = "ln -s {}.ph".format(os.path.join(path_to, fname))
        print(cmd)
        os.system(cmd)
    #        print( subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout.read())
