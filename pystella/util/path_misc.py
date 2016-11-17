import os
from os.path import isfile, join

__author__ = 'bakl'


def get_model_names(path, model_ext):
    files = [f for f in os.listdir(path) if isfile(join(path, f)) and f.endswith(model_ext)]
    names = []
    for f in files:
        names.append(os.path.splitext(f)[0])
    return names
