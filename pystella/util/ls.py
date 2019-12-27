import os
import sys
from pystella.model.stella import Stella, stella_extensions


def long(models, path='.', pars=('R', 'M', 'Mni', 'E'), sep=' |', is_verbose=False):
    """
    Print model list with parameters
        TODO  filter by parameters
    :param is_verbose:
    :param sep:
    :param models:  list of model names
    :param path: working directory. Default: ./
    :param pars: the list of printed parameters
    :return:
    """

    units = {'R': 'Rsun', 'M': 'Msun', 'Mni': 'Msun', 'E': 'FOE'}
    is_dict = type(models) is dict
    if is_dict:
        mnanes = list(models.keys())
    else:
        mnanes = models

    # print header
    s1 = '{:>40s} '.format('Name')
    s2 = '{:40s} '.format('-' * 40)
    for k in pars:
        s1 += '{} {:>3s}({:4s})'.format(sep, k, units[k])
        s2 += '{}  {:8s}'.format(sep, '-' * 8)
    s1 += '{}  {}'.format(sep, 'comment')
    s2 += '{}  {}'.format(sep, '-' * 8)
    print(s1)
    print(s2)
    # print("| %40s |  %8s |  %6s | %6s |  %s" % ('Name', 'R', 'M', 'E', 'comment'))
    # print("| %40s |  %8s |  %6s | %6s |  %s" % ('-' * 40, '-' * 8, '-' * 6, '-' * 6, '-' * 8))
    for mdl in mnanes:
        stella = Stella(mdl, path=path)
        exts = models[mdl] if is_dict else ''
        if stella.is_tt:
            info = stella.get_tt().Info
            try:
                s = '{:>40s} '.format(info.Name)
                for k in pars:
                    v = getattr(info, k)
                    s += '{}  {:8.3f}'.format(sep, v)
                s += '{}  {}'.format(sep, exts)
                print(s)
            except KeyError as ex:
                print("| %40s |  %8s |  %6s | %6s | %s  | KeyError: %s" % (info.Name, '', '', '', exts, ex))
            except:
                print("| %40s |  %8s |  %6s | %6s | %s  | Unexpected error: %s"
                      % (info.Name, '', '', '', exts, sys.exc_info()[0]))
        else:
            if is_verbose:
                print("| %40s |  %26s |  %s" % (stella.name, ' ', exts))


def summary(models, path='.', pars=('R', 'M', 'Mni', 'E'), sep='  ', is_verbose=False):
    """
    Print list main param
    :param models: list of model names
    :param path: working directory. Default: ./
    :param pars: the parameters for summary. Default: pars=('R', 'M', 'Mni', 'E')
    :param sep: string separation. Default: space
    :param is_verbose:
    :return:
    """
    data = {k: [] for k in pars}
    for mdl in models:
        stella = Stella(mdl, path=path)
        if stella.is_tt:
            tt_info = stella.get_tt().Info
            try:
                for k in data.keys():
                    v = getattr(tt_info, k)
                    if v not in data[k]:
                        data[k].append(v)
            except KeyError as ex:
                print(" KeyError: %s. %s" % (tt_info.Name, ex))
        else:
            if is_verbose:
                print("{0} No tt.".format(stella))
    if len(list(data.values())[0]) > 0:
        print('Summary in {}'.format(path))
        for k in data.keys():
            print("{:5s} {}".format(k, sep.join(map(str, sorted(data[k])))))
    else:
        if is_verbose:
            print('No tt-files in {}'.format(path))


def short(path='./', pattern='*', exts=stella_extensions):
    """Print information list about models in the path
    :param path: working directory
    :param pattern:  pathname, which must be a string containing a path specification.
    :param exts: set of extensions
    :return: None
   """
    from os import listdir
    from os.path import isfile, join
    import fnmatch

    files = [f for f in listdir(path) if isfile(join(path, f)) and fnmatch.fnmatch(f, pattern)]
    models = {}
    for f in files:
        name, ext = os.path.splitext(f)
        ext = ext.replace('.', '')
        if ext in exts:
            if name in models.keys():
                models[name] += ' ' + ext
            else:
                models[name] = ext
    return models


def info(path, cond=lambda i: True):
    """Print information list about models in the path
    :param path: working directory
    :param cond: condition function, like lambda i: 30 < i.M < 1000 and i.R > 100
    :return: None
   """
    from os import listdir
    from os.path import isfile, join

    files = [f for f in listdir(path) if isfile(join(path, f)) and f.endswith('.tt')]
    for f in files:
        # print 'Read: %s' % f
        name, ext = os.path.splitext(f)
        stella = Stella(name, path=path)
        ttinfo = stella.get_tt().Info
        #         print(info.Data)
        if cond(ttinfo):
            # #         if 30 < info.R < 40:
            ttinfo.show()
