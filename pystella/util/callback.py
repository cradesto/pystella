import os

import sys

__author__ = 'bakl'

ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
plugin_path = os.path.join(ROOT_DIRECTORY, 'plugin')


class CallBack(object):
    def __init__(self, func, path='./', load=1, a=None):
        self._path = path
        self._func = None
        self._funcName = func
        self._args = a
        if func is not None and load == 1:
            self._func = self.find_func

    @property
    def path(self):
        return self._path

    @property
    def funcName(self):
        return self._funcName

    @property
    def FuncFile(self):
        return self.funcName + '.py'

    @property
    def FuncFileFull(self):
        return os.path.join(self.path, self.FuncFile)

    def arg_totext(self, idx):
        if self._args is None:
            return ''
        if len(self._args) > idx:
            return self._args[idx]
        return ''

    def add_args(self, *args):
        self._args = args

    @property
    def find_func(self):
        possibles = globals().copy()
        possibles.update(locals())
        method = possibles.get(self.funcName)

        if not method and os.path.isfile(self.FuncFileFull):  # find in files
            sys.path.append(self.path)
            py_mod = __import__(self.funcName, fromlist=['run', 'plot'])
            if hasattr(py_mod, 'run'):
                method = getattr(__import__(self.funcName), 'run')
            elif hasattr(py_mod, 'plot'):
                method = getattr(__import__(self.funcName), 'plot')

        if not method:
            raise Exception("Method %s not implemented" % self._funcName)

        return method

    def plot(self, ax, *args):
        a = []
        if self._args is not None:
            a.extend(self._args)
        if len(args) > 0:
            a.extend(args)
        # a.append(args)
        self._func(ax, a)

    def run(self):
        self._func(self._args)