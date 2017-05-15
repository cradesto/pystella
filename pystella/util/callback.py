import os

import sys

__author__ = 'bakl'

ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
plugin_path = os.path.join(ROOT_DIRECTORY, 'plugin')


def lc_wrapper(param, path=plugin_path):
    a = param.split(':')
    func = a.pop(0)
    c = CallBack(func, path=path, args=a, load=1)
    return c


class CallBack(object):
    def __init__(self, fname, path='./', args=None, method=None, load=1):
        self._func = None
        self._path = path
        self._fname = fname
        self._args = args
        if fname is not None and load == 1:
            if method is None:
                self._func = self.find_func
            else:
                self._func = self.find_method(method)

    @property
    def Path(self):
        return self._path

    @property
    def FileName(self):
        return self._fname

    @property
    def FuncFile(self):
        return self.FileName + '.py'

    @property
    def FuncFileFull(self):
        return os.path.join(self.Path, self.FuncFile)

    @property
    def Func(self):
        return self._func

    def get_arg(self, idx):
        if self._args is None:
            return None
        if len(self._args) > idx:
            return self._args[idx]
        return None

    def set_arg(self, idx, val):
        if idx < 0:
            raise Exception("Index should be more 0 [idx = %d]" % idx)
        if self._args is None and idx > 0:
            raise Exception("Index should be 0, if self._args is None [idx = %d]" % idx)
        if len(self._args) > idx:
            raise Exception("Index should be less then len(self._args)=%d [idx = %d]" % (len(self._args), idx))

        if self._args is None and idx == 0:
            self._args = []
            self._args[idx] = val
        elif len(self._args) > idx:
            self._args[idx] = val
        elif len(self._args) == idx:
            self._args.append(val)
        return self

    def add_arg(self, val):
        if self._args is None:
            idx = 0
        else:
            idx = len(self._args)
        self.set_arg(idx, val)

    def put_args(self, args):
        self._args = args

    def arg_totext(self, idx):
        res = self.get_arg(idx)
        if res is None:
            return ''
        return res

    @property
    def find_func(self):
        possibles = globals().copy()
        possibles.update(locals())
        method = possibles.get(self.FileName)

        if not method and os.path.isfile(self.FuncFileFull):  # find in files
            sys.path.append(self.Path)
            py_mod = __import__(self.FileName, fromlist=['run', 'plot', 'load'])
            if hasattr(py_mod, 'run'):
                method = getattr(__import__(self.FileName), 'run')
            if hasattr(py_mod, 'load'):
                method_load = getattr(__import__(self.FileName), 'load')
            elif hasattr(py_mod, 'plot'):
                method = getattr(__import__(self.FileName), 'plot')

        if not method:
            raise Exception("Method %s not implemented" % self._fname)

        return method

    def find_method(self, method_name):
        method = None
        if os.path.isfile(self.FuncFileFull):  # find in files
            sys.path.append(self.Path)
            py_mod = __import__(self.FileName, fromlist=[method_name])
            if hasattr(py_mod, method_name):
                method = getattr(__import__(self.FileName), method_name)

        if not method:
            raise Exception("Method %s not implemented" % self._fname)

        return method

    @property
    def find_load(self):
        possibles = globals().copy()
        possibles.update(locals())
        method = possibles.get(self.FileName)

        if not method and os.path.isfile(self.FuncFileFull):  # find in files
            sys.path.append(self.Path)
            py_mod = __import__(self.FileName, fromlist=['load'])
            if hasattr(py_mod, 'load'):
                method_load = getattr(__import__(self.FileName), 'load')

        if not method:
            raise Exception("Method load in %s not implemented" % self._fname)

        return method

    def plot(self, ax, dic=None):

        if dic is None:
            dic = {}

        if self._args is not None:
            dic['args'] = self._args[:]
        # if len(args) > 0:
        #     dic..extend(args)
        # a.append(args)
        self._func(ax, dic)

    def run(self):
        self._func(self._args)

    def load(self, dic=None):
        if dic is None:
            dic = {}
        if self._args is not None:
            dic['args'] = self._args[:]
        return self._func(dic)


class CallBackArray(CallBack):
    def __init__(self, calls, fname=None):
        super(CallBackArray, self).__init__(fname)
        self._calls = calls

    def get_arg(self, idx):
        return super(CallBackArray, self).get_arg(idx)

    def set_arg(self, idx, val):
        super(CallBackArray, self).set_arg(idx, val)
        return self

    def add_arg(self, val):
        super(CallBackArray, self).add_arg(val)

    def put_args(self, args):
        super(CallBackArray, self).put_args(args)

    def plot(self, ax, dic=None):
        if dic is None:
            dic = {}

        if self._args is not None:
            dic['args'] = self._args[:]

        for c in self._calls:
            c.plot(ax, dic)

    def run(self):
        for c in self._calls:
            c.run(self._args)
