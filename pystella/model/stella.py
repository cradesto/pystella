import os

__author__ = 'bakl'

stella_extensions = ('tt', 'swd', 'lbol', 'res', 'dat', 'ph', "mrt", 'eve', 'rho', 'xni', 'flx')


class Stella:
    def __init__(self, name, path='./', info=False):
        """Creates a Stella model instance.  Required parameters:  name."""
        self.name = name
        self.path = path  # path to model files
        if info:
            self.show_info()

    def __str__(self):
        return "%s, path: %s" % (self.name, self.path)

    def __repr__(self):
        return "%s, path: %s" % (self.name, self.path)
        # return "%s" % self.name

    def show_info(self):
        # ext = ('tt', 'ph', 'res', 'swd')
        for e in stella_extensions:
            fname = os.path.join(self.path, self.name + '.' + e)
            if os.path.isfile(fname):
                print("Exist %s-file: %s" % (e, fname))
            else:
                print("No %s-file: %s" % (e, fname))

    @property
    def Name(self):
        """
        Alias for self.name
        :return: name
        """
        return self.name

    @property
    def is_any_data(self):
        ext = ('tt', 'ph', 'res', 'swd')
        return any(map(os.path.isfile, [os.path.join(self.path, self.name + '.' + e) for e in ext]))

    @property
    def is_ph_data(self):
        fname = os.path.join(self.path, self.name + '.ph')
        return os.path.isfile(fname)

    @property
    def is_swd_data(self):
        fname = os.path.join(self.path, self.name + '.swd')
        return os.path.isfile(fname)

    @property
    def is_res_data(self):
        fname = os.path.join(self.path, self.name + '.res')
        return os.path.isfile(fname)

    @property
    def is_tt_data(self):
        fname = os.path.join(self.path, self.name + '.tt')
        return os.path.isfile(fname)

    @property
    def is_flx_data(self):
        fname = os.path.join(self.path, self.name + '.flx')
        return os.path.isfile(fname)

    def get_eve(self, name, path=None):
        from pystella.model import sn_eve
        if path is None:
            path = self.path
        eve = sn_eve.load_rho(os.path.join(path, name + '.rho'))
        return eve

    def get_res(self):
        from pystella.model.sn_res import StellaRes
        return StellaRes(self.name, self.path)

    def get_tt(self):
        from pystella.model.sn_tt import StellaTt
        return StellaTt(self.name, self.path)

    def get_swd(self):
        from pystella.model.sn_swd import StellaShockWaveDetail
        swd = StellaShockWaveDetail(self.name, self.path)
        return swd

    def get_flx(self):
        import pystella.model.sn_flx as flx
        res = flx.flx_reader(os.path.join(self.path, self.name + '.flx'))
        return res

    def get_ph(self, t_diff=1.005, t_beg=0.0, t_end=float('inf'), is_nfrus=True):
        import pystella.model.sn_ph as ph
        res = ph.read(self.name, self.path, t_diff=t_diff, t_beg=t_beg, t_end=t_end, is_nfrus=is_nfrus)
        return res
    #
    # # old
    # def read_tt_data(self):
    #     return self.get_tt().read()


def stella_ls(path, pattern='*', exts=stella_extensions):
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


def show_info(path, cond=lambda i: True):
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
        info = stella.get_tt().Info
        #         print(info.Data)
        if cond(info):
            # #         if 30 < info.R < 40:
            info.show()
