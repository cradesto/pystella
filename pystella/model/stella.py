import os

from pystella.rf.reddening import ReddeningLaw, LawFitz

__author__ = 'bakl'

stella_extensions = ('tt', 'swd', 'lbol', 'res', 'dat', 'ph', "mrt", 'eve', 'rho', 'xni', 'flx')


class Stella:
    def __init__(self, name, path='./', info=False):
        """Creates a Stella model instance.  Required parameters:  name."""
        self.name = name
        self.path = path  # path to model files
        if info:
            self.info()

    def __str__(self):
        return "%s, path: %s" % (self.name, self.path)

    def __repr__(self):
        return "%s, path: %s" % (self.name, self.path)
        # return "%s" % self.name

    @property
    def Name(self):
        """
        Alias for self.name
        :return: name
        """
        return self.name

    @property
    def Path(self):
        """
        Alias for self.path
        :return: path
        """
        return self.path

    def is_any_data(self, ext=('tt', 'ph', 'res', 'swd')):
        return any(map(os.path.isfile, [os.path.join(self.path, self.name + '.' + e) for e in ext]))

    @property
    def is_ph(self):
        fname = os.path.join(self.path, self.name + '.ph')
        return os.path.isfile(fname)

    @property
    def is_swd(self):
        fname = os.path.join(self.path, self.name + '.swd')
        return os.path.isfile(fname)

    @property
    def is_res(self):
        fname = os.path.join(self.path, self.name + '.res')
        return os.path.isfile(fname)

    @property
    def is_tt(self):
        fname = os.path.join(self.path, self.name + '.tt')
        return os.path.isfile(fname)

    @property
    def is_flx(self):
        return self.is_any_data(ext=['flx'])
        # fname = os.path.join(self.path, self.name + '.flx')
        # return os.path.isfile(fname)

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

    def curves(self, bands, z=0., distance=10., ebv=0., Rv=None, law=LawFitz, mode=ReddeningLaw.SMC, **kwargs):
        """
            Compute magnitude in bands for the 'name' model.
        :param bands: photometric bands
        :param z: redshift, default 0
        :param distance: distance to star in parsec, default 10 pc
        :param ebv: color excess, default: 0
        :param Rv:  R_V
        :param law: the law of extinction curve
        :param mode: type of extinction curves (MW, LMC, SMC)
        :return: light curves for the photometric bands
        """
        from pystella.rf.light_curve_func import series_spec_reddening
        from pystella.rf.rad_func import pc_to_cm

        t_beg = kwargs.get("t_beg", 0.)
        t_end = kwargs.get("t_end", float('inf'))
        t_diff = kwargs.get("t_diff", 1.01)
        magnification = kwargs.get("magnification", 1.)

        if len(bands) == 0:
            raise ValueError("You have not set any bands for model: " + str(self))
        if not self.is_ph:
            self.info()
            raise ValueError("Error: No spectral data for: " + str(self))

        # Get SED(time)
        serial_spec = self.get_ph(t_diff=t_diff, t_beg=t_beg, t_end=t_end)
        # reddening
        if ebv > 0:
            ss = serial_spec.copy(wl_ab=law.LAMBDA_LIM)
            serial_spec = series_spec_reddening(ss, ebv=ebv, Rv=Rv, law=law, mode=mode)
        # light curves
        curves = serial_spec.flux_to_curves(bands, z=z, d=pc_to_cm(distance), magnification=magnification)

        return curves

    def info(self):
        """Print information list about models in the path
        :return: None
       """
        print('Stella model: {} path: {}'.format(self.Name, self.path))
        for e in stella_extensions:
            fname = os.path.join(self.path, self.name + '.' + e)
            if os.path.isfile(fname):
                print("Exist %s-file: %s" % (e, fname))
            else:
                print("No %s-file: %s" % (e, fname))
        if self.is_tt:
            info = self.get_tt().Info
            info.show()
        if self.is_res:
            info = self.get_res().Info
            info.show()


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
