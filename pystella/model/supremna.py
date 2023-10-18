import os

from pystella.rf.reddening import ReddeningLaw, LawFitz

__author__ = 'bakl'

supremna_extensions = ('swd', 'res', 'dat', 'ioh')


class Supremna:
    def __init__(self, name, path='./', info=False):
        """Creates a Supremna model instance.  Required parameters:  name."""
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

    def is_any_data(self, ext=('ioh', 'res', 'swd')):
        return any(map(os.path.isfile, [os.path.join(self.path, self.name + '.' + e) for e in ext]))

    @property
    def is_ioh(self):
        fname = os.path.join(self.path, self.name + '.ioh')
        return os.path.isfile(fname)

    @property
    def is_swd(self):
        fname = os.path.join(self.path, self.name + '.swd')
        return os.path.isfile(fname)

    @property
    def is_res(self):
        fname = os.path.join(self.path, self.name + '.res')
        return os.path.isfile(fname)

    def get_eve(self, name=None, path=None, is_hyd_abn=False, **kwargs):
        from pystella.model import sn_eve
        if name is None:
            name = self.Name
        if path is None:
            path = self.Path
        if is_hyd_abn:
            eve = sn_eve.load_hyd_abn(name, path, **kwargs)
        else:
            eve = sn_eve.load_rho(os.path.join(path, name + '.rho'))
        return eve

    def get_res(self):
        from pystella.model.sn_res import SupremnaRes
        return SupremnaRes(self.name, self.path)

    def get_swd(self):
        from pystella.model.supr_swd import SupremnaShockWaveDetail
        swd = SupremnaShockWaveDetail(self.name, self.path)
        return swd

    def get_ion(self):
        from pystella.model import SupremnaIonHistory
        return SupremnaIonHistory(self.name, self.path)

    # def get_ioh(self):
    #     from pystella.model.sn_tau import SupremnaIonHistoryDetail
    #     ioh = SupremnaIonHistoryDetail(self.name, self.path)
    #     return ioh
    
    def info(self):
        """Print information list about models in the path
        :return: None
       """
        print('Supremna model: {} path: {}'.format(self.Name, self.path))
        for e in Supremna_extensions:
            fname = os.path.join(self.path, self.name + '.' + e)
            if os.path.isfile(fname):
                print("Exist %s-file: %s" % (e, fname))
            else:
                print("No %s-file: %s" % (e, fname))
        if self.is_res:
            info = self.get_res().Info
            info.show()

