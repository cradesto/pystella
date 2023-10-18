import logging
import numpy as np
import os
from pystella.model.sn_eve import PreSN
from pystella.util.phys_var import phys

logger = logging.getLogger(__name__)
__author__ = 'bakl'


class PreSupremna(PreSN):
    """
    A class that holds data of preSupremna
    """
    sTe = 'Te'
    sTi = 'Ti'
    presupr_hydro = (PreSN.sM, PreSN.sR, sTe, sTi, PreSN.sRho, PreSN.sV)

    def __init__(self, name, nzon, elements=PreSN.stl_elements):
        """Creates a PreSupremna model instance.  Required parameters:  name, nzon"""
        super().__init__(name, nzon, elements)
        self._data_hyd = np.empty(nzon, dtype=PreSupremna.dtype_hyd())

    @property
    def T(self):
        """Temperature as Te"""
        return self.Te
    
    @property
    def Te(self):
        """Electron temperature"""
        return self.hyd(PreSupremna.sTe)
    
    @property
    def Ti(self):
        """Ion temperature"""
        return self.hyd(PreSupremna.sTi)

    @staticmethod
    def dtype_hyd():
        dt = np.dtype({'names': PreSupremna.presupr_hydro, 'formats': np.repeat('f8', len(PreSupremna.presupr_hydro))})
        return dt
    
    # ==============================================
    @staticmethod
    def load_rho(fname, path: str = None):
        import linecache
        if path is not None:
            fname = os.path.join(path, fname)
        if not os.path.isfile(fname):
            logger.error(' No rho-data for %s' % fname)
            raise ValueError(' No rho-data for %s' % fname)
            # return None
        logger.info(' Load rho-data from  %s' % fname)

        # try to get column names from header
        header = linecache.getline(fname, 1).split()
        logger.info('header(len= {}): {} '.format(len(header), ' '.join(header)))

        # if len(header) < 22 or len(header) >= 25:  # default header
        #     header = "zone mass lgR lgTpe lgTpi lgRho u Ni56 H He C N O Ne Na  Mg  Al  Si  S  Ar  Ca  Fe  Ni Ni56 Ecr".split()

        usecols, col_names = [], []
        for i, e in enumerate(header):
            if not e in col_names:
                col_names.append(e)
                usecols.append(i)

        logger.info(' col_names=  %s' % ' '.join(col_names))
        # logger.info(' usecols=  {}'.format(usecols))
        dt = np.dtype({'names': col_names, 'formats': np.repeat('f8', len(col_names))})

        data = np.loadtxt(fname, comments='#', skiprows=5, dtype=dt, usecols=usecols)

        nz = len(data['lgR'])
        ###
        name = os.path.basename(os.path.splitext(fname)[0])
        col_map = {PreSupremna.sR: 'lgR', PreSupremna.sM: 'mass', PreSupremna.sTe: 'lgTpe', PreSupremna.sTi: 'lgTpi', 
                PreSupremna.sRho: 'lgRho', PreSupremna.sV: 'u'}
        presupr = PreSupremna(name, nz)
        presupr.set_hyd('V', np.zeros(nz))
        for k, v in col_map.items():
            presupr.set_hyd(k, data[v], is_exp=v.startswith('lg'))

        # CGS
        presupr.set_hyd('M', presupr.m * phys.M_sun)

        for ename in presupr.Elements:
            presupr.set_chem(ename, data[ename], is_exp=True)

        return presupr    