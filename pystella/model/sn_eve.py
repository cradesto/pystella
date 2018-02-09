import logging
import numpy as np
import os

try:
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    is_matplotlib = True
except:
    is_matplotlib = False

from pystella.util.phys_var import phys

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

__author__ = 'bakl'

eve_elements = "Ni56 H He C N O Ne Na  Mg  Al  Si  S  Ar  Ca  Fe  Ni".split()
eve_colors = dict(Ni56="red", H="blue", He="cyan", C="darkorange", N="coral",
                  O="violet", Ne="green", Na="sandybrown",
                  Mg="skyblue", Si="olive", Al="lime",
                  S="indigo", Ar="brown", Ca="purple",
                  Fe='maroon', Ni='magenta')
eve_lntypes = dict((k, '--') for k, v in eve_colors.items())  # no y-shift
eve_lntypes['H'] = '-'
eve_lntypes['He'] = '-'
eve_lntypes['O'] = '-'
eve_lntypes['C'] = '-'
eve_lntypes['Ni56'] = '-'


class PreSN(object):
    """
    A class that holds data of presupernova
    """
    sRho = 'Rho'
    sM = 'M'
    sT = 'T'
    sR = 'R'
    sV = 'V'
    presn_hydro = (sM, sR, sT, sRho, sV)
    presn_elements = "H He C N O Ne Na  Mg  Al  Si  S  Ar  Ca  Fe  Ni Ni56".split()

    def __init__(self, name, nzon):
        """Creates a PreSN model instance.  Required parameters:  name, nzon"""
        self._name = name
        self._nzon = nzon
        self._data_hyd = np.empty(nzon, dtype=PreSN.dtype_hyd())
        self._data_chem = np.empty(nzon, dtype=PreSN.dtype_chem())
        self._params = {}
        self._loads = []

    @staticmethod
    def dtype_hyd():
        dt = np.dtype({'names': PreSN.presn_hydro, 'formats': np.repeat('f8', len(PreSN.presn_hydro))})
        return dt

    @staticmethod
    def dtype_chem():
        dt = np.dtype({'names': PreSN.presn_elements, 'formats': np.repeat('f8', len(PreSN.presn_elements))})
        return dt

    def show_info(self):
        print("-" * 20)
        print(" Name: %s    nzon: %d" % (self.Name, self.nzon))
        print(" m_tot: {:5.3f} r_cen: {:12.6e}".format(self.m_tot, self.r_cen))

    @property
    def Name(self):
        return self._name

    @property
    def nzon(self):
        """ time start"""
        return self._nzon

    @property
    def time_start(self):
        """ time start"""
        return self.par('time_start', 0.)

    @property
    def r_cen(self):
        """Center radius"""
        p = 'r_cen'
        if self.is_set(PreSN.sR):
            d = self.hyd(PreSN.sR)[0] / 2.  # todo check Rcen
        else:
            d = 0.
        return self.par(p, d)

    @property
    def m_tot(self):
        """Total mass"""
        p = 'm_tot'
        if self.is_set(PreSN.sM):
            d = self.hyd(PreSN.sM)[0]
        else:
            d = 0.
        return self.par(p, d)

    @property
    def rho_cen(self):
        """ Center density"""
        p = 'rho_cen'
        if self.is_set(PreSN.sRho):
            d = self.hyd(PreSN.sRho)[0]
        else:
            d = 0.
        return self.par(p, d)

    @property
    def lg_rho(self):
        """logarithmic density"""
        return np.log10(self.rho)

    @property
    def rho(self):
        """density"""
        return self.hyd(PreSN.sRho)

    @property
    def lg_r(self):
        """logarithmic radius"""
        return np.log10(self.r)

    @property
    def r(self):
        """logarithmic radius"""
        return self.hyd(PreSN.sR)

    @property
    def m(self):
        """Mass"""
        return self.hyd(PreSN.sM)

    @property
    def lgT(self):
        """Log T"""
        return np.log10(self.T)

    @property
    def T(self):
        """Temperature"""
        return self.hyd(PreSN.sT)

    @property
    def V(self):
        """Velocity"""
        return self.hyd(PreSN.sV)

    @property
    def lgV(self):
        """Log Velocity"""
        return np.log10(self.V)

    def hyd(self, v):
        """Hydro data"""
        if v not in self._loads:
            raise ValueError("There is no information about the parameter [%s]. You should set it." % v)
        return self._data_hyd[v]

    @property
    def chem(self):
        """Full hydro data"""
        return self._data_chem

    def par(self, name, d=None):
        return self._params.get(name, d)

    def set_par(self, name, v):
        self._params[name] = v

    def lg_el(self, el):
        return np.log10(self.el(el))

    def el(self, el):
        if el not in PreSN.presn_elements:
            raise ValueError("There is no  element [%s] in elements" % el)
        if el not in self._loads:
            raise ValueError("There is no information about the element [%s]. You should set it." % el)
        return self._data_chem[el]

    def write_hyd(self, fname):
        """
        Code readheger.trf:
          BM1=cutmass; -- core Mass
          r(0)=Rcen;
          dum=0.;
          write(12,'(1p,e12.3,i6,2e13.5)') timeStart, NzoneHyd, BM1, Rcen;
          do km=1,NzoneHyd;
            write(12,'(1x,i4,1p,e12.4,e18.10,3e15.7,2e12.4)')
               km, dum, rHyd(km), rhoHyd(km), TpHyd(km), uHyd(km), aMr(km), dum;
          enddo;
        :return:
        """
        dum = np.zeros(self.nzon)
        logger.info(' Write hyd-data to %s' % fname)
        zones = range(1, self._nzon + 1)
        with open(fname, 'w') as f:
            f.write('{:12.3e} {:6d} {:13.5e} {:13.5e} {:13.5e}\n'
                    .format(self.time_start, self.nzon, self.m_tot / phys.M_sun, self.r_cen, self.rho_cen))
            # a = '#No. Mr  dM  R  dR  Rho PRE T   V'.split()
            # f.write('           '.join(a)+'\n')
            # for _ in zip(zones, self.m/phys.M_sun, dum, self.r, dum, self.rho, dum, self.T, self.V):
            #     f.write(' %4d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e \n' % _)
            # 'evehyd.trf: idum,dum,Radius(j),RHOeve(j),TMPR(j),VELOC(j), dum,dum; '
            a = '#No. M  R  Rho T   V   M  dum '.split()
            f.write('  ' + '             '.join(a) + '\n')
            for _ in zip(zones, self.m / phys.M_sun, self.r, self.rho, self.T, self.V, self.m / phys.M_sun, dum):
                f.write(' %4d %15.7e %15.7e %15.7e %15.7e %15.7e  %15.7e  %8.1e\n' % _)
        return os.path.isfile(fname)

    def write_abn(self, fname, is_header=False, is_dum=False):
        """
        Write data to file in abn format.
        See code readheger.trf:
         _do km=1,NzoneHyd;
            write(13,'(i4,1p,19e10.3)')km,dum,dum,dum,
        -- No.  Mr    X(H He C N O Ne Na Mg Al Si S Ar Ca Fe Co Ni 56Ni)
              bh(km), bhe(km),bc(km),
              bn(km),bo(km),bne(km),bna(km),bmg(km),bal(km),bsi(km),
              bs(km),bar(km),bca(km),bfe(km),
              bni58(km),bni56(km); -- with Ni58 separated
         _od;
        :return:
        """
        dum = 0.
        logger.info(' Write abn-data to %s' % fname)
        with open(fname, 'w') as f:
            # f.write('%d\n' % self.nzon_abn)
            if is_header:
                if is_dum:
                    s = '%4s  %10s %10s %10s' % ('# zn', ' ', ' ', ' ')
                else:
                    s = '%4s' % '# zn'
                for ename in PreSN.presn_elements:
                    s += ' %10s' % ename
                f.write('%s\n' % s)
            for i in range(self.nzon):
                if is_dum:
                    s = '%4d  %10.3e %10.3e %10.3e' % (i + 1, dum, dum, dum)
                else:
                    s = '%4d' % (i + 1)
                for ename in PreSN.presn_elements:
                    s += ' %10.3e' % self.el(ename)[i]
                f.write('%s\n' % s)
        return os.path.isfile(fname)

    def plot_chem(self, x='m', elements=eve_elements, ax=None, xlim=None, ylim=None, **kwargs):
        if not is_matplotlib:
            return
        # elements = kwargs.get('elements', eve_elements)
        lntypes = kwargs.get('lntypes', eve_lntypes)
        colors = kwargs.get('colors', eve_colors)
        loc = kwargs.get('leg_loc', 8)
        leg_ncol = kwargs.get('leg_ncol', 4)
        lw = kwargs.get('lw', 2)
        marker = kwargs.get('marker', None)
        markersize = kwargs.get('markersize', 4)
        # def plot_chem(self, x='m', ax=None, elements=None, xlim=None, ylim=None,
        #               leg_loc=3, leg_ncol=4, lw=2, lntypes=None, is_save=False):
        #     if elements is None:
        #         elements = eve_elements
        #     if lntypes is None:
        #         lntypes = eve_lntypes

        is_new_plot = ax is None
        # setup figure
        if is_new_plot:
            plt.matplotlib.rcParams.update({'font.size': 14})
            fig = plt.figure(num=None, figsize=(8, 8), dpi=100, facecolor='w', edgecolor='k')

            gs1 = gridspec.GridSpec(1, 1)
            # gs1.update(wspace=0.1, hspace=0.1, top=0.97, left=0.12, right=0.98)
            gs1.update(wspace=0.1, hspace=0.1, top=0.97, left=0.12, right=0.87)
            ax = fig.add_subplot(gs1[0, 0])

            if is_new_plot:
                if x == 'r':
                    ax.set_xlabel(r'R [cm]')
                elif x == 'm':
                    ax.set_xlabel(r'M [$M_\odot$]')
                else:
                    ax.set_xscale('log')
                    ax.set_xlabel(r'R [cm]')

        is_x_lim = xlim is not None
        is_y_lim = ylim is not None

        if x == 'r':
            x = self.r
        elif x == 'm':
            x = self.m / phys.M_sun
        else:
            x = self.r

        y_min = []
        y_max = []
        for el in elements:
            if self.is_set(el):
                # y = self.lg_el(el)
                y = self.el(el)
                # x = y[np.nonzero(y)]
                # y = y[np.nonzero(y)]
                # y[y<=0] == 1e-15
                ax.plot(x, y, label='{0}'.format(el), color=colors[el], ls=lntypes[el], linewidth=lw
                            , marker=marker, markersize=markersize)
                # ax.semilogy(x, y, label='{0}'.format(el), color=colors[el], ls=lntypes[el], linewidth=lw
                #             , marker=marker, markersize=markersize)

                if not is_y_lim:
                    y_min.append(np.min(y))
                    y_max.append(np.max(y))

        if not is_y_lim and len(y_min) > 0:
            ylim = [np.min(y_min), np.max(y_min)]

        if not is_x_lim:
            xlim = np.min(x), np.max(x)

        if is_x_lim or not is_new_plot:
            ax.set_xlim(xlim)

        if is_y_lim or not is_new_plot:
            ax.set_ylim(ylim)
        ax.set_yscale('log')
        if is_new_plot:
            ax.set_ylabel(r'$X_i$')
            ax.legend(prop={'size': 9}, loc=loc, ncol=leg_ncol, fancybox=False, frameon=False, markerscale=0)
            # ax.legend(prop={'size': 9}, loc=3, ncol=4, fancybox=True, shadow=True)
            # plt.grid()
            # plt.show()
        return ax

    def plot_rho(self, x='m', ax=None, xlim=None, ylim=None, **kwargs):
        if not is_matplotlib:
            return
        lw = kwargs.get('lw', 2)
        ls = kwargs.get('ls', '-')
        label = kwargs.get('label', '')
        color = kwargs.get('color', 'black')
        marker = kwargs.get('marker', None)
        markersize = kwargs.get('markersize', 4)

        is_new_plot = ax is None
        # setup figure
        if is_new_plot:
            plt.matplotlib.rcParams.update({'font.size': 14})
            fig = plt.figure(num=None, figsize=(9, 5), dpi=100, facecolor='w', edgecolor='k')

            gs1 = gridspec.GridSpec(1, 1)
            gs1.update(wspace=0.1, hspace=0.1, top=None, left=0.13, right=0.98)
            ax = fig.add_subplot(gs1[0, 0])
            ax.set_ylabel(r'$\rho, [g/cm^3]$ ')

            if x == 'r':
                ax.set_xlabel(r'R [cm]')
            elif x == 'm':
                ax.set_xlabel(r'M [$M_\odot$]')
            else:
                ax.set_xscale('log')
                ax.set_xlabel(r'R [cm]')

        is_x_lim = xlim is not None
        is_y_lim = ylim is not None

        if x == 'r':
            xi = self.r
        elif x == 'm':
            xi = self.m / phys.M_sun
        else:
            xi = self.r

        y = self.rho
        ax.semilogy(xi, y, color=color, ls=ls, linewidth=lw, marker=marker, markersize=markersize, label=label)

        if is_new_plot:
            if not is_x_lim and len(xi) > 0:
                xlim = [np.min(xi), np.max(xi)]
                ax.set_xlim(xlim)
            if not is_y_lim and len(y) > 0:
                ylim = [np.min(y), np.max(y)]
                ax.set_ylim(ylim)

        return ax

    def is_set(self, name):
        return name in self._loads

    def set_hyd(self, name, vec, is_exp=False):
        if len(vec) != self.nzon:
            raise ValueError("The length of vector [%d] should be %d" % (len(vec), self.nzon))
        if name not in self._loads:
            self._loads.append(name)

        if is_exp:
            self._data_hyd[name] = 10. ** vec
        else:
            self._data_hyd[name] = vec

    def set_chem(self, name, vec, is_exp=False):
        if len(vec) != self.nzon:
            raise ValueError("The length of vector [%d] should be %d" % (len(vec), self.nzon))
        if name not in self._loads:
            self._loads.append(name)
        if is_exp:
            self._data_chem[name] = 10. ** vec
        else:
            self._data_chem[name] = vec

    def zone_reduce(self, by=sM, diff=1.01, start=0, end=None, mode='g'):
        """

        :param by: 'Rho' 'M' 'T' 'R' 'V', default: 'M'
        :param diff: geom progression, default: 1.01
        :param start:
        :param end:
        :param mode:
        :return:
        """
        from pystella.util.math import shrink, portion_index

        x = self.hyd(by)

        def where(a):
            return shrink(a, diff=diff, mode=mode)

        idxs = portion_index(x, where, start=start, end=end, isByEl=False)

        newPreSN = PreSN(self.Name, len(idxs))
        # hyd reshape
        for v in PreSN.presn_hydro:
            old = self.hyd(v)
            new = old[idxs]
            newPreSN.set_hyd(v, new)

        # abn reshape
        for el in PreSN.presn_elements:
            old = self.el(el)
            new = old[idxs]
            newPreSN.set_chem(el, new)

        return newPreSN

    def bad_zone_reduce(self, diff=1.05, start=0, end=None, mode='g'):
        from pystella.util.math import shrink
        x = self.m
        if end is None:
            end = len(x)
        idxs = np.arange(len(x))

        xx = x[start:end]
        idx = shrink(xx, diff=diff, mode=mode)
        idxs = np.concatenate((np.arange(start), idx, np.arange(end, len(x))))

        if start > 0:
            idxs = idxs[:start-1]
        else:
            idxs = []

        newPreSN = PreSN(self.Name, len(idxs))
        # hyd reshape
        for v in PreSN.presn_hydro:
            old = self.hyd(v)
            new = old[idxs]
            newPreSN.set_hyd(v, new)

        # abn reshape
        for el in PreSN.presn_elements:
            old = self.el(el)
            new = old[idxs]
            newPreSN.set_chem(el, new)

        return newPreSN

    def reshape(self, nz=300, nstart=0, nend=None, xmode='rlog', kind='np'):
        """
        Reshape parameters of envelope from nstart to nend to nz-zones
        :param nz: new zones
        :param nstart: zone number to start reshaping. Default: 0 (first zone)
        :param nend: zone number to end reshaping. Default: None,  (equal last zone)
        :param xmode: [lin OR rlog] - linear OR reversed log10
        :param kind: [np OR interp1d(..kind)], kind is  ('linear', 'nearest', 'zero', 'slinear', 'quadratic, 'cubic')
        :return: new preSN with reshaping zones
        """
        from scipy.interpolate import interp1d
        nznew = nstart + nz
        newPreSN = PreSN(self.Name, nznew)
        if nend is None:
            nend = self.nzon

        def rlogspace(s, e, n):
            r = np.exp(np.linspace(np.log(s), np.log(e), n))
            r = (e - r + s)
            return r[::-1]

        def interp(x, v, start=nstart, end=nend):
            res = []
            if start > 0:
                res = v[:start]  # save points before start
            xi = x[start:end]
            yi = v[start:end]
            if xmode == 'lin':
                xx = np.linspace(xi[0], xi[-1], nz)  # new x-points
            elif xmode == 'rlog':
                xx = rlogspace(xi[0], xi[-1], nz)  # new x-points
            else:
                xx = np.linspace(xi[0], xi[-1], nz)  # new x-points

            if kind == 'np':
                yy = np.interp(xx, xi, yi)
            else:
                interp_linear = interp1d(xi, yi, kind=kind)
                yy = interp_linear(xx)

            res = np.append(res, yy)
            return res

        # hyd reshape
        m = self.m
        for v in PreSN.presn_hydro:
            old = self.hyd(v)
            new = interp(m, old)
            newPreSN.set_hyd(v, new)

        # abn reshape
        for el in PreSN.presn_elements:
            old = self.el(el)
            new = interp(m, old)
            newPreSN.set_chem(el, new)

        return newPreSN


# ==============================================
def load_rho(fname, path=None):
    if path is not None:
        fname = os.path.join(path, fname)
    if not os.path.isfile(fname):
        logger.error(' No rho-data for %s' % fname)
        raise ValueError(' No rho-data for %s' % fname)
        # return None
    logger.info(' Load rho-data from  %s' % fname)
    col_names = "zone mass lgR lgTp lgRho u Ni56 H He C N O Ne Na  Mg  Al  Si  S  Ar  Ca  Fe  Ni"
    dt = np.dtype({'names': col_names.split(), 'formats': np.repeat('f8', len(col_names))})

    data = np.loadtxt(fname, comments='#', skiprows=2, dtype=dt)

    nz = len(data['lgR'])
    ###
    name = os.path.basename(os.path.splitext(fname)[0])
    col_map = {PreSN.sR: 'lgR', PreSN.sM: 'mass', PreSN.sT: 'lgTp', PreSN.sRho: 'lgRho', PreSN.sV: 'u'}
    presn = PreSN(name, nz)
    presn.set_hyd('V', np.zeros(nz))
    for k, v in col_map.items():
        presn.set_hyd(k, data[v], is_exp=v.startswith('lg'))

    # CGS
    presn.set_hyd('M', presn.m * phys.M_sun)

    for ename in PreSN.presn_elements:
        presn.set_chem(ename, data[ename], is_exp=True)

    return presn


def load_hyd_abn(name, path='.', is_dum=True):
    """
    Code readheger.trf:
      BM1=cutmass; -- core Mass
      r(0)=Rcen;
      dum=0.;
      write(12,'(1p,e12.3,i6,2e13.5)') timeStart, NzoneHyd, BM1, Rcen;
      do km=1,NzoneHyd;
        write(12,'(1x,i4,1p,e12.4,e18.10,3e15.7,2e12.4)')
           km, dum, rHyd(km), rhoHyd(km), TpHyd(km), uHyd(km), aMr(km), dum;
      enddo;
    :return: PreSN
    """

    # hydro
    ext_hyd = '.hyd'
    hyd_file = os.path.join(path, name + ext_hyd)
    if not os.path.isfile(hyd_file):
        logger.error(' No file for %s' % hyd_file)
        return None

    logger.info(' Load hyd-data from  %s' % hyd_file)

    # read table data
    if is_dum:
        col_names = "zone dm R Rho T V M dum2".split()
    else:
        col_names = "zone dm R Rho T V M".split()
    dt = np.dtype({'names': col_names,
                   'formats': ['i4'] + np.repeat('f8', len(col_names) - 1).tolist()})

    data_hyd = np.loadtxt(hyd_file, comments='#', skiprows=1, dtype=dt)
    nz = len(data_hyd['R'])
    data_hyd[PreSN.sM] = data_hyd[PreSN.sM] * phys.M_sun
    col_map = {PreSN.sR, PreSN.sM, PreSN.sT, PreSN.sRho, PreSN.sV}
    presn = PreSN(name, nz)
    for v in col_map:
        presn.set_hyd(v, data_hyd[v], is_exp=v.startswith('lg'))

    with open(hyd_file, 'r') as f:
        line = f.readline()
    if len(line) > 0:
        a = [float(x) for x in line.split()]
        if len(a) == 5:
            time_start, nzon, m_tot, r_cen, rho_cen = a
            presn.set_par('time_start', time_start)
            presn.set_par('m_tot', m_tot * phys.M_sun)
            presn.set_par('r_cen', r_cen)
            presn.set_par('rho_cen', rho_cen)
        elif len(a) == 4:
            time_start, nzon, m_tot, r_cen = a
            presn.set_par('time_start', time_start)
            presn.set_par('m_tot', m_tot * phys.M_sun)
            presn.set_par('r_cen', r_cen)
        elif len(a) == 2:
            time_start, nzon = a
            presn.set_par('time_start', time_start)

    # chemical composition
    ext_abn = '.abn'
    abn_file = os.path.join(path, name + ext_abn)
    if not os.path.isfile(abn_file):
        logger.error(' No file for %s' % abn_file)
        return None

    logger.info(' Load abn-data from  %s' % abn_file)
    abn_elements = 'H He C N O Ne Na  Mg  Al  Si  S  Ar  Ca  Fe  Ni Ni56'.split()
    abn_elements_iso = 'H He C N O Ne Na  Mg  Al  Si  S  Ar  Ca  Fe  Ni Ni56 Fe52 Cr48'.split()
    if is_dum:
        col_names = ("zone dum1 dum2 dum3 " + ' '.join(abn_elements)).split()
    else:
        col_names = ("zone " + ' '.join(abn_elements)).split()
    dt = np.dtype({'names': col_names, 'formats': np.repeat('f8', len(col_names))})
    data_chem = np.loadtxt(abn_file, comments='#', dtype=dt)

    for ename in abn_elements:
        presn.set_chem(ename, data_chem[ename])

    return presn
