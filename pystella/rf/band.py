import os
from os.path import dirname

import numpy as np

from pystella.rf.rad_func import MagAB2Flux, Flux2MagAB
from pystella.util.phys_var import phys

__author__ = 'bakl'


# see bands: http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse&gname=GALEX
# see personal page Brad Tucker: http://www.mso.anu.edu.au/~brad/filters.html


class Band(object):
    IsLoad = False
    Cache = dict()
    Alias = None
    FileFilters = 'filters.ini'
    FileSettings = 'settings.ini'
    NameBol = 'bol'
    NameUBVRI = 'ubvri'
    NameBolQuasi = 'bolq'
    NameZp = 'zp'
    NameJy = 'Jy'
    DirRoot = os.path.join(dirname(dirname(dirname(os.path.realpath(__file__)))), 'data/bands')

    def __init__(self, name=None, fname=None, zp=None, jy=None, is_load=False):
        """Creates a band instance.  Required parameters:  name and file."""
        self._fwhm = None
        self.name = name
        self._fname = fname  # location of the filter response
        self._zp = zp  # zero points for mag
        self._jy = jy  # zero points for flux [Jansky] # 1 Jy = 1.51e7 photons sec^-1 m^-2 (dlambda/lambda)^-1
        self.sync_zp()
        self.__freq = None  # frequencies of response [Hz]
        self.__wl = None  # wavelength of response [cm]
        self._resp = None  # response
        self._is_load = False
        self._wl_eff = None  # the effective wavelength
        self._norm = None
        self._normWl = None
        if is_load:  # and os.path.isfile(self.fname):
            self.load()

    @property
    def is_load(self):
        return self._is_load

    @property
    def fname(self):
        return self._fname

    @property
    def zp(self):
        if self._zp is None:
            return 0.
        return self._zp

    @property
    def is_zp(self):
        return self._zp is not None

    @property
    def Jy(self):
        if self._jy is None:
            return 0.
        return self._jy

    @property
    def is_Jy(self):
        return self._jy is not None

    @property
    def freq(self):
        return self.__freq

    @property
    def resp_fr(self):
        return self._resp[::-1]

    @property
    def resp_wl(self):
        return self._resp

    @property
    def wl(self):
        return self.__wl

    @property
    def wl2args(self):
        return self.__wl * phys.cm_to_angs

    @property
    def wlrange(self):
        return np.min(self.wl), np.max(self.wl)

    @property
    def wl_eff(self):
        """The effective wavelength"""
        if self._wl_eff is None:
            from scipy.integrate import simps
            wl = self.wl
            resp = np.array(self.resp_wl)
            num = simps(resp * wl, wl)
            den = simps(resp, wl)
            self._wl_eff = num / den
        return self._wl_eff

    @property
    def fwhm(self):
        """The FWHM"""
        if self._fwhm is None:
            from scipy.interpolate import interp1d
            x = self.wl
            y = np.array(self.resp_wl)
            imax = np.argmax(y)

            def val(xx, yy, v=0.5):
                # f = interp1d(xx, yy, kind='quadratic')
                # xx = np.exp(np.linspace(np.log(min(x)), np.log(max(x)), 100))
                # yy = f(xx)
                # return np.interp(v, yy, xx)
                return np.interp(v, yy, xx)

            wl_l = val(x[:imax], y[:imax])
            wl_r = val(x[:imax:-1], y[:imax:-1])
            self._fwhm = wl_r - wl_l
        return self._fwhm

    @property
    def wl_eff_angs(self):
        """The effective wavelength"""
        return self.wl_eff * phys.cm_to_angs

    @property
    def freq_eff(self):
        """The effective freq"""
        return phys.c / self.wl_eff

    @property
    def Name(self):
        return self.name

    @property
    def Norm(self):
        if self._norm is None:
            from scipy.integrate import simps
            nu_b = np.array(self.freq)
            resp_b = np.array(self.resp_fr)
            self._norm = simps(resp_b / nu_b, nu_b)
        return self._norm

    @property
    def NormWl(self):
        if self._normWl is None:
            # from scipy.integrate import simps
            # x = np.array(self.wl)
            # y = np.array(self.resp_wl) / (phys.c * phys.cm_to_angs * phys.h)
            # res = simps(y*x, x, even='avg')
            self._normWl = self.response(self.wl, np.ones(len(self.wl)), kind='spline')
        return self._normWl

    @property
    def wave_range(self):
        if self.wl is not None:
            return np.min(self.wl), np.max(self.wl)
        else:
            return None

    @wl.setter
    def wl(self, wl):
        self.__wl = wl
        self.__freq = phys.c / self.__wl[::-1]

    @resp_wl.setter
    def resp_wl(self, v):
        if np.any(v > 1.):
            raise ValueError("Band response  must be <=1 : " + str(self))
        self._resp = v

    def __str__(self):
        return "%s" % self.name

    def __repr__(self):
        return "%s" % self.name

    def load(self):
        """Reads the waves and response from file. Read zero point"""
        if self._fname is not None:
            f = open(self._fname)
            try:
                lines = [str.strip(line) for line in f.readlines() if not line.startswith('#')]
                wl = np.array([float(line.split()[0]) for line in lines if len(line) > 0])

                self._resp = np.array([float(line.split()[1]) for line in lines if len(line) > 0])
                self.wl = wl * phys.angs_to_cm

                self._is_load = True

            # except Exception:
            #     print"Error in band file: %s.  Exception:  %s" % (self.file, sys.exc_info()[0])
            #     sys.exit("Error while parse band-file: %s in %s" % self.file)
            finally:
                f.close()

    def sync_zp(self):
        """Compute Jy or zp if one of them is None"""
        if self.is_Jy and self.is_zp:
            return
        if not self.is_Jy and not self.is_zp:
            return
        if self.is_zp:
            self._jy = MagAB2Flux(self.zp) * 1e23
        if self.is_Jy:
            self._zp = Flux2MagAB(self.Jy * 1e-23)

    def clone(self, name):
        """
            Clone this band with new name
        :param name: the name of new cloned band
        :return:  new Band object
        """
        b = Band(name, self.fname)
        b._zp = self._zp
        b._jy = self._jy
        if self.is_load:
            b.wl = self.wl
            b.resp_wl = self.resp_wl
            b._is_load = self.is_load
        return b

    @classmethod
    def load_settings(cls, is_force=False):
        if Band.IsLoad and not is_force:
            return True

        from configparser import ConfigParser
        parser = ConfigParser()
        parser.optionxform = str
        fini = os.path.join(Band.DirRoot, Band.FileSettings)
        parser.read(fini)
        if 'alias' in parser.sections():
            Band.Alias = {k: v for k, v in parser.items('alias')}
        Band.IsLoad = True
        return True

    @staticmethod
    def response_nu(nu, flux, b, is_freq_norm=True):
        from scipy import integrate
        """
        Compute response flux using provided spectral band.
        see: http://www.astro.ljmu.ac.uk/~ikb/research/mags-fluxes/
        :param nu: the frequencies of SED
        :param flux: the flux of SED
        :param b:  photometric band
        :return:
        """
        from pystella.util.math import log_interp1d

        # nu_s = self.Freq
        # sort
        sorti = np.argsort(nu)
        nu_s = nu[sorti]
        flux_s = flux[sorti]

        nu_b = np.array(b.freq)
        resp_b = np.array(b.resp_fr)

        if np.min(nu_s) > nu_b[0] or np.max(nu_s) < nu_b[-1]:
            # decrease wave length range of the band
            f = (np.min(nu_s) < nu_b) & (nu_b < np.max(nu_s))

            # f = map(lambda x: min(nu_s) < x < max(nu_s), nu_b)
            nu_b = nu_b[f]
            if len(nu_b) < 3:
                raise ValueError("The filter bandwidth [{}] should be included in the bandwidth of the SED. "
                                 "There are only {} points in the filter wavelength range."
                                 .format(b.Name, len(nu_b)))
            resp_b = resp_b[f]

        log_interp = log_interp1d(nu_s, flux_s)
        flux_interp = log_interp(nu_b)

        if is_freq_norm:
            a = integrate.simps(flux_interp * resp_b / nu_b, nu_b)
        else:
            a = integrate.simps(flux_interp * resp_b, nu_b)
        return a

    def response_freq(self, nu, flux):
        return Band.response_nu(nu, flux, self)

    def response(self, wl, flux, z=0., kind='spline', mode_int='simps', is_out2zero=True, is_photons=True):
        """Compute the response of this filter over the flux(wl)

        kind : str or int, optional (see scipy.interpolate)
                Specifies the kind of interpolation as a string
                ('linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic',
                'previous', 'next', where 'zero', 'slinear', 'quadratic' and 'cubic'
                refer to a spline interpolation of zeroth, first, second or third
                order; 'previous' and 'next' simply return the previous or next value
                of the point) or as an integer specifying the order of the spline
                interpolator to use.
                Default is 'linear'
        """
        import scipy.integrate
        import scipy.interpolate
        from pystella.util.math import log_interp1d

        if len(wl) == 0 or len(flux) == 0:
            raise ValueError("It should be len(wl)>0 and len(flux) > 0. "
                             "Now len(wl)={},len(flux)={}".format(len(wl), len(flux)))
        if len(wl) != len(flux):
            raise ValueError("It should be len(wl) == len(flux). "
                             "Now len(wl)={},len(flux)={}".format(len(wl), len(flux)))

        if z > 0:
            wl_z = wl * (1. + z)
        elif z < 0:
            wl_z = wl / (1. + z)
        else:
            wl_z = wl

        x_min, x_max = self.wlrange

        if (x_min < wl_z[0] or x_max > wl_z[-1]) and not is_out2zero:
            return None

        i_min = 0
        try:
            i_min = np.nonzero(np.greater(wl_z - x_min, 0))[0][0]
        except:
            pass

        i_max = len(wl_z) - 1
        try:
            i_max = np.nonzero(np.greater(wl_z - x_max, 0))[0][0]
        except:
            pass

        if i_min >= 5:
            i_min -= 5
        else:
            i_min = 0

        if i_max <= len(wl_z) - 6:
            i_max += 5
        else:
            i_max = len(wl_z) - 1

        trim_flux = flux[i_min:i_max + 1]
        trim_wl = wl_z[i_min:i_max + 1]

        if kind == "spline":
            interp = scipy.interpolate.splrep(self.wl, self.resp_wl, k=1, s=0)
            resp_interp = scipy.interpolate.splev(trim_wl, interp)
        elif kind == 'log':
            interp = log_interp1d(self.wl, self.resp_wl)
            resp_interp = interp(trim_wl)
        elif kind == 'line':
            resp_interp = np.interp(trim_wl, self.wl, self.resp_wl, 0, 0)  # One-dimensional linear interpolation.
        # else:
        #     raise ValueError('No such mode: ' + mode)
        else:
            interp = scipy.interpolate.interp1d(self.wl, self.resp_wl, kind=kind)
            resp_interp = interp(trim_wl)

        if is_out2zero:
            resp_interp = np.where(np.less(trim_wl, x_min), 0, resp_interp)
            resp_interp = np.where(np.greater(trim_wl, x_max), 0, resp_interp)

        integrand = resp_interp * trim_flux

        if is_photons:
            integrand = integrand * trim_wl / (phys.c * phys.cm_to_angs * phys.h)

        if mode_int == 'simps':
            res = scipy.integrate.simps(integrand, x=trim_wl, even='avg')
        elif mode_int == 'trapz':
            res = scipy.integrate.trapz(integrand, x=trim_wl)
        else:
            res = (trim_wl[-1] - trim_wl[0]) / (len(trim_wl) - 1) * sum(integrand)

        return res


class BandUni(Band):
    def __init__(self, name='bol', wlrange=(1e1, 5e4), length=300):
        """Creates a band with uniform responce.
        :param name:  default 'Uniform'.
        :param wlrange:  the wavelength range, default (1e1, 5e4) [A]
        :param length: numbers of bins default 100
        """
        super().__init__(name)
        wl = np.exp(np.linspace(np.log(wlrange[0] * phys.angs_to_cm),
                                np.log(wlrange[1] * phys.angs_to_cm), length))
        self.wl = wl  # wavelength of response [cm]
        self.resp_wl = np.ones(len(wl))  # response

    @property
    def Norm(self, mode_int='simps'):
        import scipy.integrate  # import simps as integralfunc
        x = np.array(self.wl)
        y = np.array(self.resp_wl)
        # resp_b = np.ones(len(nu_b))
        if mode_int == 'simps':
            res = scipy.integrate.simps(y, x=x, even='avg')
        elif mode_int == 'trapz':
            res = scipy.integrate.trapz(y, x=x)
        else:
            raise ValueError('Unknown type for integration: mode_int= ' + mode_int)

        return res

        # return 1.

    @classmethod
    def get_bol(cls):
        return BandUni(name=Band.NameBol, wlrange=(1e0, 42e3), length=300)

    @classmethod
    def get_qbol(cls, Ushort=3250., Rlong=8900., length=300):
        return BandUni(name=Band.NameBolQuasi, wlrange=(Ushort, Rlong), length=length)


class BandJoin(Band):
    def __init__(self, bnames, name=None, length=300, is_norm=True, is_sum=False, is_load=False):
        """Creates a band from the join response for input bands.
        :param bnames:  input bands
        :param name:  default from bnames
        :param length: numbers of bins default 300
        :param is_norm: normalize to 1  default False
        :param is_sum: if True response is SUM, else ENVELOPE default False
        """
        self._is_sum = is_sum
        if name is None:
            name = ':'.join(bnames)
        super().__init__(name)

        self._bnames = bnames
        self._length = length
        self._is_norm = is_norm

        # self.wl = None  # wavelength of response [cm]
        # self.resp_wl = None  # response
        if is_load:
            self.load()

    def load(self):
        from scipy import interpolate

        wl = []
        # resps = []
        bands_int = {}

        for bn in self._bnames:
            b = band_by_name(bn)
            wl.extend(b.wl)
            # resps.extend(b.resp_wl)
            bands_int[bn] = interpolate.interp1d(b.wl, b.resp_wl, bounds_error=False, fill_value=(0., 0.))

        wl_min, wl_max = min(wl), max(wl)
        # print('wl_min, wl_max = ', wl_min, wl_max)
        wl_int = np.exp(np.linspace(np.log(wl_min), np.log(wl_max), self._length))
        resp_int = np.zeros_like(wl_int)
        # wl_int

        if self._is_sum:  # Response is sum for all bands
            for bn in self._bnames:
                f = bands_int[bn]
                resp_int += [f(w) for w in wl_int]
        else:  # Response is max among bands
            resps = []
            for bn in self._bnames:
                f = bands_int[bn]
                resps.append([f(w) for w in wl_int])
            resp_int = np.max(np.array(resps).transpose(), axis=1)
            # for i, w in enumerate(wl_int):
            #     resp_int[i] = np.max([bands_int[bn](w)])

        if self._is_norm:
            resp_int = resp_int / max(resp_int)  # normalization

        self.wl = wl_int  # wavelength of response [cm]
        self.resp_wl = resp_int  # response

    # @property
    # def Norm(self, mode_int='simps'):
    #     import scipy.integrate  # import simps as integralfunc
    #     x = np.array(self.wl)
    #     y = np.array(self.resp_wl)
    #     # resp_b = np.ones(len(nu_b))
    #     if mode_int == 'simps':
    #         res = scipy.integrate.simps(y, x=x, even='avg')
    #     elif mode_int == 'trapz':
    #         res = scipy.integrate.trapz(y, x=x)
    #     else:
    #         raise ValueError('Unknown type for integration: mode_int= ' + mode_int)
    #
    #     return res

    @classmethod
    def get_ubvri(cls, length=300):
        return BandJoin(name=Band.NameUBVRI, bnames=('U', 'B', 'V', 'R', 'I'), length=length,
                        is_norm=True, is_sum=False)


ROOT_DIRECTORY = dirname(dirname(dirname(os.path.abspath(__file__))))


#
# def get_full_path(fname):
#     return os.path.join(ROOT_DIRECTORY, fname)

# @lru_cache
def colors(bname=None, default='magenta'):
    c = {'U': "blue", 'B': "cyan", 'V': "darkgreen", 'R': "red", 'I': "magenta",
         'J': "blueviolet", 'H': "plum", 'K': "saddlebrown",
         'BesU': "blue", 'BesB': "cyan", 'BesV': "darkgreen", 'BesR': "red", 'BesI': "magenta",
         'SwiftU': "blue", 'SwiftB': "cyan", 'SwiftV': "darkgreen",
         'KaitU': "blue", 'KaitB': "cyan", 'KaitV': "darkgreen", 'KaitR': "red", 'KaitI': "magenta",
         'massJ': "green", 'massH': "cyan", 'massK': "black",
         'LcoJ': "green", 'LcoH': "cyan", 'LcoK': "black",
         'UVM2': "skyblue", 'UVW1': "orange", 'UVW2': "mediumpurple",
         'UVM2o': "deepskyblue", 'UVW1o': "darkviolet", 'UVW2o': "darkblue",
         'USNO40i': "blue", 'USNO40g': "cyan", 'USNO40z': "darkgreen", 'USNO40r': "darkviolet", 'USNO40u': "magenta",
         'F105W': "magenta", 'F435W': "skyblue", 'F606W': "cyan", 'F125W': "g",
         'F140W': "orange", 'F160W': "r", 'F814W': "blue",
         'Kepler': "magenta", 'Lum': "orchid",
         'g': "olive", 'r': "red", 'i': "plum", 'u': "darkslateblue", 'z': "chocolate",
         'Sdssg': "olive", 'Sdssr': "pink", 'Sdssi': "magenta", 'Sdssu': "blue", 'Sdssz': "chocolate",
         'PS1g': "olive", 'PS1r': "red", 'PS1i': "magenta", 'PS1u': "blue", 'PS1z': "chocolate", 'PS1y': "cyan",
         'PS1w': "orange", 'y': 'y', 'Y': 'y', 'w': 'tomato',
         'AtlasC': "cyan", 'AtlasO': "orange",
         'UBVRI': 'chocolate', 'GaiaG': 'g',
         'L_ubvri': 'sandybrown', 'L_bol': 'saddlebrown', 'XEUV<325': 'sienna', 'IR>890': 'peru',
         Band.NameBol: 'black', Band.NameUBVRI: 'dimgrey', Band.NameBolQuasi: 'dimgray'}
    # for Subaru HCS: colors
    for b in list('grizY'):
        c['HSC' + b] = c[b]

    if bname is None:
        return c
    else:
        if bname in c:
            return c[bname]
        return default


def lntypes(bname=None, default='-'):
    ln = {"U": "-", 'B': "-", 'V': "-", 'R': "-", 'I': "-", 'UVM2': "-.", 'UVW1': "-.", 'UVW2': "-.", 'F125W': ":",
          'F160W': "-.", 'F140W': "--", 'F105W': "-.", 'F435W': "--", 'F606W': "-.", 'F814W': "--", 'u': "--",
          'g': "--", 'r': "--", 'i': "--", 'z': "--", 'Y': '--', 'GaiaG': '--',
          Band.NameBol: '-', Band.NameUBVRI: '-.', Band.NameBolQuasi: ':'}
    # for Subaru HCS: colors
    for b in list('grizY'):
        ln['HSC' + b] = ln[b]

    if bname is None:
        return ln
    else:
        return ln[bname]
        # if bname in lntypes:
        #     return lntypes[bname]
        # return default


def band_load_names(path=Band.DirRoot):
    """Find directories with filter.dat """
    from configparser import ConfigParser

    def ini_to_bands(fini):
        parser = ConfigParser()
        parser.optionxform = str
        parser.read(fini)
        res = {}
        for bname in parser.sections():
            fn = os.path.join(dirname(fini), parser.get(bname, 'file'))
            zp = None
            if parser.has_option(bname, Band.NameZp):
                zp = float(parser.get(bname, Band.NameZp))
            jy = None
            if parser.has_option(bname, Band.NameJy):
                jy = float(parser.get(bname, Band.NameJy))

            b = Band(name=bname, fname=fn, zp=zp, jy=jy)
            res[bname] = b
        return res

    bol = BandUni.get_bol()
    qbol = BandUni.get_qbol()
    ubvri = BandJoin.get_ubvri()
    bands = {Band.NameBol: bol, Band.NameUBVRI: ubvri, Band.NameBolQuasi: qbol}
    for d, dir_names, file_names in os.walk(path):
        if Band.FileFilters in file_names:  # check that there are info-file in this directory
            fname = os.path.join(d, Band.FileFilters)
            bs = ini_to_bands(fname)
            if len(bs) > 0:
                bands.update(bs)

    return bands


def get_names():
    if len(Band.Cache) == 0:
        Band.Cache = band_load_names()

    return list(Band.Cache.keys())


def band_get_aliases():
    return Band.Alias


def is_exist_names(name):
    return name in get_names()


def is_exist_alias(name):
    if band_get_aliases() is None:
        return False
    return name in band_get_aliases().keys()


def is_exist(name):
    return is_exist_names(name) or is_exist_alias(name)


def band_get_names_alias():
    res = get_names()
    res.extend(band_get_aliases().keys())
    return res


def band_by_name(name):
    """
        Get rf by name, for example "U"
    :param name: for example "U" or "B"
    :return: class Band with waves and responses
    """
    if name in Band.Cache:
        b = Band.Cache[name]
        if not b.is_load:
            b.load()
        return Band.Cache[name]

    if Band.Alias is None:
        Band.load_settings()

    if is_exist_alias(name):
        for aname, oname in Band.Alias.items():
            if is_exist(oname):
                bo = Band.Cache[oname]
                ba = bo.clone(aname)
                if not ba.is_load:
                    ba.load()
                Band.Cache[aname] = ba
            else:
                raise ValueError("Errors in your aliases: for alias [%s] no such band [%s]. See settings in  %s"
                                 % (aname, oname, Band.FileSettings))

    if is_exist(name):
        # if name == Band.BolName:  # bolometric
        #     # todo check this filter: values less then g-u etc.
        #     b = BandUni(name=Band.BolName, wlrange=(1e0, 42e3), length=300)
        # else:
        b = Band.Cache[name]
        if not b.is_load:
            b.load()
        # Band.Cache[name] = b
        return b
    else:
        print("  Error: no band: %s " % name)
        return None


def print_bands(ncol=5):
    bands = sorted(get_names())
    print("Available bands:")
    # print "   Available bands: \n   %s" % '-'.join(sorted(bands))
    if bands is not None and len(bands) > 0:
        s = '  '
        c = bands[0][0]
        for b in bands:
            if b[0] != c:
                print(s)
                c = b[0]
                s = "   %-7s " % b
            else:
                s += " %-7s " % b
        if s != '':
            print(s)
    else:
        print("You have not load them yet. Try: band.Band.load_settings() ")

    alias = band_get_aliases()
    print("\nAvailable aliases of bands: ")
    if alias is not None and len(alias) > 0:
        col = 0
        s = ''
        for k, v in alias.items():
            col += 1
            s += " %3s => %-7s " % (k, v)
            if col == ncol:
                print(s)
                col = 0
                s = ''
        if s != '':
            print(s)
    else:
        print("     No aliases or you have not load them yet. Try: band.Band.load_settings() ")
