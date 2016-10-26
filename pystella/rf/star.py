from scipy import interpolate
import numpy as np
# from scipy.integrate import simps as integralfunc
from scipy import integrate

from pystella.rf import band
from pystella.rf.rad_func import Flux2MagAB
from pystella.util.phys_var import phys

__author__ = 'bakl'


class Star:
    def __init__(self, name, spec=None, is_flux_eq_luminosity=False):
        """Creates a Star with Spectrum instance.  Required parameters:  name."""
        self.name = name
        self._sp = spec
        self.is_flux_eq_luminosity = is_flux_eq_luminosity
        self.radius_ph = None
        self._z = None
        self._magnification = 1.
        self.distance = None
        self.Tcol = {}
        self.zeta = {}

    def set_radius_ph(self, radius):
        self.radius_ph = radius

    def set_distance(self, distance):
        """
        Set distance to the star [cm]
        :param distance:
        """
        self.distance = distance

    def set_redshift(self, z):  # shift spectrum to rest frame
        self._z = z

    def set_magnification(self, m):  # shift spectrum to rest frame
        self._magnification = m

    def set_Tcol(self, Tcol, bset):
        self.Tcol[bset] = Tcol

    def get_Tcol(self, bset):
        if bset in self.Tcol:
            return self.Tcol[bset]
        return None

    def set_zeta(self, zeta, bset):
        self.zeta[bset] = zeta

    def get_zeta(self, bset):
        if bset in self.zeta:
            return self.zeta[bset]
        return None

    @property
    def z(self):
        return self._z

    @property
    def IsRedshift(self):
        return self.z is not None and self.z > 0.

    @property
    def IsRadius(self):
        return self.radius_ph is not None

    @property
    def IsDistance(self):
        return self.distance is not None

    @property
    def IsRadiusDist(self):
        return self.IsRadius and self.IsDistance

    @property
    def Freq(self):
        if self._sp is None:
            raise ValueError("Spectrum has not been defined. ")
        if self.IsRedshift:
            return self._sp.Freq / (1. + self.z)  # redshift the flux
        else:
            return self._sp.Freq

    @property
    def Wl(self):
        if self._sp is None:
            raise ValueError("Spectrum has not been defined. ")
        return phys.c / self.Freq

    @property
    def Flux(self):
        if self._sp is None:
            raise ValueError("Spectrum has not been defined. ")
        flux = self._sp.Flux * self._magnification
        if self.IsRedshift:
            return Star.flux_to_redshift(flux, self.z)
        else:
            return flux

    @property
    def Flux_wl(self):
        return self.Flux * self.Freq ** 2 / phys.c  # flux [erg/cm^2/cm) ]

    @property
    def Luminosity(self):
        if self.is_flux_eq_luminosity:
            return self.Flux

        if self.radius_ph is None:
            raise ValueError("Photospheric radius has not been defined. ")

        return 4. * np.pi * self.radius_ph ** 2 * self.Flux

    @property
    def FluxObs(self):
        if self.IsRadiusDist:
            return self.Luminosity / (4 * np.pi * self.distance ** 2)
        elif self.IsDistance:
            return self.Flux / (4 * np.pi * self.distance ** 2)
            # return self.Flux / (4 * np.pi * self.distance ** 2)
        else:
            return self.Flux

    @property
    def FluxWlObs(self):
        if self.IsRadiusDist:
            if self.is_flux_eq_luminosity:
                return self.Flux_wl / (4 * np.pi * self.distance ** 2)
            else:
                return self.Flux_wl * (self.radius_ph / self.distance) ** 2
        elif self.IsDistance:
            return self.Flux_wl / (4 * np.pi * self.distance ** 2)
        else:
            return self.Flux_wl

    @property
    def FluxAB(self):
        return -2.5 * np.log10(self.FluxObs) + phys.ZP_AB

    def _response_lmb(self, band, is_b_spline=True):
        """
        Compute response flux using provided spectral band
        :param band:  photometric band
        :param is_b_spline:  the method of interpolation
        :return: :raise ValueError:
        """
        wl = self.Wl
        if min(wl) > band.wl[0] or max(wl) < band.wl[-1]:
            raise ValueError("Spectrum must be wider then band: " + str(band))

        flux = self.FluxWlObs / phys.cm_to_angs  # to flux [erg/cm^2/A) ]
        wl_s = wl * phys.cm_to_angs
        wl_b = band.wl * phys.cm_to_angs

        if is_b_spline:
            tck = interpolate.splrep(wl_s, flux, s=0)
            flux_spline = interpolate.splev(wl_b, tck, der=0)
        else:
            flux_spline = np.interp(wl_b, wl_s, flux, 0, 0)  # One-dimensional linear interpolation.

        a = integrate.integralfunc(flux_spline * band.resp_wl * wl_b, wl_b) / (phys.c * phys.cm_to_angs) / phys.h
        return a

    def _response_nu(self, b, is_b_spline=True):
        """
        Compute response flux using provided spectral band.
        see: http://www.astro.ljmu.ac.uk/~ikb/research/mags-fluxes/
        :param b:  photometric band
        :param is_b_spline:  the method of interpolation
        :return:
        """
        nu_s = self.Freq
        flux = self.FluxObs  # / phys.cm_to_angs  # to flux [erg/cm^2/A) ]

        # sort
        sorti = np.argsort(nu_s)
        nu_s = nu_s[sorti]
        flux = flux[sorti]

        nu_b = np.array(b.freq)
        resp_b = np.array(b.resp_fr)

        if min(nu_s) > nu_b[0] or max(nu_s) < nu_b[-1]:
            # decrease wave length range of the band
            f = (min(nu_s) < nu_b) & (nu_b < max(nu_s))
            # f = map(lambda x: min(nu_s) < x < max(nu_s), nu_b)
            nu_b = nu_b[f]
            resp_b = resp_b[f]
            # raise ValueError("Spectrum must be wider then band: " + str(band))

        if is_b_spline:
            tck = interpolate.splrep(nu_s, flux, s=0)
            flux_spline = interpolate.splev(nu_b, tck, der=0)
            if np.any(flux_spline < 0.):
                flux_spline = np.interp(nu_b, nu_s, flux, 0, 0)  # One-dimensional linear interpolation.
        else:
            flux_spline = np.interp(nu_b, nu_s, flux, 0, 0)  # One-dimensional linear interpolation.

        a = integrate.simps(flux_spline * resp_b / nu_b, nu_b)
        return a

    def flux_to_mag_not_checked(self, band):
        conv = self._response_lmb(band)
        # conv = self._response_nu(band)
        if conv <= 0:
            return None
        else:
            mag = -2.5 * np.log10(conv) + band.zp
            return mag

    def flux_to_magAB(self, b):
        response = self._response_nu(b)
        if response <= 0:
            raise ValueError("Spectrum should be more 0: %f" % response)
            # return None
        else:
            # mag = -2.5 * np.log10(conv) + phys.ZP_AB - band.zp
            mag = Flux2MagAB(response / b.Norm) - b.zp
            return mag

    def k_cor(self, band_r, band_o, z=0.):
        """
        Compute K-correction for observed and rest-frame bands.

        Args:
          band_r: Rest-frame band.
          band_o: Observed band.
          z:     redshift

        Returns:
          * K: K-correction
          * If failed return None
        """
        # todo make k-correction with b-splinesec

        if z > 0:
            self.set_redshift(z)

        z_o = z
        if self.IsRedshift:
            z_o = self.z

        self.set_redshift(0.)
        resp_0 = self._response_lmb(band_r, is_b_spline=False)

        self.set_redshift(z_o)
        resp_z = self._response_lmb(band_o, is_b_spline=False)

        if resp_0 < 0 or resp_z <= 0:
            return None
        else:
            kcor = -2.5 * np.log10(resp_z / resp_0 / (1 + z_o)) + band_r.zp - band_o.zp
            return kcor

    @staticmethod
    def flux_to_redshift(flux, z):
        if z <= 0.:
            return flux
        flux_z = flux * (1.+z)
        # flux_z = flux / (1.+z)
        # flux_z = flux
        # flux_z = np.interp(freq / (1. + z), freq[::-1], flux[::-1])
        return flux_z
