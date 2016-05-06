import numpy as np

from pystella.rf import band
from pystella.rf.lc import LightCurve
from pystella.rf.rad_func import Lum2MagBol
from pystella.util.phys_var import phys
import matplotlib.pyplot as plt


class Popov:
    def __init__(self, name, R, M, Mni, E):
        """Creates a Popov's light curves model.
        Required parameters:  name, radius in R_sun, mass in M_sun."""
        self.name = name
        self._kappa = 0.4  # electron scattering for pure H, [cm^2/g]
        # self._kappa = 0.34  # Thompson opacity for solar composition, [cm^2/g]
        self._T_ion = 5000  # H recombination temperature [K]
        # self._v_sc = None
        self._T0 = None
        self._rho0 = None

        self._R0 = R * phys.R_sun

        self._Etot = E
        self._Mtot = M * phys.M_sun
        self._Mni = Mni * phys.M_sun

    @property
    def R0(self):
        return self._R0

    @property
    def Etot(self):
        return self._Etot

    @property
    def Eth0(self):
        """The total thermal energy of the envelope just after
        propagating of the shock wave"""
        return self._Etot / 2.

    @property
    def Mtot(self):
        return self._Mtot

    @property
    def Mni(self):
        return self._Mni

    @property
    def Etot_ni(self):
        """The total energy for Ni decay"""
        e = 6.22e49 * self._Mni / phys.M_sun
        return e

    @property
    def Etot_co(self):
        """The total energy for Co decay"""
        e = 1.26e50 * self._Mni / phys.M_sun
        return e

    @property
    def Etot_decay(self):
        """The total energy from decays Ni and  Co to Fe"""
        e = self.Etot_ni + self.Etot_co
        return e

    @property
    def kappa(self):
        return self._kappa

    @property
    def Tion(self):
        return self._T_ion

    @property
    def Lambda(self):
        numerator = 54.*np.sqrt(30.)*phys.sigma_SB*self.kappa**2*self.Mtot**1.5*self.Tion**4
        denominator = np.pi**5*phys.c**2*np.sqrt(self.Etot)*self.R0
        return numerator/denominator

    @property
    def v_sc(self):
        return np.sqrt(10. / 3. * self.Etot / self.Mtot)

    @property
    def t_d(self):
        """Diffusion time"""
        # t1 = 9.*self._kappa*self.Mtot / (4. * np.pi ** 3 * phys.c * self.R0)
        pho0 = self.Mtot / (4./3. * np.pi * self.R0**3)
        t = 3.*self._kappa*pho0*self.R0**2 / phys.c
        return t

    @property
    def t_e(self):
        """Expansion time"""
        t = self.R0 / self.v_sc
        return t

    @property
    def t_a(self):
        """Characteristic time of changing SN magnitude in models without WCR (see Arnett,1980)"""
        t = np.sqrt(2.*self.t_d*self.t_e)
        return t

    @property
    def t_i(self):
        """At t < t_i the surface temperature of the envelope is higher than Tion"""
        t = self.t_a / np.sqrt(self.Lambda)
        return t

    @property
    def t_max(self):
        """The time of the maximum of the bolometric luminocity"""
        t = 0.75 * self.t_i * self.t_a ** 2 + 0.25 * self.t_i ** 3
        t **= 1. / 3.
        return t

    @property
    def t_p(self):
        """The time of the duration of the plateau"""
        t = 4. ** (1. / 3.) * self.t_max
        return t

    def R(self, time):
        t = time * phys.d2s
        return self._R0 + self.v_sc * t

    def v(self, x):
        return self.v_sc * x

    def rho(self, time):
        return self._rho0 * (self._R0/self.R(time))**3

    def Rph(self, time):
        t = time * phys.d2s
        tmp = 1./3.*self.v_sc**2 * (self.t_i*t*(3.+(self.t_i/self.t_a)**2) - (t**2/self.t_a)**2)
        return np.sqrt(tmp)

    def Lbol(self, time):
        L = np.array(map(self.lum_bol, time))
        return L

    def lum_bol(self, time, is_ni=False):
        """Return bolometric luminosity.
        Parameters: t [day]"""
        t = time * phys.d2s

        if t < self.t_i:
            lum = self.Eth0 / self.t_d * np.exp(-(t / self.t_a) ** 2)
            return lum

        lumWCR = 8. * np.pi * phys.sigma_SB * self.Tion ** 4 * self.Rph(time) ** 2
        if t < self.t_max:
            return lumWCR

        lumNiCo = self.e_rate_ni(time)
        if t > self.t_p:
            return lumNiCo

        if is_ni:
            # lumNiCo = self.e_rate_ni(time)
            lum = np.max((lumWCR, lumNiCo))
        else:
            lum = lumWCR
        return lum

    #  Ni & Co
    def e_rate_ni(self, time):
        """The total rate energy production [erg/s], see Nadyozhin D.K., 1994
         http://adsabs.harvard.edu/abs/1994ApJS...92..527N
         Parameters: time [day]"""
        e = (6.45e43*np.exp(-time/8.8) + 1.45e43*np.exp(-time/111.3)) * self.Mni / phys.M_sun  #
        return e

    def MagBol(self, time):
        mag = Lum2MagBol(self.Lbol(time))
        return mag

    def LCBol(self, time):  # todo check this code
        b = band.BandUni()
        mags = self.MagBol(time)
        lc = LightCurve(b, time, mags)
        return lc

    def plot_Lbol(self, t, is_plot_Lum=False):
        """t [day]"""
        # lum
        if is_plot_Lum:
            mags = self.Lbol(t)
            mags_ni = self.e_rate_ni(t)
        else:  # mag
            mags = self.MagBol(t)
            mags_ni = Lum2MagBol(self.e_rate_ni(t))

        fig = plt.figure()
        ax = fig.add_axes((0.1, 0.3, 0.8, 0.65))

        ax.plot(t, mags, color='blue', label='L bol', lw=2.5)
        ax.plot(t, mags_ni, color='red', ls='-.', label='Decay of Ni56 & Co', lw=2.)

        times = {'Diffusion time': self.t_d,
                 'Expansion time': self.t_e,
                 'T surf > Tion': self.t_i,
                 'Max bol time': self.t_max,
                 'Plateau duration time': self.t_p}

        for k, v in times.items():
            if is_plot_Lum:
                xy = (v / phys.d2s, self.lum_bol(v))
            else:
                xy = (v / phys.d2s, Lum2MagBol(self.lum_bol(v)))
            ax.annotate(k, xy=xy, xytext=(-10, 20), ha='left',
                        textcoords='offset points',
                        arrowprops=dict(arrowstyle='->', shrinkA=0))

        if not is_plot_Lum:
            ax.invert_yaxis()
        ax.set_xlabel('Time [day]')
        ax.set_ylabel('Absolute Magnitude')

        ax.legend()
        txt = '{0:10} {1:.4e} R_sun\n'.format('R:', self.R0 / phys.R_sun) + \
              '{0:10} {1:.4e} M_sun\n'.format('Mtot:', self.Mtot / phys.M_sun) + \
              '{0:10} {1:.4e} M_sun\n'.format('Mni:', self.Mni / phys.M_sun) + \
              '{0:10} {1} ergs'.format('Etot:', self.Etot)
        fig.text(0.17, 0.07, txt, family='monospace')
        # plt.title('R=%5.1f Mtot=%4.2f Mni=%f Etot=%e ' % (self.R0, self.Mtot, self.Mni, self.Etot))
        return ax
