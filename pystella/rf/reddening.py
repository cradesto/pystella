import numpy as np

__author__ = 'bakl'


class ReddeningLaw(object):
    MW = "MW"
    LMC = "LMC"
    SMC = "SMC"

    Names = {MW: "Milky Way", LMC: "Large Magellanic Cloud", SMC: "Small Magellanic Cloud"}
    Laws = list(Names.keys())

    @staticmethod
    def Almd(wl, ebv, law=MW):
        """
        Compute Absorption A_lamda
        :param wl: [A]
        :param ebv:  E(B-V)
        :param law:
        :return:
        """
        if law not in ReddeningLaw.Laws:
            raise ValueError('Such law [{}] is not supported. There are {}'.
                             format(law, ', '.join(ReddeningLaw.Laws)))
        wl = np.array(wl) / 1e4  # A -> microns
        return LawPei.Almd(wl, ebv, law)

    @staticmethod
    def xi(wl, law=MW):
        """
        Compute extinction curve
        :param wl: [A]
        :param law:
        :return:
        """
        if law not in ReddeningLaw.Laws:
            raise ValueError('Such law [{}] is not supported. There are {}'.
                             format(law, ', '.join(ReddeningLaw.Laws)))
        wl = np.array(wl) / 1e4  # A -> microns
        return LawPei.xi(wl, law)


class LawPei(ReddeningLaw):
    """Interstellar dust model by Pei (1992)
    See http://adsabs.harvard.edu/abs/1992ApJ...395..130P
    """
    MW = "MW"
    LMC = "LMC"
    SMC = "SMC"

    Rv = {MW: 3.08, LMC: 3.16, SMC: 2.93}

    ai = {MW: (165., 14., 0.045, 0.002, 0.002, 0.012),
          LMC: (175., 19., 0.023, 0.005, 0.006, 0.020),
          SMC: (185., 27., 0.005, 0.010, 0.012, 0.030)
          }

    li = {MW: (0.047, 0.08, 0.22, 9.7, 18.0, 25.0),
          LMC: (0.046, 0.08, 0.22, 9.7, 18.0, 25.0),
          SMC: (0.042, 0.08, 0.22, 9.7, 18.0, 25.0)
          }

    bi = {MW: (90., 4.0, -1.95, -1.95, -1.80, 0.00),
          LMC: (90., 5.5, -1.95, -1.95, -1.80, 0.00),
          SMC: (90., 5.5, -1.95, -1.95, -1.80, 0.00)
          }

    ni = {MW: (2, 6.5, 2, 2, 2, 2),
          LMC: (2, 4.5, 2, 2, 2, 2),
          SMC: (2, 4.0, 2, 2, 2, 2)
          }

    @staticmethod
    def xi(wl, law=MW):
        res = [LawPei.func_xi(l, LawPei.ai[law], LawPei.li[law], LawPei.bi[law], LawPei.ni[law], )
               for l in wl]
        # for l in wl:
        #     res += ai[i] / (bi[i] + (wl/li)**ni[i] + (li/wl)**ni[i])
        return np.array(res)

    @staticmethod
    def func_xi(wl, ai, li, bi, ni):
        res = 0.
        for i in range(6):
            res += ai[i] / (bi[i] + (wl/li[i]) ** ni[i] + (li[i]/wl) ** ni[i])
        return res

    @staticmethod
    def Almd(wl, ebv, law=MW):
        Av = LawPei.Rv[law] * ebv
        xi = LawPei.xi(wl, law)
        return xi * Av


class LawFitz(ReddeningLaw):
    """Interstellar dust model by Fitzpatrick.
    See  Fitzpatrick (1999) http://iopscience.iop.org/article/10.1086/316293
    https://github.com/sczesla/PyAstronomy/blob/master/src/pyasl/asl/unred.py
    https://github.com/Morisset/PyNeb_devel/blob/master/pyneb/extinction/red_corr.py
    https://github.com/LSSTDESC/Recipes/blob/master/python/utensils/extinction.py
    """

    # see Gordon (2003) http://arxiv.org/abs/astro-ph/0305257
    Rv = {ReddeningLaw.MW: 3.1, ReddeningLaw.LMC: 3.41, ReddeningLaw.SMC: 2.74}
    Names = {k: "{}: Rv={}".format(v, Rv[k]) for k, v in ReddeningLaw.Names}
    # Names = {ReddeningLaw.MW: "Rv:3.1", LMC: "Large Magellanic Cloud", SMC: "Small Magellanic Cloud"}

    @staticmethod
    def Almd(wl, ebv, law=ReddeningLaw.MW, model='f99'):
        """
        Origin: https://github.com/LSSTDESC/Recipes/blob/master/python/utensils/extinction.py
        :param wl: array of wave lengths [A]
        :param ebv:  E(B-V)
        :param law: R_V
        :param model: f99 or fm07. Default: f99
        :return: A_lambda
        """
        from scipy.interpolate import interp1d

        x = 1e4 / wl
        # r_v parameter. Default: the standard Milky Way average of 3.1.
        r_v = LawFitz.Rv[law]

        if np.any(x < 0.167) or np.any(x > 11.):
            raise ValueError('Wavelength(s) must be between 910 A and 6 um')
        if model == 'fm07' and abs(r_v - 3.1) > 0.001:
            raise ValueError('fm07 model not implemented for r_v != 3.1')

        k = np.zeros_like(x)
        uv_region = (x >= 1.e4 / 2700.)
        oir_region = ~uv_region

        # UV region
        y = x[uv_region]
        if model == 'f99':
            x0, gamma = 4.596, 0.99
            c3, c4, c5 = 3.23, 0.41, 5.9
            c2 = -0.824 + 4.717 / r_v
            c1 = 2.030 - 3.007 * c2
            d = y ** 2 / ((y ** 2 - x0 ** 2) ** 2 + y ** 2 * gamma ** 2)
            f = np.zeros_like(y)
            valid = (y >= c5)
            f[valid] = 0.5392 * (y[valid] - c5) ** 2 + 0.05644 * (y[valid] - c5) ** 3
            k_uv = c1 + c2 * y + c3 * d + c4 * f
        if model == 'fm07':
            x0, gamma = 4.592, 0.922
            c1, c2, c3, c4, c5 = -0.175, 0.807, 2.991, 0.319, 6.097
            D = y ** 2 / ((y ** 2 - x0 ** 2) ** 2 + y ** 2 * gamma ** 2)
            k_uv = np.zeros_like(y)
            valid = (y <= c5)
            k_uv[valid] = c1 + c2 * y[valid] + c3 * D[valid]
            valid = (y > c5)
            k_uv[valid] = c1 + c2 * y[valid] + c3 * D[valid] + c4 * (y[valid] - c5) ** 2
        k[uv_region] = k_uv

        # Calculate values for UV spline points to anchor OIR fit
        x_uv_spline = 1.e4 / np.array([2700., 2600.])
        d = (x_uv_spline ** 2 /
             ((x_uv_spline ** 2 - x0 ** 2) ** 2 + x_uv_spline ** 2 * gamma ** 2))
        k_uv_spline = c1 + c2 * x_uv_spline + c3 * d

        # Optical / IR region
        y = x[oir_region]
        if model == 'f99':
            anchors_x = 1.e4 / np.array([np.inf, 26500., 12200., 6000., 5470.,
                                         4670., 4110.])

            # The OIR anchors are from IDL astrolib, not F99.
            anchors_extinction = np.array(
                [0.,
                 0.26469 * r_v / 3.1,  # IR
                 0.82925 * r_v / 3.1,  # IR
                 -0.422809 + 1.00270 * r_v + 2.13572e-04 * r_v ** 2,  # optical
                 -5.13540e-02 + 1.00216 * r_v - 7.35778e-05 * r_v ** 2,
                 0.700127 + 1.00184 * r_v - 3.32598e-05 * r_v ** 2,
                 (1.19456 + 1.01707 * r_v - 5.46959e-03 * r_v ** 2 +
                  7.97809e-04 * r_v ** 3 - 4.45636e-05 * r_v ** 4)]
            )

            anchors_x = np.append(anchors_x, x_uv_spline)
            anchors_k = np.append(anchors_extinction - r_v, k_uv_spline)

        if model == 'fm07':
            anchors_x_ir = np.array([0., 0.25, 0.50, 0.75, 1.])
            anchors_k_ir = (-0.83 + 0.63 * r_v) * anchors_x_ir ** 1.84 - r_v
            anchors_x_opt = np.array([5530., 4000., 3300.])
            anchors_k_opt = np.array([0., 1.322, 2.055])

            anchors_x = np.append(anchors_x_ir, anchors_x_opt)
            anchors_k = np.append(anchors_k_ir, anchors_k_opt)

            anchors_x = np.append(anchors_x, x_uv_spline)
            anchors_k = np.append(anchors_k, k_uv_spline)

        # Note that interp1d requires that the input abscissa is monotonically
        # _increasing_. This is opposite the usual ordering of a spectrum, but
        # fortunately the _output_ abscissa does not have the same requirement.
        oir_spline = interp1d(anchors_x, anchors_k, kind='cubic')
        k[oir_region] = oir_spline(y)

        return (k + r_v) * ebv
