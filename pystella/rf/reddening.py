import numpy as np

__author__ = 'bakl'


class ReddeningLaw(object):
    MW = "MW"
    LMC = "LMC"
    SMC = "SMC"

    Names = {MW: "Milky Way", LMC: "Large Magellanic Cloud", SMC: "Small Magellanic Cloud"}
    Laws = list(Names.keys())

    @staticmethod
    def Rlmd(wl, law=MW):
        if law not in ReddeningLaw.Laws:
            raise ValueError('Such law [{}] is not supported. There are {}'.
                             format(law, ', '.join(ReddeningLaw.Laws)))
        return LawPei.Rlmd(wl, law)

    @staticmethod
    def ksi(wl, law=MW):
        """

        :param wl: [A
        :param law:
        :return:
        """
        if law not in ReddeningLaw.Laws:
            raise ValueError('Such law [{}] is not supported. There are {}'.
                             format(law, ', '.join(ReddeningLaw.Laws)))
        wl = np.array(wl) / 1e4  # A -> microns
        return LawPei.ksi(wl, law)


class LawPei(ReddeningLaw):
    """Interstellar dust model
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
    def ksi(wl, law=MW):
        res = [LawPei.func_ksi(l, LawPei.ai[law], LawPei.li[law], LawPei.bi[law], LawPei.ni[law], )
               for l in wl]
        # for l in wl:
        #     res += ai[i] / (bi[i] + (wl/li)**ni[i] + (li/wl)**ni[i])
        return np.array(res)

    @staticmethod
    def func_ksi(wl, ai, li, bi, ni):
        res = 0.
        for i in range(6):
            res += ai[i] / (bi[i] + (wl/li[i]) ** ni[i] + (li[i]/wl) ** ni[i])
        return res

    @staticmethod
    def Rlmd(wl, law=MW):
        # todo something...
        pass
