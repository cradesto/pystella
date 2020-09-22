import re
import unittest
from os.path import dirname, abspath, join

from pystella.model.stella import Stella

__author__ = 'bakl'


class TestStellaRes(unittest.TestCase):
    def setUp(self):
        name = 'cat_R1000_M15_Ni007_E15'
        path = join(dirname(abspath(__file__)), 'data', 'stella')
        stella = Stella(name, path=path)
        self.res = stella.get_res()

    def test_info_parse(self):
        res_header = """
                              <===== HYDRODYNAMIC RUN OF MODEL cat_R1000_M15_Ni007_E15.prf                                                     =====>

                              MASS(SOLAR)=15.000       RADIUS(SOLAR)= 1000.000
                              EXPLOSION ENERGY(10**50 ERG)= 1.50000E+01

                              <=====                                 =====>


  INPUT PARAMETERS
 EPS   =         0.00300          Rce   =     1.00000E-02 SOL.Rad.
 HMIN  =     1.00000E-11 S        AMht  =     1.00000E-02 SOL.MASS
 HMAX  =     5.00000E+04 S        Tburst=     1.00000E-01 S
 THRMAT=     1.00000E-30          Ebstht=     1.50000E+01 1e50 ergs
 METH  =               3          CONV  =               F
 JSTART=               0          EDTM  =               T
 MAXORD=               4          CHNCND=               T
 NSTA  =           -4500          FitTau=     3.00000E+02
 NSTB  =              -1          TauTol=     1.30000E+00
 NOUT  =             100          IOUT  =              -1
 TcurA =         0.00000          Rvis   =        1.00000
 TcurB =       200.00000          BQ    =     1.00000E+00
 NSTMAX=          360000          DRT   =     1.00000E-01
 XMNI  =     7.00000E-02 SOL.MASS NRT   =               1
 XNifor=     1.16561E-01
 MNicor=     1.16999E-01 SOL.MASS SCAT  =               T

               """
        pattern = filter(lambda x: len(x) > 0, res_header.splitlines())
        pattern = map(str.strip, pattern)
        p = r"(.*?)\s*=\s+([-+]?\d*\.\d+|\d+)"
        for line in pattern:
            res = re.findall(p, line)
            if len(res) > 0:
                for k, v in res:
                    print("key: %s  v: %f " % (k, float(v)))

    def test_info_parse(self):
        info = self.res.Info.parse()

        tmp = 1000.
        self.assertEquals(info.R, tmp, "Radius [%f] should be %f" % (info.R, tmp))

        tmp = 15.
        self.assertEquals(info.M, tmp, "Mass [%f] should be %f" % (info.M, tmp))

        tmp = 15.
        self.assertEquals(info.E, tmp, "Ebstht [%f] should be %f" % (info.E, tmp))

    def test_res_times(self):
        res = []
        for i, a in enumerate(self.res.blocks()):
            res.append(a)
            print(i, a)
        self.assertEqual(len(res), 141)

    def test_res_find_block(self):
        time = 99.
        i, (t, start, end) = self.res.find_block(time=time)
        print("{}: t= {} start= {} end= {}".format(i, t, start, end))
        self.assertGreater(t, time)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
