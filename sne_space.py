#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import logging
import json
import os
import sys

import pandas
# from pandas.io import json
import matplotlib.pyplot as plt


from pystella.rf import band
from pystella.rf import light_curve_plot as lcp
from pystella.rf.lc import SetLightCurve, LightCurve
from pystella.rf.spectrum import SeriesSpectrum, Spectrum
from pystella.util.phys_var import phys

__author__ = 'bakl'

logging.basicConfig()

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class SneSpace:
    def __init__(self):
        """Creates a SneSpace model instance.  Required parameters:  name."""
        self._name = None
        self._fname = None
        self._data = None
        self._serial_spec = None

    @property
    def data(self):
        return self._data

    @property
    def Name(self):
        if self._name is not None:
            return self._name
        if self.data is None:
            return None
        self._name = next(iter(self.data.keys()))
        return self._name

    @property
    def photometry(self):
        return self.property('photometry')

    @property
    def spectra(self):
        return self.property('spectra')

    @property
    def lumdist(self):
        l = self.property('lumdist')
        if l is not None:
            d = float(l[0]['value'])
            u = l[0]['u_value'].lower()
            if u == 'mpc':
                d *= 1.e6
            return d
        return None

    @property
    def comovingdist(self):
        l = self.property('comovingdist')
        if l is not None:
            dd = 0.
            for node in l:
                d = float(node['value'])
                u = node['u_value'].lower()
                if u == 'mpc':
                    d *= 1.e6
                dd += d
            return dd / len(l)
        return None

    def property(self, name):
        if self.data is None:
            return None
        return self._data[self.Name][name]

    def __str__(self):
        return "%s, path: %s" % (self.Name, self._fname)

    def __repr__(self):
        return "%s, path: %s" % (self.Name, self._fname)

    def load(self, fname):
        try:
            with open(fname, 'r') as f:
                self._data = json.loads(f.read())
        except ValueError as e:
            logger.error("Could not parse {0} due to {1}.".format(fname, e))
            return False
        else:
            self._fname = fname
        return True

    def to_curves(self):
        ph = pandas.DataFrame.from_dict(self.photometry)

        def add(d, k, v):
            if k in d:
                d[k].append(v)
            else:
                d[k] = [v]

        times = {}
        mags = {}
        err = {}
        for b, t, m, e in zip(ph.band, ph.time, ph.magnitude, ph.e_magnitude):
            add(times, b, t)
            add(mags, b, m)
            add(err, b, e)

        bands = list(times.keys())

        band.Band.load_settings()
        curves = SetLightCurve(self.Name)
        for bname in bands:
            if band.band_is_exist(bname):
                tt = list(map(float, times[bname]))
                mm = list(map(float, mags[bname]))
                ee = list(map(float, err[bname]))
                lc = LightCurve(bname, tt, mm, ee)
                curves.add(lc)

        return curves

    def to_spectra(self):
        series = SeriesSpectrum(self.Name)

        def to_pystella_spec(name, tbl):
            lmd = []
            flux = []
            for x in tbl:
                l, f = x
                lmd.append(float(l)*phys.angs_to_cm)
                flux.append(float(f))
            if len(lmd) > 0:
                sp = Spectrum.from_lambda(name, lmd, flux)
                return sp
            return None

        for i, item in enumerate(self.spectra):
            tbl = item
            # tbl = pandas.DataFrame.from_dict(item)
            if 'time' in tbl:
                t = float(tbl['time'])
                name = tbl['filename']
                spec = to_pystella_spec(name, tbl['data'])
                if spec is not None:
                    series.add(t, spec)

        return series


def get_parser():
    parser = argparse.ArgumentParser(description='Process Stella Shock Wave Details.')
    parser.add_argument('-i', '--input',
                        required=False,
                        dest="fname",
                        help="File name with json-data")
    parser.add_argument('--lumnorm',
                        required=False,
                        type=float,
                        default=1e40,
                        dest="lumnorm",
                        help="Luminously normalization, example: 1e40")
    parser.add_argument('-s', '--save',
                        action='store_const',
                        const=True,
                        dest="is_save",
                        help="To save plot to pdf-file. Default: False".format('2:10:50'))
    parser.add_argument('-w', '--write',
                        action='store_const',
                        const=True,
                        dest="is_write",
                        help="To write the data to txt-file.")
    return parser


def main():
    # fname = "~/Sn/Release/svn_kepler/stella/branches/lucy/run/res/sncurve/sn1999em/SN1999em.json"
    parser = get_parser()
    args, unknownargs = parser.parse_known_args()

    if len(unknownargs) > 0:
        fname = unknownargs[0]
    else:
        if args.fname and os.path.isfile(args.fname):
            fname = args.fname
        else:
            # logger.error(" No data. Use key '-i' ")
            parser.print_help()
            sys.exit(2)
    fname = os.path.expanduser(fname)

    sn = SneSpace()
    print("Load {0}".format(fname))
    sn.load(fname)
    print("{} is loaded.".format(sn.Name))

    curves = sn.to_curves()
    print("Curves {} is computed. There are {} bands.".format(sn.Name, '-'.join(curves.BandNames)))
    serial_spec = sn.to_spectra()
    print("Spectra {} is computed: {} times points.".format(sn.Name, serial_spec.Length))
    # curves_spectra = serial_spec.flux_to_curves(['V'], d=sn.lumdist*phys.pc)
    mu = 1.
    # mu = (10./sn.lumdist)**2
    curves_spectra = serial_spec.flux_to_curves(curves.BandNames, d=None, magnification=mu)

    # plotting
    ax = lcp.curves_plot(curves, title='Photometry', is_line=False)
    lcp.curves_plot(curves_spectra, title='From spectra', is_line=False)
    # lcp.curves_plot(curves_spectra, ax=ax, is_line=True)
    plt.show()



if __name__ == '__main__':
    main()
