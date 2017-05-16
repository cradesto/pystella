#!/usr/bin/python3
# -*- coding: utf-8 -*-
import os
import readline
import subprocess
from cmd import Cmd
from os.path import dirname
import numpy as np

ROOT_DIRECTORY = dirname(os.path.abspath(__file__))

# hist_file = os.path.join(ROOT_DIRECTORY, '.pystella_history')
hist_file = os.path.expanduser('~/.pystella_history')
hist_file_size = 1000


class HistConsole(Cmd):
    @property
    def hfile(self):
        local = '.pystella_history'
        if os.path.exists(local):
            fname = local
        else:
            fname = hist_file
        return fname

    def preloop(self):
        if readline and os.path.exists(self.hfile):
            readline.read_history_file(self.hfile)

    def postloop(self):
        if readline:
            readline.set_history_length(hist_file_size)
            readline.write_history_file(self.hfile)

    @staticmethod
    def do_h(args):
        """Show history
        :type args: number of commands
        """
        try:
            lim = int(args)
        except:
            lim = 10

        hlen = readline.get_current_history_length()
        for i in np.arange(max(0, hlen-lim), hlen):
            print("{0}: {1}".format(i, readline.get_history_item(i)))


class MyPrompt(HistConsole):

    @staticmethod
    def do_snespace(args):
        """Plot SN using json-file from https://sne. For detailed help type 'snespace -h'.
        """
        script = os.path.join(ROOT_DIRECTORY, 'sne_space.py')
        subprocess.call("{0} {1}".format(script, args), shell=True)

    @staticmethod
    def do_snec(args):
        """Convert SNEC models to  Stella models. For detailed help type 'snec -h'.
        """
        script = os.path.join(ROOT_DIRECTORY, 'snec.py')
        subprocess.call("{0} {1}".format(script, args), shell=True)

    @staticmethod
    def do_ls(args):
        """Show Stella models. For detailed help type 'ls -h'.
        """
        script = os.path.join(ROOT_DIRECTORY, 'ls.py')
        subprocess.call("{0} {1}".format(script, args), shell=True)

    @staticmethod
    def do_ubv(args):
        """Plot Stella Light Curves. For detailed help type 'ubv -h'.
        """
        if len(args) == 0:
            name = 'Please provide a ph-file.'
        else:
            name = args
        print("ubv %s" % name)
        script = os.path.join(ROOT_DIRECTORY, 'ubv.py')
        subprocess.call("{0} {1}".format(script, args), shell=True)

    @staticmethod
    def do_spec(args):
        """Plot Stella Light Curves. For detailed help type 'spec -h'.
        """
        if len(args) == 0:
            name = 'Please provide a ph-file.'
        else:
            name = args
        print("spec %s" % name)
        script = os.path.join(ROOT_DIRECTORY, 'plot_spec.py')
        subprocess.call("{0} {1}".format(script, args), shell=True)

    @staticmethod
    def do_eve(args):
        """Plot Stella Light Curves. For detailed help type 'eve -h'.
        """
        if len(args) == 0:
            name = 'No data. Please provide a rho-file.'
        else:
            name = args
        print("eve %s" % name)
        script = os.path.join(ROOT_DIRECTORY, 'eve.py')
        subprocess.call("{0} {1}".format(script, args), shell=True)

    @staticmethod
    def do_swd(args):
        """Plot  Stella Shock Wave Details. For detailed help type 'swd -h'.
        """
        if len(args) == 0:
            name = 'No data. Please provide a swd-file.'
        else:
            name = args
        print("swd %s" % name)
        script = os.path.join(ROOT_DIRECTORY, 'swd.py')
        subprocess.call("{0} {1}".format(script, args), shell=True)

    @staticmethod
    def do_bands(args):
        """Plot passbands. For detailed help type 'bands -h'.
        """
        if len(args) == 0:
            name = 'All the passbands.'
        else:
            name = args
        print("bands %s" % name)
        script = os.path.join(ROOT_DIRECTORY, 'plot_bands.py')
        subprocess.call("{0} {1}".format(script, args), shell=True)

    def do_quit(self, args):
        """Quits the program."""
        self.postloop()
        print("Do Svidanya!")
        raise SystemExit

    def do_q(self, args):
        """Quits the program."""
        self.do_quit(args)

    def do_ipython(self, args):
        """Run ipython shell"""

        import math
        import matplotlib.pyplot as plt
        from matplotlib import gridspec

        import pystella.util.callback as cb
        from pystella import velocity as vel
        from pystella.rf import band
        from pystella.model.stella import Stella
        from pystella.rf import light_curve_func as lcf
        from pystella.rf import light_curve_plot as lcp
        from pystella.util.path_misc import get_model_names
        from pystella.util.phys_var import cosmology_D_by_z

        from IPython import embed
        embed(banner1="Hit Ctrl-D to exit interpreter and continue pystella",
              exit_msg="Back in pystella, moving along...")


if __name__ == '__main__':
    prompt = MyPrompt()
    prompt.prompt = 'pystella> '
    prompt.cmdloop('Let''s study a supernova ... \n type "help" for help ')
