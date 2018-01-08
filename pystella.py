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
        MyPrompt.call_cmd(script, args)

    @staticmethod
    def call_cmd(script, args, is_proc=False, is_shell=True):
        if is_proc:
            with subprocess.Popen([script, args.strip()], stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as p:
                for line in p.stdout:
                    print(line, end='')
        else:
            command = "{0} {1}".format(script, args).strip()
            subprocess.call(command, shell=is_shell)

    @staticmethod
    def do_snec(args):
        """Convert SNEC models to  Stella models. For detailed help type 'snec -h'.
        """
        script = os.path.join(ROOT_DIRECTORY, 'snec.py')
        MyPrompt.call_cmd(script, args)

    @staticmethod
    def do_ls(args):
        """Show Stella models. For detailed help type 'ls -h'.
        """
        script = os.path.join(ROOT_DIRECTORY, 'ls.py')
        MyPrompt.call_cmd(script, args)

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
        MyPrompt.call_cmd(script, args)

    @staticmethod
    def do_fit(args):
        """Fit with Stella  Obs Light Curves. For detailed help type 'fit -h'.
        """
        if len(args) == 0:
            name = 'Please provide a ph-file.'
        else:
            name = args
        print("fit %s" % name)
        script = os.path.join(ROOT_DIRECTORY, 'fit.py')
        MyPrompt.call_cmd(script, args)

    @staticmethod
    def do_obs(args):
        """Plot Observed Light Curves. For detailed help type 'obs -h'.
        """
        if len(args) == 0:
            name = 'Please provide a obs. data. You may use "lcobs:fname:marker:dt:dm" '
        else:
            name = args
        print("obs %s" % name)
        script = os.path.join(ROOT_DIRECTORY, 'obs.py')
        MyPrompt.call_cmd(script, args)

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
        MyPrompt.call_cmd(script, args)

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
        MyPrompt.call_cmd(script, args)

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
        MyPrompt.call_cmd(script, args)

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
        MyPrompt.call_cmd(script, args)

    def do_quit(self, args):
        """Quits the program."""
        self.postloop()
        print("Do Svidanya!")
        raise SystemExit

    def do_q(self, args):
        """Quits the program."""
        self.do_quit(args)

    @staticmethod
    def do_ipython(args):
        """Run ipython shell
        Usage:
          It was already done:  import pystella as ps
          so you may:
          > s = ps.Stella('model3xni')
          > curves = s.curves(bands=('U','B','R','V'))
          > ax = ps.light_curve_plot.curves_plot(curves)
          > plt.show()
          OR
          > ax.get_figure().savefig('zz.pdf')
          > ...
        """

        import math
        import matplotlib.pyplot as plt
        from matplotlib import gridspec
        import numpy as np
        #
        import pystella as ps
        # from pystella import callback as cb
        # from pystella import velocity as vel
        # from pystella import band
        # from pystella import Stella
        # from pystella import light_curve_func as lcf
        # from pystella import light_curve_plot as lcp

        from IPython import embed
        embed(banner1="Hit Ctrl-D to exit interpreter and continue pystella",
              exit_msg="Back in pystella, moving along...")


if __name__ == '__main__':
    prompt = MyPrompt()
    prompt.prompt = '\x01\u001b[35m\x02' + 'pystella>' + '\x01\u001b[0m\x02'
    # prompt.prompt = '\x1b[1;35;47m' + 'pystella>' + '\x1b[0m'
    # prompt.prompt = "pystella>"
    prompt.cmdloop('Let''s study a supernova ... \n type "help" for help ')
