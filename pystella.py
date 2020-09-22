#!/usr/bin/env python3
# #!/usr/bin/python3

import os
import readline
import subprocess
import shlex
from cmd import Cmd
from os.path import dirname
import argparse

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
    options = {}

    @staticmethod
    def insert_options(line):
        res = line
        if len(MyPrompt.options) > 0:
            for k, v in MyPrompt.options.items():
                s = '$'+k
                if s in res:
                    res = res.replace(s, v)
        return res

    def bad_insert_options(line):
        words = line.split()
        keys = [w[1:] for w in words if w.startswith('$')]
        res = line
        if len(keys) > 0:
            for k in keys:
                v = MyPrompt.options.get(k, False)
                if v:
                    res = res.replace('$'+k, v)
                else:
                    print(f'You should set up the var: {k}')
        return res

    @staticmethod
    def do_set(args, sep='=', sep_group=';'):
        """Set global options via name = value.
        To insert option to command use prefix: $.
        For example:
           pystella>set ebv=0.07
           pystella>set mdl= cR500M20Ni06_3eps
           pystella>ubv -i $mdl -e $ebv
           >>/home/bakl/Sn/Release/python/pystella/ubv.py -i cR500M20Ni06_3eps -e 0.07
        Run 'set' without arguments to see the current options: pystella>set
        You may set  many options separating them via semicolon ;
           pystella>set ebv=0.07; mdl= cR500M20Ni06_3eps
        To remove name: pystella>set name= None
        """
        if not args:
            print(MyPrompt.options)
            return

        if sep not in args:
            print(f'Format:  key {sep} value')
            return
        if sep_group not in args:
            args += sep_group  # add 1 sep_group

        # print('Before: ', MyPrompt.options)
        for group in args.split(sep_group):
            group = group.strip()
            if sep not in group:
                print(f'No {sep} in "{group}"')
                continue
            k, v = group.split(sep)
            k, v = k.strip(), v.strip()
            print(k, v)
            if v == 'None':
                MyPrompt.options.pop(k, None)  # remove key
            else:
                MyPrompt.options[k] = v.strip()
        print(MyPrompt.options)

    @staticmethod
    def call_cmd(script, args, is_proc=False, is_shell=True):
        # print(args, is_proc)
        if is_proc:
            print(f'Detach thread. script= {script}')
            print(f'   args.strip()= {args.strip()}')
            with subprocess.Popen([script, args.strip()], stdout=subprocess.PIPE, bufsize=1,
                                  shell=True,
                                  universal_newlines=True) as p:
                stdout, stderr = p.communicate()
                for line in stdout:
                    print(line, end='')
        else:
            command = "{0} {1}".format(script, args).strip()
            if '$' in command:
                # print(command)
                command = MyPrompt.insert_options(command)
                print(f">>> {command}")
            subprocess.call(command, shell=is_shell)

    @staticmethod
    def do_snespace(args):
        """Plot SN using json-file from https://sne. For detailed help type 'snespace -h'.
        """
        script = os.path.join(ROOT_DIRECTORY, 'sne_space.py')
        MyPrompt.call_cmd(script, args)

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
            name = 'Please provide a model-file.'
        else:
            name = args
        print(f"ubv {name}")
        script = os.path.join(ROOT_DIRECTORY, 'ubv.py')
        MyPrompt.call_cmd(script, args)

    @staticmethod
    def do_fit(args):
        """Fit with Stella  Obs Light Curves. For detailed help type 'fit -h'.
        """
        if len(args) == 0:
            name = 'Please provide a model-file AND/OR some options'
        else:
            name = args
        print("fit %s" % name)
        script = os.path.join(ROOT_DIRECTORY, 'fit.py')
        MyPrompt.call_cmd(script, args)

    # @staticmethod
    # def do_run(args):
    #     """Run Stella. For detailed help type 'fit -h'.
    #     """
    #     if len(args) == 0:
    #         name = 'Please provide a model-file AND/OR some options'
    #     else:
    #         name = args
    #     print("run %s" % name)
    #     script = os.path.join(ROOT_DIRECTORY, 'run_stella.py')
    #     MyPrompt.call_cmd(script, args)

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
        script = os.path.join(ROOT_DIRECTORY, 'spec.py')
        MyPrompt.call_cmd(script, args)

    @staticmethod
    def do_tau(args):
        """Plot the optical depth of the model. For detailed help type 'tau -h'.
        """
        if len(args) == 0:
            name = 'Please provide a tau-file.'
        else:
            name = args
        print("tau %s" % name)
        script = os.path.join(ROOT_DIRECTORY, 'tau.py')
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
        script = os.path.join(ROOT_DIRECTORY, 'bands.py')
        MyPrompt.call_cmd(script, args)

    @staticmethod
    def do_run(args):
        """Run  the Stella simulations. For detailed help type 'run -h'.
        """
        script = os.path.join(ROOT_DIRECTORY, 'run_stella.py')
        MyPrompt.call_cmd(script, args)

    @staticmethod
    def do_zeta(args):
        """Plot Tcolor-Zeta diagrams. For detailed help type 'zeta -h'.
        """
        if len(args) == 0:
            name = 'Please provide some options'
        else:
            name = args
        print("zeta %s" % name)
        script = os.path.join(ROOT_DIRECTORY, 'zeta.py')
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
