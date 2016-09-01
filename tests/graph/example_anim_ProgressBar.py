#!/usr/bin/env python3

#  see https://github.com/matplotlib/matplotlib/issues/6985

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab, cm
from matplotlib.animation import FuncAnimation
import sys
import time

barWidth = 40
n_steps = 100


class ProgressBar:
    """
    Provides a simple progress bar class
    """

    def __init__(self, nsteps, width=barWidth):
        self._start = self._stop = time.time()
        self._nsteps = nsteps
        self._width = width
        self._status = ""

    def update(self, step):
        """
        This function produces and updates a progress bar.
        It was stolen and modified from
        http://stackoverflow.com/a/15860757/1552338
        """
        self._start = self._stop
        self._stop = time.time()
        self._status = self._stop - self._start

        progress = float(step) / float(self._nsteps - 1)
        if progress >= 1:
            progress = 1
            self._status = "Complete...\r\n"
        block = int(round(self._width * progress))
        text = "\rProcessing: [{}] {:.1%} {:.3}".format("#" * block
                                                        + "-" * (self._width
                                                                 - block),
                                                        progress,
                                                        self._status)
        sys.stdout.write(text)
        sys.stdout.flush()


# Create x, y data
delta = 0.5

extent = (-3, 4, -4, 3)

x = np.arange(-3.0, 4.001, delta)
y = np.arange(-4.0, 3.001, delta)
X, Y = np.meshgrid(x, y)

# Boost the upper limit to avoid truncation errors.
levels = np.arange(-2.0, 1.601, 0.4)

fig, ax = plt.subplots(1, 1)  # create our figure and axes
progress = ProgressBar(n_steps)  # initialize the progress bar


def init():
    """
    This is where our contour is created
    """
    sigma = np.random.uniform(0, 1, 4)
    Z1 = mlab.bivariate_normal(X, Y, sigma[0], sigma[1], 0.0, 0.0)
    Z2 = mlab.bivariate_normal(X, Y, sigma[2], sigma[3], 1, 1)
    Z = (Z1 - Z2) * 10

    norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=-abs(Z).max())
    cmap = cm.PRGn

    contf = plt.contourf(X, Y, Z, levels,
                         cmap=cm.get_cmap(cmap, len(levels) - 1),
                         norm=norm)

    return contf


def update(frame_number, contf):
    """
    This is where our contour is updated
    """
    sigma = np.random.uniform(0, 1, 4)
    Z1 = mlab.bivariate_normal(X, Y, sigma[0], sigma[1], 0.0, 0.0)
    Z2 = mlab.bivariate_normal(X, Y, sigma[2], sigma[3], 1, 1)
    Z = (Z1 - Z2) * 10

    norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=-abs(Z).max())
    cmap = cm.PRGn

    contf.set_array(Z)
    contf.set_cmap(cm.get_cmap(cmap, len(levels) - 1))
    contf.set_norm(norm)

    progress.update(frame_number)

# Construct the animation, using the update function as the animation
# director.
contf = init()
anim = FuncAnimation(fig, update, interval=n_steps, fargs=(contf, ))
anim.save("AnimContourf.mp4")
plt.close(fig)
sys.stdout.write("\n")
