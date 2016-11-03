import sys
import time

barWidth = 40


class ProgressBar:
    """
    Provides a simple progress bar class
    """

    def __init__(self, nsteps, width=barWidth):
        self._start = self._stop = time.time()
        self._nsteps = nsteps
        self._width = width
        self._status = ""

    def update(self, step, text=''):
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
        text = "\rProcessing {}: [{}] {:.1%} {:.3}".format(text, "#" * block + "-" * (self._width - block),
                                                           progress,
                                                           self._status)
        sys.stdout.write(text)
        sys.stdout.flush()
