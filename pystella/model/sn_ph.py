import os
import numpy as np


def read(name, path='./', t_diff=1.005, t_beg=0.0, t_end=float('inf'), is_nfrus=True):
    from pystella.rf.spectrum import SeriesSpectrum, Spectrum

    # read first line with frequencies
    fname = os.path.join(path, name + '.ph')
    f = open(fname, 'r')
    try:
        header1 = f.readline()
    finally:
        f.close()

    freqs = [float(x) for x in header1.split()]
    freqs = np.array(freqs).reshape(-1)
    freqs = [10. ** nu for nu in freqs]
    # freqs = np.exp(math.log(10) * freqs)

    data = np.loadtxt(fname, comments='!', skiprows=1)

    times = np.array(data[:, 0])
    is_times = np.zeros(len(times), dtype=bool)
    k = 1
    for i in range(len(times)):
        if times[i] < t_beg:
            k = i
            continue
        if times[i] > t_end:
            break
        if np.abs(times[i] / times[k]) > t_diff:
            is_times[k] = True
            k = i
    is_times[0] = True  # times[0] > 0.
    is_times[-1] = True

    series = SeriesSpectrum(name)
    for i in range(len(times)):
        if is_times[i]:
            t = times[i]
            if is_nfrus:
                nfrus = int(data[i, 1])  # exact number of used (saved) freqs
                freqs = freqs[:nfrus]
                fl = np.array(data[i, 3:nfrus + 3])
            else:
                fl = np.array(data[i, 3:])
            fl[fl < 0] = 0.
            fl = np.exp(np.log(10) * fl)
            s = Spectrum(name, freq=freqs, flux=fl, is_sort_wl=True)
            series.add(t, s)

    series.set_freq(freqs)
    # series.set_times(times_thin)
    return series
