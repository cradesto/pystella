__author__ = 'bakl'

import matplotlib.pyplot as plt

import rf.band as band


def plot_griuz():
    plt.title('griuz filter response')
    bands = dict(g='g', i='r+', r='r', u='o', z='*')
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl,  b.resp, v, label=k)

#    b = band.band_by_name('g')
#    plt.plot(b.wl,  np.log10(b.resp), 'r')
    plt.legend()
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def plot_UBVRI():
    plt.title('UBVRI filter response')
    bands = dict(U='b', B='c', V='g', R='r', I='p')
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl,  b.resp, v, label=k)

#    b = band.band_by_name('g')
#    plt.plot(b.wl,  np.log10(b.resp), 'r')
    plt.legend()
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def compute_mag():
    pass


def main():
    compute_mag()
    plot_UBVRI()
    plot_griuz()

if __name__ == '__main__':
    main()