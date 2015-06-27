__author__ = 'bakl'

import matplotlib.pyplot as plt

import rf.band as band


def plot_griuz():
    plt.title('griuz filter response')
    bands = dict(g='g', i='r+', r='r', u='o', z='*')
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl,  b.resp, v, label=k)

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
    plt.legend()
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def plot_SWIFT():
    plt.title('SWIFT filter response')
    bands = dict(UVM2='m', UVW1='r', UVW2='b', U_UVOT='k', B_UVOT='c', V_UVOT='g')
    for k, v in bands.items():
        b = band.band_by_name(k)
        plt.plot(b.wl,  b.resp, v, label=k)
    plt.legend()
    plt.ylabel('Amplitude Response')
    plt.xlabel('Wave [A]')
    plt.grid()
    plt.show()


def main():
    plot_UBVRI()
    plot_griuz()
    plot_SWIFT()

if __name__ == '__main__':
    main()