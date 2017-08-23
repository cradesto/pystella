import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle

from matplotlib import gridspec
# from matplotlib import colors as mcolors
# import seaborn as sns

from pystella.model import sn_swd
from pystella.rf import band

__author__ = 'bakl'

lc_colors = band.bands_colors()

lc_lntypes = dict(U="-", B="-", V="-", R="-", I="-",
                  UVM2="-.", UVW1="-.", UVW2="-.",
                  F125W=":", F160W="-.", F140W="--", F105W="-.", F435W="--", F606W="-.", F814W="--",
                  u="--", g="--", r="--", i="--", z="--",
                  bol='-')

lines = ["-", "--", "-.", ":"]

markers = {u'D': u'diamond', 6: u'caretup', u's': u'square', u'x': u'x',
           5: u'caretright', u'^': u'triangle_up', u'd': u'thin_diamond', u'h': u'hexagon1',
           u'+': u'plus', u'*': u'star', u'o': u'circle', u'p': u'pentagon', u'3': u'tri_left',
           u'H': u'hexagon2', u'v': u'triangle_down', u'8': u'octagon', u'<': u'triangle_left'}
markers = list(markers.keys())


def lbl(b, band_shift):
    shift = band_shift[b]
    l = b
    if shift == int(shift):
        shift = int(shift)
    if shift > 0:
        l += '+' + str(shift)
    elif shift < 0:
        l += '-' + str(abs(shift))
    return l


def plot_ubv_models(ax, models_dic, bands, **kwargs):
    # bshift=None, xlim=None, ylim=None, colors=lc_colors, is_time_points=False):
    xlim = kwargs.get('xlim', None)
    ylim = kwargs.get('ylim', None)
    bshift = kwargs.get('bshift', None)
    ls1 = kwargs.get('ls1', "-")
    ls_multi = kwargs.get('ls_multi', ":")
    lw = kwargs.get('lw', 2)
    markersize = kwargs.get('markersize', 6)
    is_time_points = kwargs.get('is_time_points', False)

    is_compute_x_lim = xlim is None
    is_compute_y_lim = ylim is None

    t_points = [0.2, 1, 2, 3, 4, 5, 10, 20, 40, 80, 150]
    colors = band.bands_colors()
    band_shift = dict((k, 0) for k, v in colors.items())  # no y-shift
    if bshift is not None:
        for k, v in bshift.items():
            band_shift[k] = v

    mi = 0
    x_max = []
    y_mid = []
    lc_min = {}
    for mname, mdic in models_dic.items():
        mi += 1
        for bname in bands:
            lc = mdic[bname]
            x = lc.Time
            y = lc.Mag + band_shift[bname]
            bcolor = colors[bname]

            if len(models_dic) == 1:
                ax.plot(x, y, label='%s  %s' % (lbl(bname, band_shift), mname), color=bcolor, ls=ls1, linewidth=lw)
            elif len(models_dic) <= len(lines):
                ax.plot(x, y, label='%s  %s' % (lbl(bname, band_shift), mname), color=bcolor, ls=lines[mi-1],
                        linewidth=lw)
            else:
                ax.plot(x, y, marker=markers[mi % (len(markers) - 1)], label='%s  %s' % (lbl(bname, band_shift), mname),
                        markersize=markersize, color=bcolor, ls=ls_multi, linewidth=lw)

            if is_time_points:
                integers = [np.abs(x - t).argmin() for t in t_points]  # set time points
                for (X, Y) in zip(x[integers], y[integers]):
                    ax.annotate('{:.0f}'.format(X), xy=(X, Y), xytext=(-10, 20), ha='right',
                                textcoords='offset points', color=bcolor,
                                arrowprops=dict(arrowstyle='->', shrinkA=0))
            idx = np.argmin(y)
            lc_min[bname] = (x[idx], y[idx])
            if is_compute_x_lim:
                x_max.append(np.max(x))
            if is_compute_y_lim:
                y_mid.append(np.min(y))

    if is_compute_x_lim:
        xlim = [-10, np.max(x_max) + 10.]
    if is_compute_y_lim:
        ylim = [np.min(y_mid) + 7., np.min(y_mid) - 2.]

    ax.set_xlim(xlim)
    ax.invert_yaxis()
    ax.set_ylim(ylim)

    return lc_min


def plot_models_band(ax, models_dic, bname, **kwargs):
    # , xlim=None, ylim=None,colors=lc_colors, is_time_points=False):
    xlim = kwargs.get('xlim', None)
    ylim = kwargs.get('ylim', None)
    sep = kwargs.get('sep', 'm')  # line separator

    is_compute_x_lim = xlim is None
    is_compute_y_lim = ylim is None

    lw = 1.5
    mi = 0
    x_min = []
    x_max = []
    y_mid = []
    lc_min = {}

    dashes = get_dashes(len(models_dic))

    dashes_cycler = cycle(dashes)

    for mname, mdic in models_dic.items():
        mi += 1
        lc = mdic[bname]
        x = lc.Time
        y = lc.Mag
        # bcolor = colors[mi % len(colors)]

        if len(models_dic) == 1:
            ax.plot(x, y, label='%s  %s' % (bname, mname), ls="-", linewidth=lw)
        else:
            if sep == 'm':
                ax.plot(x, y, marker=markers[mi % (len(markers) - 1)], label='%s  %s' % (bname, mname),
                        markersize=4, ls=":", linewidth=lw)
            else:
                ax.plot(x, y, label='%s  %s' % (bname, mname), dashes=next(dashes_cycler), linewidth=lw)
                # ax.plot(x, y, label='%s  %s' % (bname, mname), ls=lines[mi % (len(lines) - 1)], linewidth=lw)

        idx = np.argmin(y)
        lc_min[bname] = (x[idx], y[idx])
        x_min.append(np.min(x))
        x_max.append(np.max(x))
        if is_compute_y_lim:
            y_mid.append(np.min(y))

    # x-axe
    if is_compute_x_lim:
        xlim = [-10, np.max(x_max) + 10.]
    elif xlim[0] == float('-inf'):
        xlim[0] = np.min(x_min)
    elif xlim[1] == float('inf'):
        xlim[1] = np.max(x_max)
    ax.set_xlim(xlim)

    # y-axe
    ax.invert_yaxis()
    if is_compute_y_lim:
        ylim = [np.min(y_mid) + 7., np.min(y_mid) - 2.]
    ax.set_ylim(ylim)

    return lc_min


def get_dashes(nums):
    dashes = []
    for i in range(nums):
        if i < 5:
            dashes.append((4, 1 + int(i / 2)))
        elif i < 10:
            dashes.append((i - 3, 1 + int(i / 2), 2, 1 + i))
        else:
            dashes.append((i - 8, 1, 2, 1 + int(i / 2), 2, 1 + i))
    return dashes


def plot_models_curves(ax, models_curves, band_shift=None, xlim=None, ylim=None, lc_types=None, colors=lc_colors,
                       lw=2.):
    is_compute_x_lim = xlim is None
    is_compute_y_lim = ylim is None

    if lc_types is None:
        lc_types = dict((name, '-') for name in models_curves.keys())  # solid line

    mi, ib = 0, 0
    x_max = []
    y_mid = []
    for mname, curves in models_curves.items():
        mi += 1
        bands = curves.BandNames
        if band_shift is None:
            mshift = dict((bname, 0.) for bname in bands)  # no y-shift
        for bname in bands:
            ib += 1
            x = curves.TimeDef
            y = curves[bname].Mag + mshift[bname]
            ax.plot(x, y, label='%s  %s' % (lbl(bname, mshift), mname),
                    color=colors[bname], ls=lc_types[mname], linewidth=lw)
            if is_compute_x_lim:
                x_max.append(np.max(x))
            if is_compute_y_lim:
                y_mid.append(np.min(y))

    if is_compute_x_lim:
        xlim = [-10, np.max(x_max) + 10.]
    if is_compute_y_lim:
        ylim = [np.min(y_mid) + 7., np.min(y_mid) - 2.]

    ax.set_xlim(xlim)
    ax.invert_yaxis()
    ax.set_ylim(ylim)


def plot_models_curves_fixed_bands(ax, models_curves, bands, band_shift=None, xlim=None, ylim=None, lc_types=None,
                                   colors=lc_colors,
                                   lw=2.):
    is_compute_x_lim = xlim is None
    is_compute_y_lim = ylim is None

    if band_shift is None:
        band_shift = dict((bname, 0) for bname in bands)  # no y-shift

    if lc_types is None:
        lc_types = dict((name, '-') for name in models_curves.keys())  # solid line

    mi, ib = 0, 0
    x_max = []
    y_mid = []
    for mname, curves in models_curves.items():
        mi += 1
        for bname in bands:
            ib += 1
            x = curves.TimeDef
            y = curves[bname].Mag
            ax.plot(x, y, label='%s  %s' % (lbl(bname, band_shift), mname),
                    color=colors[bname], ls=lc_types[mname], linewidth=lw)
            if is_compute_x_lim:
                x_max.append(np.max(x))
            if is_compute_y_lim:
                y_mid.append(np.min(y))

    if is_compute_x_lim:
        xlim = [-10, np.max(x_max) + 10.]
    if is_compute_y_lim:
        ylim = [np.min(y_mid) + 7., np.min(y_mid) - 2.]

    ax.set_xlim(xlim)
    ax.invert_yaxis()
    ax.set_ylim(ylim)


def plot_bands(dict_mags, bands, title='', fname='', distance=10., xlim=(-10, 200), ylim=(-12, -19),
               is_time_points=True):
    fig = plt.figure()
    ax = fig.add_axes((0.1, 0.3, 0.8, 0.65))
    # ax.set_title(''.join(bands) + ' filter response')

    # colors = band.bands_colors()
    # colors = dict(U="blue", B="cyan", V="black", R="red", I="magenta",
    #               J="blue", H="cyan", K="black",
    #               UVM2="green", UVW1="red", UVW2="blue",
    #               g="black", r="red", i="magenta", u="blue", z="magenta")
    # band_shift = dict(U=6.9, B=3.7, V=0, R=-2.4, I=-4.7,
    #                   UVM2=11.3, UVW1=10, UVW2=13.6,
    #                   u=3.5, g=2.5, r=-1.2, i=-3.7, z=-4.2)
    band_shift = dict((k, 0) for k, v in lc_colors.items())  # no y-shift

    t_points = [2, 5, 10, 20, 40, 80, 150]

    dm = 5 * np.log10(distance) - 5  # distance module
    # dm = 0
    ylim += dm
    is_auto_lim = True
    if is_auto_lim:
        ylim = [0, 0]

    x = dict_mags['time']
    is_first = True
    for n in bands:
        y = dict_mags[n]
        y += dm + band_shift[n]
        ax.plot(x, y, label=lbl(n, band_shift), color=lc_colors[n], ls=lc_lntypes[n], linewidth=2.0)
        # ax.plot(x, y, label=lbl(n, band_shift), color=colors[n], ls=lntypes[n], linewidth=2.0, marker='s')
        if is_time_points and is_first:
            is_first = False
            integers = [np.abs(x - t).argmin() for t in t_points]  # set time points
            for (X, Y) in zip(x[integers], y[integers]):
                plt.annotate('{:.0f}'.format(X), xy=(X, Y), xytext=(10, -30), ha='right',
                             textcoords='offset points',
                             arrowprops=dict(arrowstyle='->', shrinkA=0))

        if is_auto_lim:
            if ylim[0] < max(y[len(y) / 2:]) or ylim[0] == 0:
                ylim[0] = max(y[len(y) / 2:])
            if ylim[1] > min(y) or ylim[1] == 0:
                ylim[1] = min(y)

    ylim = np.add(ylim, [1, -1])
    ax.invert_yaxis()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.legend()
    ax.set_ylabel('Magnitude')
    ax.set_xlabel('Time [days]')
    # ax.set_title(title)
    fig.text(0.17, 0.07, title, family='monospace')
    ax.grid()
    if fname != '':
        plt.savefig("ubv_%s.png" % fname, format='png')
    # ax.show()
    # plt.close()
    return ax


def curves_plot(curves, ax=None, xlim=None, ylim=None, title=None, fname='', **kwargs):
    ls = kwargs.get('ls', {lc.Band.Name: '-' for lc in curves})
    if isinstance(ls, str):
        c = ls.strip()
        ls = {lc.Band.Name: c for lc in curves}
    is_legend = kwargs.get('is_legend', True)
    is_line = kwargs.get('is_line', True)
    if 'lt' in kwargs:
        is_line = False
    lt = kwargs.get('lt', {lc.Band.Name: 'o' for lc in curves})
    colors = kwargs.get('colors', lc_colors)
    linewidth = kwargs.get('linewidth', 2.0)
    markersize = kwargs.get('markersize', 5)
    rect = kwargs.get('rect', (0.1, 0.3, 0.8, 0.65))
    fontsize = kwargs.get('fontsize', 18)
    figsize = kwargs.get('figsize', (20, 10))

    is_new_fig = ax is None
    if is_new_fig:
        plt.matplotlib.rcParams.update({'font.size': 14})
        fig = plt.figure(figsize=figsize)
        # fig = plt.figure(num=None, figsize=(7, 11), dpi=100, facecolor='w', edgecolor='k')
        # ax = fig.add_axes()
        ax = fig.add_axes(rect)
        for item in ([ax.title, ax.xaxis.label,
                      ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fontsize)
            # ax = fig.add_axes((0.1, 0.3, 0.8, 0.65))

    # plt.title(''.join(bands) + ' filter response')
    is_xlim = False
    is_ylim = False
    if xlim is None:
        is_xlim = True
        xlim = [float('inf'), float('-inf')]
    if ylim is None:
        is_ylim = True
        ylim = [float('-inf'), float('inf')]

    for lc in curves:
        x = lc.Time
        y = lc.Mag
        bname = lc.Band.Name
        if is_line:
            ax.plot(x, y, label='{0} {1}'.format(bname, curves.Name),
                    color=colors[bname], ls=ls[bname], linewidth=linewidth)
        else:
            if lc.IsErr:
                yyerr = abs(lc.MagErr)
                ax.errorbar(x, y, label='{0} {1}'.format(bname, fname), yerr=yyerr, fmt=lt[bname],
                            color=colors[bname], ls='')
            else:
                # ax.plot(x, y, label='{0} {1}'.format(bname, fname), color=bcolors[bname], ls='',
                #         marker=marker, markersize=markersize)
                ax.plot(x, y, label='{0} {1}'.format(bname, curves.Name),
                        color=colors[bname], ls='', marker=lt[bname], markersize=markersize)

        if is_xlim:
            xlim[0] = min(xlim[0], np.min(x))
            xlim[1] = max(xlim[1], np.max(x))
        if is_ylim:
            ylim[0] = max(ylim[0], np.max(y))
            ylim[1] = min(ylim[1], np.min(y))

    if is_ylim:
        ylim = [ylim[1] + 10, ylim[1] - 2]
    # ylim = np.add(ylim, [1, -1])
    ax.invert_yaxis()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if is_legend:
        ax.legend()
    ax.set_ylabel('Magnitude')
    ax.set_xlabel('Time [days]')
    ax.grid()
    if title is not None:
        plt.title(title)
        # ax.text(0.17, 0.07, title, family='monospace')
    if fname != '':
        # plt.savefig("ubv_%s.png" % fname, format='png')
        plt.savefig(fname)
    # if is_new_fig:
    #     plt.show()
    # plt.close()
    return ax


def plot_shock_details(swd, times, **kwargs):
    is_legend = kwargs.get('is_legend', True)
    rnorm = kwargs.get('rnorm', 'lgr')
    vnorm = kwargs.get('vnorm', 1e8)
    lumnorm = kwargs.get('lumnorm', 1e40)
    font_size = kwargs.get('font_size', 12)
    is_grid = kwargs.get('is_grid', False)

    xlim = None
    ylim = None
    nrow = len(times)
    ncol = 2
    fig = plt.figure(figsize=(12, nrow * 4))
    plt.matplotlib.rcParams.update({'font.size': font_size})
    # fig = plt.figure(num=None, figsize=(12, len(times) * 4), dpi=100, facecolor='w', edgecolor='k')
    # gs1 = gridspec.GridSpec(len(times), 2)

    for i, t in enumerate(times):
        b = swd.block_nearest(t)
        # plot radius
        ax = fig.add_subplot(nrow, ncol, ncol * i + 1)
        # plot swd(radius)
        sn_swd.plot_swd(ax, b, is_xlabel=(i == len(times) - 1), vnorm=vnorm, lumnorm=lumnorm,
                        rnorm=rnorm, is_legend=False, is_yrlabel=False, text_posy=0.92,
                        is_grid=is_grid)
        x = ax.get_xlim()
        if xlim is None:
            xlim = x
        else:
            xlim = (min(x[0], xlim[0]), max(x[1], xlim[1]))
        y = ax.get_ylim()
        if ylim is None:
            ylim = y
        else:
            ylim = (min(y[0], ylim[0]), max(y[1], ylim[1]))

    # for i, t in list(reversed(list(enumerate(times)))):
    # ylim = None
    for i, t in enumerate(times):
        # plot mass
        b = swd.block_nearest(t)
        ax2 = fig.add_subplot(nrow, ncol, ncol * i + 2)
        sn_swd.plot_swd(ax2, b, is_xlabel=(i == len(times) - 1), vnorm=vnorm, lumnorm=lumnorm,
                        rnorm='m', is_legend=i == 0, is_yllabel=False, text_posy=0.92, is_day=False)
        # Set limits
        ax = fig.add_subplot(nrow, ncol, ncol * i + 1)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax2.set_ylim(ylim)

    fig.subplots_adjust(wspace=0, hspace=0)
    return fig


def plot_shock_details_old(swd, times, **kwargs):
    is_legend = kwargs.get('is_legend', True)
    rnorm = kwargs.get('rnorm', 'm')
    vnorm = kwargs.get('vnorm', 1e8)
    lumnorm = kwargs.get('lumnorm', 1e40)
    font_size = kwargs.get('font_size', 12)

    fig = plt.figure(num=None, figsize=(12, len(times) * 4), dpi=100, facecolor='w', edgecolor='k')
    gs1 = gridspec.GridSpec(len(times), 2)
    plt.matplotlib.rcParams.update({'font.size': font_size})

    i = 0
    for t in times:
        b = swd.block_nearest(t)
        ax = fig.add_subplot(gs1[i, 0])
        sn_swd.plot_swd(ax, b, is_xlabel=i == len(times) - 1, vnorm=vnorm, lumnorm=lumnorm, rnorm=rnorm,
                        is_legend=is_legend)
        ax2 = fig.add_subplot(gs1[i, 1])
        sn_swd.plot_swd(ax2, b, is_xlabel=i == len(times) - 1, vnorm=vnorm, lumnorm=lumnorm,
                        is_legend=is_legend, is_ylabel=False)
        i += 1
    plt.show()
    return fig
