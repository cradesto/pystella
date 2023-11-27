import os
import sys
from itertools import cycle
from collections import OrderedDict

import numpy as np
from pystella.rf import band
from pystella.rf import MagBol2Lum

try:
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
except ImportError as ex:
    # import traceback
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fn = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fn, exc_tb.tb_lineno, ex)
    print('  Probably, you should install module: {}'.format('matplotlib'))
    #    print(ex)
    plt = None
    gridspec = None
    pass

__author__ = 'bakl'

lc_colors = band.colors()

lc_lntypes = band.lntypes()

linestyles = ('-', '--', '-.', ':')
linestyles_extend = list( OrderedDict(   # ref https://stackoverflow.com/a/54804349
    [
        ('solid', (0, ())),
        ('dashed', (0, (5, 5))),
        ('dashdotted', (0, (3, 5, 1, 5))),
        ('densely dotted', (0, (1, 1))),
        ('densely dashed', (0, (5, 1))),
        ('dotted', (0, (1, 5))),

        ('loosely dashed', (0, (5, 10))),

        ('loosely dotted', (0, (1, 10))),
        ('loosely dashdotted', (0, (3, 10, 1, 10))),
        ('densely dashdotted', (0, (3, 1, 1, 1))),

        ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
        ('dashdotdotted', (0, (3, 5, 1, 5, 1, 5))),
        ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]).values())

markers = {u'D': u'diamond', 6: u'caretup', u's': u'square', u'x': u'x',
           5: u'caretright', u'^': u'triangle_up', u'd': u'thin_diamond', u'h': u'hexagon1',
           u'+': u'plus', u'*': u'star', u'o': u'circle', u'p': u'pentagon', u'3': u'tri_left',
           u'H': u'hexagon2', u'v': u'triangle_down', u'8': u'octagon', u'<': u'triangle_left'}
markers = list(markers.keys())


# def lbl(b, band_shift):
#     shift = band_shift[b]
#     s = b
#     if shift == int(shift):
#         shift = int(shift)
#     if shift > 0:
#         s += '+' + str(shift)
#     elif shift < 0:
#         s += '-' + str(abs(shift))
#     return s
def lbl(b, band_shift, length=0):
    shift = band_shift[b]
    if shift == int(shift):
        shift = int(shift)
    # if shift == 0:
    #     return b

    s = b
    if shift > 0:
        s += '+'
        s += str(abs(shift))
    elif shift < 0:
        s += '-'
        s += str(abs(shift))

    if length > 0:
        s = ("{0:<" + str(length) + "s}").format(s)
    return s


def lbl_length(bshifts):
    return max((len(lbl(b, bshifts)) for b in bshifts.keys()))


def plot_lum_models(ax, models_dic, bands, **kwargs):
    # bshift=None, xlim=None, ylim=None, colors=lc_colors, is_time_points=False):
    global linestyles
    xlim = kwargs.get('xlim', None)
    ylim = kwargs.get('ylim', None)
    bshift = kwargs.get('bshift', None)
    ls1 = kwargs.get('ls1', "-")
    ls_multi = kwargs.get('ls_multi', ":")
    lw = kwargs.get('lw', 2)
    markersize = kwargs.get('markersize', 6)
    is_time_points = kwargs.get('is_time_points', False)
    is_dashes = kwargs.get('is_dashes', False)
    line_styles = kwargs.get('linestyles', linestyles)
    #    linestyles = kwargs.get('linestyles', ['-'])

    is_compute_x_lim = xlim is None
    is_compute_y_lim = ylim is None

    t_points = [0.2, 1, 2, 3, 4, 5, 10, 20, 40, 80, 150]
    colors = band.colors()
    band_shift = dict((k, 0) for k, v in colors.items())  # no y-shift
    if bshift is not None:
        for k, v in bshift.items():
            band_shift[k] = v

    lbl_len = lbl_length(band_shift)

    mi = 0
    x_max = []
    y_mid = []
    lc_min = {}
    line_cycle = cycle(line_styles)
    dashes = get_dashes(len(bands) + 1, scale=2)

    for mname, mdic in models_dic.items():
        mi += 1
        ls = next(line_cycle)
        for ib, bname in enumerate(bands):
            b = band.band_by_name(bname)
            lc = mdic[bname]
            x = lc.Time
            y = b.mag2lum(lc.Mag)
            bcolor = colors[bname]
            dash = dashes[ib]
            if len(models_dic) == 1:
                if is_dashes:
                    ax.plot(x, y, label='%s  %s' % (lbl(bname, band_shift, lbl_len), mname), color=bcolor, ls=ls1,
                            linewidth=lw, dashes=dash)
                else:
                    ax.plot(x, y, label='%s  %s' % (lbl(bname, band_shift, lbl_len), mname), color=bcolor, ls=ls,
                            linewidth=lw)
            # ax.plot(x, y, label='%s  %s' % (lbl(bname, band_shift), mname), color=bcolor, ls=ls1, linewidth=lw)
            elif len(models_dic) <= len(line_styles):
                ax.plot(x, y, label='%s  %s' % (lbl(bname, band_shift, lbl_len), mname), color=bcolor, ls=ls,
                        linewidth=lw)
            else:
                ax.plot(x, y, marker=markers[mi % (len(markers) - 1)],
                        label='%s  %s' % (lbl(bname, band_shift, lbl_len), mname),
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
                y_mid.append(np.max(y))

    if is_compute_x_lim:
        xlim = [-10, np.max(x_max) + 10.]
    if is_compute_y_lim:
        ylim = [np.max(y_mid)*1e-5, np.max(y_mid) * 5.]

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    return lc_min


def plot_ubv_models(ax, models_dic, bands, **kwargs):
    # bshift=None, xlim=None, ylim=None, colors=lc_colors, is_time_points=False):
    global linestyles
    xlim = kwargs.get('xlim', None)
    ylim = kwargs.get('ylim', None)
    bshift = kwargs.get('bshift', None)
    ls1 = kwargs.get('ls1', "-")
    ls_multi = kwargs.get('ls_multi', ":")
    lw = kwargs.get('lw', 2)
    markersize = kwargs.get('markersize', 6)
    is_time_points = kwargs.get('is_time_points', False)
    is_dashes = kwargs.get('is_dashes', False)
    line_styles = kwargs.get('linestyles', linestyles)
    #    linestyles = kwargs.get('linestyles', ['-'])

    is_compute_x_lim = xlim is None
    is_compute_y_lim = ylim is None

    t_points = [0.2, 1, 2, 3, 4, 5, 10, 20, 40, 80, 150]
    colors = band.colors()
    band_shift = dict((k, 0) for k, v in colors.items())  # no y-shift
    if bshift is not None:
        for k, v in bshift.items():
            band_shift[k] = v

    lbl_len = lbl_length(band_shift)

    mi = 0
    x_max = []
    y_mid = []
    lc_min = {}
    line_cycle = cycle(line_styles)
    dashes = get_dashes(len(bands) + 1, scale=2)

    for mname, mdic in models_dic.items():
        mi += 1
        ls = next(line_cycle)
        for ib, bname in enumerate(bands):
            lc = mdic[bname]
            x = lc.Time
            y = lc.Mag + band_shift[bname]
            bcolor = colors[bname]
            dash = dashes[ib]
            if len(models_dic) == 1:
                if is_dashes:
                    ax.plot(x, y, label='%s  %s' % (lbl(bname, band_shift, lbl_len), mname), color=bcolor, ls=ls1,
                            linewidth=lw, dashes=dash)
                else:
                    ax.plot(x, y, label='%s  %s' % (lbl(bname, band_shift, lbl_len), mname), color=bcolor, ls=ls,
                            linewidth=lw)
            # ax.plot(x, y, label='%s  %s' % (lbl(bname, band_shift), mname), color=bcolor, ls=ls1, linewidth=lw)
            elif len(models_dic) <= len(line_styles):
                ax.plot(x, y, label='%s  %s' % (lbl(bname, band_shift, lbl_len), mname), color=bcolor, ls=ls,
                        linewidth=lw)
            else:
                ax.plot(x, y, marker=markers[mi % (len(markers) - 1)],
                        label='%s  %s' % (lbl(bname, band_shift, lbl_len), mname),
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


def get_dashes(nums, scale=1):
    dashes = []

    for i in range(nums):
        if i < 5:
            dashes.append((4, scale * (1 + int(i / 2))))
        elif i < 10:
            dashes.append((i - 3, scale * (1 + int(i / 2)), 2, 1 + i))
        else:
            dashes.append((i - 8, 1, 2, scale * (1 + int(i / 2)), 2, 1 + i))
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
    bshifts = band_shift
    for mname, curves in models_curves.items():
        mi += 1
        bands = curves.BandNames
        if band_shift is None:
            bshifts = {bname: 0. for bname in bands}  # no y-shift
        for bname in bands:
            ib += 1
            mshift = bshifts[bname]
            x = curves.TimeCommon
            y = curves[bname].Mag + mshift
            ax.plot(x, y, label='%s  %s' % (lbl(bname, bshifts), mname),
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
    # ax.invert_yaxis()
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
            x = curves.TimeCommon
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
        fig.savefig("ubv_%s.png" % fname, format='png')
    # ax.show()
    # plt.close()
    return ax


def curves_plot(curves, ax=None, xlim=None, ylim=None, title=None, fname=None, **kwargs):
    """
    Plot curves.
    If err = -1, it's upper limit, if err = -2 it's lower limit
    :param curves:
    :param ax: Axis. If ax is None, it would be created.
    :param xlim:
    :param ylim:
    :param title:
    :param fname:
    :param kwargs:
        linewidth = kwargs.get('linewidth', 2.0)
        markersize = kwargs.get('markersize', 5)
        fontsize = kwargs.get('fontsize', 18)
        figsize = kwargs.get('figsize', (20, 10))
        legncol = kwargs.get('legncol', 1)
        legloc = kwargs.get('legloc', 1)
        alpha = kwargs.get('alpha', 1.)
        is_legend = kwargs.get('is_legend', True)
        is_line = kwargs.get('is_line', True)
        is_fill = kwargs.get('is_fill', False)
        if 'marker' in kwargs:
            is_line = False
        marker = kwargs.get('marker', 'o')
        if not isinstance(marker, (list, dict, tuple)):
            marker = {lc.Band.Name: marker for lc in curves}
        colors = like {'B': 'blue', 'V': 'green}
        if not isinstance(colors, (list, dict, tuple)):
            colors = {lc.Band.Name: colors for lc in curves}
    :return: ax
    """
    ls = kwargs.get('ls', {lc.Band.Name: '-' for lc in curves})
    if isinstance(ls, str):
        c = ls.strip()
        ls = {lc.Band.Name: c for lc in curves}
    is_legend = kwargs.get('is_legend', True)
    is_line = kwargs.get('is_line', True)
    is_fill = kwargs.get('is_fill', False)
    if 'marker' in kwargs:
        is_line = False
    marker = kwargs.get('marker', 'o')
    if not isinstance(marker, (list, dict, tuple)):
        marker = {lc.Band.Name: marker for lc in curves}
    colors = kwargs.get('colors', None)
    if colors is not None and not isinstance(colors, (list, dict, tuple)):
        colors = {lc.Band.Name: colors for lc in curves}
    linewidth = kwargs.get('linewidth', 2.0)
    markersize = kwargs.get('markersize', 5)
    # rect = kwargs.get('rect', (0.1, 0.2, 0.8, 0.65))
    fontsize = kwargs.get('fontsize', 18)
    figsize = kwargs.get('figsize', (20, 10))
    legncol = kwargs.get('legncol', 1)
    legloc = kwargs.get('legloc', 1)
    alpha = kwargs.get('alpha', 1.)
    flabel = kwargs.get('flabel', None)
    label = kwargs.get('label', None)
    length_lo_up_lims = kwargs.get('length_lo_up_lims', 0.5)

    is_new_fig = ax is None
    if is_new_fig:
        plt.matplotlib.rcParams.update({'font.size': 14})
        fig = plt.figure(figsize=figsize)
        # fig = plt.figure(num=None, figsize=(7, 11), dpi=100, facecolor='w', edgecolor='k')
        # ax = fig.add_axes(rect)
        ax = fig.add_subplot(1, 1, 1)

        for item in ([ax.title, ax.xaxis.label,
                      ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fontsize)
            # ax = fig.add_axes((0.1, 0.3, 0.8, 0.65))

    # ax.set_title(''.join(bands) + ' filter response')
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
        lbl = '{0} {1}'.format(bname, curves.Name.replace("_", ""))
        if flabel is not None:
            lbl = flabel(bname)
        elif label is not None:
            lbl = label.format(bname)

        if colors is not None:
            color = colors[bname]
        else:
            color = band.colors(bname)

        if is_line:
            ax.plot(x, y, label=lbl, color=color, ls=ls[bname], linewidth=linewidth)
        else:
            if lc.IsErr:
                y_el = np.copy(lc.MagErr)
                y_eu = np.copy(lc.MagErr)
                lolims = np.array(y_el == -2, dtype=bool)
                uplims = np.array(y_eu == -1, dtype=bool)
                y_el[lolims] = length_lo_up_lims
                y_eu[uplims] = length_lo_up_lims
                ax.errorbar(x, y, label=lbl, yerr=[y_el, y_eu], fmt=marker[bname],
                            lolims=lolims, uplims=uplims, xlolims=lolims, xuplims=uplims,
                            color=color, ls='', markersize=markersize, )
            else:
                # ax.plot(x, y, label='{0} {1}'.format(bname, fname), color=bcolors[bname], ls='',
                #         marker=marker, markersize=markersize)
                ax.plot(x, y, label=lbl, color=color, ls='', marker=marker[bname], markersize=markersize)
        if is_fill and lc.IsErr:
            yy_err = abs(lc.MagErr)
            # ax.fill(np.concatenate([x, x[::-1]]), np.concatenate([y - yyerr, (y + yyerr)[::-1]]),
            #         alpha=.3, fc=color, ec='None', label=label)  # '95% confidence interval')
            ax.fill_between(x, y - yy_err, y + yy_err, facecolor=color, alpha=alpha)

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
        ax.legend(ncol=legncol, loc=legloc)
    ax.set_ylabel('Magnitude')
    ax.set_xlabel('Time [days]')
    # ax.grid()
    if title is not None:
        ax.get_figure().title(title)
        # ax.text(0.17, 0.07, title, family='monospace')
    if fname is not None:
        print('Save plot to {}'.format(fname))
        ax.get_figure().savefig(fname)
    return ax


def lc_plot(lc, ax=None, xlim=None, ylim=None, title=None, fname=None, **kwargs):
    ls = kwargs.get('ls', {lc.Band.Name: '-'})
    if isinstance(ls, str):
        c = ls.strip()
        ls = {lc.Band.Name: c}
    is_legend = kwargs.get('is_legend', True)
    is_line = kwargs.get('is_line', True)
    if 'lt' in kwargs:
        is_line = False
    lt = kwargs.get('lt', {lc.Band.Name: 'o'})
    if isinstance(lt, str):
        c = lt.strip()
        lt = {lc.Band.Name: c}
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

    x = lc.Time
    y = lc.Mag
    bname = lc.Band.Name
    if is_line:
        ax.plot(x, y, label='{0} {1}'.format(bname, lc.Name),
                color=colors[bname], ls=ls[bname], linewidth=linewidth)
    else:
        if lc.IsErr:
            yyerr = abs(lc.Err)
            ax.errorbar(x, y, label='{0} {1}'.format(bname, fname), yerr=yyerr, fmt=lt[bname],
                        color=colors[bname], ls='', markersize=markersize)
        else:
            # ax.plot(x, y, label='{0} {1}'.format(bname, fname), color=bcolors[bname], ls='',
            #         marker=marker, markersize=markersize)
            ax.plot(x, y, label='{0} {1}'.format(bname, lc.Name),
                    color=colors[bname], ls='', marker=lt[bname], markersize=markersize)

    if is_xlim:
        xlim[0] = np.min(x)
        xlim[1] = np.max(x)
    if is_ylim:
        ylim[0] = np.max(y)
        ylim[1] = np.min(y)

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
    if fname is not None:
        # plt.savefig("ubv_%s.png" % fname, format='png')
        plt.savefig(fname)
    return ax


def ticks_on(ax, minor=3, major=6):
    ax.minorticks_on()
    ax.tick_params(direction='in', which='minor', length=minor)
    ax.tick_params(direction='in', which='major', length=major)
    return ax


def plot_shock_details(swd, times, **kwargs):
    from pystella.model import sn_swd

    is_legend = kwargs.get('is_legend', False)
    # ylim_par = kwargs.get('ylim_par', (0.001, 11))
    font_size = kwargs.get('font_size', 12)
    # is_grid = kwargs.get('is_grid', False)
    is_adjust = kwargs.get('is_adjust', True)
    is_axes = kwargs.get('is_axes', False)
    dic_axes = kwargs.get('dic_axes', None)
    is_ax_old = False
    xlim = None
    ylim_rho = None
    nrow = len(times)
    ncol = 2

    axes1 = []
    if dic_axes is None:
        dic_axes = {'r': [], 'm': []}
        fig = plt.figure(figsize=(12, nrow * 4))
        # plt.minorticks_on()
        # fig = plt.figure(num=None, figsize=(12, len(times) * 4), dpi=100, facecolor='w', edgecolor='k')
        # gs1 = gridspec.GridSpec(len(times), 2)
    else:
        fig = (dic_axes['r'][0]['rho']).get_figure()
        is_ax_old = True
        is_adjust = False
        # is_legend = False
        kwargs['is_day'] = False

    plt.matplotlib.rcParams.update({'font.size': font_size})
    # plot radius column
    for i, t in enumerate(times):
        if is_ax_old:
            axrho, axpar = dic_axes['r'][i]['rho'], dic_axes['r'][i]['par']
        else:
            axrho = fig.add_subplot(nrow, ncol, ncol * i + 1, label='radius {}'.format(i))
            axpar = None
        axes1.append(axrho)
        legmask = sn_swd.LEGEND_MASK_None
        if is_legend and i == 0:
            legmask = sn_swd.LEGEND_MASK_Rho
        # plot swd(radius)
        b = swd.block_nearest(t)
        axrho, axpar = sn_swd.plot_swd((axrho, axpar), b, name=swd.Name, is_xlabel=(i == len(times) - 1),
                                       legmask=legmask, is_yrlabel=False, text_posy=0.88,
                                       **kwargs)
        if not is_ax_old:
            x = axrho.get_xlim()
            if xlim is None:
                xlim = x
            else:
                xlim = (min(x[0], xlim[0]), max(x[1], xlim[1]))
            y = axrho.get_ylim()
            if ylim_rho is None:
                ylim_rho = y
            else:
                ylim_rho = (min(y[0], ylim_rho[0]), max(y[1], ylim_rho[1]))
            # axpar.tick_params(direction='in', which='both', length=4)
            ticks_on(axrho)
            ticks_on(axpar)
            dic_axes['r'].append({'itime': i, 't': t, 'rho': axrho, 'par': axpar})

    if 'axeX' in kwargs:
        kwargs.pop('axeX')
    axes2 = []
    # Plot mass column
    for i, t in enumerate(times):
        if is_ax_old:
            # ax2 = dic_axes['m'][i]['rho']
            axrho, axpar = dic_axes['m'][i]['rho'], dic_axes['m'][i]['par']
        else:
            axrho = fig.add_subplot(nrow, ncol, ncol * i + 2, label='mass {}'.format(i))
            axrho.tick_params(direction='in', which='minor', length=3)
            axpar = None
        axes2.append(axrho)
        legmask = sn_swd.LEGEND_MASK_None

        if is_legend and i == 0:
            legmask = sn_swd.LEGEND_MASK_Vars
        b = swd.block_nearest(t)
        axrho, axpar = sn_swd.plot_swd((axrho, axpar), b, name=swd.Name, is_xlabel=(i == len(times) - 1),
                                       axeX='m', legmask=legmask, is_yllabel=False, text_posy=0.88,
                                       **kwargs)
        if not is_ax_old:
            dic_axes['m'].append({'itime': i, 't': t, 'rho': axrho, 'par': axpar})
            ticks_on(axrho)
            axpar.tick_params(direction='in', which='major', length=5)
            ticks_on(axpar)

    # Set limits
    for i, ax in enumerate(axes1):
        ax.set_xlim(xlim)
        ax.set_ylim(ylim_rho)
        # remove labels between subplots
        if not i == len(times) - 1:
            plt.setp(ax.get_xticklabels(), visible=False)

    for i, ax2 in enumerate(axes2):
        ax2.set_ylim(ylim_rho)
        # remove labels between subplots
        if not i == len(times) - 1:
            plt.setp(ax2.get_xticklabels(), visible=False)

    if is_adjust:
        fig.subplots_adjust(wspace=0., hspace=0.)

    # print(len(axes1), len(axes2))
    if is_axes:
        return fig, dic_axes
    return fig


def plot_swd_chem(dic_axes, argsrho, path_relativ, alpha=0.25):
    """Add chemical data to the plot_shock_details plot
    :type dic_axes: dict, dic_axes['r' or 'm'].append({'itime': i, 't': t, 'rho': axrho, 'par': axpar})
    :param argsrho: rho-file + Elements
    :param path_relativ:
    :param alpha: the transparanc
    """
    from pystella.model import sn_eve
    from pystella.util.phys_var import phys

    colors = sn_eve.eve_colors
    s_elements = None
    if '+' in argsrho:
        frho, s_elements = argsrho.split('+')
    else:
        frho = argsrho
    frho = os.path.join(path_relativ, frho)
    pathrho, namerho = os.path.split(frho)

    eve = sn_eve.load_rho(os.path.join(pathrho, namerho))
    elements = eve.Elements
    if s_elements is not None:
        elements = s_elements.split(':')
    print('Add chem [{}] from {}'.format(':'.join(elements), frho))
    x = eve.m / phys.M_sun
    if min(x) > 0.1:
        print('Chem is shifted at {} Msun'.format(-min(x)))
        x = x - min(x)

    def two_y(x, lims, up=0.85):
        return (lims[1]*up - lims[0]) * x + lims[0]

    for i, d in enumerate(dic_axes['m']):
        ax = d['par']
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        # print(ylim)
        xx = np.linspace(xlim[1]*0.8, xlim[0], len(elements))
        yy_tot = np.zeros_like(x)  # + ylim[1]
        yy_prev = np.zeros_like(x)  # + ylim[0]
        for ie, el in enumerate(elements):
            # for el in reversed(elements):
            y = eve.el(el)
            # yy += y
            # yy = y
            # yyy = yy
            yy_tot += y
            # yyy = np.log10(yy) + ylim[1]-1
            ax.fill_between(x, y1=two_y(yy_tot, ylim), y2=two_y(yy_prev, ylim), color=colors[el], alpha=alpha)
            # ax.plot(x, yyy, color=colors[el], alpha=alpha)
            # Print the element name
            # xp = np.average(x, weights=y)
            # yp = np.average(yyy, weights=y) * 0.9
            xp = xx[ie]
            ax.text(xp, ylim[1] * 0.9, el, color=colors[el], fontsize=11)
            yy_prev = yy_tot.copy()


def plot_swd_tau(dic_axes, stella, times, bnames=('B',), tau_ph=2. / 3., is_obs_time=False, **kwargs):
    """Add photospheric data to the plot_shock_details plot
    :type dic_axes: object
    :param stella:
    :param times:
    :param bnames:
    :param tau_ph: the photosphere location. Default: 2/3
    :param is_obs_time:  If True to compute Obs.Time as ProperTime - R(N-1)/c and use them. Default: False
    :param tnorm: the T normalization. Default: None = log10(T)
    :param vnorm: the V normalization. Default: 1e8
    :param alpha: the transparent. Default: 0.5
    """
    from pystella.rf.band import band_by_name

    tnorm = kwargs.get('tnorm', None)
    vnorm = kwargs.get('vnorm', 1e8)
    alpha = kwargs.get('alpha', 0.5)
    markersize = kwargs.get('markersize', 6)
    marker = kwargs.get('marker', 'o')

    if not stella.is_tau:
        print('There is no tau-file for model {} in path: {}'.format(stella.Name, stella.Path))
        return

    print('Add tau [{}] for {} at tau_ph= {:.3f}'.format(':'.join(bnames), stella.Name, tau_ph))
    pars_data = ['T', 'V', 'R']
    tau = stella.get_tau().load(is_info=False)
    tau_data = tau.params_ph(pars=pars_data, moments=times, tau_ph=tau_ph, is_obs_time=is_obs_time)

    # Extract phot data
    data = {bn: {p: [] for p in pars_data} for bn in bnames}

    for bname in bnames:
        b = band_by_name(bname)
        fr_eff = b.freq_eff
        for p in pars_data:
            for i, (t, freq, y) in enumerate(tau_data[p]):
                s = '{:9.4f} '.format(t)
                idx = (np.abs(freq - fr_eff)).argmin()
                s += ' {:10e}'.format(y[idx])
                data[bname][p].append(y[idx])

    # Plot
    for ii, d in enumerate(dic_axes['r']):
        ax = d['par']
        # t = d['t']
        # print('{:9s} {}'.format('t_real',  ' '.join([f'{p}_{b:10s}' for b in bnames])))
        # s = '{:9.4f} '.format(t)
        for i, bname in enumerate(bnames):
            r_ph = np.array(data[bname]['R'])
            color = band.colors(bname)
            # print(bname)
            xr = r_ph[ii]
            ax.text(xr, 0.5 + i, bname, fontsize=12, color=color)
            ax.axvline(x=xr, ymin=0., ymax=0.99, linestyle='--', color=color, alpha=alpha)
            # Temperature
            if tnorm is None:
                yt = np.log10(data[bname]['T'][ii])
                ax.plot(xr, yt, color='green', ls='', marker=marker, markersize=markersize, alpha=alpha)
            else:
                yt = data[bname]['T'][ii] / tnorm
                ax.plot(xr, yt, color='green', ls='', marker=marker, markersize=markersize, alpha=alpha)
            # Velocity
            yv = data[bname]['V'][ii] / vnorm
            ax.plot(xr, yv, color='blue', ls='', marker=marker, markersize=markersize, alpha=alpha)
        #     print(f't={t:9.3f}  R= {xr:e}  V= {yv:e} T= {yt:e}')
        #     s += f'R= {xr:e}  V= {yv:e} T= {yt:e}'
        # print(s)

    # Print
    for p in pars_data:
        print('{:9s} {}'.format(
            't_real', ' '.join([f'{p}({b:4s}:{band_by_name(b).wl_eff_angs:.0f})' for b in bnames])))
        # print(p)
        for ii, d in enumerate(dic_axes['r']):
            t = d['t']
            s = '{:9.4f} '.format(t)
            for i, bname in enumerate(bnames):
                v = np.array(data[bname][p][ii])
                s += f'{v:12e} '
            print(s)


########################################
# from https://stackoverflow.com/questions/7358118/matplotlib-black-white-colormap-with-dashes-dots-etc


def setAxLinesBW(ax, color='black', markersize=3):
    """
    Take each Line2D in the axes, ax, and convert the line style to be
    suitable for black and white viewing.
    """
    # https://matplotlib.org/gallery/lines_bars_and_markers/linestyles.html
    # dashList = [(5, 2), (2, 5), (4, 10), (3, 3, 2, 2), (5, 2, 20, 2)]

    COLORMAP = {
        'blue': {'marker': 's', 'dash': (5, 2, 5, 2, 5, 10)},
        'darkgreen': {'marker': 'x', 'dash': (5, 5)},
        'red': {'marker': '*', 'dash': (5, 3, 1, 3)},
        'cyan': {'marker': 'd', 'dash': (1, 3, 1, 1)},
        'magenta': {'marker': 'o', 'dash': (2, 5)},
        'yellow': {'marker': '<', 'dash': (5, 3, 1, 2, 1, 10)},
        'k': {'marker': 'P', 'dash': (5, 2)}  # (3, 3, 2, 2)  [1,2,1,10]}
    }

    lines_to_adjust = ax.get_lines()
    # try:
    #     lines_to_adjust += ax.get_legend().get_lines()
    # except AttributeError:
    #     pass

    for line in lines_to_adjust:
        origColor = line.get_color()
        origLineType = line.get_linestyle()
        line.set_color(color)
        if origLineType is not None:
            line.set_dashes(COLORMAP[origColor]['dash'])
        else:
            line.set_marker(COLORMAP[origColor]['marker'])
            line.set_markersize(markersize)
            line.set_linestyle('')


def setFigLinesBW(fig):
    """
    Take each axes in the figure, and for each line in the axes, make the
    line viewable in black and white.
    """
    for ax in fig.get_axes():
        setAxLinesBW(ax)
