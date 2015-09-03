#!/usr/bin/env python

from time import sleep
import argparse

import logging
from scipy.ndimage.filters import median_filter
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from sherpa import ui
from Chandra.Time import DateTime


def get_opt():
    parser = argparse.ArgumentParser(description='Plot pixel values in real time')
    parser.add_argument('--pix-filename',
                        default='pixel_values.dat',
                        help='Input pixel values filename')

    parser.add_argument('--start',
                        help='Start time (default=2000:001)')

    parser.add_argument('--fit-scaling',
                        action='store_true',
                        help='Fit the dark and scale factor (default=False)')
    parser.add_argument('--plot-fits',
                        action='store_true')
    parser.add_argument('--n-brightest',
                        default=64,
                        type=int,
                        help='Plot the N brightest (either 64 or n**2)')

    args = parser.parse_args()
    return args

opt = get_opt()


def scale_model(pars, t_ccd):
    scale, y_const, x_const = pars
    return 1 / np.exp(np.log(scale) / 4.0 * (t_ccd - x_const))


def fit_pix_values(t_ccd, esec):
    logger = logging.getLogger("sherpa")
    logger.setLevel(logging.WARN)
    data_id = 1
    ui.set_method('simplex')
    ui.set_stat('cash')
    ui.load_user_model(scale_model, 'model')
    ui.add_user_pars('model', ['scale', 'y_const', 'x_const'])
    ui.set_model(data_id, 'model')
    ui.load_arrays(data_id,
                   np.array(t_ccd),
                   np.array(esec)
                   )
    model.scale.val = 0.70
    ui.fit(data_id)
    return ui.get_fit_results()


def print_info_block(fits, dat):
    print("*************************************************")
    print("Time = {}".format(DateTime(dat[-1]['time']).date))
    print("CCD temperature = {}".format(dat[-1]['TEMPCD']))
    print("Slot = {}\n".format(dat[-1]['SLOT']))
    print("Fit values:\n")
    mini_table = []
    for pix_id in sorted(fits):
        fit = fits[pix_id]
        t_sf = scale_model(fit.parvals, dat[-1]['TEMPCD'])
        m_sf = scale_model(fit.parvals, -19)
        minus_19_val = dat[-1][pix_id] * m_sf / t_sf
        mini_table.append([pix_id, dat[-1][pix_id], minus_19_val, fit.parvals[0]])
    mini_table = Table(rows=mini_table,
                       names=['PixId', 'Val', 'Val(-19)', 'Scale'])
    mini_table['Val(-19)'].format='.2f'
    mini_table['Scale'].format='.4f'
    print mini_table
    print "*************************************************"


plt.close(1)
plt.close(2)
plt.ion()

if opt.n_brightest == 64:
    N = 8
    n_fig = 1
else:
    N = np.int(np.sqrt(opt.n_brightest))
    n_fig = 1

figs = []
axess = []
for num in range(n_fig):
    fig, axes = plt.subplots(N, N, sharex=True, sharey=True,
                             num=num, figsize=(8, 8))
    figs.append(fig)
    axess.append(axes)
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    axes[0][0].set_xticklabels([])
    axes[0][0].set_yticklabels([])

colnames = ['r{}_c{}'.format(r, c)
            for r in range(8)
            for c in range(8)]

start = DateTime(opt.start or '2000:001')

while True:
    dat = Table.read(opt.pix_filename, format='ascii.basic', guess=False,
                     fast_reader=True)
    dat = dat[dat['time'] > start.secs]
    dat['dt'] = dat['time'] - dat['time'][0]

    for colname in colnames:
        dat[colname] = median_filter(dat[colname], 5)
    maxes = [np.max(dat[colname]) for colname in colnames]

    i_brightest = np.argsort(maxes)[-opt.n_brightest:]
    cols = []
    for i, colname in enumerate(colnames):
        if i in i_brightest:
            cols.append(dat[colname])

    for im in range(n_fig):
        fig = figs[im]
        axes = axess[im]
        i_col = 0
        fits = {}
        for r in range(N):
            for c in range(N):
                ax = axes[r][c]
                x = dat['dt']
                y = cols[i_col]
                i_col += 1
                if ax.lines:
                    l0 = ax.lines[0]
                    l0.set_data(x, y)
                    ax.relim()
                    ax.autoscale_view()
                else:
                    ax.plot(x, y)
                ax.annotate("{}".format(y.name),
                            xy=(0.5, 0.5), xycoords="axes fraction",
                            ha='center', va='center',
                            color='lightgrey')
                if opt.fit_scaling:
                    nonzero = y != 0
                    fit = fit_pix_values(dat['TEMPCD'][nonzero],
                                         y[nonzero])
                    fits[y.name] = fit
                    if opt.plot_fits:
                        plt.figure("fitfig_{}".format(y.name))
                        ui.plot_fit(replot=True, overplot=True)

        print_info_block(fits, dat)
        plt.draw()
        fig.canvas.draw()
        fig.canvas.flush_events()

    sleep(5)
