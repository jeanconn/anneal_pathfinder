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

T_CCD_REF = -19 # Reference temperature for dark current values in degC
def dark_scale_model(pars, t_ccd):
    """
    dark_t_ref : dark current of a pixel at the reference temperature
    scale : dark current model scale factor
    returns : dark_t_ref scaled to the observed temperatures t_ccd
    """
    scale, dark_t_ref = pars
    scaled_dark_t_ref = dark_t_ref / np.exp(np.log(scale) / 4.0 * (t_ccd - T_CCD_REF))
    return scaled_dark_t_ref


def fit_pix_values(t_ccd, esec, id=1):
    logger = logging.getLogger("sherpa")
    logger.setLevel(logging.WARN)
    data_id = id
    ui.set_method('simplex')
    ui.set_stat('cash')
    ui.load_user_model(dark_scale_model, 'model')
    ui.add_user_pars('model', ['scale', 'dark_t_ref'])
    ui.set_model(data_id, 'model')
    ui.load_arrays(data_id,
                   np.array(t_ccd),
                   np.array(esec)
                   )
    model.scale.val = 0.70
    ui.freeze(model.scale)
    # Fit first for dark_t_ref
    ui.fit(data_id)
    ui.thaw(model.scale)
    # And then fit for scale
    # (though dark_t_ref is not frozen)
    ui.fit(data_id)
    return ui.get_fit_results()


def print_info_block(fits, last_dat):
    print("*************************************************")
    print("Time = {}".format(DateTime(last_dat['time']).date))
    print("CCD temperature = {}".format(last_dat['TEMPCD']))
    print("Slot = {}\n".format(last_dat['SLOT']))
    print("Fit values:\n")
    mini_table = []
    for pix_id in sorted(fits):
        fit = fits[pix_id]
        if fit is None:
            mini_table.append([pix_id, last_dat[pix_id], np.nan, np.nan])
            continue
        t_sf = dark_scale_model(fit.parvals, last_dat['TEMPCD'])
        m_sf = dark_scale_model(fit.parvals, -19)
        minus_19_val = last_dat[pix_id] * m_sf / t_sf
        mini_table.append([pix_id, last_dat[pix_id], minus_19_val, fit.parvals[0]])
    mini_table = Table(rows=mini_table,
                       names=['PixId', 'Val', 'Val(-19)', 'Scale'])
    mini_table['Val(-19)'].format = '.2f'
    mini_table['Scale'].format = '.4f'
    print mini_table
    print "*************************************************"


plt.close(1)
plt.close("fitplots")
plt.ion()

if opt.n_brightest == 64:
    N = 8
else:
    N = np.int(np.sqrt(opt.n_brightest))

fig, axes = plt.subplots(N, N, sharex=True, sharey=True,
                         num=1, figsize=(8, 8))
fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
axes[0][0].set_xticklabels([])
axes[0][0].set_yticklabels([])

if opt.plot_fits:
    fitfig, fitaxes = plt.subplots(N, N, sharex=True, sharey=True,
                                   num="fitplots", figsize=(8, 8))
    fitfig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    fitaxes[0][0].set_xticklabels([])
    fitaxes[0][0].set_yticklabels([])



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
            ax.texts = []
            ax.annotate("{}".format(y.name),
                        xy=(0.5, 0.5), xycoords="axes fraction",
                        ha='center', va='center',
                        color='lightgrey')
            if opt.fit_scaling:
                # Use only non-zero pixel data for fits
                nonzero = y != 0
                # Only fit if more than 5 degC spread in t_ccd
                t_ccd = dat['TEMPCD'][nonzero]
                if np.max(t_ccd) - np.min(t_ccd) < 5:
                    fits[y.name] = None
                    continue
                fit = fit_pix_values(t_ccd,
                                     y[nonzero],
                                     id=i_col)
                fits[y.name] = fit
                if opt.plot_fits:
                    fitax = fitaxes[r][c]
                    fitax.clear()
                    fitax.plot(t_ccd, y[nonzero], '.',
                               markersize=2.5, color='red')
                    mp = ui.get_model_plot(i_col)
                    fitax.plot(mp.x, mp.y, 'k')
                    fitax.texts = []
                    fitax.annotate("{}".format(y.name),
                                   xy=(0.5, 0.5), xycoords="axes fraction",
                                   ha='center', va='center',
                                   color='lightgrey')

    print_info_block(fits, dat[-1])
    plt.draw()
    fig.canvas.draw()
    fig.canvas.flush_events()

    sleep(5)
