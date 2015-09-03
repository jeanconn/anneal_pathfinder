#!/usr/bin/env python

from time import sleep
import argparse

from scipy.ndimage.filters import median_filter
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
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

    parser.add_argument('--n-brightest',
                        default=64,
                        type=int,
                        help='Plot the N brightest (either 64 or n**2)')

    args = parser.parse_args()
    return args

opt = get_opt()

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

        plt.draw()
        fig.canvas.draw()
        fig.canvas.flush_events()

    sleep(5)
