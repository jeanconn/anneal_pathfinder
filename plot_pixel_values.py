#!/usr/bin/env python

from time import sleep

import matplotlib.pyplot as plt
from astropy.table import Table

min_time = 557518477.65
pix_filename = 'pixel_values.dat'

plt.close(1)
plt.close(2)

plt.ion()

figs = []
axess = []
for num in (0, 1):
    fig, axes = plt.subplots(8, 8, sharex=True, sharey=True,
                             num=num, figsize=(8, 8))
    figs.append(fig)
    axess.append(axes)
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    axes[0][0].set_xticklabels([])
    axes[0][0].set_yticklabels([])

while True:
    dat = Table.read(pix_filename, format='ascii.basic', guess=False,
                     fast_reader=True)
    dat = dat[dat['time'] > min_time]
    dat['dt'] = dat['time'] - dat['time'][0]

    for im in (0, 1):
        fig = figs[im]
        axes = axess[im]
        for r in range(8):
            for c in range(8):
                ax = axes[r][c]
                x = dat['dt']
                y = dat['im{}_r{}_c{}'.format(im, r, c)]
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
