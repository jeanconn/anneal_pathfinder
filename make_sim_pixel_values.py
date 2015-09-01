#!/usr/bin/env python

import numpy as np
from time import sleep
import argparse

from Chandra.Time import DateTime
from mica.archive.aca_dark import dark_cal
from Ska.Numpy import interpolate

now = DateTime()
t_readout = 0.016 * 7  # seconds
pix_filename = 'pixel_values.dat'


def get_opt():
    parser = argparse.ArgumentParser(description='Plot pixel values in real time')
    parser.add_argument('--pix-filename',
                        default='pixel_values.dat',
                        help='Input pixel values filename')

    parser.add_argument('--delay',
                        type=float,
                        default=0.0,
                        help='Delay between outputs (secs, default=0)')

    args = parser.parse_args()
    return args


def get_dark_image():
    dark_props = dark_cal.get_dark_cal_props(now, include_image=True)
    image = dark_props['image']
    t_ccd0 = dark_props['ccd_temp']
    return image, t_ccd0


def get_t_ccd(time):
    d_time = np.array([0.0, 600, 900, 1500, 1800, 2400, 2700, 3000, 3600])
    t_ccd = np.array([-14.0, -14, -10, -10, 0, 0, 10, 10, -15])
    out = interpolate(t_ccd, d_time, [time - now.secs], method='linear')
    return out[0]


opt = get_opt()

if 'image' not in globals():
    image, t_ccd0 = get_dark_image()
    r0, c0 = 300, 300
    r1, c1 = 800, 800

pixels = np.concatenate([image[r0:r0+8, c0:c0+8].flatten(),
                         image[r1:r1+8, c1:c1+8].flatten()])

colnames = ['time', 't_ccd'] + ['im{}_r{}_c{}'.format(im, r, c)
                                for im in 0, 1
                                for r in range(8)
                                for c in range(8)]

with open(opt.pix_filename, 'w') as fh:
    fh.write('# Pixel values in e-/sec\n')
    fh.write(' '.join(colnames) + '\n')

for time in np.arange(now.secs, now.secs + 3600, 4.1):
    print(time - now.secs)
    t_ccd = get_t_ccd(time)
    scale = dark_cal.dark_temp_scale(t_ccd0, t_ccd_ref=t_ccd)
    pix_readout_electrons = pixels * scale * t_readout

    # Add gaussian count noise
    count_noise = np.sqrt(pix_readout_electrons)
    pix_readout_electrons += np.random.normal(loc=0,
                                              scale=count_noise,
                                              size=128)

    pix_readout_electrons = np.trunc(pix_readout_electrons / 5) * 5
    pix_readout_e_per_sec = pix_readout_electrons / t_readout

    vals = [time, t_ccd] + pix_readout_electrons.tolist()
    vals = ['{:.2f}'.format(val) for val in vals]
    with open(opt.pix_filename, 'a') as fh:
        fh.write(' '.join(vals) + '\n')

    if opt.delay:
        sleep(opt.delay)
