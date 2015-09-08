#!/usr/bin/env python

import numpy as np
from time import sleep
import argparse

from Chandra.Time import DateTime
from mica.archive.aca_dark import dark_cal
from Ska.Numpy import interpolate

now = DateTime()
t_readout = 0.016 * 7  # seconds


def get_opt():
    parser = argparse.ArgumentParser(description='Plot pixel values in real time')
    parser.add_argument('--pix-filename1',
                        default='pixel_values1.dat',
                        help='Input pixel values filename')
    parser.add_argument('--pix-filename2',
                        default='pixel_values2.dat',
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

colnames = ['time', 'TEMPCD', 'SLOT', 'INTEG'] + ['r{}_c{}'.format(r, c)
                                                  for r in range(8)
                                                  for c in range(8)]

for filename in [opt.pix_filename1, opt.pix_filename2]:
    with open(filename, 'w') as fh:
        fh.write('# Pixel values in e-/sec\n')
        fh.write(' '.join(colnames) + '\n')


for time in np.arange(now.secs, now.secs + 3600, 4.1):
    print(time - now.secs)
    t_ccd = get_t_ccd(time)
    scale = dark_cal.dark_temp_scale(t_ccd0, t_ccd_ref=t_ccd)
    pix_readout_dn = pixels * scale * t_readout

    # Add gaussian count noise
    count_noise = np.sqrt(pix_readout_dn)
    pix_readout_dn += np.random.normal(loc=0,
                                       scale=count_noise,
                                       size=128)

    pix_readout_dn = np.trunc(pix_readout_dn / 5) * 5

    first_image = pix_readout_dn.tolist()[0:64]
    second_image = pix_readout_dn.tolist()[64:]

    for image, pix_filename, slot in zip(
        [first_image, second_image],
        [opt.pix_filename1, opt.pix_filename2],
        [6, 7]):
        vals = [time, t_ccd, slot, t_readout] + image
        vals = ['{:.2f}'.format(val) for val in vals]
        with open(pix_filename, 'a') as fh:
            fh.write(' '.join(vals) + '\n')

    if opt.delay:
        sleep(opt.delay)
