#!/usr/bin/env python
import os
from time import sleep
import argparse

import logging
import pyyaks.logger
from scipy.ndimage.filters import median_filter
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from sherpa import ui
from Chandra.Time import DateTime
from Ska.Matplotlib import plot_cxctime, cxctime2plotdate, set_time_ticks


def get_opt():
    parser = argparse.ArgumentParser(description='Plot pixel values in real time')
    parser.add_argument('--pix-filename',
                        default='acaslot7.dat',
                        help='Input pixel values filename')
    parser.add_argument('--logfile',
                        help='Output log filename')
    parser.add_argument('--start',
                        default="2015:252:13:35:00.000",
                        help='Start time (default=2000:001)')
    parser.add_argument('--stop',
                        default="2015:252:14:50:00.000",
                        help='Stop time (default=2099:001)')
    parser.add_argument('--plot-fit-curves',
                        action='store_true',
                        help="Plot dark current vs t_ccd curves and fits")
    parser.add_argument('--n-brightest',
                        default=1,
                        type=int,
                        help='Plot the N brightest (must be n**2)')

    args = parser.parse_args()
    return args

opt = get_opt()
if opt.logfile is None:
    root, ext = os.path.splitext(opt.pix_filename)
    opt.logfile = "{}.log".format(root)

pix_log = pyyaks.logger.get_logger(name='pix_log',
                                   filename=opt.logfile,
                                   filemode='a',
                                   level=logging.INFO)


T_CCD_REF = -19 # Reference temperature for dark current values in degC
def dark_scale_model(pars, t_ccd):
    """
    dark_t_ref : dark current of a pixel at the reference temperature
    scale : dark current model scale factor
    returns : dark_t_ref scaled to the observed temperatures t_ccd
    """
    scale, dark_t_ref = pars
    scaled_dark_t_ref = dark_t_ref * np.exp(np.log(scale) / 4.0 * (T_CCD_REF - t_ccd))
    return scaled_dark_t_ref


def fit_pix_values(t_ccd, esec, id=1):
    logger = logging.getLogger("sherpa")
    logger.setLevel(logging.WARN)
    data_id = id
    ui.clean()
    ui.set_method('simplex')
    ui.load_user_model(dark_scale_model, 'model')
    ui.add_user_pars('model', ['scale', 'dark_t_ref'])
    ui.set_model(data_id, 'model')
    ui.load_arrays(data_id,
                   np.array(t_ccd),
                   np.array(esec),
                   0.1*np.ones(len(t_ccd)),
                   )
    model.scale.val = 0.70
    model.dark_t_ref.val = 500
    ui.freeze(model.scale)
    # If more than 5 degrees in the temperature range,
    # thaw and fit for model.scale.  Else just use/return
    # the fit of dark_t_ref
    ui.fit(data_id)
    ui.thaw(model.scale)
    ui.fit(data_id)
    return ui.get_fit_results(), ui.get_model(data_id)


def print_info_block(fits, last_dat):
    pix_log.info("*************************************************")
    pix_log.info("Time = {}".format(DateTime(last_dat['time']).date))
    pix_log.info("CCD temperature = {}".format(last_dat['TEMPCD']))
    pix_log.info("Slot = {}\n".format(last_dat['SLOT']))
    pix_log.info("Fit values:\n")
    mini_table = []
    other_t_ccd = [-10, -5, 0, 5, 10]
    for pix_id in sorted(fits):
        fitinfo = fits[pix_id]
        if fitinfo is None:
            continue
        m = fitinfo['modpars']
        dc = dark_scale_model((m.scale.val, m.dark_t_ref.val), last_dat['TEMPCD'])
        ref_dc = dark_scale_model((m.scale.val, m.dark_t_ref.val), -19)
        scale_factor = ref_dc / dc
        rec_esec = last_dat[pix_id] * GAIN / last_dat['INTEG']
        minus_19_esec =  rec_esec * scale_factor
        new_rec = [pix_id, rec_esec, minus_19_esec, m.scale.val, dc / ref_dc]
        for t_ccd in other_t_ccd:
            dc_temp = dark_scale_model((m.scale.val, m.dark_t_ref.val), t_ccd)
            new_rec.append(dc_temp / ref_dc)
        mini_table.append(new_rec)
    if not len(mini_table):
        return
    colnames = ['PixId', 'e-/sec', 'e-/sec(-19)', 'Scale', 'r({:.1f})'.format(last_dat['TEMPCD'])]
    for t_ccd in other_t_ccd:
        colnames.append("r({})".format(t_ccd))
    mini_table = Table(rows=mini_table,
                       names=colnames)
    mini_table['e-/sec'].format = '.2f'
    mini_table['e-/sec(-19)'].format = '.2f'
    mini_table['Scale'].format = '.4f'
    for col in mini_table.colnames:
        if col.startswith('r('):
            mini_table[col].format = '.2f'
    pix_log.info(mini_table)
    pix_log.info("*************************************************")


def plot_two(fig_id, x, y, x2, y2,
             linestyle='-', linestyle2='-',
             color='blue', color2='magenta',
             ylim=None, ylim2=None,
             xlabel='', ylabel='', ylabel2='', title='',
             figsize=(7, 3.5),
             markersize2=None, marker2=None,
             ):
    """Plot two quantities with a date x-axis"""
    xt = cxctime2plotdate(x)
    fig = plt.figure(fig_id, figsize=figsize)
    fig.clf()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot_date(xt, y, fmt='-', linestyle=linestyle, color=color)
    ax.set_xlim(min(xt), max(xt))
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid()

    ax2 = ax.twinx()

    xt2 = cxctime2plotdate(x2)
    ax2.plot_date(xt2, y2, fmt='-', linestyle=linestyle2, color=color2,
                  markersize=markersize2, marker=marker2)
    ax2.set_xlim(min(xt), max(xt))
    if ylim2:
        ax2.set_ylim(*ylim2)
    ax2.set_ylabel(ylabel2, color=color2)
    ax2.xaxis.set_visible(False)

    set_time_ticks(ax)
    [label.set_rotation(30) for label in ax.xaxis.get_ticklabels()]
    [label.set_color(color2) for label in ax2.yaxis.get_ticklabels()]

    fig.subplots_adjust(bottom=0.22)

    return {'fig': fig, 'ax': ax, 'ax2': ax2}


plt.close(1)
plt.close("separate")
plt.close("overplots")
plt.ion()

best_pix = 'r2_c2'

N = np.int(np.sqrt(opt.n_brightest))

onefig = plt.figure(figsize=(5.5,2.75), num='overplots')
#dualfig, axes = plt.subplots(1, 2, figsize=(5,2.5), num='separate')
#fig, axes = plt.subplots(N, N, sharex=True, sharey=True,
#                         num=1, figsize=(8, 8))
#fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
#if N > 1:
#    axes[0][0].set_xticklabels([])
#    axes[0][0].set_yticklabels([])

#if opt.plot_fit_curves:
#    fitfig, fitaxes = plt.subplots(N, N, sharex=True, sharey=True,
#                                   num="fitplots", figsize=(8, 8))
#    fitfig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
#    fitaxes[0][0].set_xticklabels([])
#    fitaxes[0][0].set_yticklabels([])


colnames = ['r{}_c{}'.format(r, c)
            for r in range(8)
            for c in range(8)]

start = DateTime(opt.start or '2000:001')
stop = DateTime(opt.stop or '2099:001')
GAIN = 5.0


dat = Table.read(opt.pix_filename, format='ascii.basic', guess=False,
                 fast_reader=True)
# Filter:
#  first record (usually bad in splat output)
#  unknown/unset temperature data (TEMPCD = -99)
#  records with very large INTEG time (bad decom?)
#  records before set start time
dat = dat[1:]
dat = dat[dat['TEMPCD'] != -99]
dat = dat[dat['INTEG'] < 2.0]
dat = dat[dat['time'] > start.secs]
dat = dat[dat['time'] < stop.secs]
dat['dt'] = dat['time'] - dat['time'][0]
integ = dat['INTEG']

#ccdfig = plt.figure("ccdplot", figsize=(1,1))
#ccdax = axes[0] #plt.gca()
#if ccdax.lines:
#    ccdline = ccdax.lines[0]
#    ccdline.set_data(cxctime2plotdate(dat['time']), dat['TEMPCD'])
#    ccdax.relim()
#    ccdax.autoscale_view()
#else:
#    plot_cxctime(dat['time'], dat['TEMPCD'], 'b.', markersize=4, ax=ccdax)
#    ccdax.grid('on')
#
#
#for colname in colnames:
#    dat[colname] = median_filter(dat[colname], 5)
#maxes = [np.max(dat[colname]) for colname in colnames]
#
#i_brightest = np.argsort(maxes)[-opt.n_brightest:]
#cols = []
#for i, colname in enumerate(colnames):
#    if i in i_brightest:
#        cols.append(dat[colname])

i_col = 0
fits = {}
#ax = axes[1]
x = dat['dt']
y = dat[best_pix] * GAIN / integ
#ax.yaxis.set_label_position("right")
#ax.yaxis.tick_right()
#ax.grid()
#plot_cxctime(dat['time'], y, 'k', ax=ax)
#ax.set_ylabel('Dark Current (e-/sec)')
#ccdax.set_ylabel('CCD Temp (DegC)')
#ax.texts = []
#ax.annotate("{}".format(y.name),
#            xy=(0.5, 0.5), xycoords="axes fraction",
#            ha='center', va='center',
#            color='lightgrey')
t_ccd = dat['TEMPCD']
fit, modpars = fit_pix_values(t_ccd,
                              y,
                              id=i_col)
fits[y.name] = {'fit': fit,
                'modpars': modpars}
fitmod = ui.get_model_plot(i_col)
#plot_cxctime(DateTime(dat['time'][0]).secs + x, fitmod.y, color='red', ax=ax)
#plt.tight_layout()
plot2 = plot_two(fig_id='overplots',
                 x2=dat['time'], y2=y,
                 x=dat['time'], y=t_ccd,
                 color='blue', color2='black')
plot_cxctime(DateTime(dat['time'][0]).secs + x, fitmod.y,
             color='red',
             ax=plot2['ax2'],
             linewidth=6, alpha=.5,
             label="Dark current model"
             )
plot2['ax2'].legend(loc='upper left', fontsize=10)
plot2['ax'].yaxis.label.set_color('blue')
plot2['ax'].tick_params(axis='y', colors='blue')
plot2['ax2'].set_ylabel('Dark Current (e-/sec)')
plot2['ax'].set_ylabel('CCD Temp (DegC)')
plot2['ax'].tick_params(axis='x', which='major', labelsize=10)
#plot2['fig'].suptitle('Dark Current, CCD Temp, and Fit vs Time')
onefig.tight_layout()

print_info_block(fits, dat[-1])
plt.draw()
onefig.savefig("annealing.png")

#plot2['ax'].xaxis.set_ticklabels([])
#onefig.tight_layout()
#plt.draw()
#onefig.savefig("annealing_notimes.png")
