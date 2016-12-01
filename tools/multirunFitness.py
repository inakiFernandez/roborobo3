# -*- coding: utf-8 -*-
"""
script plot several runs: quartiles and median through time
"""
import numpy as np
import sys
from pylab import *
import brewer2mpl


def read_logfile(fname):
    d = []
    fh = open(fname, 'r')
    for line in fh:
        l = []
        data = line.split()
        for o in data:
            l.append(float(o))
        d.append(l)
    fh.close()
    return d


def perc(data_l):
    data = np.asarray(data_l)
    median = np.zeros(data.shape[0])
    perc_25 = np.zeros(data.shape[0])
    perc_75 = np.zeros(data.shape[0])
    for i in range(0, len(median)):
        median[i] = np.median(data[i, :])
        perc_25[i] = np.percentile(data[i, :], 25)
        perc_75[i] = np.percentile(data[i, :], 75)
        # perc_25[i] = np.percentile(data[i, :], 5)
        # perc_75[i] = np.percentile(data[i, :], 95)
    return median, perc_25, perc_75


def plot_mean_curve(data, color, axis, label):
    mean = np.mean(data, 1)
    axis.plot(mean, lw=2, label=label, color=color)
    axis.grid(axis='y', color="0.9", linestyle='-', linewidth=1)
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.spines['left'].set_visible(False)
    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()
    axis.tick_params(axis='x', direction='out')
    axis.tick_params(axis='y', length=0)
    for spine in axis.spines.values():
        spine.set_position(('outward', 5))
    axis.set_axisbelow(True)


def plot_one_curve(data, color, axis, label, quartiles=False):
    med, perc_25, perc_75 = perc(data)
    print med
    if quartiles:
        axis.fill_between(np.arange(0, len(med)), perc_25, perc_75,
                          alpha=0.25, linewidth=0, color=color)
    axis.plot(med, lw=2, label=label, color=color)
    axis.grid(axis='y', color="0.9", linestyle='-', linewidth=1)
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.spines['left'].set_visible(False)
    axis.get_xaxis().tick_bottom()
    #axis.get_yaxis().tick_left()
    axis.tick_params(axis='x', direction='out')
    axis.tick_params(axis='y', length=0)
    for spine in axis.spines.values():
        spine.set_position(('outward', 5))
    axis.set_axisbelow(True)

def taskIntervals(horSize,interv=25):
    #vertical coordinates for taskswitch
    isT1 = True
    xcoords =  np.arange(0,horSize,interv)
    plt.locator_params(axis='y',nbins=20)
    for xc in xcoords:
        plt.axvline(x=xc,color='gray')
        if(isT1):
            plt.axvspan(xc, xc+interv, facecolor='grey', alpha=0.07)
        else:
            plt.axvspan(xc, xc+interv, facecolor='yellow', alpha=0.07)
        isT1 = not isT1



if __name__ == "__main__":
    # args: file name
    if len(sys.argv) != 2:
        sys.exit("Error: wrong number of arguments\n" +
                 "Usage: python multirunFitness filename")

    dat = read_logfile(sys.argv[1])

    bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
    colors = bmap.mpl_colors
    axis = subplot2grid((1, 1), (0, 0))

    print len(dat[-1])
    plot_one_curve(dat, colors[1], axis, "Collecting", True)
    show()
