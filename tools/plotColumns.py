# -*- coding: utf-8 -*-
"""
script plot a curve per column
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


def plot_curve(data, color, axis, label, quartiles=False):

    axis.plot(data, lw=1.5, label=label, color=color)

    axis.grid(axis='y', color="0.9", linestyle='-', linewidth=1)
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.spines['left'].set_visible(False)
    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()
    axis.tick_params(axis='x', direction='out')
    axis.tick_params(axis='y', length=0)
    axis.set_ylabel(label)
    for spine in axis.spines.values():
        spine.set_position(('outward', 5))
    axis.set_axisbelow(True)

if __name__ == "__main__":
    # args: file name
    if len(sys.argv) != 2:
        sys.exit("Error: wrong number of arguments\n" +
                 "Usage: python multirunFitness filename")

    dat = read_logfile(sys.argv[1] + '-0.log')
    dat2 = read_logfile(sys.argv[1] + '-1.log')
    allData = [x + y for x, y in zip(dat, dat2)]
    bmap = brewer2mpl.get_map('Set2', 'qualitative', 8)
    colors = bmap.mpl_colors

    nb_graphs = len(allData[0]) + 1
    labels = ["Fitness", "D2Item", "D2Robot"] * 2 + ['SumFitness']

    figure(1)
    axis = []
    for i in xrange(nb_graphs):
        axis.append(subplot2grid((nb_graphs, 1), (i, 0)))
    idx_graph = 0
    for c in zip(*allData):
        plot_curve(c, colors[idx_graph + 1], axis[idx_graph],
                   labels[idx_graph], True)
        idx_graph += 1

    sumFitness = [x + y for x, y in zip(zip(*dat)[0], zip(*dat2)[0])]

    plot_curve(sumFitness, colors[idx_graph + 1], axis[idx_graph],
               labels[idx_graph], True)
    '''
    figure(2)
    axis2 = []
    for i in xrange(nb_graphs - 1):
        axis2.append(subplot2grid((nb_graphs, 1), (i, 0)))
    idx_graph = 0
    for c in zip(*allData):
        plot_curve(c[35000:40000], colors[idx_graph + 1], axis2[idx_graph],
                   labels[idx_graph], True)
        idx_graph += 1

    plt.show()
    '''
