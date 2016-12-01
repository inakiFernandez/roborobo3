# -*- coding: utf-8 -*-
"""
script plot a curve per column
"""
import numpy as np
import sys
from pylab import *
import brewer2mpl
#from pprint import pprint


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


def plot_boxplot(stats, colors, axis, labels, sig=False):

    bp = axis.boxplot(stats)

    for i in range(0, len(bp['boxes'])):
        bp['boxes'][i].set_color(colors[i])
        # we have two whiskers!
        #bp['whiskers'][i*2].set_color(colors[i])
        #bp['whiskers'][i*2 + 1].set_color(colors[i])
        #bp['whiskers'][i*2].set_linewidth(2)
        #bp['whiskers'][i*2 + 1].set_linewidth(2)
        # top and bottom fliers
        # (set allows us to set many parameters at once)
        #bp['fliers'][i * 2].set(markerfacecolor=colors[i],
                                #marker='o', alpha=0.75, markersize=6,
                                #markeredgecolor='none')
        #bp['fliers'][i * 2 + 1].set(markerfacecolor=colors[i],
                                   # marker='o', alpha=0.75, markersize=6,
                                   # markeredgecolor='none')
        #bp['medians'][i].set_color('black')
        #bp['medians'][i].set_linewidth(3)
        # and 4 caps to remove
        #for c in bp['caps']:
        #    c.set_linewidth(0)

    for i in range(len(bp['boxes'])):
        box = bp['boxes'][i]
        box.set_linewidth(0)
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
            boxCoords = zip(boxX, boxY)
            boxPolygon = Polygon(boxCoords, facecolor=colors[i], linewidth=0)
            axis.add_patch(boxPolygon)

    axis.set_xticklabels(labels)
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.spines['left'].set_visible(False)
    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()
    axis.tick_params(axis='x', direction='out')
    axis.tick_params(axis='y', length=0)
    axis.grid(axis='y', color="0.9", linestyle='-', linewidth=1)
    axis.set_axisbelow(True)

    '''
    # stat test
    if sig :
        for i in xrange(len(stats)) :
            for j in xrange(i+1,len(stats)) :
                y_max = max(concatenate((stats[i], stats[j])))
                y_min = min(concatenate((stats[i], stats[j])))
                z,p = stat_test(stats[i], stats[j])

                axis.annotate("", xy=(i+1, y_max), xycoords='data',
                             xytext=(j+1, y_max), textcoords='data',
                             arrowprops=dict(arrowstyle="-", ec='#aaaaaa',
                                             connectionstyle="bar,
                                             fraction=0.1"))
                axis.text((j-i)/2.0 + j, y_max + abs(y_max - y_min)*0.1,
                          stars(p*2.0), horizontalalignment='center',
                         verticalalignment='center')
    '''


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
    if len(sys.argv) != 5:
        sys.exit("Error: wrong number of arguments\n" +
                 "Usage: python statsOnlineER filename filename2 \
                 backPerc targetFitPerc")

    dat = read_logfile(sys.argv[1])
    dat2 = read_logfile(sys.argv[2])
    maxFit = np.max(dat + dat2)
    targetFit = float(sys.argv[4]) * maxFit
    lastEpoch = len(dat)
    backEpoch = int(np.floor(float(sys.argv[3]) * lastEpoch))

    avgCum = []
    fixedBudgetFit = []
    timeToTarget = []
    accOverTgt = []
    for c in zip(*dat):
        avgCum.append(np.mean(c[backEpoch:]))
        fixedBudgetFit.append(c[backEpoch])
        ttT = lastEpoch
        found = False
        aOT = 0.0
        for i in xrange(len(c)):
            if c[i] >= targetFit:
                if not(found):
                    ttT = i
                    found = True
                aOT += c[i] - targetFit

        timeToTarget.append(ttT)
        accOverTgt.append(aOT)
    stats1 = [avgCum, fixedBudgetFit, timeToTarget, accOverTgt]

    avgCum = []
    fixedBudgetFit = []
    timeToTarget = []
    accOverTgt = []
    for c in zip(*dat2):
        avgCum.append(np.mean(c[backEpoch:]))
        fixedBudgetFit.append(c[backEpoch])
        ttT = lastEpoch
        found = False
        aOT = 0.0
        for i in xrange(len(c)):
            if c[i] >= targetFit:
                if not(found):
                    ttT = i
                    found = True
                aOT += c[i] - targetFit

        timeToTarget.append(ttT)
        accOverTgt.append(aOT)

    stats2 = [avgCum, fixedBudgetFit, timeToTarget, accOverTgt]
    nb_graphs = 4
    figure(num=None, figsize=(10, 5), dpi=100)
    clf()
    axis = []
    row = 0
    column = 0
    for i in xrange(nb_graphs):
        axis.append(subplot2grid((nb_graphs/2, 2), (row, column)))
        column += 1
        if column > 1:
            column = 0
            row += 1

    bmap = brewer2mpl.get_map('Set2', 'qualitative', 8)
    colors = bmap.mpl_colors
    labels = ["Avg. Accum. Fitness", "Fixed Budget Fitness",
              "Time to Target", "Fitness Over Target"]

    for i in xrange(nb_graphs):
        stats = []
        stats.append(stats1[i])
        stats.append(stats2[i])
        plot_boxplot(stats, colors, axis[i], '', sig=False)


'''
    figure(1)
    for c in zip(*allData):
        plot_curve(c, colors, axis[idx_graph], labels[idx_graph], True)
'''
