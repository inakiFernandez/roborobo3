"""
fnames += [
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/multiAverageAll.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/multiAvgDepth.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/multiBestAll.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/multiDepthBest.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/multiLinksBest.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/multiNodesBest.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/multiOtherAverageAll.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/multiOtherBestAll.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/multiSpecies.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/multiSpeciesSize.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/noMultiAverageAll.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/noMultiAvgDepth.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/noMultiBestAll.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/noMultiDepthBest.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/noMultiLinksBest.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/noMultiNodesBest.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/noMultiOtherAverageAll.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/noMultiOtherBestAll.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/noMultiSpecies.log",
"/home/fernandi/git/roborobo3/learnFunction/NEAT/expMultiP50TogetherNumberTaskProbMutW.6/noMultiSpeciesSize.log",
]
"""
import sys
#sys.path.insert(0, '../tools/')
import matplotlib 
matplotlib.use('svg')
import multirunFitness
import argparse
import matplotlib.pyplot as plt
import os, time
import numpy as np
import brewer2mpl
import importlib
############################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='To png or not')
    
    parser.add_argument('outfile', help='output filename for png')
    parser.add_argument('nbExp', type=int,
                        help='Number of variants of the experiments')
    parser.add_argument('--png', action='store_true',
                    help='output to png file?')
    parser.add_argument('itPerTask', type=int,
                        help='Number of iterations per task')                
    parser.add_argument('runId', type=int,
                        help='ID of the current run')                
    args = parser.parse_args()

    isToPng = args.png

    outfile = args.outfile
    nbExp =  args.nbExp
    runId = args.runId
    filenames = importlib.import_module("tmpDataFilesExp" + str(runId)) 
    # import tmpDataFilesExp as filenames    
    fnames = filenames.fnames

    labels = ["".join(n.split('/')[-1].split('.')[:-1]) for n in fnames]  #file name without .log
    nameExp = fnames[0].split('/')[-2]
    #if isToPng:
    #   plt.ioff()

    pixw, pixh = 2100, 1205
    dpi = 100
    winch = pixw/dpi # width in inches is pixel width / dpi
    hinch = pixh/dpi
    bgcolor="gainsboro"
    figurePlt = plt.figure(nameExp,figsize=(winch, hinch), dpi=dpi,facecolor=bgcolor)
    
    data = []
    xlength = len(multirunFitness.read_logfile(fnames[0]))

    interval = args.itPerTask # 50

    doIntervals = (xlength >= interval)
    bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
    colors = bmap.mpl_colors
    axis = plt.subplot2grid((1, 1), (0, 0),facecolor=bgcolor) # (2, 1), (0, 0))
    print(fnames)
    for i in range(len(fnames)):
        if((i % (len(fnames) / nbExp)) in [1,3,7,8]): #depends on alphabetical order of files
            multirunFitness.plot_one_curve(multirunFitness.read_logfile(fnames[i]), 
                                           colors[(i)%len(colors)], 
                                           axis, labels[i], True)
    
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    axis.legend(loc='upper left', fontsize = 20,ncol=4)
    axis.set_ylim(ymin=-0.01)
    if doIntervals:
        multirunFitness.taskIntervals(xlength,interval)
    """
    #titleAxis = 'Random ANN: multi vs no multi. Same Samples. Mut(n0.1,l0.2,tog0.5,s0.1,indiv1.0)'
    titleAxis = 'Best and Ave fitness NEAT 2 random functions w/wo multisynapses'
    axis.set_title(titleAxis)
    axis.set_ylim([0,1.0])

    axis2 = subplot2grid((3, 1), (1, 0)) #axis.twinx()

    for i in range(len(fnames)):
        if(not((i % (len(fnames) / nbExp)) in [1,3,7,8,9,10,0])): #depends on alphabetical order of files
            plot_one_curve(read_logfile(fnames[i]), colors[(i+1)%len(colors)], axis2, labels[i], True)

    axis2.legend(loc='upper left')
    taskIntervals(xlength,interval)

    axis3 = subplot2grid((3, 1), (2, 0)) #axis.twinx()
    for i in range(len(fnames)):
        if((i % (len(fnames) / nbExp)) in [9,10]): #depends on alphabetical order of files
            plot_one_curve(read_logfile(fnames[i]),
                           colors[(i+1)%len(colors)], axis3, labels[i], True)

    axis3.legend(loc='upper right',ncol=5)
    taskIntervals(xlength,interval)
    """
    #manager = plt.get_current_fig_manager()
    #manager.window.showMaximized() # TODO set fixed resolution (not depending on screen)

    if isToPng:
        #plt.show() # TODO programmatically close window
        print("Saved to ", outfile)
        time.sleep(2)
        figurePlt.savefig(outfile, dpi=dpi)
        plt.close(figurePlt)
        #plt.close('all')
    else:
        manager.window.showMaximized()
        plt.show()

    fnames = []

