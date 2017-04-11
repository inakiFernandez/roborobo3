# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 10:25:48 2017

@author: fernandi
"""
import matplotlib 
matplotlib.use('svg')
import matplotlib.pyplot as plt
import os, time
import numpy as np
import brewer2mpl
import importlib
#import multirunFitness
import argparse
import glob

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='To png or not')
    
    parser.add_argument('itemsPerColorFile', help='input filename items per color')
    parser.add_argument('itemsGatheredAndPossibleFile', help='input filename items gathered and possible')
    parser.add_argument('numberChangesFile', help='input filename number of changes')
    parser.add_argument('outfile', help='output basename for png')
    #parser.add_argument('nbExp', type=int, help='Number of variants of the experiments')
    
    parser.add_argument('--png', action='store_true', help='output to png file?')
    #parser.add_argument('itPerTask', type=int, help='Number of iterations per task')                
    #parser.add_argument('runId', type=int, help='ID of the current run')                
    args = parser.parse_args()

    isToPng = args.png
    #nbExp =  args.nbExp
    outfile = args.outfile    
    #runId = args.runId
    dpi = 100
    
    datapath = os.path.dirname(os.path.abspath(args.itemsPerColorFile))
    dirFiles = glob.glob(datapath + "/items-run*.log") # os.listdir(datapath)    
    #dirFiles.sort(key=lambda f: int(f.split('Distance')[1].split('.')[0]))    
    
    print(dirFiles)
    #quit()
    
    nColors = 8
    colorArr = []
    for j in range(nColors):
        absVal = (256/nColors) * j;
        colorArr += [(absVal/256.0, (255 - absVal)/256.0, 0)]
    
    #change colors to rainbow [], .append
    colorArr = [(0.0,0.0,1.0),(0.0,145.0/255.0,1.0),(0.0,1.0,218.0/255.0),(0.0,1.0,72.0/255.0),
                (72.0/255.0,1.0,0.0),(218.0/255.0,1.0,0.0),(1.0,145.0/255,0.0),(1.0,0.0,0.0)]
    
    fname =  '/home/fernandi/git/roborobo3/logs/items.log'
    fname = args.itemsPerColorFile
    
    with open(fname) as f:
        content = f.readlines()
    
    content = [x.strip() for x in content]
    
    contentSplit = [x.split(" ") for x in content][1:]
    
    contentIncr = []
    for it in contentSplit:
        tmp = []
        aggreg = 0.0
        for val in it:
            aggreg += float(val)
            tmp += [aggreg]
        contentIncr += [tmp]
    
    contentIncr= [list(x) for x in zip(*contentIncr)] # contentSplit)]
    
    bgcolor="gainsboro"
    fig1 = plt.figure(1)
    axis1 = plt.subplot2grid((1, 1), (0, 0), facecolor=bgcolor)
    for i,lColor in enumerate(contentIncr):
        #axis1.plot(lColor, color = colorArr[i])
        if i == 0:
            axis1.fill_between(np.arange(0, len(lColor)), np.zeros(len(lColor)) ,  lColor,
                               alpha=0.55, linewidth=0, color=colorArr[i])
        else:
            axis1.fill_between(np.arange(0, len(lColor)), contentIncr[i-1],  lColor,
                               alpha=0.55, linewidth=0, color=colorArr[i])
        #axis1.plot(lColor, color = colorArr[i], linewidth = 0.2)           
                
    ###############################################################
    nbItems = 100 # 40 # 80 # 200 # 30
    nbSteps = 800 #500
    fname = '/home/fernandi/git/roborobo3/logs/itemsIter.log'
    fname = args.itemsGatheredAndPossibleFile
    with open(fname) as f:
        content = f.readlines()
    
    content = [x.strip() for x in content]
    
    contentSplit = [x.split(" ") for x in content]
                    
    #contentIncr= [list(x) for x in  zip(*contentSplit)]                         
    
    valuesItems = [ [(float(x[0]) / float(x[1])),  # / float(nbItems), 
                     float(x[1]) / float(nbItems)] if(x[1] != '0') else [0.0, 0.0] for x in contentSplit]
    
    valuesItems= [list(x) for x in zip(*valuesItems)]
    
    valuesItemsAggGeneration = [[],[]]
    for i,x in enumerate(valuesItems):    
        nbGen = int(len(x) / nbSteps)
        aggPerIter = [np.average(x[nbSteps * i:nbSteps * (i + 1)]) for i in range(nbGen)]
        valuesItemsAggGeneration[i] = aggPerIter
    
    bgcolor="gainsboro"
    fig2 = plt.figure(2)
    axis2 = plt.subplot2grid((1, 1), (0, 0), facecolor=bgcolor)
    
    axis2.plot(valuesItems[1], color = "blue", linewidth = 0.2)    
    axis2.plot(valuesItems[0], color = "orange", linewidth = 0.2)    
    #axis2.fill_between(np.arange(0, len(valuesItems[0])), 
    #                   np.zeros(len(valuesItems[0])) ,  valuesItems[1],
    #                   alpha=0.35, linewidth=0, color="blue")
    #axis2.fill_between(np.arange(0, len(valuesItems[0])), 
    #                   np.zeros(len(valuesItems[0])) ,  valuesItems[0],
    #                   alpha=0.95, linewidth=0, color="red")
    
    bgcolor="gainsboro"
    fig21 = plt.figure(21)
    axis21 = plt.subplot2grid((1, 1), (0, 0), facecolor=bgcolor)
    
    axis21.plot(valuesItemsAggGeneration[1], color = "blue", linewidth = 0.7)    
    axis21.plot(valuesItemsAggGeneration[0], color = "orange", linewidth = 0.7)    
    
                       
    ##############################################################
    
    fname = '/home/fernandi/git/roborobo3/logs/colorChanges.log'
    fname = args.numberChangesFile
    with open(fname) as f:
        content = f.readlines()
    
    content = [x.strip() for x in content]
    
    contentSplit = [x.split(" ") for x in content]         
    
    contentAgg = [np.average(np.array([int(y) for y in x])) for x in contentSplit] # Only if measured per generation /500.0
    bgcolor="gainsboro"
    fig3 = plt.figure(3)
    axis3 = plt.subplot2grid((1, 1), (0, 0), facecolor=bgcolor)
    
    axis3.plot(contentAgg, color = "blue", linewidth = 0.2) 
    
    if isToPng:        
        print("Saved to ", outfile)
        time.sleep(2)
        fig1.savefig(outfile + "1.png"  , dpi=dpi)        
        plt.close(fig1)
        fig2.savefig(outfile + "2.png"  , dpi=dpi)        
        plt.close(fig2)
        fig21.savefig(outfile + "21.png"  , dpi=dpi)        
        plt.close(fig21)
        fig3.savefig(outfile + "3.png"  , dpi=dpi)        
        plt.close(fig3)
        #plt.close('all')     
