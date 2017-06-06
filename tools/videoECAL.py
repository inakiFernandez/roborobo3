# -*- coding: utf-8 -*-
"""
Created on Mon May 29 11:41:41 2017
Fitness curves for video ECAL
@author: fernandi
"""
import matplotlib
matplotlib.use('svg')
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, '../../tools/')
import multirunFitness
import numpy as np
from os.path import expanduser
home = expanduser("~")

fnameFit = home + "/git/roborobo3/logs/stdoutInfo.log"
fnameColorItems = home + "/git/roborobo3/logs/items.log"

fitness = multirunFitness.read_logfile(fnameFit)
colorItems = multirunFitness.read_logfile(fnameColorItems)

bgcolor = "white" # "gainsboro"
maxGen = len(fitness)
maxGen2 = len(colorItems)

for i in range(maxGen - 1):
    fig = plt.figure("Cooperative Foraging",figsize=(21, 12.05),dpi=100,facecolor=bgcolor)
    axis = plt.subplot2grid((1, 1), (0, 0),facecolor=bgcolor)
    multirunFitness.plot_one_curve(fitness[0:i+1], "blue", axis, "Swarm Fitness", True)
    
    plt.ylabel("#Items")
    plt.xlabel("Generations")
    plt.xticks(fontsize = 35)
    plt.yticks(fontsize = 35)
    axis.legend(loc='upper left', fontsize=40,ncol=4)
    axis.set_axis_on()    
    axis.set_ylim(ymin=-0.01,ymax=np.max(fitness))
    axis.set_xlim(xmin=-0.01,xmax=200)
    plt.ylabel("Number of Items",fontsize = 40)
    plt.xlabel("Generations",fontsize = 40)
    plt.title("Collected items over time",fontsize=50) 

    fig.savefig(home + "/git/roborobo3/logs/imgVideoECAL/fitness" + str(i) + ".png", dpi=100)
    plt.close(fig)
    

colorItemsAggreg = []
print(len(colorItems))
for vectorColors in colorItems:
    print("At iteration:")    
    print(vectorColors)
    vectorAggreg = []
    aggreg = 0.0
    for c in vectorColors:        
        aggreg = aggreg + c
        vectorAggreg.append(aggreg)
    print(c)
    print(vectorAggreg)
    colorItemsAggreg.append(vectorAggreg)

print(colorItemsAggreg[0:3])    
colorArr = [(0.0,0.0,1.0),(0.0,145.0/255.0,1.0),(0.0,1.0,218.0/255.0),(0.0,1.0,72.0/255.0),
            (72.0/255.0,1.0,0.0),(218.0/255.0,1.0,0.0),(1.0,145.0/255,0.0),(1.0,0.0,0.0)]

colorItemsAggregInv = [list(x) for x in zip(*colorItemsAggreg)]             

for i in range(maxGen - 1):
#for i in range(90,95):
    fig = plt.figure("Cooperative Foraging",figsize=(21, 12.05),dpi=100,facecolor=bgcolor)
    axis = plt.subplot2grid((1, 1), (0, 0),facecolor=bgcolor) 
    for j in range(len(colorArr)):
        if j == 0:            
            #print(len(np.arange(0,  i+1)), ", ", len(np.zeros(i+1)), ", ", len(np.array(colorItemsAggregInv[j][:i+1])))            
            axis.fill_between(np.arange(0, i+1), np.zeros(i+1) , colorItemsAggregInv[j][:i+1],
                              alpha=0.55, linewidth=0, color=colorArr[j])
        else:
            print("From previous")
            print(len(np.arange(0,  i+1)), ", ", len(colorItemsAggregInv[j-1][:i+1]), ", ", len(colorItemsAggregInv[j][:i+1]))
            print(colorItemsAggregInv[j-1][:i+1])
            print()
            print(colorItemsAggregInv[j][:i+1])
            axis.fill_between(np.arange(0,  i+1), colorItemsAggregInv[j-1][:i+1], colorItemsAggregInv[j][:i+1],
                              alpha=0.55, linewidth=0, color=colorArr[j])
    plt.ylabel("Number of Items",fontsize = 40)
    plt.xlabel("Generations",fontsize = 40)
    plt.xticks(fontsize = 35)
    plt.yticks(fontsize = 35)
    axis.legend(loc='upper left', fontsize = 40,ncol=4)
    axis.set_ylim(ymin=-0.01,ymax=np.max(fitness))
    axis.set_xlim(xmin=-0.01,xmax=200)
    plt.title("Collected items per color over time",fontsize=50) 
    fig.savefig(home + "/git/roborobo3/logs/imgVideoECAL/fitnessPerColor" + str(i) + ".png", dpi=100)
    plt.close(fig)
    
  
        
       

        itemsPerColorRunB = [np.sum([float(y) for y in x]) for x in itemsPerColorRunB]        
        itemsPerColorAllB.append(itemsPerColorRunB)        
        
        contentIncrB = []
        for it in contentSplitB:
            tmp = []
            aggreg = 0.0
            for val in it:
                aggreg += float(val)
                tmp += [aggreg]
            contentIncrB += [tmp]
        
